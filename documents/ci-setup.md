# CI Setup — semseeker

This document describes the Continuous Integration pipelines for semseeker,
how to run them locally, and known issues that have been encountered and
resolved.

---

## Workflows

Two GitHub Actions workflows run on every push to `develop` and on every
pull request targeting `main`, `master`, or `develop`.

| Workflow | File | Purpose |
|---|---|---|
| `R-CMD-check` | `.github/workflows/R-CMD-check.yml` | Runs `R CMD check --as-cran` on macOS, Windows, Ubuntu |
| `test-coverage` | `.github/workflows/test-coverage.yml` | Runs `covr::package_coverage()` on macOS and uploads to Codecov |

---

## Platform Matrix

| Platform | Status | Notes |
|---|---|---|
| `macos-latest` | ✅ Active | Primary development platform; `brew install cairo fontconfig pkg-config` required |
| `windows-latest` | ✅ Active | No special system libs needed; path separator handled via `file.path()` |
| `ubuntu-latest` | ✅ Active | Pre-built binaries served by Posit Package Manager (P3M) for fast installs |

---

## Dependency Installation Strategy

Both workflows use `remotes::install_deps()` instead of `pak` because:

- `pak` resolves all dependencies before configuring extra repositories, so
  `polars` (from `community.r-multiverse.org`) cannot be found.
- `remotes` reads `Additional_repositories` from `DESCRIPTION` first, then
  resolves dependencies — polars installs correctly.

**Key parameter: `dependencies = NA`**

Both workflows pass `dependencies = NA` to `remotes::install_deps()`. This
installs only `Imports`, `Depends`, and `LinkingTo` — **not** `Suggests`.

Skipping Suggests is intentional:
- Suggests includes heavy Bioconductor packages (`GEOquery`, `pathfindR`,
  `ggkegg`) that are not needed for `R CMD check` or coverage measurement.
- Installing Suggests would require adding Bioconductor repositories and
  significantly increases CI time.
- All code paths that use Suggests packages are guarded with
  `requireNamespace("pkg", quietly = TRUE)` so they degrade gracefully when
  the package is absent (tests call `testthat::skip()` in that case).

---

## Known Issues

### 1. `pathfindR` → `ggkegg` Bioconductor chain (resolved)

**Symptom:** The `test-coverage` workflow failed with:

```
! dependency 'ggkegg' is not available
Error: running the tests in 'tastthat.R' failed
```

**Root cause:**

`pathfindR` is listed in semseeker's `Suggests`. `pathfindR` in turn
**requires** `ggkegg`, a Bioconductor package:

```
semseeker (Suggests) → pathfindR (Imports) → ggkegg (Bioconductor)
```

When the coverage workflow used `dependencies = TRUE`, `remotes` tried to
install all Suggests including `pathfindR`. Installing `pathfindR` required
`ggkegg`, which is not available on CRAN or `community.r-multiverse.org`.
The broken install caused tests to fail when `pathfindR` was loaded.

**Fix (commit `242a64f`):**

Changed `dependencies = TRUE` → `dependencies = NA` in
`.github/workflows/test-coverage.yml`:

```yaml
remotes::install_deps(
  dependencies = NA,   # Imports/Depends/LinkingTo only; skip Suggests
  repos = c(           # (pathfindR → ggkegg Bioconductor chain not resolvable from CRAN)
    rpolars = "https://community.r-multiverse.org",
    CRAN    = "https://packagemanager.posit.co/cran/latest"
  )
)
```

The test for `pathway_pathfindR()` already guards with
`requireNamespace("pathfindR", quietly = TRUE)` and calls
`testthat::skip()` when the package is absent, so coverage is not
materially affected.

**General rule:** Do not add Bioconductor packages to `Imports` or
`Depends`. If a Bioconductor package is needed as an optional feature,
add it to `Suggests` only and guard every call site with
`requireNamespace()`.

---

### 2. polars API breaking changes (0.22 → 0.23)

Several polars R API methods were renamed or removed in version 0.23:

| Old (≤ 0.22) | New (≥ 0.23) | Affected files |
|---|---|---|
| `df$to_data_frame()` | `as.data.frame(df)` | All R source + tests |
| `df$cast(list(...))` | `df$with_columns(pl$col(...)$cast(...))` | `annotate_position_pivots.R` |
| `df$rename(list(...))` | `df$select(pl$col(...)$alias(...))` | `create_position_pivots.R` |
| `$group_by(..., maintain_order=)` | `$group_by(..., .maintain_order=)` | `annotate_position_pivots.R`, `get_pivot_both.R` |

All occurrences have been migrated. polars is not pinned to a specific
version in `DESCRIPTION`; if a future polars release introduces further
breaking changes, search for these patterns first:

```bash
grep -r "to_data_frame\|maintain_order\|cast(list\|rename(list" R/ tests/
```

---

### 3. `arrow` undeclared dependency (resolved)

Several files called `arrow::read_parquet()` / `arrow::write_parquet()` but
`arrow` was not listed in `DESCRIPTION`. Under `R CMD check` this produces
an `ERROR`.

**Fix:** Replaced all `arrow` calls with polars equivalents:

```r
# Read
as.data.frame(polars::pl$read_parquet(path))

# Write
polars::as_polars_df(as.data.frame(df))$write_parquet(path)
```

`arrow` is no longer a dependency of semseeker.

---

### 4. `semseeker:::` prefix required in tests

Under `devtools::load_all()` all internal functions are visible globally, so
tests can call `my_internal_fn()` directly. Under `R CMD check` the package
is installed and only exported functions are visible; calling an internal
function without the `:::` prefix causes:

```
Error: could not find function "my_internal_fn"
```

**Rule:** Every call to a non-exported (internal) function in a test file
**must** use the `semseeker:::` prefix:

```r
# Wrong — works locally, fails under R CMD check:
result <- pivot_file_name_parquet("SIGNAL", "MEAN", "PROBE", "WHOLE")

# Correct:
result <- semseeker:::pivot_file_name_parquet("SIGNAL", "MEAN", "PROBE", "WHOLE")
```

---

### 5. macOS `/var` → `/private/var` symlink in path assertions

`dir_check_and_create()` calls `normalizePath()` internally, which on macOS
resolves the `/var` → `/private/var` symlink. Path equality assertions in
tests must do the same:

```r
# Wrong — fails on macOS:
testthat::expect_equal(ssEnv$result_folderData, file.path(tempFolder, "Data"))

# Correct:
testthat::expect_equal(
  ssEnv$result_folderData,
  normalizePath(file.path(tempFolder, "Data"), mustWork = FALSE)
)
```

---

## Running CI Locally (Docker)

A `Dockerfile.ci` and helper script `ci-local.sh` are provided for running
the full `R CMD check` locally without waiting for GitHub Actions.

### Prerequisites

- Docker Desktop running
- From the repo root: `chmod +x ci-local.sh`

### Commands

```bash
# First run — build image and run checks (~7 min first time, ~30 sec after)
./ci-local.sh check

# Rebuild the image from scratch (e.g. after DESCRIPTION changes)
./ci-local.sh rebuild

# Inspect the full .Rcheck output from the last run
./ci-local.sh logs

# Open an interactive shell inside the container
./ci-local.sh shell
```

### Output

Errors and warnings from the last `check` run are written to
`ci-check-output/semseeker.Rcheck/`. The most useful files:

| File | Contents |
|---|---|
| `00check.log` | Full `R CMD check` transcript |
| `tests/testthat.Rout` | testthat output (passing tests) |
| `tests/testthat.Rout.fail` | testthat output (only present when tests fail) |

### Docker layer cache

The `Dockerfile.ci` is structured so that R package installation is cached
in a separate layer from the source code:

```dockerfile
COPY DESCRIPTION .          # ← cache invalidated only when DESCRIPTION changes
RUN Rscript -e "..."        # ← heavy install step (~6 min), cached after first run
COPY . .                    # ← invalidated on every code change, but fast (~5 sec)
```

Subsequent builds after a code-only change take ~30 seconds.

---

## Adding a New Suggests Package

If you need to add an **optional** feature backed by a Bioconductor or
hard-to-install package:

1. Add the package to `Suggests:` in `DESCRIPTION` only (never `Imports:`).
2. Guard every call site:
   ```r
   if (!requireNamespace("mypkg", quietly = TRUE)) {
     message("mypkg is required for this feature. Install with: ...")
     return(invisible(NULL))
   }
   mypkg::some_function(...)
   ```
3. In the corresponding test:
   ```r
   test_that("feature works when mypkg is available", {
     if (!requireNamespace("mypkg", quietly = TRUE))
       testthat::skip("mypkg not installed")
     # ... test body
   })
   ```
4. Do **not** add the package's repository to the CI `repos` vector unless
   it is on CRAN. If it requires Bioconductor, document it here and note
   that CI will skip the test.
