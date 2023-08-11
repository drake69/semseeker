loadNamespace("future")
loadNamespace("stats")
nitem <- 4e3
nsamples <- 21

probe_features <- semseeker::PROBES
probe_features <- probe_features[!is.na(probe_features$START),c("CHR","START","PROBE")]
probe_features <- unique(probe_features)
probe_features$END <- probe_features$START

nitem <- min(nitem, nrow(probe_features))
probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]
probe_features$ABSOLUTE <- paste(probe_features$CHR, probe_features$START, sep="_")

methylation_data <- stats::rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))

row.names(methylation_data) <- probe_features$PROBE

Sample_ID <- stringi::stri_rand_strings(nsamples, 15, pattern = "[A-Za-z]")
colnames(methylation_data) <- Sample_ID
Sample_Group <- c(rep("Control",nsamples/3),rep("Case",nsamples/3),rep("Reference",nsamples/3))
mySampleSheet <- data.frame(Sample_Group, Sample_ID)

mySampleSheet$Phenotest <- stats::rnorm(nsamples, mean= 1000, sd= 567)
mySampleSheet$Group <- c(rep(TRUE,ceiling(nsamples/2)), rep(FALSE,floor(nsamples/2)))
mySampleSheet$Covariates1 <- stats::rnorm(nsamples, mean= 567, sd= 1000)
mySampleSheet$Covariates2 <- stats::rnorm(nsamples, mean= 67, sd= 100)

mySampleSheet_batch <<- list(mySampleSheet, mySampleSheet, mySampleSheet)
methylation_data_batch <<- list(methylation_data, methylation_data, methylation_data)

beta_superior_thresholds <- data.frame(stats::rnorm(nitem, mean = 1, sd=0.2))
beta_inferior_thresholds <- data.frame(stats::rnorm(nitem, mean=0.2, sd=0.2))
colnames(beta_superior_thresholds) <- "HIGH"
colnames(beta_inferior_thresholds) <- "LOW"

iqr <- data.frame(stats::rnorm(nitem, mean=0.1, sd=0.3))
row.names(beta_superior_thresholds) <- probe_features$PROBE
row.names(beta_inferior_thresholds) <- probe_features$PROBE
beta_medians <- (beta_superior_thresholds + beta_inferior_thresholds) / 2

thresholds <- data.frame("thresholds"= stats::rnorm(nitem, mean=0.5, sd= 0.5))
row.names(thresholds) <- probe_features$PROBE

#####
### copy as global

mySampleSheet <<- mySampleSheet

methylation_data <<- abs(methylation_data)
beta_medians <<- beta_medians
beta_inferior_thresholds <<- abs(beta_inferior_thresholds)
beta_superior_thresholds <<- abs(beta_superior_thresholds)
beta_superior_thresholds[beta_superior_thresholds<beta_inferior_thresholds,1] <- beta_inferior_thresholds[beta_superior_thresholds<beta_inferior_thresholds,1] + 0.1
thresholds <<- abs(thresholds)
beta_thresholds <- data.frame("beta_median_values"=beta_medians,
                              "beta_inferior_thresholds"=beta_inferior_thresholds,
                              "beta_superior_thresholds"=beta_superior_thresholds,
                              "iqr" = iqr)
colnames(beta_thresholds) <- c("beta_median_values","beta_inferior_thresholds","beta_superior_thresholds","iqr")
sliding_window_size <<- 11
bonferroni_threshold <<- 0.1
batch_id <<- 1
iqrTimes <<- 3
# multisession
# sequential
parallel_strategy <<- "sequential"

# TODO
# recover session stored
#
tmp <- tempdir()
tempFolders <<- paste(tmp,"/semseeker/",stringi::stri_rand_strings(50, 7, pattern = "[A-Za-z0-9]"),sep="")
