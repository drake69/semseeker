loadNamespace("future")
loadNamespace("stats")
nprobes <- 6e3
nsamples <- 240

# intersample_mean <- 5
# intersample_sd <- 1
# perc_epimutation <- 0.05

probe_features <- semseeker::PROBES
probe_features <- probe_features[!is.na(probe_features$START),c("CHR","START","PROBE")]
probe_features <- unique(probe_features)
probe_features$END <- probe_features$START

nprobes <- min(nprobes, nrow(probe_features))
probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nprobes),]
probe_features$ABSOLUTE <- paste(probe_features$CHR, probe_features$START, sep="_")

# test <- stats::rnorm(nsamples,mean = intersample_mean, sd = intersample_sd)
# mean(test)
# plot(density(as.numeric(test)))

set.seed(474693)
nCpG <- nprobes * 0.9
nTypeI <- 90
nTypeII <- 10
Meth <- as.data.frame(matrix(rbeta(nCpG * nsamples, nTypeI, nTypeII), nrow = nCpG))

nCpG <- nprobes * 0.1
nTypeI2 <- 4
nTypeII2 <- 6
Unmeth <- as.data.frame(matrix(rbeta(nCpG * nsamples, nTypeI2, nTypeII2), nrow = nCpG))

signal_data <- rbind(Meth, Unmeth)
# signal_data <- rep(stats::rnorm(nprobes,mean = intersample_mean, sd = intersample_sd),nsamples)

# plot(density(as.numeric(signal_data)))

q1 <-stats::quantile(signal_data, probs = 0.25, na.rm = TRUE)
q3 <-stats::quantile(signal_data, probs = 0.75, na.rm = TRUE)

# noise_position <- unique(floor(runif(nsamples*nprobes*perc_epimutation*2, min = 1, max = nsamples*nprobes)))
# noise_position <- noise_position[1:(nsamples*nprobes*perc_epimutation)]
# min(noise_position)
# max(noise_position)
# median(noise_position)
# noise <- c(runif(nsamples*nprobes*perc_epimutation/2,min = q3,  max= q3 + 2*(abs(q3)) ),runif(nsamples*nprobes*perc_epimutation/2,min = q1-2*(abs(q1)),max=q1))
# signal_data[noise_position] <- noise
# min(noise)
# max(noise)
# signal_data <- signal_data + abs(min(noise))
# min(signal_data)
# signal_data <- as.data.frame(matrix(data = signal_data,nrow = nprobes,ncol = nsamples, byrow = TRUE))
# plot(density(as.numeric(signal_data[1,])))

row.names(signal_data) <- probe_features$PROBE

Sample_ID <- stringi::stri_rand_strings(nsamples, 15, pattern = "[A-Za-z]")
colnames(signal_data) <- Sample_ID
Sample_Group <- c(rep("Control",nsamples/3),rep("Case",nsamples/3),rep("Reference",nsamples/3))
mySampleSheet <- data.frame(Sample_Group, Sample_ID)

mySampleSheet$Phenotest <- stats::rnorm(nsamples, mean= 1000, sd= 567)
mySampleSheet$Group <- c(rep(TRUE,ceiling(nsamples/2)), rep(FALSE,floor(nsamples/2)))
mySampleSheet$Covariates1 <- stats::rnorm(nsamples, mean= 567, sd= 1000)
mySampleSheet$Covariates2 <- stats::rnorm(nsamples, mean= 67, sd= 100)

mySampleSheet_batch <<- list(mySampleSheet, mySampleSheet, mySampleSheet)
signal_data_batch <<- list(signal_data, signal_data, signal_data)

q1 <- apply(signal_data, 1, function (x) { stats::quantile(x, probs = 0.25, na.rm = TRUE) })
q3 <- apply(signal_data, 1, function (x) { stats::quantile(x, probs = 0.75, na.rm = TRUE) })
iqr <- data.frame(q3 - q1)

signal_superior_thresholds <- data.frame("HIGH" = q3 + 3 * iqr)
signal_inferior_thresholds <- data.frame("LOW" = q1 - 3 * iqr)
colnames(signal_inferior_thresholds) <- "LOW"
colnames(signal_superior_thresholds) <- "HIGH"

# sum(as.numeric(signal_data[,1] > signal_superior_thresholds))
# sum(as.numeric(signal_data[,1] < signal_inferior_thresholds))

row.names(signal_superior_thresholds) <- probe_features$PROBE
row.names(signal_inferior_thresholds) <- probe_features$PROBE
signal_medians <- apply(signal_data, 1, median)

signal_thresholds <- data.frame("signal_median_values"=signal_medians,
  "signal_inferior_thresholds"=signal_inferior_thresholds,
  "signal_superior_thresholds"=signal_superior_thresholds,
  "iqr" = iqr,
  "q1"=q1,
  "q3"=q3)

colnames(signal_thresholds) <- c("signal_median_values","signal_inferior_thresholds","signal_superior_thresholds","iqr","q1","q3")

#####
### copy as global
mySampleSheet <<- mySampleSheet
signal_data <<- (signal_data)
signal_medians <<- signal_medians
signal_inferior_thresholds <<- (signal_inferior_thresholds)
signal_superior_thresholds <<- (signal_superior_thresholds)
signal_thresholds <<- signal_thresholds
sliding_window_size <<- 11
bonferroni_threshold <<- 0.1
batch_id <<- 1
iqrTimes <<- 3
# perc_epimutation <<- perc_epimutation
nsamples <<- nsamples
nprobes <<- nprobes
# multisession
# sequential
parallel_strategy <<- "multicore"
probe_features <<- probe_features
markers <<- c("MUTATIONS","DELTAQ")


# TODO
# recover session stored
#
tmp <- tempdir()
tempFolders <<- paste(tmp,"/semseeker/",stringi::stri_rand_strings(50, 7, pattern = "[A-Za-z0-9]"),sep="")
