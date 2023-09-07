loadNamespace("future")
loadNamespace("stats")
nprobes <- 4e3
nsamples <- 240

intersample_mean <- 5
intersample_sd <- 1
perc_epimutation <- 0.05

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

methylation_data <- rep(stats::rnorm(nprobes,mean = intersample_mean, sd = intersample_sd),nsamples)
plot(density(as.numeric(methylation_data)))

q1 <-stats::quantile(methylation_data, probs = 0.25, na.rm = TRUE)
q3 <-stats::quantile(methylation_data, probs = 0.75, na.rm = TRUE)

noise_position <- unique(floor(runif(nsamples*nprobes*perc_epimutation*2, min = 1, max = nsamples*nprobes)))
noise_position <- noise_position[1:(nsamples*nprobes*perc_epimutation)]
min(noise_position)
max(noise_position)
median(noise_position)

noise <- c(runif(nsamples*nprobes*perc_epimutation/2,min = q3,  max= q3 + 2*(abs(q3)) ),runif(nsamples*nprobes*perc_epimutation/2,min = q1-2*(abs(q1)),max=q1))

methylation_data[noise_position] <- noise
# min(noise)
# max(noise)
methylation_data <- methylation_data + abs(min(noise))
# min(methylation_data)

methylation_data <- as.data.frame(matrix(data = methylation_data,nrow = nprobes,ncol = nsamples, byrow = TRUE))
plot(density(as.numeric(methylation_data[1,])))


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

q1 <- apply(methylation_data, 1, function (x) { stats::quantile(x, probs = 0.25, na.rm = TRUE) })
q3 <- apply(methylation_data, 1, function (x) { stats::quantile(x, probs = 0.75, na.rm = TRUE) })
iqr <- data.frame(q3 - q1)

beta_superior_thresholds <- data.frame("HIGH" = q3 + 3 * iqr)
beta_inferior_thresholds <- data.frame("LOW" = q1 - 3 * iqr)
colnames(beta_inferior_thresholds) <- "LOW"
colnames(beta_superior_thresholds) <- "HIGH"

# sum(as.numeric(methylation_data[,1] > beta_superior_thresholds))
# sum(as.numeric(methylation_data[,1] < beta_inferior_thresholds))

row.names(beta_superior_thresholds) <- probe_features$PROBE
row.names(beta_inferior_thresholds) <- probe_features$PROBE
beta_medians <- apply(methylation_data, 1, median)

beta_thresholds <- data.frame("beta_median_values"=beta_medians,
  "beta_inferior_thresholds"=beta_inferior_thresholds,
  "beta_superior_thresholds"=beta_superior_thresholds,
  "iqr" = iqr,
  "q1"=q1,
  "q3"=q3)

colnames(beta_thresholds) <- c("beta_median_values","beta_inferior_thresholds","beta_superior_thresholds","iqr","q1","q3")

#####
### copy as global
mySampleSheet <<- mySampleSheet
methylation_data <<- (methylation_data)
beta_medians <<- beta_medians
beta_inferior_thresholds <<- (beta_inferior_thresholds)
beta_superior_thresholds <<- (beta_superior_thresholds)
beta_thresholds <<- beta_thresholds
sliding_window_size <<- 11
bonferroni_threshold <<- 0.1
batch_id <<- 1
iqrTimes <<- 3
perc_epimutation <<- perc_epimutation
nsamples <<- nsamples
nprobes <<- nprobes
# multisession
# sequential
parallel_strategy <<- "sequential"
probe_features <<- probe_features
markers <<- c("MUTATIONS","DELTAQ")


# TODO
# recover session stored
#
tmp <- tempdir()
tempFolders <<- paste(tmp,"/semseeker/",stringi::stri_rand_strings(50, 7, pattern = "[A-Za-z0-9]"),sep="")
