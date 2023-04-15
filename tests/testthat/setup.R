
nitem <- 4e4
nsamples <- 21

probes <- semseeker::PROBES
probe_features <- probes[!is.na(probes$START),c("CHR","START","PROBE")]
probe_features <- unique(probe_features)
probe_features$END <- probe_features$START

nitem <- min(nitem, nrow(probe_features))
probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]
probe_features$ABSOLUTE <- paste(probe_features$CHR, probe_features$START, sep="_")

methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))

row.names(methylation_data) <- probe_features$PROBE

Sample_ID <- stringi::stri_rand_strings(nsamples, 15, pattern = "[A-Za-z]")
colnames(methylation_data) <- Sample_ID
Sample_Group <- c(rep("Control",nsamples/3),rep("Case",nsamples/3),rep("Reference",nsamples/3))
mySampleSheet <- data.frame(Sample_Group, Sample_ID)

mySampleSheet$Phenotest <- rnorm(nsamples, mean= 1000, sd= 567)
mySampleSheet$Group <- c(rep(TRUE,ceiling(nsamples/2)), rep(FALSE,floor(nsamples/2)))
mySampleSheet$Covariates1 <- rnorm(nsamples, mean= 567, sd= 1000)
mySampleSheet$Covariates2 <- rnorm(nsamples, mean= 67, sd= 100)

mySampleSheet_batch <<- list(mySampleSheet, mySampleSheet, mySampleSheet)
methylation_data_batch <<- list(methylation_data, methylation_data, methylation_data)

beta_superior_thresholds <- data.frame(rnorm(nitem, mean = 1, sd=0.2))
beta_inferior_thresholds <- data.frame(rnorm(nitem, mean=0.2, sd=0.2))
row.names(beta_superior_thresholds) <- probe_features$PROBE
row.names(beta_inferior_thresholds) <- probe_features$PROBE
beta_medians <- (beta_superior_thresholds + beta_inferior_thresholds) / 2

tresholds <- data.frame("tresholds"= rnorm(nitem, mean=0.5, sd= 0.5))
row.names(tresholds) <- probe_features$PROBE

#####
### copy as global

mySampleSheet <<- mySampleSheet

methylation_data <<- methylation_data
beta_medians <<- beta_medians
beta_inferior_thresholds <<- beta_inferior_thresholds
beta_superior_thresholds <<- beta_superior_thresholds
tresholds <<- tresholds

sliding_window_size <<- 11
bonferroni_threshold <<- 0.1
batch_id <<- 1
iqrTimes <<- 3

parallel_strategy <<- "sequential"

