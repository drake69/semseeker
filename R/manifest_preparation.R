#
#  temp <- tempfile()
#  utils::download.file("http://webdata.illumina.com.s3-website-us-east-1.amazonaws.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip", temp)
#
#  library(readr)
#  # infinium_methylationepic_v_1_0_b5_manifest_file <- readr::read_csv("~/Downloads/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip = 7)
#  infinium_methylationepic_v_1_0_b5_manifest_file <- read_csv("~/Desktop/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip = 7)
#  # PROBES <- utils::read.table(unz(temp, "infinium-methylationepic-v-1-0-b5-manifest-file.csv"))
#  # PROBES <- utils::read.table("~/Desktop/infinium-methylationepic-v-1-0-b5-manifest-file.csv")
#
#  PROBES <- data.frame(CHR = infinium_methylationepic_v_1_0_b5_manifest_file$CHR, START= infinium_methylationepic_v_1_0_b5_manifest_file$MAPINFO,
#                              END= infinium_methylationepic_v_1_0_b5_manifest_file$MAPINFO, PROBE = infinium_methylationepic_v_1_0_b5_manifest_file$IlmnID )
#
#  PROBES <- PROBES[with(PROBES, order(PROBE)), ]
#  row.names(PROBES) <- PROBES$PROBE
#
#  usethis::use_data(PROBES, overwrite = TRUE)
#
#  unlink(temp)
#  # download manifest if not existent
#
#  Regulatory_Feature_Name <- stringr::str_split(infinium_methylationepic_v_1_0_b5_manifest_file$Regulatory_Feature_Name, ";", simplify = TRUE)
#  Regulatory_Feature_Group <- stringr::str_split(infinium_methylationepic_v_1_0_b5_manifest_file$Regulatory_Feature_Group, ";", simplify = TRUE)
#
#
#  UCSC_RefGene_Accession <- stringr::str_split(infinium_methylationepic_v_1_0_b5_manifest_file$UCSC_RefGene_Accession, ";", simplify = TRUE)
#  UCSC_RefGene_Name <- stringr::str_split(infinium_methylationepic_v_1_0_b5_manifest_file$UCSC_RefGene_Name, ";", simplify = TRUE)
#  UCSC_RefGene_Group <- stringr::str_split(infinium_methylationepic_v_1_0_b5_manifest_file$UCSC_RefGene_Group, ";", simplify = TRUE)
#
#  Relation_to_UCSC_CpG_Island <- stringr::str_split(infinium_methylationepic_v_1_0_b5_manifest_file$Relation_to_UCSC_CpG_Island, ";", simplify = TRUE)
#  UCSC_CpG_Islands_Name <- stringr::str_split(infinium_methylationepic_v_1_0_b5_manifest_file$UCSC_CpG_Islands_Name, ";", simplify = TRUE)
#
#  PROBES_ID <- infinium_methylationepic_v_1_0_b5_manifest_file$IlmnID
#  CHR <- infinium_methylationepic_v_1_0_b5_manifest_file$CHR
#  START <- infinium_methylationepic_v_1_0_b5_manifest_file$MAPINFO
#
#  rm(infinium_methylationepic_v_1_0_b5_manifest_file)
#
#  keys <- data.frame( PROBE="",GENE="", ACCESSION="", GROUP="", CHR="", START="", POSITION ="")
#  PROBES_Gene_Body <- keys[keys$GROUP=="Body",]
#  PROBES_Gene_TSS200 <- keys[keys$GROUP=="TSS200",]
#  PROBES_Gene_TSS1500 <- keys[keys$GROUP=="TSS1500",]
#  PROBES_Gene_5UTR <- keys[keys$GROUP=="5'UTR",]
#  PROBES_Gene_3UTR <- keys[keys$GROUP=="3'UTR",]
#  PROBES_Gene_ExonBnd <- keys[keys$GROUP=="ExonBnd",]
#  PROBES_Gene_1stExon <- keys[keys$GROUP=="1stExon",]
#  dim_matrix <- dim(UCSC_RefGene_Accession)[2]
#  for (i in 1:dim_matrix ) {
#
#
#    temp_keys <- data.frame( PROBE=PROBES_ID,GENE=UCSC_RefGene_Name[,i], ACCESSION=UCSC_RefGene_Accession[,i], GROUP = UCSC_RefGene_Group[,i], CHR = CHR, START = START, POSITION = i)
#    temp_keys <- (temp_keys[!is.na(temp_keys$GENE),])
#
#    temp_keys_Gene_Body <- temp_keys[temp_keys$GROUP=="Body",]
#    temp_keys_Gene_TSS200 <- temp_keys[temp_keys$GROUP=="TSS200",]
#    temp_keys_Gene_TSS1500 <- temp_keys[temp_keys$GROUP=="TSS1500",]
#    temp_keys_Gene_5UTR <- temp_keys[temp_keys$GROUP=="5'UTR",]
#    temp_keys_Gene_3UTR <- temp_keys[temp_keys$GROUP=="3'UTR",]
#    temp_keys_Gene_ExonBnd <- temp_keys[temp_keys$GROUP=="ExonBnd",]
#    temp_keys_Gene_1stExon <- temp_keys[temp_keys$GROUP=="1stExon",]
#
#    PROBES_Gene_Body <- rbind(PROBES_Gene_Body,temp_keys_Gene_Body )
#    PROBES_Gene_TSS200 <- rbind(PROBES_Gene_TSS200,temp_keys_Gene_TSS200 )
#    PROBES_Gene_TSS1500 <- rbind(PROBES_Gene_TSS1500,temp_keys_Gene_TSS1500 )
#    PROBES_Gene_5UTR <- rbind(PROBES_Gene_5UTR,temp_keys_Gene_5UTR )
#    PROBES_Gene_3UTR <- rbind(PROBES_Gene_3UTR,temp_keys_Gene_3UTR )
#    PROBES_Gene_ExonBnd <- rbind(PROBES_Gene_ExonBnd,temp_keys_Gene_ExonBnd )
#    PROBES_Gene_1stExon <- rbind(PROBES_Gene_1stExon,temp_keys_Gene_1stExon )
#
#  }
#
#  RefGene_Group <- unique(as.vector(UCSC_RefGene_Group))
#
#  rm(UCSC_RefGene_Group)
#  rm(UCSC_RefGene_Accession)
#  rm(UCSC_RefGene_Name)
#
#  usethis::use_data(PROBES_Gene_Body, overwrite = TRUE)
#  usethis::use_data(PROBES_Gene_1stExon, overwrite = TRUE)
#  usethis::use_data(PROBES_Gene_3UTR, overwrite = TRUE)
#  usethis::use_data(PROBES_Gene_5UTR, overwrite = TRUE)
#  usethis::use_data(PROBES_Gene_ExonBnd, overwrite = TRUE)
#  usethis::use_data(PROBES_Gene_TSS1500, overwrite = TRUE)
#  usethis::use_data(PROBES_Gene_TSS200, overwrite = TRUE)
#
#
#  keys <- data.frame( PROBE="",ISLAND="", RELATION_TO_CPGISLAND="", CHR="", START="", POSITION ="")
#  PROBES_Island_Island <- keys[keys$RELATION_TO_CPGISLAND=="Island",]
#  PROBES_Island_N_Shore <- keys[keys$RELATION_TO_CPGISLAND=="N_Shore",]
#  PROBES_Island_S_Shore <- keys[keys$RELATION_TO_CPGISLAND=="S_Shore",]
#  PROBES_Island_S_Shelf <- keys[keys$RELATION_TO_CPGISLAND=="S_Shelf",]
#  PROBES_Island_N_Shelf <- keys[keys$RELATION_TO_CPGISLAND=="N_Shelf",]
#
#  temp_keys <- data.frame( PROBE=PROBES_ID, ISLAND= UCSC_CpG_Islands_Name , RELATION_TO_CPGISLAND = Relation_to_UCSC_CpG_Island, CHR = CHR, START = START, POSITION = 1)
#  temp_keys <- (temp_keys[!is.na(temp_keys$ISLAND),])
#
#  temp_keys_Island_Island <- temp_keys[temp_keys$RELATION_TO_CPGISLAND=="Island",]
#  PROBES_Island_Island <- rbind(PROBES_Island_Island,temp_keys_Island_Island )
#
#  temp_keys_Island_N_Shore <- temp_keys[temp_keys$RELATION_TO_CPGISLAND=="N_Shore",]
#  PROBES_Island_N_Shore <- rbind(PROBES_Island_N_Shore,temp_keys_Island_N_Shore )
#
#  temp_keys_Island_S_Shore <- temp_keys[temp_keys$RELATION_TO_CPGISLAND=="S_Shore",]
#  PROBES_Island_S_Shore <- rbind(PROBES_Island_S_Shore,temp_keys_Island_S_Shore )
#
#
#  temp_keys_Island_N_Shelf <- temp_keys[temp_keys$RELATION_TO_CPGISLAND=="N_Shelf",]
#  PROBES_Island_N_Shelf <- rbind(PROBES_Island_N_Shelf,temp_keys_Island_N_Shelf )
#
#  temp_keys_Island_S_Shelf <- temp_keys[temp_keys$RELATION_TO_CPGISLAND=="S_Shelf",]
#  PROBES_Island_S_Shelf <- rbind(PROBES_Island_S_Shelf,temp_keys_Island_S_Shelf )
#
#  rm(UCSC_CpG_Islands_Name)
#  rm(Relation_to_UCSC_CpG_Island)
#
#  usethis::use_data(PROBES_Island_Island, overwrite = TRUE)
#  usethis::use_data(PROBES_Island_N_Shore, overwrite = TRUE)
#  usethis::use_data(PROBES_Island_S_Shore, overwrite = TRUE)
#  usethis::use_data(PROBES_Island_N_Shelf, overwrite = TRUE)
#  usethis::use_data(PROBES_Island_S_Shelf, overwrite = TRUE)
#
#
#
# unique(as.factor(Regulatory_Feature_Group))
#
# keys <- data.frame( PROBE="",REGULATORY_FEATURE="", REGULATORY_FEATURE_GROUP="", CHR="", START="", POSITION ="")
# PROBES_Regulatory_Gene_Associated <- keys[keys$REGULATORY_FEATURE_GROUP=="Gene_Associated",]
# PROBES_Regulatory_Gene_Associated_Cell_type_specific<- keys[keys$REGULATORY_FEATURE_GROUP=="Gene_Associated_Cell_type_specific",]
# PROBES_Regulatory_NonGene_Associated<- keys[keys$REGULATORY_FEATURE_GROUP=="NonGene_Associated",]
# PROBES_Regulatory_NonGene_Associated_Cell_type_specific<- keys[keys$REGULATORY_FEATURE_GROUP=="NonGene_Associated_Cell_type_specific",]
# PROBES_Regulatory_Promoter_Associated<- keys[keys$REGULATORY_FEATURE_GROUP=="Promoter_Associated",]
# PROBES_Regulatory_Promoter_Associated_Cell_type_specific<- keys[keys$REGULATORY_FEATURE_GROUP=="Promoter_Associated_Cell_type_specific",]
# PROBES_Regulatory_Unclassified<- keys[keys$REGULATORY_FEATURE_GROUP=="Unclassified",]
# PROBES_Regulatory_Unclassified_Cell_type_specific<- keys[keys$REGULATORY_FEATURE_GROUP=="Unclassified_Cell_type_specific",]
#
# temp_keys <- data.frame( PROBE="",REGULATORY_FEATURE=Regulatory_Feature_Name, REGULATORY_FEATURE_GROUP=Regulatory_Feature_Group, CHR=CHR, START=START, POSITION =1)
# temp_keys <- (temp_keys[!is.na(temp_keys$REGULATORY_FEATURE),])
#
# temp_keys_Regulatory_Gene_Associated <- temp_keys[temp_keys$REGULATORY_FEATURE_GROUP=="Regulatory_Gene_Associated",]
# PROBES_Regulatory_Gene_Associated <- rbind(PROBES_Regulatory_Gene_Associated,temp_keys_Regulatory_Gene_Associated )
#
# da finire ....
#


# library(readr)
#
# EPIC_MANIFEST_DMR <- read_delim("~/Desktop/EPIC_MANIFEST-DMR.csv", ";", escape_double = FALSE, trim_ws = TRUE)
# PROBES_DMR <- EPIC_MANIFEST_DMR
# PROBES_DMR_1 <- data.frame(CHR = PROBES_DMR$CHR, START= PROBES_DMR$MAPINFO, END= PROBES_DMR$MAPINFO, PROBE = PROBES_DMR$IlmnID , DMR = PROBES_DMR$DMR, GROUP=1, POSITION=1)
#
# usethis::use_data(PROBES_DMR_1, overwrite = TRUE)

