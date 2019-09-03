


reqd.packages <- c("dplyr","tidyr","mirbase.db")

source("https://raw.githubusercontent.com/srkoppolu/PackageLoader/master/PackageLoad2.R")
PackageLoad2(reqd.packages)


# ## Read the mirBase annotation dataset
# mirBase.data <- fread(fname.mirbase, sep = "\t", header = T, stringsAsFactors = F)

## Isolate the miRNAs for homo sapiens
speciess <- mirbaseID2SPECIES
speciess.mk <- mappedkeys(speciess)
mir.species <- as.data.frame(speciess[speciess.mk])
mir.select <- mir.species[mir.species$organism %in% "hsa", ]

## Get the mirbase IDs for the human miRNAs
acc <- mirbaseMATURE
mir.acc <- as.data.frame(acc[mir.select$mirna_id])

## Get the miRNA sequences for human miRNAs
sequ <- mirbaseSEQUENCE
mir.sequ <- as.data.frame(sequ[mir.select$mirna_id])

## Merge the mirbase IDs and sequences for human miRNAs.
mir.info.raw <- merge(mir.acc[,c("mirna_id","mature_acc","mature_name","mature_from","mature_to","similarity")], 
                      mir.sequ, 
                      by = "mirna_id",
                      allow.cartesian=TRUE)

## Merge the mirbase IDs and sequences for human miRNAs.
mir.info.raw <- merge(mir.acc, mir.sequ, by = "mirna_id")

mir.info.2 <- mir.info.raw[ ,c("sequence","mature_from","mature_to")]
mir.info.2$mature_seq <- apply(as.matrix(mir.info.2), 1, 
                               FUN = function(x){return(substr(as.character(x[1]),x[2],x[3]))})

mir.info <- mir.info.raw[ ,c("mature_acc","mirna_id","mature_name")]
mir.info$mature_seq <- mir.info.2$mature_seq

## Create an alias for all mirna annotations and remove the duplicates
mirbase.data <- melt(mir.info, id.vars = c("mature_acc","mature_seq"))
mirbase.data$variable <- NULL
mirbase.data <- unique(mirbase.data)
mirbase.data <- mirbase.data %>%
  group_by(mature_acc) %>%
  mutate(mirna_Aliases = paste0(value, collapse = "|")) 
mirbase.data$value <- NULL
mirbase.data <- unique(mirbase.data)

## Remove the other files
rm(mir.info, mir.info.2, mir.info.raw)

## Clear Screen
cat("\014")
