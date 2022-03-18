# Merge datatables prior to analysis - this script should generate the final table(s) for analysis and the output to be shared in a datasharing repository (journal specific)

### merge with PCRID ############
test <- rbind(sampleid_size_pcrid_concentration, pcrid_qubit_19)
data <- merge(data, test, by = "SampleID_size")

#### match asvtable data with lab meta data (not all samples were sequenced) ####
keep <- colnames(asvs)
dplyr::setdiff(data$PCRID, colnames(asvs))
tasvs <- as.data.frame(t(asvs))
tasvs <- tasvs %>% rownames_to_column(var = "PCRID")
missingsamples <- anti_join(data, tasvs) # Rows in the data that do not have a match in the asv table. 266 samples are not in the ASV table and at least some of them were sequenced but apparently they were removed during the bioinformatics. It seems library 38 and 44 is largely missing - check this.
labdata <- data %>% filter(PCRID %in% keep)
nomatchsamples <- anti_join(tasvs, labdata) # rows in the asv table that don't have a match in the labdata - many bird droppings samples and some samples from lib 28
keep <- labdata$PCRID

asvtable <- asvs[ ,colnames(asvs) %in% keep]
rowSums(asvtable)

otus <-
  asvtable[apply(asvtable[,-1], 1, function(x)
    ! all(x == 0)),] # remove rows that contain only zeros (OTUs that are not present in the subsetted samples)

keep <- rownames(otus)

taxonomy <- taxonomy %>% filter(occurrenceId %in% keep)
min(colSums(otus))