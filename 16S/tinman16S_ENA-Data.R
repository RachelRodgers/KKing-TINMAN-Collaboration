# tinman16S_ENA-DATA.R

library("here")
library("tidyverse")

sampleDataFinal <- readRDS(here("../data/RDataObjects/16S/sampleDataFinal.RDS"))

# Data for Sample Registration
sampleReg <- sampleDataFinal %>% 
  dplyr::rename(sample_alias = sample_name) %>% 
  mutate(sample_title = paste0("Subject_", unique_id, "_", timepoint_condition),
         sample_description = sample_title,
         collection_date = "NA",
         geographic_location = "USA",
         tax_id = "9606",
         scientific_name = "Homo sapiens") %>% 
  select(tax_id,
         scientific_name,
         sample_alias,
         sample_title,
         sample_description,
         collection_date,
         geographic_location)

# save
write.table(sampleReg, file = here("../documents/16S_ENA_sampleRegInfo.txt"),
            append = FALSE, quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")

# Data for Reads Submission
seqReg <- read.delim("../data/tinman16S_sequencingFileNames.txt",
                       header = FALSE, col.names = "sequencingFileName") %>% 
  mutate("type" = ifelse(grepl("\\.1\\.", sequencingFileName),
                         yes = "forward_file_name",
                         no = "reverse_file_name"),
         "sample" = str_remove_all(string = sequencingFileName,
                                   pattern = c("^King_|.1.fq.bz2$|.2.fq.bz2$"))) %>% 
  pivot_wider(names_from = type, values_from = sequencingFileName) %>% 
  mutate("study" = "PRJEB72348",
         "instrument_model" = "Illumina MiSeq",
         "library_name" = "",
         "library_source" = "METAGENOMIC",
         "library_selection" = "PCR",
         "library_strategy" = "AMPLICON",
         "library_layout" = "PAIRED",
         "forward_file_md5" = "",
         "reverse_file_md5" = "") %>% 
  select(sample, study, instrument_model, library_name, library_source,
         library_selection, library_strategy, library_layout, forward_file_name,
         forward_file_md5, reverse_file_name, reverse_file_md5)

write.table(seqReg, file = here("../documents/16S_ENA_sequencingFileRegInfo.txt"),
            append = FALSE, quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")
