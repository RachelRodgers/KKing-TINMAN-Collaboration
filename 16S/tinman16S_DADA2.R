# tinman16S_DADA2.R

library("dada2")
library("msa")
library("phangorn")
library("phyloseq")
library("ggpubr")
library("tidyverse")

#----- Set Paths -----#

inputDir <- "../data/raw_data/"
filteredDir <- "../data/filtered"

taxonomyDBPath <- "/mnt/pathogen1/rrodgers/databases/dada2_taxonomy"

# RDP v18
rdpDB <- file.path(taxonomyDBPath, "rdp_train_set_18.fa.gz")
rdpDBSpecies <- file.path(taxonomyDBPath, "rdp_species_assignment_18.fa.gz")

#----- Analyze Read Quality Data -----#

# Only the R1 files are needed:
rawR1Files <- list.files(path = inputDir, pattern = ".1.fq") # 76 (removed empty blank)

# Raw Read Quality
rawQualPlot <- plotQualityProfile(paste0(inputDir, rawR1Files), 
                                  aggregate = TRUE)
rawQualPlot # Good quality maintained throughout. Will trim last 10bp to 240bp

# Filter & Trim
filterAndTrim(fwd = paste0(inputDir, rawR1Files), filt = filteredDir, 
              trimLeft = 10, maxEE = 2, 
              truncQ = 11, maxN = 0, rm.phix = TRUE, compress = TRUE, 
              truncLen = 240, verbose = TRUE)

#----- Dereplication -----#

# Get list of filtered files and assign sample names to them
filteredFiles <- list.files(filteredDir, pattern = "fq", full.names = TRUE)
r1SampleNames <- map_chr(basename(filteredFiles),
                         ~ str_remove_all(string = .x, 
                                          pattern = c("^King_|.1.fq$")))

names(filteredFiles) <- r1SampleNames

# Infer Error Rates
set.seed(98669693)
errF <- learnErrors(filteredFiles, multithread = TRUE)
errorPlot <- plotErrors(errF, nominalQ = TRUE)
errorPlot

# Create a list that will hold dereplication objects for each sample
singles <- vector("list", length(r1SampleNames))
names(singles) <- r1SampleNames

for(sample in r1SampleNames) {
  derepF <- derepFastq(filteredFiles[[sample]])
  singles[[sample]] <- dada(derepF, err = errF, multithread = TRUE)
}

#----- Construct Sequence Table -----#

sequenceTable <- makeSequenceTable(singles)
sequenceTableNoChimeras <- removeBimeraDenovo(sequenceTable, multithread = TRUE)

#----- Assign Taxonomy -----#

taxaRankNamesFull <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxaRankNamesTrunc <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# RDP
taxaRDP <- assignTaxonomy(sequenceTableNoChimeras, rdpDB, multithread = TRUE)
colnames(taxaRDP) <- taxaRankNamesTrunc
taxaRDPPlus <- addSpecies(taxaRDP, rdpDBSpecies)

#----- Construct Phylogenetic Tree -----#

# Get the sequences from the sequence table
seqs <- getSequences(sequenceTableNoChimeras)
names(seqs) <- seqs
# Multiple sequence alignment
mult <- msa(seqs, method = "ClustalW", type = "dna", order = "input")
# Convert MSA to phyDAT format
phangAlign <- as.phyDat(mult, type = "dna", order = "input")
# Compute pairwise distances on phangAlign
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm)
# Compute likelihood of tree
fit <- pml(tree = treeNJ, data = phangAlign)
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", 
                    control = pml.control(trace = 0))

#----- Build Phyloseq Objects -----#

physeqDir <- "../data/RDataObjects/physeqObjects"
dir.create(physeqDir)

# RDP
ps0.rdp <- phyloseq(otu_table(sequenceTableNoChimeras, taxa_are_rows = FALSE),
                    tax_table(taxaRDPPlus), phy_tree(fitGTR$tree))
saveRDS(ps0.rdp, file.path(physeqDir, "ps0.rdp18_single.RDS"))

#----- Save Data -----#

writeLines(capture.output(sessionInfo()), 
           "tinman16S_DADA2_session_info.txt")
Sys.Date()
getwd()
sessionInfo()
