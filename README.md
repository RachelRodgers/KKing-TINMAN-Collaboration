# KKing-TINMAN-Collaboration

Repository of scripts associated with the study "Antibiotic-associated neutropenia is marked by decreased predominance of intestinal <em>Lachnospiraceae</em> in pediatric patients."

The raw data files for this study are available through the [European Nucleotide Archive Accession Number PRJEB72348] (https://www.ebi.ac.uk/ena/browser/text-search?query=PRJEB72348).

Associated data files needed to run the scripts are available in the documents subdirectory.

## Bacterial 16S Analyses
All code associated with the analyses of the 16S data can be found in the 16S subdirectory. Processing and analyses is segmented into several scripts and should be run in the following order:

**1.) timan16S_Format-Metadata.Rmd** - Scripts for formatting the original metadata into an R-friendly format.
**2.) tinman16S_DADA2.R** - Generation and annotation of an amplicon sequence variant (ASV) table using the quality controlled forward reads only. Results are stored in a phyloseq object for downstream analyses.
**3.) tinman_16S_Analysis-SILVA.Rmd** - Primary ecological analyses.
**4.) tinman_16S_Analysis-of-DA-Taxa.Rmd** - Additional analyses of potential biomarkers.
**5.) tinman_16S_ForestPlots.Rmd** - Code used to generate Forest Plots.
**6.) tinman16S_ENA-Data.R** - Code used to organize and format data for submission to the ENA.
