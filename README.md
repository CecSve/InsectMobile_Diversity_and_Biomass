# InsectMobile: landscape-level drivers of flying insect diversity and biomass 
WORK IN PROGRESS
## 
This repository (R project) contains all scripts necessary to run the ecological/statistical analyses from the study "Landscape-level drivers of flying insect diversity and biomass" (Svenningsen et al. XXXX, DOI: [add link]).

The data were collected in June/July 2018 and June 2019 as part of the citizen science project InsectMobile ("Insektmobilen") at the Natural History Museum of Denmark and iDiv, Germany.

* **Data**: data (e.g. the proportional land cover and land use data, flying insect biomass, diversity estimates etc. for each buffer zone) used in the analyses is deposited in XX: [add link]

---
**NOTE ON BIOINFORMATICS**

Two sequencing platforms were used to generate the data: HiSeq 4000 and NovaSeq 6000. The NovaSeq processing assigns quality scores differently from the HiSeq platform, where NovaSeq [simplify the error rates](https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/novaseq-hiseq-q30-app-note-770-2017-010.pdf) by binning the 40 possible quality scores into just 4 categories which [vastly reduces the amount of information dada2 can work off of to infer errors in the data](https://github.com/benjjneb/dada2/issues/791). This discrepancy between platforms were not dealt with during the dada2 processing as we were unaware of the problem when we ran the bioinformatics. It seems to affect how many rare species/sequences are detected/retained in that fewer rare species are detected from NovaSeq data if the error rate step is not updated in the dada2 pipeline to accomodate the fewer quality scores. If anyone wish to use the data for analysis, we encourage users to find a way to deal with the NovaSeq processing in dada2 (for example by testing the four options mentioned [here](https://github.com/ErnakovichLab/dada2_ernakovichlab#learn-the-error-rates)) and apply the best suited fix for this data by evaluating the quality score plots.

---

## Description of the sub directories for data processing ##

* **reports**: a step-wise list of scripts (01_, 02_ etc.) used for processing and analysing the data

## Statistical analyses and modelling ##
Landscape-level effects on flying insect diversity and biomass was modelled with linear mixed-effects models at the buffer size with the most pronounced effect size.
