# README

This repository contains the dataflow and scripts I used to process the COI metabarcode reads in the paper Porter et al., 2019 available from BioRxiv https://doi.org/10.1101/693499 .

## Supporting Files

1. The denoised ESV fasta files for the BE and F230 amplicons are in denoised_ESVs.fastas.zip

2. The original taxonomic assignment matrix from the SCVUC bioinformatic pipeline is LV2016_1.tar.gz

## Data analysis outline

1. Raw reads were processed using the SCVUC v2.0 pipeline available from https://github.com/EcoBiomics-Zoobiome/SCVUC_COI_metabarcode_pipeline . 

2. The resulting taxonomic assignments for each marker were concatenated (LV2016_1.tar.gz) and headers edited for analysis in R (LV2016_2.tar.gz) using GRDIname.map and the Perl script reviseGRDInames2.plx .  

3. Figures were pepared in R as follows:
  * Fig 1 was generated with Fig1_Richness_FigS5_Venn_ESV.R . 
  * Fig 2 was generated with Fig2_NMDS.R from LV2016_2.csv and legends were edited in Inkscape.
  * Fig 3 was generated with Fig3_FigS8_FigS9_FigS10_indicator_heattree_ESV.R .  
  
4. Supplementary figures were prepared in R as follows:
  * Fig S2 was generated with FigS2_phyla.R . 
  * Fig S3 was generated with FigS3_PropConfID.R . 
  * Fig S4 was generated with FigS4_Rarefaction.R . 
  * Fig S5 was generated with Fig1_Richness_FigS5_Venn_ESV.R . 
  * Fig S6 was generated with FigS6_RichnessDNAextn.R .
  * Fig S7 was generated with FigS7_heatmap_indic_ESVs.R
  * Fig S8, S9, S10 were generated with Fig3_FigS8_FigS9_FigS10_indicator_heattree_ESV.R .

## Acknowledgements

I would like to acknowledge funding from the Canadian government from the Genomics Research and Development Initiative (GRDI) Ecobiomics project.

Last updated: July 8, 2019
