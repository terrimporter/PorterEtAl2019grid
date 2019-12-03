# README

This repository contains the dataflow and scripts I used to process the COI metabarcode reads in the paper Porter et al., 2019 available from BioRxiv https://doi.org/10.1101/693499 .

If you use these data files or scripts in future work, please cite:
Porter, T.M., Morris, D.M., Basiliko, N., Hajibabaei, M., Doucet, D., Bowman, S., Emilson, E.J.S., Emilson, C.E., Chartrand, D., Wainio-Keizer, K., Séguin, A., Venier, L.  2019.  Variations in terrestrial arthropod DNA metabarcoding methods recovers robust beta diversity but variable richness and site indicators.  Scientific Reports, 9: 18218. https://www.nature.com/articles/s41598-019-54532-0

## Infiles

1. The denoised ESV fasta files for the BE and F230 amplicons are in denoised_ESVs.fastas.zip

2. The original taxonomic assignment matrix from the SCVUC bioinformatic pipeline is LV2016_1.tar.gz

3. Sample x filename map GRDIname.map

4. Taxonomic assignment matrix for ESVs with edited headers LV2016_2.tar.gz.  Taxonomic assignment matrix for OTUs with edited headers LV2016_4.tar.gz.

5. Metadata including latitude and longitude in metadata.csv

## Data analysis outline

1. Raw reads were processed using the SCVUC v2.0 pipeline available from https://github.com/EcoBiomics-Zoobiome/SCVUC_COI_metabarcode_pipeline . 

2. The resulting taxonomic assignments for each marker were concatenated (LV2016_1.tar.gz) and headers edited for analysis in R (LV2016_2.tar.gz) using GRDIname.map and the Perl script reviseGRDInames2.plx .  

3. Figures were pepared in R as follows:
  * Fig 2 was generated with Fig2_Richness_FigS4_Venn_ESV_100919.R . 
  * Fig 3 was generated with Fig3_BetaDiv_100819.R from LV2016_2.csv and metadata.csv .
  * Fig 4 was generated with Fig4_FigS7_FigS8_FigS9_indicator_heattree_ESV_100819.R .  
  
4. Supplementary figures were prepared in R as follows:
  * Fig S1 was generated with FigS1_phyla.R . 
  * Fig S2 was generated with FigS2_PropConfID.R . 
  * Fig S3 was generated with FigS3_Rarefaction.R . 
  * Fig S4 was generated with Fig2_Richness_FigS4_Venn_ESV_100919.R . 
  * Fig S5 was generated with FigS5_RichnessDNAextn.R .
  * Fig S6 was generated with FigS6_heatmap_indic_ESVs_100919.R
  * Fig S7, S8, S9 were generated with Fig4_FigS7_FigS8_FigS9_indicator_heattree_ESV_100819.R .

## Acknowledgements

I would like to acknowledge funding from the Canadian government from the Genomics Research and Development Initiative (GRDI) Ecobiomics project.

Last updated: December 2, 2019
