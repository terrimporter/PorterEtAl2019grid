# README

This repository contains the dataflow and scripts I used to process the COI metabarcode reads in the paper Porter et al., 2019 (in prep).

## Data analysis outline

1. Raw reads were processed using the SCVUC v2.0 pipeline available from https://github.com/EcoBiomics-Zoobiome/SCVUC_COI_metabarcode_pipeline 
2. The resulting taxonomic assignments for each marker were concatenated (LV2016_1.tar.gz) and headers edited for analysis in R (LV2016_2.tar.gz) using GRDIname.map and the Perl script reviseGRDInames2.plx .
3. Figures were pepared in R as follows:
  * Fig 1 was generated with Fig1_PropConfID.R from LV2016_2.csv
  * Fig 2 was generated with Fig2_Richness.R from LV2016_2.csv
  * Fig 3 was generated with Fig3_BetaDiv.R and legends were edited in Inkscape.  This script uses metadata.csv and LV2016_2.csv .
  * Fig 4 was generated with Fig4_Fig5_Indicator.R from LV2016_2.csv
  * Fig 5 was generated with Fig4_Fig5_Indicator.R from LV2016_2.csv
4. Supplementary figures were prepared in R as follows:
  * S1 Fig was generated with FigS1_TaxSummary.R from LV2016_2.csv
  * S2 Fig was generated with FigS2_OrderSummary.R from LV2016_2.csv
  * S3 Fig was generated with FigS3_Rarefaction.R from LV2016_2.csv
  * S4 Fig was generated with Fig4_Fig5_Indicator.R from LV2016_2.csv
  * S5 Fig was generated with Fig4_Fig5_Indicator.R from LV2016_2.csv
  * S6 Fig was generated with Fig4_Fig5_Indicator.R from LV2016.csv

## Acknowledgements

I would like to acknowledge funding from the Canadian government from the Genomics Research and Development Initiative (GRDI) Ecobiomics project.

Last updated: February 12, 2019
