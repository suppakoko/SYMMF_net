# SYMMF-net - SYnergistic Meta-Microbe Feature Network
SYMMF-net is a network analysis tool for Microbiome 


## Files for analysis

1. ERP004264_PD.table -- example file (Parkinson's disease)
2. NMF_analysis.R -- R script for NMF analysis
3. SYMMF_net_analysis.py -- Python script for SYMMF-net analysis

## Running Script example

1. Running NMF_analysis.R script

  ```
  $ Rscript NMF_analysis.R -m ERP004264_PD.table
  ```
  
2. Copy SYMMF_net_analysis.py to same folder as result of NMF_analsis.R

3. Running SYMMF_net_analysis.py script
  
  ```
  $ python SYMMF_net_analysis.py -p 0.05
  ```

## Development environment of R
  ```
  R ver. = 3.4.0
  required R packages
  Bioconductor biobase
  pheatmap
  RColorBrewer
  gplots
  NMF
  dendextend
  dplyr
  tidyr
  reshape2
  ggplot2
  argparse
  ```
### Output files after executing NMF_analysis.R
  
  ```
  Clustering_purity_value.txt -- matrix of clustering purity results (m = 2 ~ 25)
  raw_coefmat_[K#].txt -- coef matrix of each [K#]
  raw_bmap_mat_[K#].txt -- basis matrix of each [K#]
  group_coefmat_[K#].txt -- coef matrix of each [K#] group by each sample group
  extracted_microb_feat_by_[K#].txt -- microbe list of each [K#]
  disease_k_connect_by_[K#].tiff -- heatmap image of each group_coefmat_[K#].txt
  basismatrix.pdf -- heatmap image of all basismap (m = 2 ~ 25)
  coefmatrix.pdf -- heatmap image of all coefmap (m = 2 ~ 25)
  micro_count_mat.txt -- input data matrix afeter normalization and filtering
  micro_count_mat_pre.txt -- input data matrix before normalization and filtering
  microbe_list_original.txt -- microbe list of input data from original matrix
  microbe_list_trimmed.txt -- microbe list of input data from trimmed matrix
  sample_lable_list_trimmed.txt -- sample list of input data from trimmed matrix
  MMF_microbe_list.txt -- microbe list and each contribution value from basismatrix
  ```
  
## Development environment of Python
  ```
  Python = 3.7.5
  pandas = 0.25.3
  numpy =  1.17.4
  seaborn = 0.9.0
  sklearn = 0.22
  scipy = 1.3.2
  matplotlib = 3.1.1
  argparse = 1.1
  ```
  
### Output files after executing SYMMF_net_analysis.py
  
  ```
  mi_[Sample].data -- all of AUC value of each MMFs
  kk_[Sample].data -- all of AUC value of each microbe
  options_record.txt -- analysis option record
  [Sample]_pairs_list.sif -- All microbe list from SYMMFs
  AUC_value_list.txt -- All of AUC values all microbe and MMFs
  MMF_high_AUC_and_list.data -- Sample specific MMF list and each AUC and highest micro AUC each MMF
  microbe_AUC_calculation.xlsx -- AUC calculation result of each microbe
  MMF_AUC_calculation.xlsx -- AUC calculation result of each MMF
  SYMMF_summary_table.xlsx -- SYMMF information table
  MMF_microbe_AUC_data_summary.xlsx -- all information of each MMF
  NMF_result_summary.txt -- MMF information table
  MMF_microbe_list_ordered.txt -- All microbe list of each MMF order by contribution value
  SYMMF_feature_selected_microbe_list.txt -- microbe list from selected SYMMF of SYMMF-net
  ```
  
### matrix of SYMMF
 
 ```
 SYMMF_feature_bmap_selected.txt -- W_symmf matrix
 SYMMF_feature_coef_selected.txt -- H_symmf matrix
 SYMMF_feature_selected_microbe_count_mat.txt
 cor_SYMMF_feature_bmap_selected.txt -- calculation result of correlation of each microbe from SYMMF_feature_bmap_selected.txt
 ```
 
### Output for network generation
 
 ```
 whole_network_cor[value]_p[value]_auc[value]_ECDF[value].sif -- cytoscape network file
 whole_network_nodes_size_p[value]_auc[value]_ECDF[value].attrs -- network node info
 whole_network_edges_cor[value]_p[value]_auc[value]_ECDF[value].attrs -- network edges info
 whole_network_node_piechart_count_ratio_cor[value]_p[value]_auc[value]_ECDF[value].table -- network node info of pie-chart count
 whole_network_node_piechart_contribution_ratio_cor[value]_p[value]_auc[value]_ECDF[value].table -- network node info of pie-chart contribution
 ```

### R script for sub-analysis

  ```
  SYMMF_feature_cor_gen.R -- Rscript for Pearson's correlation calculation
  AUC_comparison_rscript.R -- Rscript for AUC comparison boxplot
  ```
  
### Output Figures

  ```
  NMF_microbe_AUC_comparison.tiff -- boxplot for AUC comparison between NMF and microbe
  purification_comparison_plot.png -- Clustering purity comparison plot
  all_SYMMF_MMF_scatter_plot.jpg -- scatter plot of all MMFs
  SYMMF_MMF_scatter_plot_ECDF[value].jpg -- scatter plot of selected MMFs ECDF cut-off
  MMF_AUC_save_kde_ECDF[value].png -- scatter plot of selected MMFs ECDF cut-off kde
  MMF_AUC_save_hex_ECDF[value].png -- scatter plot of selected MMFs ECDF cut-off hex
  SYMMFs_coefmap_mean_reord_save_ECDF[value].png -- heatmap of coefmap selected SYMMFs group by mean value
  SYMMFs_coefmap_save_ECDF[value].png -- heatmap of coefmap selected SYMMFs
  SYMMFs_basismap_reord_save_ECDF[value].png -- heatmap of basismap selected SYMMFs
  ```
  
## Contributor

  ```
  Keunwan Park (keunwan@kist.re.kr)
  Young-Joon Ko (yjko@kist.re.kr)
  ```
# Docker image
  ```
  docker pull suppak/r340_symmf
  ```
