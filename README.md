# Tmuris_BulkRNA
Code employed for BulkRNA analysis

- KallistoRun details the commands used for running Kallisto on the trimmed FastQ files

- PrincipalComponentAnalysis details the creation of the PCA using the Kallisto data, as well as plotting of the scree plots and PCA plots. 

- Sleuth_Running details the preparation and eventual running of Sleuth analysis using the Kallisto data. 

- Sleuth_LTR_Upset_Heatmap details the code on creating the Upset plot ans the heatmaps from normalized counts organized by Sleuth

- GProfiler GMT Genertion details the use of Omer Faruk Bay's Gene Ontology annotation to generate a custom GMT file to use for online GProfiler analysis. 

- Sleuth_Wald_Volcano_Gprofiler_Loop details the code on creating the Volcano plots of various comparisons, using the Wald test results from Sleuth. It also shows the plotting of the gene enrichment analysis from GProfiler's output. 

- Clust_Preparation details the, well, preparation of the data and the command employed to run BaselAbujamous' clust for co-expression analysis. 
