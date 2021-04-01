**Meta-analysis of heterogeneous Down syndrome transcriptomic data**

Access the result of the meta-analysis through this web application:

https://ilariodetoma.shinyapps.io/shiny_meta-analysis/

This repository contains the Rmarkdown files with all the code to reproduce the analysis.

The *meta_analysis_DEs.Rmd* file contains the code of the differential expression analyses for each dataset used in the meta-analysis.
For the RNA-seq data, the actual code was loaded after having performed the analysis remotely. You can find the code of the analysis performed with DESeq2 in the separate _.R_ files and for the mapping in the _script_star_ file, both presents in the  folder _otherscripts_ .
Due to the size limitation of the online repository the input data to run this script were not committed. Therefore, in order to reproduce this script, researchers should download the raw-data from the GEO and the SRA databases and name the data in the same way as shown in the script.

The _meta_analysis.Rmd_ file contains the actual analysis for reproducing the figures in the manuscript and can be completely reproduced with the input files provided upon installation of the required packages.

All input and output files are also included in the manuscript, together with the figures of the paper.

The *Rmd* files were also run and knitted locally and the resulting *html* files committed to this online repository.

Files used as input by the *meta_analysis.Rmd* script:
* **201911_de_list_ds.RData** r data file with the results of all the differentially expression analyses performed in the *meta_analysis_DEs.Rmd* script
* gene pair combinations in the __combinazioni__ folder
* *coords2.rda* with the node coordinates for plotting
* **summary_metanalysis_DS.csv** csv file with the list of the comparisons of all the studies and their information
* **ds_metanalysis.csv** csv file with the list of DS genes found by the study of Vilardell et al.
* **hESC_tads_dixon** csv file with the TAD borders of hESC from Dixon et al.
* **chrsizehg18.txt** txt file with the chromosome sizes for human genome release hg18


Output files of the *meta_analysis.Rmd* script:

* **result_comorbidity_top500.csv** csv file with the output of the `disgenet2r` package (gene/disease associations)
* **coDEgenes_microcategories.csv** gene statistics with tissue micro-categories
* **coDEgenes_macrocategories.csv** gene statistics with tissue macro-categories

