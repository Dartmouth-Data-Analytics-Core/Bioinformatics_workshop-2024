Setting up an R project

We will be using R-Studio to explore and analyze genomics data on day 3, therefore we ask that you have R and R-Studio installed prior to attending the workshop. You will need to be running at least version 3.6 to ensure all of the packages needed will run smoothly. The latest versions for R and R-Studio can be found here and here.

Next you will need to set up a new project in R-Studio for this workshop. Projects in R are like containers for various jobs that you will perform, the history of the project will be loaded when you open a new project. By using a project you can install all of the tools ahead of time and they will be there for you when you need them during the workshop. In R-Studio under File select New directory and then select New Project and name the project something you will remember (bioinfo_workshop).

Now that you have your project loaded, run the following code to install all of the packages we will be using during the workshop. For those of you that haven't used RStudio before we have made a video showing the successful installation of the R packages you will need using the commands below.

I bumbled the description of the code chunks with the nested loops in the video, so here is a better description for those that are interested:

There are two loops and each loop starts with an if statement. The first loop states "if biomaRt is not installed enter this loop" and the second one "if BioCManager is not installed enter this loop", when the condition is fulfilled (the package is not installed) the loop is entered and the function install.packages is used to install the package. Each loop is exited once the packages are installed and the package is loaded with the 'library' function to make the functions contained in the package available during the current session.

```R
if (!any(rownames(installed.packages()) == "biomaRt")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
library(biomaRt)

if (!any(rownames(installed.packages()) == "IRanges")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("IRanges")
}
library(IRanges)

if (!any(rownames(installed.packages()) == "GenomicRanges")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("GenomicRanges")
}
library(GenomicRanges)

if (!any(rownames(installed.packages()) == "Gviz")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("Gviz")
}
library(Gviz)

if (!any(rownames(installed.packages()) == "org.Hs.eg.db")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

if (!any(rownames(installed.packages()) == "EnsDb.Hsapiens.v86")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("EnsDb.Hsapiens.v86")
}
library(EnsDb.Hsapiens.v86)

if (!any(rownames(installed.packages()) == "GenomicFeatures")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("GenomicFeatures")
}
library(GenomicFeatures)

if (!any(rownames(installed.packages()) == "VariantAnnotation")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("VariantAnnotation")
}
library(VariantAnnotation)

if (!any(rownames(installed.packages()) == "TxDb.Hsapiens.UCSC.hg38.knownGene")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
}
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

if (!any(rownames(installed.packages()) == "TxDb.Mmusculus.UCSC.mm10.knownGene")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
}
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

if (!any(rownames(installed.packages()) == "BSgenome")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("BSgenome")
}
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

if (!any(rownames(installed.packages()) == "ChIPseeker")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("ChIPseeker")
}
library(ChIPseeker)

BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

BiocManager::install("BSgenome.Mmusculus.UCSC.mm10.masked")
```

After all of the packages have been loaded run the following command 

```R
sessionInfo()
```

The output should look like this:


```r
> sessionInfo()
R version 4.1.2 (2021-11-01)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Monterey 12.6

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
 [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ChIPseeker_1.28.3                         TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
 [3] TxDb.Hsapiens.UCSC.hg38.knownGene_3.13.0  VariantAnnotation_1.38.0                 
 [5] Rsamtools_2.8.0                           Biostrings_2.60.2                        
 [7] XVector_0.32.0                            SummarizedExperiment_1.22.0              
 [9] MatrixGenerics_1.4.3                      matrixStats_0.63.0                       
[11] EnsDb.Hsapiens.v86_2.99.0                 ensembldb_2.16.4                         
[13] AnnotationFilter_1.16.0                   GenomicFeatures_1.44.2                   
[15] org.Hs.eg.db_3.13.0                       AnnotationDbi_1.54.1                     
[17] Biobase_2.52.0                            Gviz_1.36.2                              
[19] GenomicRanges_1.44.0                      GenomeInfoDb_1.28.4                      
[21] IRanges_2.26.0                            S4Vectors_0.30.2                         
[23] BiocGenerics_0.38.0                       biomaRt_2.48.3                           

loaded via a namespace (and not attached):
  [1] shadowtext_0.1.2                        backports_1.4.1                        
  [3] fastmatch_1.1-3                         Hmisc_4.7-0                            
  [5] BiocFileCache_2.0.0                     plyr_1.8.8                             
  [7] igraph_1.3.5                            lazyeval_0.2.2                         
  [9] splines_4.1.2                           BiocParallel_1.26.2                    
 [11] ggplot2_3.4.0                           digest_0.6.30                          
 [13] yulab.utils_0.0.5                       htmltools_0.5.3                        
 [15] GOSemSim_2.18.1                         viridis_0.6.2                          
 [17] GO.db_3.13.0                            fansi_1.0.3                            
 [19] magrittr_2.0.3                          checkmate_2.1.0                        
 [21] memoise_2.0.1                           BSgenome_1.60.0                        
 [23] cluster_2.1.4                           graphlayouts_0.8.3                     
 [25] enrichplot_1.12.3                       prettyunits_1.1.1                      
 [27] jpeg_0.1-9                              colorspace_2.0-3                       
 [29] blob_1.2.3                              rappdirs_0.3.3                         
 [31] ggrepel_0.9.2                           xfun_0.35                              
 [33] dplyr_1.0.10                            jsonlite_1.8.3                         
 [35] crayon_1.5.2                            RCurl_1.98-1.9                         
 [37] scatterpie_0.1.8                        TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
 [39] ape_5.6-2                               survival_3.4-0                         
 [41] glue_1.6.2                              polyclip_1.10-4                        
 [43] gtable_0.3.1                            zlibbioc_1.38.0                        
 [45] DelayedArray_0.18.0                     scales_1.2.1                           
 [47] DOSE_3.18.3                             DBI_1.1.3                              
 [49] Rcpp_1.0.9                              plotrix_3.8-2                          
 [51] viridisLite_0.4.1                       progress_1.2.2                         
 [53] htmlTable_2.4.1                         tidytree_0.4.1                         
 [55] gridGraphics_0.5-1                      foreign_0.8-83                         
 [57] bit_4.0.5                               Formula_1.2-4                          
 [59] htmlwidgets_1.5.4                       httr_1.4.4                             
 [61] fgsea_1.18.0                            gplots_3.1.3                           
 [63] RColorBrewer_1.1-3                      ellipsis_0.3.2                         
 [65] pkgconfig_2.0.3                         XML_3.99-0.12                          
 [67] farver_2.1.1                            nnet_7.3-18                            
 [69] dbplyr_2.2.1                            deldir_1.0-6                           
 [71] utf8_1.2.2                              ggplotify_0.1.0                        
 [73] reshape2_1.4.4                          tidyselect_1.2.0                       
 [75] rlang_1.0.6                             munsell_0.5.0                          
 [77] tools_4.1.2                             cachem_1.0.6                           
 [79] cli_3.4.1                               generics_0.1.3                         
 [81] RSQLite_2.2.18                          stringr_1.4.1                          
 [83] fastmap_1.1.0                           yaml_2.3.6                             
 [85] ggtree_3.0.4                            knitr_1.41                             
 [87] bit64_4.0.5                             tidygraph_1.2.2                        
 [89] caTools_1.18.2                          purrr_0.3.5                            
 [91] KEGGREST_1.32.0                         ggraph_2.1.0                           
 [93] nlme_3.1-160                            aplot_0.1.8                            
 [95] DO.db_2.9                               xml2_1.3.3                             
 [97] compiler_4.1.2                          rstudioapi_0.14                        
 [99] filelock_1.0.2                          curl_4.3.3                             
[101] png_0.1-7                               treeio_1.16.2                          
[103] tibble_3.1.8                            tweenr_2.0.2                           
[105] stringi_1.7.8                           lattice_0.20-45                        
[107] ProtGenerics_1.24.0                     Matrix_1.5-1                           
[109] vctrs_0.5.1                             pillar_1.8.1                           
[111] lifecycle_1.0.3                         BiocManager_1.30.19                    
[113] cowplot_1.1.1                           data.table_1.14.6                      
[115] bitops_1.0-7                            patchwork_1.1.2                        
[117] qvalue_2.24.0                           rtracklayer_1.52.1                     
[119] R6_2.5.1                                BiocIO_1.2.0                           
[121] latticeExtra_0.6-30                     KernSmooth_2.23-20                     
[123] gridExtra_2.3                           dichromat_2.0-0.1                      
[125] gtools_3.9.3                            boot_1.3-28                            
[127] MASS_7.3-58.1                           assertthat_0.2.1                       
[129] rjson_0.2.21                            GenomicAlignments_1.28.0               
[131] GenomeInfoDbData_1.2.6                  hms_1.1.2                              
[133] ggfun_0.0.8                             rpart_4.1.19                           
[135] tidyr_1.2.1                             biovizBase_1.40.0                      
[137] ggforce_0.4.1                           base64enc_0.1-3                        
[139] interp_1.1-3                            restfulr_0.0.15  
```

