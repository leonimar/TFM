
                          ### TEST COMPLETE PIPELINE ###
              ######################################################
              ######################################################
package_name <- "here"
  if (!requireNamespace(package_name, quietly = TRUE))
      install.packages(package_name)
  
library(here)

## CASE OF STUDY: FANCONI ANEMIA KEGG PATHWAY ID : hsa03460

path_id <- "03460"

 ## 1.INSTALL PACKAGES AND LIBRARIES

    source(here("scripts","00_librariesPackages.R"))

 ## 2.DISEASE MAP CREATION AND ANNOTATION

    source(here("scripts","01_DiseaseMapConstruction_Annotation.R"))

 ## 3.HPO ENRICHMENT AND DESCRIPTIVE ANALYSIS

    source(here("scripts","02_HPOenrichment_DescriptiveAnalysis.R"))

 ## 4.DATA PREPROCESSEMENT PRE-HIPATHIA

    source(here("scripts","03_DataPreparation_preHiPathia.R"))
    
 ## 5. HIPATHIA  VALIDATION
  
     # 5.1 HiPathia method of GSEA16334 dataset over FA KEGG pathway  
      
        source(here("scripts","04_HiPathiaMethod_4FAkeggPathway.R"))
    
     # 5.2 HiPathia method of GSEA16334 dataset over Disease MAP (Extended pathway)
              
        source(here("scripts","05_HiPathiaMethod_4DiseaseMap.R"))