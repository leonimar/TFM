# INSTRUCTIONS FOR THE CORRECT PERFORMANCE OF THE PIEPLINE

The pipeline must be excuted under R base 3.5.1 or higher

In order to correctly charge all the scripts and paths it must we performed in the PIPELINE_TFM_Marina.Rproj
All the scripts needed to repropduce the results presented in the TFM
"Mechanistic models for drug repositioining in rare diseases" are located
in "/scripts", they MUST be executed in their numerated order: 00-01-02-03-04-05

00_librariesPackages.R --> Loads all the packages and libraries needed in this project.

01_DiseaseMapConstruction_Annotation.R --> Builds up the Disease Map from a given KEGGpathwayID(XXXX), i.e. For hsa03460 => KEGGpathID = "03460" 
                                           and two tabular files containing extra genes and extra interactions named after the KEGGpathID, i.e. For hsa03460 => vertex_03460.xls; edges_03460.xls
                                           Annotation of the Disease Map with DrugBank and HPO

02_HPOenrichment_DescriptiveAnalysis.R --> Enrichment analysis of the HPO terms associated to the genes of the Disease Map ; Enrichment of the HPO associated to the Disease-network

03_DataPreparation_preHiPathia.R --> Preprocessment of the genesets chosen:  1.Normalization with RMA of the raw data 2. Tranformation Probe to gene (preprocessment for HiPathia method)

04_HiPathiaMethod_4FAkeggPathway.R --> HiPathia validation Part I: Performance of HiPathia method on previously preproccessed geneset GSEA16334 over Fanconi Anemia KEGG Pathway

05_HiPathiaMethod_4DiseaseMap.R --> HiPathia validation Part II: Performance of HiPathia method on previously preproccessed geneset GSEA16334 over Disease Map (Extended FA Disease pathway)


"/data"

Contains all the tables and data sets needed to carry out the analyses

https://drive.google.com/open?id=12y5PyQXmlwzaXKgnZqGZcpIfX70lMW-O


"/library"

Contains information needed to buld up the new metaginfo for the inclusion of the Disease Map in HiPathia method.


"/OUTPUT"

All the tables and results will be exported to this folder.





