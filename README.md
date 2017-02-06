# Conserved noncoding sequence transcription factor binding prediction


Our algorithm, as illustrated in Fig. 1,  rst builds an initial gene regulatory network based on the integration of multiple heterogeneous, transcriptional-regulation related datasets (Fig. 1(A)). Data integration has proven necessary in the context of network inference for higher organisms11,14. A general challenge of data integration is the limited availability of di erent datasets. For example, transcription factor binding information can be pow- erful for establishing directed regulatory networks. However, binding information is available only for a limited number of transcription factors for most organisms10.  erefore, in order to produce an initial network, GRACE uses genome-wide datasets  rst, followed by a network re nement based on additional, more sparsely available, datasets. First, GRACE constructs an expression based gene regulatory network. To this end, GRACE implements a random forest regression model similar to the one used by GENIE320, which is a state of the art gene expres- sion based network inference algorithm. However, GRACE’s highly scalable random forest regression model performs several times faster than the one used in GENIE3 (see supplement for details and speed comparison). Subsequently, an empirical cumulative distribution over all link predictions is constructed and only the top 5% of all expression based link predictions are kept. Finally, these top predictions are further  ltered with available transcription factor binding within conserved non-coding promoter sequences in order to obtain a direct binding based gene regulatory network (see methods).

To construct the initial A. thaliana developmental gene regulatory network, we integrated three types of datasets. 
First, we incorporated conserved non-coding sequences within 2000 bp promoter regions of 17610 A. thaliana genes. 
Conserved non-coding promoter sequences were shown to be reliable predictors of regulatory elements controlling gene expression. 
Second, we added DNA binding predictions within these sequences for 120 transcription factors as provided by Van de Velde et al. 
In addition, we predicted binding within these promoter sequences for curated experimental DNA binding motifs of an additional set of 270 transcription factors. 
Therefore, we mapped the curated binding elements to all conserved non-coding promoter elements within the 17610 A. thaliana genes using the biocon- ductor TFBSTools package (p-value threshold p<0.001). 
As a result, we obtained a regulatory blueprint of 390 regulators and 17610 targets with 219000 link predictions. 


ird, we added an expression atlas of A. thaliana development57 comprising RNA samples from 83 tissues and developmental stages.  e expression data was used to derive a condition speci c co-expression network.  e expression dataset had already been normalized using Robust Multichip Averaging (RMA)57. Subsequently, we averaged tissue and developmental stage speci c exper- imental replicates. Finally, a variance based  ltering (using the gene lter R package) was applied to remove genes that exhibited little variation across all tissues and developmental stages.


General usage 
Step 1) Initial gene regulatory network inference (Figure 1 A) is based on fast random forest regression followed by DNA binding prediction map filtering (GRACE_initial_GRN_inference.R) 

Specific parameters and datasets 
th.grnCutoff <- 0.999 // initial cutoff for random forest regression on gene expression network (should result in not more than 20000 - 30000 links) 
n.cpus <- 2 // number of available cpus for parallel random forest regression 
m.expression <- ... // load expression matrix for grn inference - matrix - rows (genes) x cols (conditions) 
v.tfs <- ... // character vector of transcription factor genes 
df.dna_binding <- // load DNA binding set (dataframe with 2 columns) 

Step 2) Configure GRACE_load_datasets.R to set paths to the location of all needed datasets.


Step 3) Configure GRACE_pipeline_template.R, as described within the script, to run GRACE on the initial network using cofunctional network information and regulatory as well as cofunctional evidence data as training sets (Figure 1 B-C)

Specific parameters and datasets 
n.cpus <- 2 
beta <- 1 # equal emphasis on precision (atrm) and recall (bp coreg pairs) - introducing (novel) combined f-measure (physical + functional evidence, singularity handling, minimum size) 
b.normalize_precision <- TRUE 
b.jaccard = TRUE 
n.sets <- 20 
lambda.gridSearch <- c(0.01,seq(0.5,2.5,0.5)) 
th.percentage.hyperparameters = 0.99 # grid search 
max.call = 200 # simulated annealing 
n.models <- 100 # 100 models 
n.sample_size <- 0.632 # as in traditional bootrapping approaches 

Tutorials
In addition, we have prepared two standalone tutorials (for A. thaliana and D. melanogaster) in order to reproduce results 
GRACE_tutorial_athaliana.R - represents a less optimized version of GRACE - can be used to recompute the prioritized predictions (unpack A. thaliana dataset in GRACE folder) 
GRACE_tutorial_dmelanogaster.R - based on the current version of GRACE - recompute the prioritized predictions (unpack D. melanogaster dataset in GRACE folder) 

By default, the precomputed GRACE models are loaded for further processing, otherwise they can be recomputed. 

References
Datasets used within the D. melanogaster model
Marbach D, Roy S, Ay F, Meyer PE, Candeias R, Kahveci T, Bristow CA, Kellis M. Predictive regulatory models in Drosophila melanogaster by integrative inference of transcriptional networks. Genome Res. 2012 Jul;22(7):1334-49. 

GRACE algorithm and results
Banf M, and Rhee S. Enhancing gene regulatory network inference through data integration with markov random fields. accepted in Nature Scientific Reports.


For help or questions please contact: 
mbanf.research(at)gmail.com

