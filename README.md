# Conserved noncoding sequence transcription factor binding prediction


To construct the initial A. thaliana developmental gene regulatory network, we integrated three types of datasets. 
First, we incorporated conserved non-coding sequences within 2000 bp promoter regions of 17610 A. thaliana genes. 
Conserved non-coding promoter sequences were shown to be reliable predictors of regulatory elements controlling gene expression. 
Second, we added DNA binding predictions within these sequences for 120 transcription factors as provided by Van de Velde et al. 
In addition, we predicted binding within these promoter sequences for curated experimental DNA binding motifs of an additional set of 270 transcription factors. 
Therefore, we mapped the curated binding elements to all conserved non-coding promoter elements within the 17610 A. thaliana genes using the biocon- ductor TFBSTools package (p-value threshold p<0.001). 
As a result, we obtained a regulatory blueprint of 390 regulators and 17610 targets with 219000 link predictions. 


ird, we added an expression atlas of A. thaliana development57 comprising RNA samples from 83 tissues and developmental stages.  e expression data was used to derive a condition speci c co-expression network.  e expression dataset had already been normalized using Robust Multichip Averaging (RMA)57. Subsequently, we averaged tissue and developmental stage speci c exper- imental replicates. Finally, a variance based  ltering (using the gene lter R package) was applied to remove genes that exhibited little variation across all tissues and developmental stages.
