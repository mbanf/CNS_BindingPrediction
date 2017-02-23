# Using conserved non-coding sequences for robust transcription factor binding prediction

Accurate discovery of cis-regulatory elements remains a challenging, suffering from many false positive predictions. A large-scale comparison of cis-regulatory module detection approaches concluded that methods considering evolutionary conservation have a stronger predictive power than methods designed to be run on a single genome. In addition, conserved non-coding sequences (CNS) in promoter regions are reliable pointers to regulatory elements controlling gene expression. Recent examples of conserved regulatory element recovery in plant science include comparative, phylogenetic footprinting approaches. The key assumption in these approaches is that mutations within functional regions of genes are likely to accumulate more slowly than those in regions without sequence-specific function and, therefore, the comparison of conserved sequences from orthologous genes can indicate segments that might direct transcription. These approaches predicted hundreds of conserved non-coding sequences upstream of A. thaliana genes. The length, specific positioning and enrichment for transcription factor binding sites suggest these conserved non-coding sequences play a functional role in transcriptional regulation. 

A general framework to identify conserved non-coding sequences in multiple species is BLSSpeller, an algorithm that supports both alignment-free and alignment-based motif discovery in the promoter sequences of related species. Putative motifs are exhaustively enumerated as words over the IUPAC alphabet and screened for conservation using the branch length score. Additionally, a confidence score is established in a genome-wide fashion. (https://academic.oup.com/bioinformatics/article/31/23/3758/209257/BLSSpeller-exhaustive-comparative-discovery-of)

<br/>
*Matching DNA binding motifs to CNS followed by matching CNS to promoters to identify genes* <br/>

![Alt text](/CNS_DNA_binding.png?raw=true "CNS_DNA_binding")

To reduce computation time of CNS based binding prediction, we build an example R script that matches binding motifs to identify binding within a set of CNS regions instead of embedding CNS into promoters, followed by scanning each gene's promoter for motif similarity. We mapped the curated binding elements to all conserved non-coding promoter elements within the 17610 A. thaliana genes using the biocon- ductor TFBSTools package (p-value threshold p<0.001). Subsequently, we identify promoter regions that contain the CNS, therefore most likely bound by the transcription factors, as well as the corresponding genes. For gene association, we incorporated conserved non-coding sequences within 2000 bp promoter regions of 17610 A. thaliana genes. Our small example uses a collection of 270 A. thaliana experimentally validated transcription factor binding motifs, in the pwm (position weight matrix) format. Further, it uses a straight forward parallelization based on the foreach/doParallel packages, to process each motif independently.


While the growing number of sequenced genomes creates new opportunities for comparative approaches to motif discovery and CNS, a general limiting factor remains the availbility of DNA binding motifs. A valuable resource of experimental is cisBP (http://cisbp.ccbr.utoronto.ca). The database also curates a series of inferred binding motifs. Inference of binding motifs can be of great value in the absence of experimental binding information (see https://www.ncbi.nlm.nih.gov/pubmed/25215497)





For help or questions please contact: 
mbanf.research(at)gmail.com

