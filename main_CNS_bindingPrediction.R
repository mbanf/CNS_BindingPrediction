
setwd("/Users/michaelbanf/Documents/postdoctoral_work/Projects/CNS_Binding/")

# the main function
extract_CNS <- function(file.promoterSeq = "Datasets/TAIR10_upstream_3000_20101028.txt", df.CNS = df.CNS , th.promoter_length = 2000, n.cores = 10, load_from_file = FALSE){
  
  message("compute binding matrix")
  
  if(!load_from_file){
    
    #source("http://bioconductor.org/biocLite.R")
    #   biocLite("ChIPpeakAnno")
    #   biocLite("biomaRt")
    #   biocLite("Biostrings")
    #   install.packages("VennDiagram")
    library(ChIPpeakAnno)
    library(biomaRt)
    library(Biostrings)
    
    library(seqLogo)
    #seqLogo(pwm)
    
    library(seqinr)
    upstream_sequences <- read.fasta(file = file.promoterSeq, seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
    upstream_sequences <- lapply(upstream_sequences, function(m) {toupper(m)})
    df.promSequences <- DNAStringSet(unlist(upstream_sequences))
    
    
    #     library(TxDb.Athaliana.BioMart.plantsmart25)
    #     seqlevels(TxDb.Athaliana.BioMart.plantsmart25) <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9)
    #     transcriptCoordsByGene.GRangesList <- transcriptsBy (TxDb.Athaliana.BioMart.plantsmart25, by = "gene")
    #     isCircular(transcriptCoordsByGene.GRangesList) <- isCircular(BSgenome.Athaliana.TAIR.TAIR9)
    #     
    
    #biocLite("BSgenome.Athaliana.TAIR.TAIR9")
    library(BSgenome.Athaliana.TAIR.TAIR9)
    genome <- BSgenome.Athaliana.TAIR.TAIR9
    
    listMarts(host="plants.ensembl.org")
    ensmart = useMart("plants_mart",dataset="athaliana_eg_gene", host="plants.ensembl.org")
    
    #ensmart = useMart("plants.ensembl.org", dataset="athaliana_eg_gene")
    #ensmart = useMart('ENSEMBL_MART_PLANT', "athaliana_eg_gene")
    #df.CNS_2014 <- readRDS("partial_datasets/df.CNS_2014.rds") 
    
    g1.r <- BED2RangedData(df.CNS, header=FALSE)
    
    annotatedData = getAnnotation(ensmart, featureType = "TSS") # transcription start site
    annotatedPeaks = annotatePeakInBatch(g1.r, AnnotationData= annotatedData)
    
    df.cns_positions <- as.data.frame(annotatedPeaks, row.names=NULL) 
    df.cns_positions <- subset(df.cns_positions, df.cns_positions$insideFeature %in% c("upstream", "overlapStart"))
    df.cns_positions$space <- gsub("chr","", df.cns_positions$space)
    df.cns_positions["cons_score"] <- NA
    df.cns_positions["cns_sequence"] <- NA
    
    pb <- txtProgressBar(min = 0, max = nrow(df.cns_positions), style = 3)
    
    for(i in 1:nrow(df.cns_positions)){
      
      setTxtProgressBar(pb, i)
      
      
      #cat("Processing contigs ", round(i/nrow(df.cns_positions) * 100, digits = 2) , "%, ", "\r"); flush.console() 
      
      chr <- as.integer(df.cns_positions$space[i])
      start <- as.integer(df.cns_positions$start[i])
      end <- as.integer(df.cns_positions$end[i])
      
      df.CNS_2014.set <- df.CNS[df.CNS$Chr == chr,]
      df.CNS_2014.set <- df.CNS_2014.set[df.CNS_2014.set$start == start,]
      df.CNS_2014.set <- df.CNS_2014.set[df.CNS_2014.set$end == end,]
      
      cons_score <- as.integer(df.CNS_2014.set$cons_score)
      
      df.cns_positions$cons_score[i] <- cons_score
      
      genome_by_chrom <- DNAString(genome[[as.numeric(chr)]]) 
      
      df.cns_positions$cns_sequence[i] <- as.character(substring(genome_by_chrom,start,end))
    }
    
    close(pb)
    
    
    saveRDS(df.cns_positions, "df.cns_positions.rds")
    
    
  }else{
    
    #df.CNS_2014 <- subset(df.CNS_2014, !is.na(df.CNS_2014$cns_sequence))
    
    ## map the 400 tfs dna binding elements into (unique) cns sets, p-value (cutoff later )
    
    #df.cns_positions <- readRDS("../GRACE_prerelease/partial_datasets/df.cns_positions.rds")
    
    df.cns_positions <- readRDS("df.cns_positions.rds")
    
  }
  
  df.cns_positions <- subset(df.cns_positions, df.cns_positions$shortestDistance <= th.promoter_length) # subset to 2 kb promoter region 
  
  
}


# load CNS datasets - format
df.CNS <- read.table("AllFootPrintsFDR0.10_scores.bed", header = FALSE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
names(df.CNS) <- c("Chr", "start", "end", "val", "cons_score")



df.cns_positions <- extract_CNS(file.promoterSeq = "TAIR10_upstream_3000_20101028.txt", df.CNS = df.CNS , th.promoter_length = 2000, n.cores = 10, load_from_file = FALSE)
 