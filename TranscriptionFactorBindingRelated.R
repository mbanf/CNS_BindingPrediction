# TODO: Add comment
# 
# Author: michaelbanf
###############################################################################



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


is.installed <- function(lib)
{
  if (!require(lib,character.only = TRUE))
  {
    install.packages(lib,dep=TRUE)
    if(!require(lib,character.only = TRUE)) stop("Package not found")
  }
}


request_libraries <- function(){
  
  is.installed("graphics")
  is.installed("Matrix")
  is.installed("ggplot2")
  is.installed("ROCR")
  is.installed("reshape2")
  is.installed("CRF")
  is.installed("caTools")
  is.installed("plyr")
  is.installed("foreach")
  is.installed("doParallel")
  is.installed("igraph")
  is.installed("GenSA")
  is.installed("ranger")
  
  #detach("package:TFBSTools", unload=TRUE)
  
  if(!require(BSgenome.Athaliana.TAIR.TAIR9)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Athaliana.TAIR.TAIR9")
  }
  
  if(!require(ChIPpeakAnno)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("ChIPpeakAnno")
  }
  
  if(!require(biomaRt)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("biomaRt")
  }
  
  if(!require(Biostrings)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("Biostrings")
  }
  
  if(!require(seqinr)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("seqinr")
  }
  
  if(!require(seqLogo)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("seqLogo")
  }

  library(BSgenome.Athaliana.TAIR.TAIR9)

  library(ChIPpeakAnno)
  library(biomaRt)
  library(Biostrings)
  library(seqLogo)
  library(seqinr)
  
}



# get type A motifs - complete promoters
match_PWM_DNA_binding <- function(file.promoterSeq = "Datasets/TAIR10_upstream_1000_20101104.fasta", v.mode = "exp_only", n.cores = 10){

  library(TFBSTools)
  #detach("package:TFBSTools", unload=TRUE)
  library(Biostrings)
  library(seqinr)
  library(plyr)
  
  #library(seqLogo)
  #seqLogo(pwm)
  #library(seqinr)
  upstream_sequences <- read.fasta(file = file.promoterSeq, seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
  upstream_sequences <- lapply(upstream_sequences, function(m) {toupper(m)})
  df.promSequences <- DNAStringSet(unlist(upstream_sequences))
  
  genomic.acgt <- readRDS("genomic.acgt.rds")
  
  # lst.pwm.data.inf <- get_cell_and_pnasPaper_and_jaspar_pwmSet(v.mode = "inferred_only")
  lst.pwm.data.exp <- get_cell_and_pnasPaper_and_jaspar_pwmSet(v.mode = v.mode)
  
  df.motifs <- lst.pwm.data.exp$df.motifs
  df.motifs["pos"] <- seq(1:nrow(df.motifs))
  
  lst.pwm.data <- lst.pwm.data.exp
  
  ###
  
  v.dna.tfs <- unique(df.motifs$TF.locus)
  
  v.tfs.pwm <- unique(lst.pwm.data$df.motifs$TF.locus)
  lst.tfs.motifsets <- vector(mode = "list", length(v.tfs.pwm))
  names(lst.tfs.motifsets) <- v.tfs.pwm
  
  for(i in 1:length(v.tfs.pwm)){
    df.motifs.i <- subset(df.motifs, df.motifs$TF.locus %in% v.tfs.pwm[i])
    lst.tfs.motifsets[[i]] <- df.motifs.i$pos 
  }
  
  library(foreach)
  library(doParallel)
  #library(Biostrings)
  
  
  #lst.crfs <- foreach(l = 1:length(v.lambda), .packages=c("Matrix", "reshape2", "CRF", "caTools", "ROCR", "ggplot2", "igraph")) %dopar% { 
  
  #length(lst.tfs.motifsets)
  
  subject = df.promSequences
  
  #strt<-Sys.time()
  #lst.v.cns_dna_binding <- foreach(i = 1:length(lst.tfs.motifsets), .packages=c("Biostrings", "TFBSTools", "seqinr")) %dopar% { 
    
    #mat.dna_binding <- matrix(1, nrow = length(v.tfs.pwm), ncol = length(v.cns.tgs), dimnames = list(v.tfs.pwm, v.cns.tgs))
  for(i in 1:length(lst.tfs.motifsets)){
    
    #v.dna_binding <- rep(1, length(v.cns.tgs))
    #names(v.dna_binding) <- v.cns.tgs
    
    #v.cns_dna_binding <- rep(1, length(v.cns_positions))
    #names(v.cns_dna_binding) <- v.cns_positions
    
    pwmList <- lapply(lst.tfs.motifsets[[i]], function(j){
      profileMatrix <- lst.pwm.data$lst.pwm.motif[[j]];
      ID <- lst.pwm.data$df.motifs$TF.locus[j];
      name <- lst.pwm.data$df.motifs$TF.name[j];
      matrixClass <- lst.pwm.data$df.motifs$TF.family[j];
      PWMatrix(ID=ID, name=name, matrixClass=matrixClass,
               strand="+", bg=genomic.acgt, profileMatrix=profileMatrix)
    })
    
    cat("Processing... ", round(i/length(lst.tfs.motifsets) * 100, digits = 2) , "% (", length(pwmList), ")\r"); flush.console() 
    
    
    strt<-Sys.time()
    n.cpus <- min(n.cores, length(pwmList))
    
    cl <- makeCluster(n.cpus)
    registerDoParallel(cl)
    
    # search genome wide promoters for binding (high strigency - 90% min score) - compute pvalues - save
    lst.results <- foreach(k = 1:length(pwmList), .packages=c("Biostrings", "TFBSTools", "seqinr")) %dopar% { 
    
      sitesetList = searchSeq(pwmList[[k]], subject, seqname="", min.score="90%", strand="+") # match the pwm against the entire 1000 kb promoter region
      
      df.results <- as(sitesetList, "data.frame")
      df.results["p.values"] <- unlist(pvalues(sitesetList, type="TFMPvalue"))

      df.results
    }
    
    stopCluster(cl)
    print(Sys.time()-strt)
    
    
    # store as TF specific binding
    tf <- ID(pwmList[[1]])      
    df.dna_binding <- ldply(lst.results, data.frame)
    saveRDS(df.dna_binding, paste("tmp/df.dna_binding_",tf,".rds", sep =""))
 
    #subset(df.dna_binding, df.dna_binding$seqnames %in% names(which(table(df.dna_binding$seqnames) > 5)))
  }
  
}



# get type B motifs
identify_code_between_the_code <- function(){
  
  library(TFBSTools)
  #detach("package:TFBSTools", unload=TRUE)
  library(Biostrings)
  library(seqinr)
  #library(seqLogo)
  #seqLogo(pwm)
  #library(seqinr)
  
  upstream_sequences <- read.fasta(file = "Datasets/TAIR10_upstream_1000_20101104.fasta", seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
  upstream_sequences <- lapply(upstream_sequences, function(m) {toupper(m)})
  df.promSequences <- DNAStringSet(unlist(upstream_sequences))
  
  genomic.acgt <- readRDS("genomic.acgt.rds")
  
  
  # lst.pwm.data.inf <- get_cell_and_pnasPaper_and_jaspar_pwmSet(v.mode = "inferred_only")
  lst.pwm.data.exp <- get_cell_and_pnasPaper_and_jaspar_pwmSet(v.mode = "exp_only")
  
  df.motifs <- lst.pwm.data.exp$df.motifs
  df.motifs["pos"] <- seq(1:nrow(df.motifs))
  
  lst.pwm.data <- lst.pwm.data.exp
  
  ###
  
  v.dna.tfs <- unique(df.motifs$TF.locus)
  
  v.tfs.pwm <- unique(lst.pwm.data$df.motifs$TF.locus)
  lst.tfs.motifsets <- vector(mode = "list", length(v.tfs.pwm))
  names(lst.tfs.motifsets) <- v.tfs.pwm
  
  for(i in 1:length(v.tfs.pwm)){
    df.motifs.i <- subset(df.motifs, df.motifs$TF.locus %in% v.tfs.pwm[i])
    lst.tfs.motifsets[[i]] <- df.motifs.i$pos 
  }


  df.dna_binding <- data.frame()
  
  for(i in 1:length(lst.tfs.motifsets)){
  
    pwmList <- lapply(lst.tfs.motifsets[[i]], function(j){
      profileMatrix <- lst.pwm.data$lst.pwm.motif[[j]];
      ID <- lst.pwm.data$df.motifs$TF.locus[j];
      name <- lst.pwm.data$df.motifs$TF.name[j];
      matrixClass <- lst.pwm.data$df.motifs$TF.family[j];
      PWMatrix(ID=ID, name=name, matrixClass=matrixClass,
               strand="+", bg=genomic.acgt, profileMatrix=profileMatrix)
    })
    
    cat("Processing... ", round(i/length(lst.tfs.motifsets) * 100, digits = 2) , "% (", length(pwmList), ")\r"); flush.console() 

    # store as TF specific binding
    tf <- ID(pwmList[[1]])      
    df.dna_binding.tmp <- readRDS(paste("tmp/df.dna_binding_",tf,".rds", sep =""))
    df.dna_binding.tmp <- subset(df.dna_binding.tmp, df.dna_binding.tmp$p.values < 0.0001)
    df.dna_binding <- rbind(df.dna_binding, df.dna_binding.tmp[,c(1,4,9,12)])
    
  }
  
  
  # get consensus motifs - construct set of fasta sequences (for repeating homo and heterogeneous motifs)
  
  # count the number of TFs per 
  v.tgs.dna_bound <- as.character(unique(df.dna_binding$seqnames))
  mat.dna_binding.repeats <- matrix(0, nrow = length(v.tfs.pwm), ncol = length(v.tgs.dna_bound), dimnames = list(v.tfs.pwm, v.tgs.dna_bound)) 
  
  for(i in 1:length(v.tfs.pwm)){
    cat("Processing... ", round(i/length(v.tfs.pwm) * 100, digits = 2) , "% \r"); flush.console() 
    df.dna_binding.i <- subset(df.dna_binding, df.dna_binding$ID == v.tfs.pwm[[i]])
    tb.tgs.i <- table(as.character(df.dna_binding.i$seqnames))
    mat.dna_binding.repeats[v.tfs.pwm[i], names(tb.tgs.i)] <- as.numeric(tb.tgs.i)
  }
    
  # identify motifs with at least 3 repeats per sequence (set a minimum threshold - incorporated distance distribution per TF)
  idx.repeats <- which(mat.dna_binding.repeats >= 4, arr.ind = TRUE)
  
  #v.tgs.repeats <- unique(idx.repeats[,2])
  #v.tgs.repeats <- v.tgs.dna_bound[v.tgs.repeats]
  #tb.tfs <- table(idx.repeats[,1])
  
  
  v.window_sizes <- c(10,30,50) #seq(10,150,10)
  mat.dna_binding.peaks <- matrix(0, nrow = length(v.tfs.pwm), ncol = length(v.tgs.dna_bound), dimnames = list(v.tfs.pwm, v.tgs.dna_bound)) 
  
  for(i in 1:nrow(idx.repeats)){
    cat("Processing... ", round(i/nrow(idx.repeats) * 100, digits = 2) , "% \r"); flush.console() 
  
    df.dna_binding.i <- subset(df.dna_binding, df.dna_binding$ID == v.tfs.pwm[idx.repeats[i,1]] & df.dna_binding$seqnames == v.tgs.dna_bound[idx.repeats[i,2]])
    df.dna_binding.i <- unique(df.dna_binding.i)
    
    df.dna_binding.i <- df.dna_binding.i[order(df.dna_binding.i$start),]
    
    # identify binding clusters - non-overlapping and not-dispersed
    v.motif_seq <- df.dna_binding.i$start
    n.seq <- nchar(df.dna_binding.i$siteSeqs) # length of the motifs 
    
    d.motif_seq <- v.motif_seq[2:length(v.motif_seq)] - v.motif_seq[1:(length(v.motif_seq)-1)]        
    n.seq <- n.seq[1:(length(v.motif_seq)-1)]
    
    # prior filter - no overlapping
    idx.overlaps <- (which(d.motif_seq < n.seq))
    
    n.peaks <- numeric(length(v.window_sizes))
    #v.spreads <- numeric(length(v.window_sizes))
    #v.pos <- numeric(length(v.window_sizes))
    
    for(w in 1:length(v.window_sizes)){
      
      # number of neighbored motifs within window size
      idx.selection <- (which(d.motif_seq < v.window_sizes[w]))
      
      # peak per run 
      n.peaks[w] <- length(idx.selection)
      # v.spreads[w] <- (v.motif_seq[max(idx.selection) + 1] - v.motif_seq[min(idx.selection)]) + n.seq[idx.selection]
      
      # center
      # v.pos[w] <- v.motif_seq[min(idx.selection)] + spread/2
      
    }
    
    # average detections 
    n.peak <- mean(n.peaks)
    # n.spread <- mean(v.spreads)
    # v.pos <- mean(v.pos)
    
    # remove all overlapping double counts
    n.peak <- n.peak - length(idx.overlaps)
    
    mat.dna_binding.peaks[idx.repeats[i,1], idx.repeats[i,2]] <- n.peak
    
    
  }
  
  saveRDS(mat.dna_binding.peaks, "mat.dna_binding.peaks.rds")

}


# load a multiset of binding motifs - last curated spring 2016
# select experimental or inferred
get_cell_and_pnasPaper_and_jaspar_pwmSet <- function(v.mode = "exp_only"){
  
  print("prepare PWM motifs and mappings")
  library(fume)
  ## create motif - gene mapping
  df.motifs <- data.frame(Motif.id = character(), TF.locus = character(), TF.name = character(), TF.family = character(), TF.subfamily = character(), src = character(), stringsAsFactors = FALSE)
  #lst.motifs <- vector(mode = "list")
  
  # Cell paper
  files <- list.files("Datasets/novelTFmotifs/Ath_CellPaper/pwms_all_motifs/")
  motif.mapping <- read.table("Datasets/novelTFmotifs/Ath_CellPaper/TF_Information_all_motifs_plus.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
  
  motif.mapping$DBID <- gsub("\\-.*", "", motif.mapping$DBID)
  motif.mapping <- subset(motif.mapping, motif.mapping$TF_Status != "N")
  idx <- which(grepl("^AT", motif.mapping$DBID) == TRUE)
  motif.mapping <- motif.mapping[idx,]
  
  if(v.mode == "inferred_only"){
    motif.mapping <- subset(motif.mapping, motif.mapping$TF_Status == "I")  
  }else if(v.mode == "exp_only"){
    motif.mapping <- subset(motif.mapping, motif.mapping$TF_Status == "D")  
  }
  
  idx <- which(gsub(".txt", "", files) %in% motif.mapping$Motif_ID)
  files <- files[idx]
  n.cell <- length(files)
  
  for(i in 1:n.cell){  
    motif.id <- substr(files[i], 1, nchar(files[i])-4) 
    #lst.motifs[[i]] <- vector(mode = "list", length = 2)
    #lst.motifs[[i]][[1]] <- motif.id
    idx <- match(motif.id, motif.mapping$Motif_ID)
    #lst.motifs[[i]][[2]] <- motif.mapping$DBID[idx]
    newrow <- data.frame(Motif.id = motif.id, TF.locus = motif.mapping$DBID[idx], TF.name = motif.mapping$TF_Name[idx], TF.family = motif.mapping$Family_Name[idx], TF.subfamily = "", src = "Cell", stringsAsFactors = FALSE)
    names(newrow) = c("Motif.ID", "TF.locus", "TF.name", "TF.family", "TF.subfamily", "src")
    df.motifs <- rbind(df.motifs, newrow)
  }
  
  ## map motifs
  lst.pwm.motif <- vector(mode = "list", length = length(files))
  for(i in 1:n.cell){ 
    #print(paste("TF Motif ", i, "of", length(files)))
    lst.pwm.motif[[i]] <- read.table(paste("Datasets/novelTFmotifs/Ath_CellPaper/pwms_all_motifs/", files[i], sep = ""), header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
    
    lst.pwm.motif[[i]] <- lst.pwm.motif[[i]][,-1]
    lst.pwm.motif[[i]] <- t(as.matrix(lst.pwm.motif[[i]]))
    names(lst.pwm.motif)[i] <- df.motifs$Motif.ID[i]
    colnames(lst.pwm.motif[[i]]) <- as.character(seq(1:ncol(lst.pwm.motif[[i]])))
  }
  
  ### --------- 
  if(v.mode == "inferred_only"){
    return(list(lst.pwm.motif=lst.pwm.motif, df.motifs=df.motifs))
  }
  
  # PNAS paper
  df.pwms.pnas <- read.table("Datasets/novelTFmotifs/pnas_paper/all_pnas_pwms.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
  
  n.pnas <- nrow(df.pwms.pnas)/4
  
  idx <- 1
  for(i in 1:n.pnas){
    motif.id <- paste("motif_pnas",i, sep = "")
    vec.rows <-  c(idx, idx + 1, idx + 2, idx + 3)
    pnas.motif <- as.matrix(df.pwms.pnas[vec.rows, seq(2,11)])
    rownames(pnas.motif) <- c("A","C","G","T")
    colnames(pnas.motif) <- as.character(seq(1:ncol(pnas.motif)))
    #lst.motifs[[i]] <- vector(mode = "list", length = 2)
    #lst.motifs[[i]][[1]] <- motif.id
    lst.pwm.motif <- lappend(lst.pwm.motif, pnas.motif)
    names(lst.pwm.motif)[[(n.cell + i)]] <- motif.id #df.pwms.pnas$TF.locus[idx]
    
    # check for identical list elements
    newrow <- data.frame(Motif.id = motif.id, TF.locus = df.pwms.pnas$TF.locus[idx], TF.name = df.pwms.pnas$TF.name[idx], TF.family = df.pwms.pnas$TF.Family[idx], TF.subfamily = df.pwms.pnas$TF.Subfamily[idx], src = "PNAS")
    names(newrow) = c("Motif.ID", "TF.locus", "TF.name", "TF.family", "TF.subfamily", "src")
    df.motifs <- rbind(df.motifs, newrow)
    idx <- idx + 4
  }
  
  
  #biocLite("MotifDb")
  library (MotifDb)
  library (MotIV)
  library (seqLogo)
  
  library(rtfbs)
  
  #require("rtfbs")
  #install.packages("../Datasets/novelTFmotifs/Jaspar2016/rtfbs_0.3.4.tar", repos = NULL, type="source")
  
  df.motifs.jaspar_2016 <- read.table("Datasets/novelTFmotifs/Jaspar2016/motif2016.txt", header = FALSE, sep = " ", stringsAsFactors = FALSE)
  df.motifs.jaspar_2016 <- df.motifs.jaspar_2016[,2:6]
  names(df.motifs.jaspar_2016) <- c("motif.id", "tf.name","tf.id", "tf.fam","species")
  
  lst.jaspar_2016 <- read.pwm("Datasets/novelTFmotifs/Jaspar2016/JASPAR_CORE_REDUNDANT_2016_plants.meme")
  lst.jaspar_2016 <- lapply(lst.jaspar_2016, exp)
  names(lst.jaspar_2016) <- df.motifs.jaspar_2016$tf.id
  
  idx.ath <- which(df.motifs.jaspar_2016$species == "ath")
  
  lst.jaspar_2016 <- lst.jaspar_2016[idx.ath]
  df.motifs.jaspar_2016 <- df.motifs.jaspar_2016[idx.ath,]
  df.motifs.jaspar_2016$tf.id <- toupper(df.motifs.jaspar_2016$tf.id)
  
  
  #   mFile <- system.file("../Datasets/novelTFmotifs/Jaspar2014/pfm_plants_jaspar2016.txt", package="seqLogo")
  #   m <- read.table(mFile)
  #   pwm <- makePWM(m)
  
  
  for(i in 1:nrow(df.motifs.jaspar_2016)){  
    if(df.motifs.jaspar_2016$species[i] == "ath"){
      newrow <- data.frame(Motif.id = df.motifs.jaspar_2016$motif.id[i], TF.locus = df.motifs.jaspar_2016$tf.id[i], 
                           TF.name = df.motifs.jaspar_2016$tf.name[i], TF.family = df.motifs.jaspar_2016$tf.fam[i],
                           TF.subfamily = "", src = "Jaspar2016", stringsAsFactors = FALSE)
      names(newrow) = c("Motif.ID", "TF.locus", "TF.name", "TF.family", "TF.subfamily", "src")
      df.motifs <- rbind(df.motifs, newrow)
    }
  }
  
  lst.jaspar_2016 <- lapply(lst.jaspar_2016, t)
  lst.pwm.motif <- c(lst.pwm.motif, lst.jaspar_2016)
  
  return(list(lst.pwm.motif=lst.pwm.motif, df.motifs=df.motifs))
}

  
  




extract_CNS2014 <- function(file.promoterSeq = "Datasets/TAIR10_upstream_3000_20101028.txt", df.CNS = df.CNS , th.promoter_length = 2000, n.cores = 10, load_from_file = FALSE){
  
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
    
    
  message("compute binding matrix")
  
  if(!load_from_file){
    
    library(TFBSTools)
    #detach("package:TFBSTools", unload=TRUE)
    library(Biostrings)
    library(seqinr)
    
    genomic.acgt <- readRDS("genomic.acgt.rds")
    
    # lst.pwm.data.inf <- get_cell_and_pnasPaper_and_jaspar_pwmSet(v.mode = "inferred_only")
    lst.pwm.data.exp <- get_cell_and_pnasPaper_and_jaspar_pwmSet(v.mode = "exp_only")
    
    df.motifs <- lst.pwm.data.exp$df.motifs
    df.motifs["pos"] <- seq(1:nrow(df.motifs))
    
    lst.pwm.data <- lst.pwm.data.exp
    
    ###
    
    v.dna.tfs <- unique(df.motifs$TF.locus)
    v.cns.tgs <- unique(df.cns_positions$feature)
    
    mat.grn.cns_novel <- matrix(0, nrow = length(v.dna.tfs), ncol = length(v.cns.tgs), dimnames = list(v.dna.tfs, v.cns.tgs))
    
    ## background distributions 
    #Background motif distributions for a custom set of PWMs can be easily calculated for all model organisms. 
    #We will illustrate this by creating a new lognormal background for two de-novo motifs in Drosophila.
    #To load in the motifs the package provides functions to read standard JASPAR and TRANSFAC formats.
    #library(PWMEnrich.Dmelanogaster.background)
    
    v.tfs.pwm <- unique(lst.pwm.data$df.motifs$TF.locus)
    lst.tfs.motifsets <- vector(mode = "list", length(v.tfs.pwm))
    names(lst.tfs.motifsets) <- v.tfs.pwm
    
    for(i in 1:length(v.tfs.pwm)){
      df.motifs.i <- subset(df.motifs, df.motifs$TF.locus %in% v.tfs.pwm[i])
      lst.tfs.motifsets[[i]] <- df.motifs.i$pos 
    }
    
    #th.p_value <- 5e-2
    #v.promoter_length <- 2000
    
    
    #library(biomaRt)
    library(Biostrings)
    #biocLite("BSgenome.Athaliana.TAIR.TAIR9")
    #     library(BSgenome.Athaliana.TAIR.TAIR9)
    #     genome <- BSgenome.Athaliana.TAIR.TAIR9
    #     #source("http://bioconductor.org/biocLite.R")
    #     #biocLite("TxDb.Athaliana.BioMart.plantsmart25")
    #     library(TxDb.Athaliana.BioMart.plantsmart25)
    #     
    #     seqlevels(TxDb.Athaliana.BioMart.plantsmart25) <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9)
    #     transcriptCoordsByGene.GRangesList <- transcriptsBy (TxDb.Athaliana.BioMart.plantsmart25, by = "gene")
    #     isCircular(transcriptCoordsByGene.GRangesList) <- isCircular(BSgenome.Athaliana.TAIR.TAIR9)
    #     #transcriptCoordsByGene.GRangesList <- renameSeqlevels( transcriptCoordsByGene.GRangesList, c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrM", "ChrC") )
    #     transcriptCoordsByGene.GRangesList <- transcriptCoordsByGene.GRangesList[names(transcriptCoordsByGene.GRangesList) %in% v.tgs.selection]
    #     
    #     df.promSequences <- getPromoterSeq(transcriptCoordsByGene.GRangesList, genome, upstream=v.promoter_length, downstream=0)
    #     
    #     
    #length(which(v.tgs_selection %in% names(df.promSequences)))
    #df.promSequences <- df.promSequences[v.tgs_selection,]
    
    
    v.cns_positions <- unique(df.cns_positions$cns_sequence)
    
    
    library(foreach)
    library(doParallel)
    #library(Biostrings)
    
    
    cl<-makeCluster(n.cores)
    registerDoParallel(cl)
    
    #lst.crfs <- foreach(l = 1:length(v.lambda), .packages=c("Matrix", "reshape2", "CRF", "caTools", "ROCR", "ggplot2", "igraph")) %dopar% { 
    
    #length(lst.tfs.motifsets)
    
    strt<-Sys.time()
    lst.v.cns_dna_binding <- foreach(i = 1:length(lst.tfs.motifsets), .packages=c("Biostrings", "TFBSTools", "seqinr")) %dopar% { 
      
      #mat.dna_binding <- matrix(1, nrow = length(v.tfs.pwm), ncol = length(v.cns.tgs), dimnames = list(v.tfs.pwm, v.cns.tgs))
      #for(i in 1:length(lst.tfs.motifsets)){
      cat("Processing... ", round(i/length(lst.tfs.motifsets) * 100, digits = 2) , "%", "\r"); flush.console() 
      #v.dna_binding <- rep(1, length(v.cns.tgs))
      #names(v.dna_binding) <- v.cns.tgs
      
      v.cns_dna_binding <- rep(1, length(v.cns_positions))
      names(v.cns_dna_binding) <- v.cns_positions
      
      pwmList <- lapply(lst.tfs.motifsets[[i]], function(j){
        profileMatrix <- lst.pwm.data$lst.pwm.motif[[j]];
        ID <- lst.pwm.data$df.motifs$TF.locus[j];
        name <- lst.pwm.data$df.motifs$TF.name[j];
        matrixClass <- lst.pwm.data$df.motifs$TF.family[j];
        PWMatrix(ID=ID, name=name, matrixClass=matrixClass,
                 strand="+", bg=genomic.acgt, profileMatrix=profileMatrix)
      })
      
      tf <- ID(pwmList[[1]])    
      
      for(k in 1:length(pwmList)){
        #for(j in 1:length(v.cns_positions)){
        #cat("Processing... ", round(j/length(v.cns_positions) * 100, digits = 2) , "%", "\r"); flush.console() 
        
        # print(paste(j, "of", length(df.promSequences)))
        #hits <- matchPWM(lst.pwm.motif[[i]], DNAString(toString(df.promSequences[[j]][[1]])) , with.score = TRUE)
        #hits <- matchPWM(pcm, DNAString(toString(df.promSequences[[j]][[1]])) , with.score = TRUE)
        #hits <- matchPWM(pcm, DNAString(toString(df.promSequences[[j]][[1]])) , with.score = TRUE)
        #tg <- df.cns_positions$feature[j]
        #subject = DNAString(toString(v.cns_positions[j])) 
        subject = DNAStringSet(v.cns_positions)
        
        sitesetList = searchSeq(pwmList[[k]], subject, seqname="", min.score="60%", strand="+")
        
        for(j in 1:length(sitesetList)){
          if(length(sitesetList[[j]]) > 0){
          #if(length(siteset) > 0){
            ## calculate the empirical p-values of the scores
            p.value.min <- min(pvalues(sitesetList[[j]], type="TFMPvalue"))
            #print(p.value.min)
            #if(p.value.min <= th.p_value){
            #mat.dna_binding[tf, tg] <-  min(p.value.min, mat.dna_binding[tf, tg])  
            #v.dna_binding[tg] <- min(p.value.min, v.dna_binding[tg])
            v.cns_dna_binding[j] <- min(p.value.min, v.cns_dna_binding[j])
            #}
          }
        }
        #}
      }
      #set <- PWMatrixList(pwmList[[1]], pwmList[[2]])
      #sitesetList = searchSeq(set, subject, seqname="seq1", min.score="60%", strand="*")
      #v.dna_binding
      #mat.dna_binding
      #v.dna_binding
      
      v.cns_dna_binding
    }
    
    # double check with coexpression data (stress and development)
    #which(mat.dna_binding[,tg] < 0.01)
    stopCluster(cl)
    print(Sys.time()-strt)
    
    names(lst.v.cns_dna_binding) <- names(lst.tfs.motifsets)
    saveRDS(lst.v.cns_dna_binding, "lst.v.cns_dna_binding.exp.rds")
    
  }else{
    lst.v.cns_dna_binding <- readRDS("lst.v.cns_dna_binding.exp.rds")
  }
  
  ###
    
  message("compute binding matrix")
  
    th.p.value <- 0.001
    
    mat.cns.pvalue <- matrix(1, nrow = length(v.dna.tfs), ncol = length(v.cns.tgs), dimnames = list(v.dna.tfs, v.cns.tgs))
    mat.cns.cscore <- matrix(0, nrow = length(v.dna.tfs), ncol = length(v.cns.tgs), dimnames = list(v.dna.tfs, v.cns.tgs))
    mat.cns.tts <- matrix(5000, nrow = length(v.dna.tfs), ncol = length(v.cns.tgs), dimnames = list(v.dna.tfs, v.cns.tgs))
    
    for(i in 1:length(lst.v.cns_dna_binding)){
      
      cat("Processing... ", round(i/length(lst.v.cns_dna_binding) * 100, digits = 2) , "%", "\r"); flush.console() 
      tf <- names(lst.v.cns_dna_binding)[i]
      cns.seq.i <- lst.v.cns_dna_binding[[i]]
      #cns.seq.i <- cns.seq.i[cns.seq.i < 1]
      
      
      cns.seq.i <- cns.seq.i[cns.seq.i <= th.p.value]
      
      if(length(cns.seq.i) > 0){
        
        
        df.cns.tgs <- subset(df.cns_positions, df.cns_positions$cns_sequence %in% names(cns.seq.i))
        
        for(j in 1:length(cns.seq.i)){
          df.cns.tgs.j <- subset(df.cns.tgs, df.cns.tgs$cns_sequence == names(cns.seq.i[j]))
          tgs.j <- df.cns.tgs.j$feature
          p.value <- cns.seq.i[j]
          #mat.cns.pvalue[tf, tgs.j]  <- pmin(mat.cns.pvalue[tf, tgs.j], p.value)
          
          #if(p.value <= th.p.value){
            cscore.j <- unique(df.cns.tgs.j$cons_score)
            mat.cns.cscore[tf, tgs.j]  <- pmax(mat.cns.cscore[tf, tgs.j], cscore.j)          
          #}
          
          
  #         if(p.value <= th.p.value){
  #           tts.j <- df.cns.tgs.j$shortestDistance
  #           mat.cns.tts[tf, tgs.j]  <- pmin(mat.cns.tts[tf, tgs.j], tts.j)  
  #         }
          #         
          
        }
      }
    }
    
    saveRDS(mat.cns.cscore, "mat.cns.cscore_p0001.rds")
    
      #saveRDS(mat.cns.pvalue, "partial_datasets/mat.cns.pvalue.rds")
#     saveRDS(mat.cns.cscore, "partial_datasets/mat.cns.cscore.rds")
#     saveRDS(mat.cns.tts,    "partial_datasets/mat.cns.tts.rds")
#     
#     
  
    
    #mat.cns.cscore <- readRDS("mat.cns.cscore.rds")
    
    mat.cns.cscore
    
} 
    
    #saveRDS(mat.dna_binding, "partial_datasets/mat.dna_binding.rds")
    
    
    
    ### run rapid for all regulators but only these targets on the development set 
    
    v.cns.tgs
    
    
    
    #length(intersect(v.tfs, rownames(mat.cns.tts)))

    
    
    
    tfs.tg <- names(which(mat.dna_binding[,tg] < 0.01))
    
    
    
    quantile(mat.genie3, 0.999)
    
    round(mat.genie3[tfs.tg, tg], digits = 5)
    
    ### combine with the vandepoele network (max value if multiple)
    
    
    
    
    
    
    
    ## map dna binding sets 
    
    
    
    
    
    perform_TFBS_mapping_CNS2014 <- function(v.tgs, v.th.bind_score = 0, n.cores = 15){
      
      print("compute scored CNS2014-GRN")
      
      lst.motifs <- get_cell_and_pnas_Paper_PWMs()
      lst.pwm.motif <- lst.motifs[[1]]
      df.motifs <- lst.motifs[[2]]
      
      if(b.CALC){
        df.CNS_2014 <- readRDS("../Datasets/df.CNS_2014.rds") 
      }else{
        df.CNS_2014 <- readRDS("Datasets/df.CNS_2014.rds") 
      }
      df.CNS_2014 <- subset(df.CNS_2014, df.CNS_2014$Target %in% v.tgs)
      # df.CNS_2014 <- readRDS(paste("Datasets/novelTFmotifs/df.CNS_2014_2000kb.rds", sep = ""))
      
      #df.CNS_2014 <- subset(df.CNS_2014, df.CNS_2014$end_dist_to_TSS < 1000)  
      vec.cns.sequences <- as.character(unique(df.CNS_2014$cns_sequence))
      n.cns.seq <- length(vec.cns.sequences)
      
      #system.time(for(i in 1:1){#length(lst.pwm.motif)){ 
      #foreach(i = 1:2) %dopar% { 
      cl<-makeCluster(n.cores)
      registerDoParallel(cl)
      lst.dna_binding <- foreach(i = 1:length(lst.pwm.motif), .packages=c("Biostrings")) %dopar% { 
        #print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
        df.grn.CNS2014.TFB_map <- data.frame(TF = character(), cns_sequence = character(), TF.Motif = character(), 
                                             bind_score = numeric(), cmp_bind_score = numeric(), BS_pos_to_TSS = numeric(), 
                                             stringsAsFactors = FALSE)    
        names(df.grn.CNS2014.TFB_map) <- c("TF", "cns_sequence", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_in_CNS")
        pcm <- round(100 * lst.pwm.motif[[i]])
        for(j in 1:n.cns.seq){
          hits <- matchPWM(pcm,  vec.cns.sequences[j], with.score = TRUE)
          nhits <- length(hits)  
          if(nhits >= 1){
            cmp_bind_score <- min(mcols(hits)$score / maxScore(pcm)) # should be >= 0.8
            motif.score <- mcols(hits)$score / 100
            if(cmp_bind_score >= 0.8){
              for(k in 1:nhits){
                if(motif.score[k] >= v.th.bind_score){
                  newrow <- data.frame(TF =  as.character(df.motifs$TF.locus[i]), 
                                       cns_sequence = vec.cns.sequences[j],
                                       TF.Motif =  names(lst.pwm.motif)[i],
                                       bind_score = motif.score[k], 
                                       cmp_bind_score = cmp_bind_score,
                                       BS_pos_in_CNS = abs(end(hits)[k] + start(hits)[k])/2,
                                       stringsAsFactors = FALSE)
                  df.grn.CNS2014.TFB_map <- rbind(df.grn.CNS2014.TFB_map, newrow)   
                }
              }
            }
          }
        }
        names(df.grn.CNS2014.TFB_map) <- c("TF", "cns_sequence", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_in_CNS")
        #saveRDS(df.grn.CNS2014.TFB_map, paste("Datasets/CNS_GRNS/tmp/df.cns2014_grn_map_",i,".rds", sep = ""))
        df.grn.CNS2014.TFB_map
      }
      stopCluster(cl)
      saveRDS(lst.dna_binding, "Datasets/novelTFmotifs/lst.dna_binding.rds")
      
      
      df.grn.CNS2014.TFB_map <- do.call("rbind", lst.dna_binding)
      
      # combine
      #   df.grn.CNS2014.TFB_map <- data.frame(TF = character(), cns_sequence = character(), TF.Motif = character(), 
      #                                        bind_score = numeric(), cmp_bind_score = numeric(),
      #                                        stringsAsFactors = FALSE)      
      #   for(i in 1:length(lst.pwm.motif)){ 
      #     print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
      #     df.grn.CNS2014.TFB_map <- rbind(df.grn.CNS2014.TFB_map, readRDS(paste("Datasets/CNS_GRNS/tmp/df.cns2014_grn_map_",i,".rds", sep = "")))
      #   }
      
      #names(df.grn.CNS2014.TFB_map) <- c("TF", "cns_sequence", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_in_CNS") 
      #saveRDS(df.grn.CNS2014.TFB_map, paste("Datasets/CNS_GRNS/df.grn.CNS2014.TFB_map.rds", sep = ""))
      #saveRDS(df.grn.CNS2014.TFB_map, paste("Datasets/CNS_GRNS/df.grn.CNS2014.TFB_map_complete.rds", sep = ""))
      
      
      
      # finalize
      #df.grn.CNS2014.TFB_map <- readRDS(paste("Datasets/CNS_GRNS/df.grn.CNS2014.TFB_map_complete.rds", sep = ""))
      df.grn.CNS2014.TFB_map <- df.grn.CNS2014.TFB_map[,c(-3,-5)]
      df.grn.CNS2014.TFB_map <- unique(df.grn.CNS2014.TFB_map)
      
      #df.CNS_2014["BS_pos_to_TSS"] <-  apply(df.CNS_2014[,c('start_dist_to_TSS','end_dist_to_TSS')], 1, function(x) mean(x) )
      #df.CNS_2014 <- df.CNS_2014[, c(-1,-2,-3,-8,-9,-10)]
      
      df.CNS_grn <- merge(df.grn.CNS2014.TFB_map, df.CNS_2014, by = "cns_sequence")
      
      
      ## compare for enrichment (with gold standard links)
      
      
      
      
      
      
      
      
      
      
    }
    
    
    
    # #     annotatedData = getAnnotation(ensmart, featureType = "TSS")
    # #     annotatedPeaks = annotatePeakInBatch(g1.r, AnnotationData= annotatedData)
    # #     #annotatedPeaks <- annotatedPeaks[order(rownames(annotatedPeaks)),]
    # #     annotatedPeaks.sset  = 	annotatedPeaks[!is.na(annotatedPeaks$distancetoFeature)  & 
    # #                                              annotatedPeaks$fromOverlappingOrNearest == "NearestStart" & 
    # #                                              annotatedPeaks$distancetoFeature < 0 &
    # #                                              abs(annotatedPeaks$distancetoFeature) < v.promoter_length,]
    # 
    # df.CNS_2014["Target"] <- NA
    # df.CNS_2014["cns_sequence"] <- NA
    # df.CNS_2014["start_dist_to_TSS"] <- NA
    # df.CNS_2014["shortest_dist_to_TSS"] <- NA
    # df.CNS_2014["feature_position"] <- NA
    # 
    # for(i in 1:nrow(annotatedPeaks.sset)){
    #   print(paste(i, " of", nrow(annotatedPeaks.sset)))
    #   nr.row <- as.numeric(gsub("\\ .*", "",  rownames(annotatedPeaks.sset)[i]))
    #   chr <- df.CNS_2014$Chr[nr.row]
    #   start <- df.CNS_2014$start[nr.row]
    #   end  <- df.CNS_2014$end[nr.row]
    #   df.CNS_2014$Target[nr.row] <- as.character(annotatedPeaks.sset$feature[i])
    #   genome_by_chrom <- DNAString(genome[[chr]]) 
    #   df.CNS_2014$cns_sequence[nr.row] <- as.character(substring(genome_by_chrom,start,end))
    #   df.CNS_2014$start_dist_to_TSS[nr.row] <- abs(as.numeric(annotatedPeaks.sset$distancetoFeature[i]))
    #   df.CNS_2014$shortest_dist_to_TSS[nr.row]   <- abs(as.numeric(annotatedPeaks.sset$shortestDistance[i]))
    #   df.CNS_2014$feature_position[nr.row]   <- as.character(annotatedPeaks.sset$insideFeature[i])
    # }
    # df.CNS_2014 <- subset(df.CNS_2014, !is.na(df.CNS_2014$cns_sequence))
    # 
    
    
    
    
    
    
    
    
    
    
    
    
    
    ##################
    
    # load grn 
    
    v.tgs.selection <- colnames(lst.grn[[1]])
    
    
    
    
    
    
    # JASPAR 2014
    # map tf name to .. 
    # 
    #source("http://bioconductor.org/biocLite.R")
    #biocLite("PWMEnrich")
    #detach("package:PWMEnrich", unload=TRUE)
    #library(PWMEnrich)
    # 
    # registerCoresPWMEnrich(10)
    # useBigMemoryPWMEnrich(TRUE)
    
    library(TFBSTools)
    #detach("package:TFBSTools", unload=TRUE)
    library(Biostrings)
    library(seqinr)
    # upstream_sequences <- read.fasta(file = "../Datasets/TAIR10_upstream_1000_20101104.txt", seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
    # #upstream_sequences <- read.fasta(file = "../Datasets/TAIR10_upstream_3000_20101028.txt", seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
    # upstream_sequences <- lapply(upstream_sequences, function(m) {toupper(m)})
    # df.promSequences <- DNAStringSet(unlist(upstream_sequences))
    
    
    # background distritbution A,C,G,T
    #genomic.acgt = makePriors(df.promSequences,1)
    #saveRDS(genomic.acgt, "genomic.acgt.rds")
    
    genomic.acgt <- readRDS("genomic.acgt.rds")
    #pwms.denovo = toPWM(motifs.denovo, prior=genomic.acgt)
    #library(BSgenome.Athaliana.TAIR.TAIR9)
    #bg.denovo = makeBackground(lst.pwm.data$lst.pwm.motif, organism=BSgenome.Athaliana.TAIR.TAIR9, type="logn", quick=TRUE)
    #bg.denovo = makeBackground(lst.pwm.data$lst.pwm.motif, type="empirical", quick=TRUE, bg.seq = df.promSequences)
    # 
    # motifEnrichment(DNAString("TGCATCAAGTGTGTAGTG"), bg.pval)
    # 
    # biocLite("PWMEnrich.Dmelanogaster.background")
    # library(PWMEnrich.Dmelanogaster.background)
    
    
    
    
    
    
    #lst.pwm.data.inf <- get_cell_and_pnasPaper_and_jaspar_pwmSet(v.mode = "inferred_only")
    lst.pwm.data.exp <- get_cell_and_pnasPaper_and_jaspar_pwmSet(v.mode = "exp_only")
    
    df.motifs <- lst.pwm.data.exp$df.motifs
    df.motifs["pos"] <- seq(1:nrow(df.motifs))
    
    lst.pwm.data <- lst.pwm.data.exp
    
    ###
    
    v.dna.tfs <- unique(df.motifs$TF.locus)
    
    mat.grn.cns_novel <- matrix(0, nrow = length(v.dna.tfs), ncol = length(df.cns_positions$feature)
                                
                                
                                ## background distributions 
                                #Background motif distributions for a custom set of PWMs can be easily calculated for all model organisms. 
                                #We will illustrate this by creating a new lognormal background for two de-novo motifs in Drosophila.
                                #To load in the motifs the package provides functions to read standard JASPAR and TRANSFAC formats.
                                #library(PWMEnrich.Dmelanogaster.background)
                                
                                v.tfs.pwm <- unique(lst.pwm.data$df.motifs$TF.locus)
                                lst.tfs.motifsets <- vector(mode = "list", length(v.tfs.pwm))
                                names(lst.tfs.motifsets) <- v.tfs.pwm
                                
                                for(i in 1:length(v.tfs.pwm)){
                                  df.motifs.i <- subset(df.motifs, df.motifs$TF.locus %in% v.tfs.pwm[i])
                                  lst.tfs.motifsets[[i]] <- df.motifs.i$pos 
                                }
                                
                                th.p_value <- 5e-2
                                v.promoter_length <- 2000
                                
                                
                                library(biomaRt)
                                library(Biostrings)
                                #biocLite("BSgenome.Athaliana.TAIR.TAIR9")
                                library(BSgenome.Athaliana.TAIR.TAIR9)
                                genome <- BSgenome.Athaliana.TAIR.TAIR9
                                #source("http://bioconductor.org/biocLite.R")
                                #biocLite("TxDb.Athaliana.BioMart.plantsmart25")
                                library(TxDb.Athaliana.BioMart.plantsmart25)
                                
                                seqlevels(TxDb.Athaliana.BioMart.plantsmart25) <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9)
                                transcriptCoordsByGene.GRangesList <- transcriptsBy (TxDb.Athaliana.BioMart.plantsmart25, by = "gene")
                                isCircular(transcriptCoordsByGene.GRangesList) <- isCircular(BSgenome.Athaliana.TAIR.TAIR9)
                                #transcriptCoordsByGene.GRangesList <- renameSeqlevels( transcriptCoordsByGene.GRangesList, c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrM", "ChrC") )
                                transcriptCoordsByGene.GRangesList <- transcriptCoordsByGene.GRangesList[names(transcriptCoordsByGene.GRangesList) %in% v.tgs.selection]
                                
                                df.promSequences <- getPromoterSeq(transcriptCoordsByGene.GRangesList, genome, upstream=v.promoter_length, downstream=0)
                                
                                
                                #length(which(v.tgs_selection %in% names(df.promSequences)))
                                #df.promSequences <- df.promSequences[v.tgs_selection,]
                                
                                v.gns_proms <- (names(df.promSequences))
                                
                                
                                
                                n.cores <- 15
                                library(foreach)
                                #library(Biostrings)
                                
                                
                                cl<-makeCluster(n.cores)
                                registerDoParallel(cl)
                                
                                #lst.crfs <- foreach(l = 1:length(v.lambda), .packages=c("Matrix", "reshape2", "CRF", "caTools", "ROCR", "ggplot2", "igraph")) %dopar% { 
                                
                                #length(lst.tfs.motifsets)
                                
                                lst.v.dna_binding <- foreach(i = 1:length(lst.tfs.motifsets), .packages=c("Biostrings", "TFBSTools", "seqinr")) %dopar% { 
                                  
                                  mat.dna_binding <- matrix(1, nrow = length(v.tfs.pwm), ncol = length(v.gns_proms), dimnames = list(v.tfs.pwm, v.gns_proms))
                                  
                                  for(i in 1:length(lst.tfs.motifsets)){
                                    
                                    cat("Processing... ", round(i/length(lst.tfs.motifsets) * 100, digits = 2) , "%", "\r"); flush.console() 
                                    #v.dna_binding <- rep(1, length(v.gns_proms))
                                    #names(v.dna_binding) <- v.gns_proms
                                    
                                    pwmList <- lapply(lst.tfs.motifsets[[i]], function(j){
                                      profileMatrix <- lst.pwm.data$lst.pwm.motif[[j]];
                                      ID <- lst.pwm.data$df.motifs$TF.locus[j];
                                      name <- lst.pwm.data$df.motifs$TF.name[j];
                                      matrixClass <- lst.pwm.data$df.motifs$TF.family[j];
                                      PWMatrix(ID=ID, name=name, matrixClass=matrixClass,
                                               strand="+", bg=genomic.acgt, profileMatrix=profileMatrix)
                                    })
                                    tf <- ID(pwmList[[1]])    
                                    for(k in 1:length(pwmList)){
                                      for(j in 1:length(df.promSequences)){
                                        # cat("Processing... ", round(j/length(df.promSequences) * 100, digits = 2) , "%", "\r"); flush.console() 
                                        # print(paste(j, "of", length(df.promSequences)))
                                        #hits <- matchPWM(lst.pwm.motif[[i]], DNAString(toString(df.promSequences[[j]][[1]])) , with.score = TRUE)
                                        #hits <- matchPWM(pcm, DNAString(toString(df.promSequences[[j]][[1]])) , with.score = TRUE)
                                        #hits <- matchPWM(pcm, DNAString(toString(df.promSequences[[j]][[1]])) , with.score = TRUE)
                                        tg <- v.gns_proms[j]
                                        subject = DNAString(toString(df.promSequences[[j]][[1]])) 
                                        siteset = searchSeq(pwmList[[k]], subject, seqname=tg, min.score="60%", strand="*")
                                        if(length(siteset) > 0){
                                          ## calculate the empirical p-values of the scores
                                          p.value.min <- min(pvalues(siteset, type="TFMPvalue"))
                                          #if(p.value.min <= th.p_value){
                                          mat.dna_binding[tf, tg] <-  min(p.value.min, mat.dna_binding[tf, tg])  
                                          #v.dna_binding[tg] <- min(p.value.min, v.dna_binding[tg])
                                          #}
                                        }
                                      }
                                    }
                                    #set <- PWMatrixList(pwmList[[1]], pwmList[[2]])
                                    #sitesetList = searchSeq(set, subject, seqname="seq1", min.score="60%", strand="*")
                                    #v.dna_binding
                                  }
                                  
                                  
                                  saveRDS(mat.dna_binding, "mat.dna_binding.rds")
                                  
                                  stopCluster(cl)
                                  
                                  names(lst.v.dna_binding) <- names(lst.tfs.motifsets)
                                  saveRDS(lst.v.dna_binding, "lst.v.dna_binding.exp.rds")
                                  
                                  #pwm = toPWM(pfm, pseudocounts=0.8)
                                  
                                  #group motifs per tf 
                                  
                                  # compute p value... 
                                  
                                  
                                  
                                  perform_TFBS_mapping_all_genes <- function(v.tgs, v.promoter_length = 1500, v.th.bind_score = 0, n.cores = 15){
                                    
                                    print("compute scored all gene GRN")
                                    
                                    library(GenomicFeatures)
                                    #source("http://bioconductor.org/biocLite.R")
                                    # biocLite("GenomicFeatures")
                                    #   biocLite("ChIPpeakAnno")
                                    #   biocLite("biomaRt")
                                    #   biocLite("Biostrings")
                                    #   install.packages("VennDiagram")
                                    #library(ChIPpeakAnno)
                                    library(biomaRt)
                                    library(Biostrings)
                                    #biocLite("BSgenome.Athaliana.TAIR.TAIR9")
                                    library(BSgenome.Athaliana.TAIR.TAIR9)
                                    genome <- BSgenome.Athaliana.TAIR.TAIR9
                                    #source("http://bioconductor.org/biocLite.R")
                                    #biocLite("TxDb.Athaliana.BioMart.plantsmart25")
                                    library(TxDb.Athaliana.BioMart.plantsmart25)
                                    
                                    seqlevels(TxDb.Athaliana.BioMart.plantsmart25) <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9)
                                    transcriptCoordsByGene.GRangesList <- transcriptsBy (TxDb.Athaliana.BioMart.plantsmart25, by = "gene")
                                    isCircular(transcriptCoordsByGene.GRangesList) <- isCircular(BSgenome.Athaliana.TAIR.TAIR9)
                                    #transcriptCoordsByGene.GRangesList <- renameSeqlevels( transcriptCoordsByGene.GRangesList, c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrM", "ChrC") )
                                    transcriptCoordsByGene.GRangesList <- transcriptCoordsByGene.GRangesList[names(transcriptCoordsByGene.GRangesList) %in% v.tgs]
                                    
                                    df.promSequences <- getPromoterSeq(transcriptCoordsByGene.GRangesList, genome, upstream=v.promoter_length, downstream=0)
                                    
                                    #install.packages("foreach")
                                    
                                    lst.motifs <- get_cell_and_pnas_Paper_PWMs()
                                    lst.pwm.motif <- lst.motifs[[1]]
                                    df.motifs <- lst.motifs[[2]]
                                    
                                    
                                    # adapt to current approach
                                    library(foreach)
                                    library(Biostrings)
                                    
                                    
                                    cl<-makeCluster(n.cores)
                                    registerDoParallel(cl)
                                    
                                    #lst.crfs <- foreach(l = 1:length(v.lambda), .packages=c("Matrix", "reshape2", "CRF", "caTools", "ROCR", "ggplot2", "igraph")) %dopar% { 
                                    
                                    lst.dna_binding <- foreach(i = 1:length(lst.pwm.motif), .packages=c("Biostrings")) %dopar% { 
                                      
                                      #for(i in 1:length(lst.pwm.motif)){
                                      print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
                                      df.grn.all_motifs <- data.frame(TF = character(), Target = character(), TF.Motif = character(), 
                                                                      bind_score = numeric(), cmp_bind_score = numeric(), BS_pos_to_TSS = numeric(),
                                                                      stringsAsFactors = FALSE)    
                                      names(df.grn.all_motifs) <- c("TF", "Target", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_to_TSS")
                                      pcm <- round(100 * lst.pwm.motif[[i]])
                                      for(j in 1:length(df.promSequences)){
                                        cat("Processing... ", round(j/length(df.promSequences) * 100, digits = 2) , "%", "\r"); flush.console() 
                                        # print(paste(j, "of", length(df.promSequences)))
                                        #hits <- matchPWM(lst.pwm.motif[[i]], DNAString(toString(df.promSequences[[j]][[1]])) , with.score = TRUE)
                                        hits <- matchPWM(pcm, DNAString(toString(df.promSequences[[j]][[1]])) , with.score = TRUE)
                                        #matchPWM(reverseComplement(pcm), DNAString(toString(df.promSequences[[j]][[1]])) , with.score = TRUE)
                                        nhits <- length(hits)  
                                        if(nhits >= 1){
                                          # cmp_bind_score <- min(mcols(hits)$score / maxScore(lst.pwm.motif[[i]])) # should be >= 0.8
                                          cmp_bind_score <- min(mcols(hits)$score / maxScore(pcm)) # should be >= 0.8
                                          motif.score <- mcols(hits)$score / 100
                                          if(cmp_bind_score >= 0.8){
                                            for(k in 1:nhits){
                                              if(motif.score[k] >= v.th.bind_score){
                                                newrow <- data.frame(TF =  as.character(df.motifs$TF.locus[i]), 
                                                                     Target = names(df.promSequences[[j]])[1],
                                                                     TF.Motif =  names(lst.pwm.motif)[i],
                                                                     bind_score = motif.score[k], 
                                                                     cmp_bind_score = cmp_bind_score,
                                                                     BS_pos_to_TSS = abs(end(hits)[k] + start(hits)[k])/2,
                                                                     stringsAsFactors = FALSE)
                                                df.grn.all_motifs <- rbind(df.grn.all_motifs, newrow)   
                                              }
                                            }
                                          }
                                        }
                                      }
                                      names(df.grn.all_motifs) <- c("TF", "Target", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_to_TSS")
                                      #saveRDS(df.grn.all_motifs, paste("Datasets/novelTFmotifs/tmp/df.grn.all_motifs_",i,".rds", sep = ""))
                                      df.grn.all_motifs
                                    }
                                    
                                    
                                    stopCluster(cl)
                                    #saveRDS(lst.dna_binding, "Datasets/novelTFmotifs/lst.dna_binding.rds")
                                    lst.dna_binding <- readRDS("Datasets/novelTFmotifs/lst.dna_binding.rds")
                                    
                                    
                                    
                                    #df.grn.dnabinding_map <- do.call("rbind", lst.dna_binding)
                                    
                                    df.grn.dnabinding_map <- data.frame(TF = character(), Target = character(), TF.Motif = character(), 
                                                                        bind_score = numeric(), cmp_bind_score = numeric(), BS_pos_to_TSS = numeric(),
                                                                        stringsAsFactors = FALSE)      
                                    for(i in 1:length(lst.pwm.motif)){ 
                                      cat("Processing... ", round(i/length(lst.pwm.motif) * 100, digits = 2) , "%", "\r") 
                                      df.grn.dnabinding_map <- rbind(df.grn.dnabinding_map, lst.dna_binding[[i]])
                                    }
                                    names(df.grn.dnabinding_map) <- c("TF", "Target", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_to_TSS") 
                                    #saveRDS(df.grn.dnabinding_map, paste("Datasets/novelTFmotifs/df.grn.dnabinding_map.rds", sep = ""))
                                    
                                    # as with cns apply a mean based filter
                                    #df.grn.dnabinding_map <- subset(df.grn.dnabinding_map, df.grn.dnabinding_map$bind_score > mean(df.grn.dnabinding_map$bind_score))
                                    df.grn.dnabinding_map <- readRDS(paste("Datasets/novelTFmotifs/df.grn.dnabinding_map.rds", sep = ""))
                                    
                                    
                                    
                                    ##
                                    
                                    l.grn <- nrow(unique(df.grn.dnabinding_map[,1:2]))
                                    df.dnabind_w_counter_grn <- data.frame(TF = character(l.grn), Target = character(l.grn), v.max_bindscore = numeric(l.grn), bind_counter = numeric(l.grn), v.mean_bindscore = numeric(l.grn), v.sd_bindscore = numeric(l.grn), stringsAsFactors = FALSE) 
                                    tfs.dna_binding <- unique(df.grn.dnabinding_map$TF)
                                    
                                    idx <- 1
                                    for(i in 1:length(tfs.dna_binding)){
                                      cat("Processing... ", round(i/length(tfs.dna_binding) * 100, digits = 2) , "%", "\r") 
                                      flush.console()
                                      df.grn.dnabinding_map.i <- subset(df.grn.dnabinding_map, df.grn.dnabinding_map$TF == tfs.dna_binding[i])
                                      tgs.dna_binding <- unique(df.grn.dnabinding_map.i$Target)
                                      for(j in 1:length(tgs.dna_binding)){
                                        df.grn.dnabinding_map.ij <- subset(df.grn.dnabinding_map.i, df.grn.dnabinding_map.i$TF == tfs.dna_binding[i] & df.grn.dnabinding_map.i$Target == tgs.dna_binding[j])
                                        if(nrow(df.grn.dnabinding_map.ij) > 1){
                                          df.dnabind_w_counter_grn$TF[idx] = tfs.dna_binding[i]
                                          df.dnabind_w_counter_grn$Target[idx] = tgs.dna_binding[j]
                                          df.dnabind_w_counter_grn$v.max_bindscore[idx] = max(df.grn.dnabinding_map.ij$bind_score)
                                          df.dnabind_w_counter_grn$bind_counter[idx] = nrow(df.grn.dnabinding_map.ij)
                                          df.dnabind_w_counter_grn$v.mean_bindscore[idx] = mean(df.grn.dnabinding_map.ij$bind_score)
                                          df.dnabind_w_counter_grn$v.sd_bindscore[idx] = sd(df.grn.dnabinding_map.ij$bind_score)
                                          idx <- idx + 1
                                        }
                                      } 
                                    }
                                    
                                    saveRDS(df.dnabind_w_counter_grn, "df.dnabind_w_counter_grn.rds")
                                    ##
                                    
                                  }
                                  
                                  
                                  
                                  #biocLite("JASPAR2014")
                                  #biocLite("TFBSTools")
                                  #biocLite("RSQLite")
                                  #library(RSQLite)
                                  #library(TFBSTools)
                                  #library(JASPAR2015)
                                  #   opts = list()
                                  #opts[["tax_group"]] = "plants"
                                  #PFMatrixList = getMatrixSet(JASPAR2014, opts)
                                  
                                  #   
                                  #   
                                  #   vec.jaspar.motifs <- c("MA0008.1", "MA0121.1", "MA0548.1", "MA0549.1", "MA0550.1", "MA0551.1",
                                  #                          "MA0552.1", "MA0553.1", "MA0554.1", "MA0555.1", "MA0556.1", "MA0557.1",
                                  #                          "MA0558.1","MA0559.1","MA0560.1", "MA0561.1", "MA0562.1", "MA0563.1",
                                  #                          "MA0564.1","MA0565.1","MA0566.1","MA0567.1","MA0568.1","MA0569.1","MA0570.1",
                                  #                          "MA0005.2","MA0571.1","MA0572.1","MA0110.2","MA0573.1","MA0574.1",
                                  #                          "MA0575.1","MA0576.1","MA0577.1","MA0578.1","MA0579.1","MA0580.1","MA0581.1",
                                  #                          "MA0582.1","MA0583.1","MA0584.1","MA0001.2","MA0585.1","MA0586.1","MA0587.1",
                                  #                          "MA0588.1","MA0589.1","MA0590.1")
                                  # 
                                  # 
                                  #   vec.TF.names <- sapply(PFMatrixList.Ath, function(x) x@name, simplify = TRUE)
                                  #   vec.TF.names[40] <- "RAV1"
                                  # 
                                  #   df.locus.set <- subset(df.gene_names, df.gene_names$primary_gene_symbol %in% vec.TF.names)
                                  
                          
                                  
                                  length(unique(df.pwms.pnas$TF.locus))
                                  length(unique(df.motifs$TF.locus))
                                  length(unique(df.motifs.jaspar_2016$tf.id))
                                  (unique(c(df.motifs.jaspar_2016$tf.name, df.pwms.pnas$TF.name, df.motifs$TF.name)))
                                  
                                  
                                  
                                  # to paper evaluation 
                                  
                                  library(seqinr)
                                  
                                  
                                  
                                  
                                  
                                  
                                  #####
                                  
                                  pcm <- matrix(c(1, 965, 0, 0, 0, 827, 930, 50, 421, 0, 999, 0, 276, 173, 0, 599, 1, 35, 0, 999, 69, 0, 0, 50, 578, 0, 0, 0, 655, 0, 69, 300), nrow = 4, ncol = 8, dimnames = list(c("A","C","G","T"), c(1,2,3,4,5,6,7,8)))
                                  #pwm <- PWM(pcm)
                                  
                                  
                                  library(foreach)
                                  library(doParallel)
                                  
                                  
                                  
                                  library(seqLogo)
                                  seqLogo(pwm)
                                  
                                  #   lst.dna_binding <- vector( mode = "list", length = nrow(df.motifs))
                                  #   names(lst.dna_binding) <- df.motifs$TF.locus
                                  
                                  
                                  # STEP I - compute prior network (complete network) 
                                  # main merit (based on lam_meeting_merit)
                                  #upstream_sequences <- read.fasta(file = "../Datasets/TAIR10_upstream_1000_20101104.txt", seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
                                  upstream_sequences <- read.fasta(file = "../Datasets/TAIR10_upstream_3000_20101028.txt", seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
                                  
                                  #upstream_sequences <- upstream_sequences[names(upstream_sequences) %in% unique(df.GRACE.ensemble.ams$TG)]
                                  upstream_sequences <- lapply(upstream_sequences, function(m) {toupper(m)})
                                  df.promSequences <- DNAStringSet(unlist(upstream_sequences))
                                  
                                  v.gns.ath <- intersect(names(df.promSequences), tgs)
                                  df.promSequences <- df.promSequences[v.gns.ath]
                                  #v.gns.ath <- names(df.promSequences)
                                  
                                  
                                  #   which(df.motifs$TF.locus == "AT1G06160")
                                  #   
                                  # STEP II - compute binding 
                                  
                                  n.cores <- 15
                                  
                                  strt<-Sys.time()
                                  cl<-makeCluster(n.cores)
                                  registerDoParallel(cl)
                                  lst.dna_binding <- foreach(i = 1:nrow(df.motifs), .packages=c("Biostrings")) %dopar% { 
                                    
                                    cat("Processing... ", round(i/nrow(df.motifs) * 100, digits = 2) , "%", "\r"); flush.console() 
                                    pwm <- lst.pwm.motif[[i]]
                                    
                                    lst.dna_binding <- vector(mode = "list", length = length(df.promSequences))
                                    names(lst.dna_binding) <- names(df.promSequences)
                                    
                                    #strt<-Sys.time()
                                    for(j in 1:length(df.promSequences)){
                                      cat("Processing... ", round(j/length(df.promSequences) * 100, digits = 2) , "%", "\r"); flush.console() 
                                      
                                      # + strand
                                      hits <- matchPWM(pwm, df.promSequences[[j]] , with.score = TRUE)
                                      
                                      nhits <- length(hits)  
                                      if(nhits >= 1){
                                        cmp_bind_score <- min(mcols(hits)$score / maxScore(pwm)) # should be >= 0.8
                                        motif.score <- mcols(hits)$score
                                        if(cmp_bind_score >= 0.8){
                                          scores <- numeric(nhits)
                                          BS_pos_to_TSS <- numeric(nhits) 
                                          for(k in 1:nhits){
                                            scores[k] = motif.score[k]
                                            BS_pos_to_TSS[k] = abs(end(hits)[k] + start(hits)[k])/2
                                          }
                                          lst.dna_binding[[j]] <- vector(mode = "list", length = 2)
                                          names(lst.dna_binding[[j]]) <- c("score", "position")
                                          lst.dna_binding[[j]][[1]] <- scores
                                          lst.dna_binding[[j]][[2]] <- BS_pos_to_TSS
                                        }
                                      }
                                      
                                      # - strand
                                      #       hits <- matchPWM(reverseComplement(pwm), df.promSequences[[j]] , with.score = TRUE)
                                      #       
                                      #       nhits <- length(hits)  
                                      #       if(nhits >= 1){
                                      #         cmp_bind_score <- min(mcols(hits)$score / maxScore(pwm)) # should be >= 0.8
                                      #         motif.score <- mcols(hits)$score
                                      #         if(cmp_bind_score >= 0.8){
                                      #           scores <- numeric(nhits)
                                      #           BS_pos_to_TSS <- numeric(nhits) 
                                      #           for(k in 1:nhits){
                                      #             scores[k] = motif.score[k]
                                      #             BS_pos_to_TSS[k] = abs(end(hits)[k] + start(hits)[k])/2
                                      #           }
                                      #           lst.dna_binding[[j]][[2]] <- vector(mode = "list", length = 2)
                                      #           names(lst.dna_binding[[j]][[2]]) <- c("score", "position")
                                      #           lst.dna_binding[[j]][[2]][[1]] <- scores
                                      #           lst.dna_binding[[j]][[2]][[2]] <- BS_pos_to_TSS
                                      #         }
                                      #       }
                                    }
                                    #print(Sys.time()-strt)
                                    
                                    lst.dna_binding <- lst.dna_binding[!sapply(lst.dna_binding, is.null)] 
                                    lst.dna_binding
                                  }
                                  
                                  stopCluster(cl)
                                  print(Sys.time()-strt)
                                  print("..finished.")
                                  
                                  
                                  names(lst.dna_binding) <- df.motifs$TF.locus
                                  #saveRDS(lst.dna_binding, "lst.dna_binding.rds")
                                  
                                  
                                  #lst.dna_binding <- lapply(lst.dna_binding, function(m) m[!sapply(m, is.null)])
                                  
                                  
                                  # STEP 2 - BP associattion network? 
                                  
                                  
                                  ## conserved only 
                                  
                                  perform_TFBS_mapping_CNS2014 <- function(v.tgs, v.th.bind_score = 0, n.cores = 15){
                                    
                                    print("compute scored CNS2014-GRN")
                                    
                                    lst.motifs <- get_cell_and_pnas_Paper_PWMs()
                                    lst.pwm.motif <- lst.motifs[[1]]
                                    df.motifs <- lst.motifs[[2]]
                                    
                                    if(b.CALC){
                                      df.CNS_2014 <- readRDS("../Datasets/df.CNS_2014.rds") 
                                    }else{
                                      df.CNS_2014 <- readRDS("Datasets/df.CNS_2014.rds") 
                                    }
                                    df.CNS_2014 <- subset(df.CNS_2014, df.CNS_2014$Target %in% v.tgs)
                                    # df.CNS_2014 <- readRDS(paste("Datasets/novelTFmotifs/df.CNS_2014_2000kb.rds", sep = ""))
                                    
                                    #df.CNS_2014 <- subset(df.CNS_2014, df.CNS_2014$end_dist_to_TSS < 1000)  
                                    vec.cns.sequences <- as.character(unique(df.CNS_2014$cns_sequence))
                                    n.cns.seq <- length(vec.cns.sequences)
                                    
                                    #system.time(for(i in 1:1){#length(lst.pwm.motif)){ 
                                    #foreach(i = 1:2) %dopar% { 
                                    cl<-makeCluster(n.cores)
                                    registerDoParallel(cl)
                                    lst.dna_binding <- foreach(i = 1:length(lst.pwm.motif), .packages=c("Biostrings")) %dopar% { 
                                      #print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
                                      df.grn.CNS2014.TFB_map <- data.frame(TF = character(), cns_sequence = character(), TF.Motif = character(), 
                                                                           bind_score = numeric(), cmp_bind_score = numeric(), BS_pos_to_TSS = numeric(), 
                                                                           stringsAsFactors = FALSE)    
                                      names(df.grn.CNS2014.TFB_map) <- c("TF", "cns_sequence", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_in_CNS")
                                      pcm <- round(100 * lst.pwm.motif[[i]])
                                      for(j in 1:n.cns.seq){
                                        hits <- matchPWM(pcm,  vec.cns.sequences[j], with.score = TRUE)
                                        nhits <- length(hits)  
                                        if(nhits >= 1){
                                          cmp_bind_score <- min(mcols(hits)$score / maxScore(pcm)) # should be >= 0.8
                                          motif.score <- mcols(hits)$score / 100
                                          if(cmp_bind_score >= 0.8){
                                            for(k in 1:nhits){
                                              if(motif.score[k] >= v.th.bind_score){
                                                newrow <- data.frame(TF =  as.character(df.motifs$TF.locus[i]), 
                                                                     cns_sequence = vec.cns.sequences[j],
                                                                     TF.Motif =  names(lst.pwm.motif)[i],
                                                                     bind_score = motif.score[k], 
                                                                     cmp_bind_score = cmp_bind_score,
                                                                     BS_pos_in_CNS = abs(end(hits)[k] + start(hits)[k])/2,
                                                                     stringsAsFactors = FALSE)
                                                df.grn.CNS2014.TFB_map <- rbind(df.grn.CNS2014.TFB_map, newrow)   
                                              }
                                            }
                                          }
                                        }
                                      }
                                      names(df.grn.CNS2014.TFB_map) <- c("TF", "cns_sequence", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_in_CNS")
                                      #saveRDS(df.grn.CNS2014.TFB_map, paste("Datasets/CNS_GRNS/tmp/df.cns2014_grn_map_",i,".rds", sep = ""))
                                      df.grn.CNS2014.TFB_map
                                    }
                                    stopCluster(cl)
                                    saveRDS(lst.dna_binding, "Datasets/novelTFmotifs/lst.dna_binding.rds")
                                    
                                    
                                    df.grn.CNS2014.TFB_map <- do.call("rbind", lst.dna_binding)
                                    
                                    # combine
                                    #   df.grn.CNS2014.TFB_map <- data.frame(TF = character(), cns_sequence = character(), TF.Motif = character(), 
                                    #                                        bind_score = numeric(), cmp_bind_score = numeric(),
                                    #                                        stringsAsFactors = FALSE)      
                                    #   for(i in 1:length(lst.pwm.motif)){ 
                                    #     print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
                                    #     df.grn.CNS2014.TFB_map <- rbind(df.grn.CNS2014.TFB_map, readRDS(paste("Datasets/CNS_GRNS/tmp/df.cns2014_grn_map_",i,".rds", sep = "")))
                                    #   }
                                    
                                    #names(df.grn.CNS2014.TFB_map) <- c("TF", "cns_sequence", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_in_CNS") 
                                    #saveRDS(df.grn.CNS2014.TFB_map, paste("Datasets/CNS_GRNS/df.grn.CNS2014.TFB_map.rds", sep = ""))
                                    #saveRDS(df.grn.CNS2014.TFB_map, paste("Datasets/CNS_GRNS/df.grn.CNS2014.TFB_map_complete.rds", sep = ""))
                                    
                                    
                                    
                                    # finalize
                                    #df.grn.CNS2014.TFB_map <- readRDS(paste("Datasets/CNS_GRNS/df.grn.CNS2014.TFB_map_complete.rds", sep = ""))
                                    df.grn.CNS2014.TFB_map <- df.grn.CNS2014.TFB_map[,c(-3,-5)]
                                    df.grn.CNS2014.TFB_map <- unique(df.grn.CNS2014.TFB_map)
                                    
                                    #df.CNS_2014["BS_pos_to_TSS"] <-  apply(df.CNS_2014[,c('start_dist_to_TSS','end_dist_to_TSS')], 1, function(x) mean(x) )
                                    #df.CNS_2014 <- df.CNS_2014[, c(-1,-2,-3,-8,-9,-10)]
                                    
                                    df.CNS_grn <- merge(df.grn.CNS2014.TFB_map, df.CNS_2014, by = "cns_sequence")
                                    
                                    
                                    #   df.CNS_grn["BS_pos_to_TSS"] <- NA
                                    #   for(i in 1:nrow(df.CNS_grn)){
                                    #     print(paste(i, "of", nrow(df.CNS_grn)))
                                    #    if(df.CNS_grn$start_dist_to_TSS[i] == df.CNS_grn$shortest_dist_to_TSS[i]){
                                    #      df.CNS_grn$BS_pos_to_TSS[i] <- df.CNS_grn$shortest_dist_to_TSS[i] + df.CNS_grn$BS_pos_in_CNS[i]
                                    #    }else{
                                    #      df.CNS_grn$BS_pos_to_TSS[i] <- df.CNS_grn$start_dist_to_TSS[i] - df.CNS_grn$BS_pos_in_CNS[i] 
                                    #    }
                                    #   }
                                    #   
                                    #saveRDS(df.CNS_grn, paste("Datasets/df.CNS_grn2014_complete.rds", sep = ""))
                                    df.CNS_grn <- unique(df.CNS_grn[,c(2,3,8,9,10)])
                                    saveRDS(df.CNS_grn, paste("Datasets/df.CNS_novel_pwm_motifs_grn.rds", sep = ""))
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    #
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    #saveRDS(lst.dna_binding.nac025, "lst.dna_binding.nac025.rds")
                                    
                                    
                                    dist.i <- unlist(lapply(lst.dna_binding[[i]], function(m) { m[[1]] }))
                                    
                                    
                                    # Cis Regulatory Modules (at a reasonable distance)
                                    
                                    tgs.i <- unique(names(which(dist.i >= quantile(dist.i, 0.99))))
                                    
                                    
                                    
                                    
                                    [[]]
                                    quantile(dist.i, 0.5)
                                    quantile(dist.i, 0.8)
                                    quantile(dist.i, 0.85)
                                    quantile(dist.i, 0.9)
                                    which(dist.i >= quantile(dist.i, 0.95))
                                    quantile(dist.i, 0.99)
                                    
                                    
                                    
                                    
                                    lst.dna_binding.nac025["AT2G16910"]
                                    lst.dna_binding.nac025["AT5G40260"] # Swwet 
                                    lst.dna_binding.nac025["AT4G20050"]
                                    
                                    
                                    
                                    
                                    ####
                                    
                                    
                                    installed.packages("rtfbs")
                                    library(rtfbs)
                                    require("rtfbs")
                                    
                                    biocLite("MotIV")
                                    library(MotIV)
                                    
                                    sort (table (values (MotifDb)$dataSource), decreasing=TRUE)
                                    
                                    sort (table (values (MotifDb)$organism), decreasing=TRUE)
                                    
                                    jaspar <- readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
                                    
                                    jaspar <- readPWMfile("../Datasets/novelTFmotifs/Jaspar2014/pfm_plants_jaspar2016.txt")
                                    
                                    df.pwms.jaspar2014 <- read.table("../Datasets/novelTFmotifs/Jaspar2014/pfm_plants_JASPAR_Ath_only.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
                                    #PFMatrixList.Ath <- PFMatrixList[vec.jaspar.motifs]
                                    
                                    # saveRDS(PFMatrixList.Ath, "Datasets/novelTFmotifs/Jaspar2014/PFMatrixList.Ath.rds")
                                    if(b.CALC){
                                      PFMatrixList.Ath <- readRDS("../Datasets/novelTFmotifs/Jaspar2014/PFMatrixList.Ath.rds")  
                                      df.jaspar.motif.info <- read.table("../Datasets/novelTFmotifs/Jaspar2014/jaspar2014AthMotifTable.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
                                      df.jaspar.motif.info <- read.table("../Datasets/novelTFmotifs/Jaspar2014/pfm_plants_jaspar2016.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
                                    }else{
                                      PFMatrixList.Ath <- readRDS("Datasets/novelTFmotifs/Jaspar2014/PFMatrixList.Ath.rds")  
                                      df.jaspar.motif.info <- read.table("Datasets/novelTFmotifs/Jaspar2014/jaspar2014AthMotifTable.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)   
                                    }
                                    
                                    
                                    for(i in 1:length(PFMatrixList.Ath)){
                                      motif.id <- PFMatrixList.Ath[[i]]@ID
                                      pfm <- PFMatrixList.Ath[[i]]@profileMatrix
                                      jaspar.motif <- apply(pfm, 2, function(x) x / sum(x))
                                      colnames(jaspar.motif) <- as.character(seq(1:ncol(jaspar.motif)))
                                      lst.pwm.motif <- lappend(lst.pwm.motif, jaspar.motif)
                                      names(lst.pwm.motif)[[(n.cell + n.pnas + i)]] <- motif.id 
                                      df.jaspar.motif.info.sset <- subset(df.jaspar.motif.info, df.jaspar.motif.info$Motif.ID == motif.id)
                                      # check for identical list elements
                                      newrow <- data.frame(Motif.id = motif.id, TF.locus = df.jaspar.motif.info.sset$TF.locus, TF.name = df.jaspar.motif.info.sset$TF.name,
                                                           TF.family = df.jaspar.motif.info.sset$TF.family, TF.subfamily = "", src = "JASPAR2014")
                                      names(newrow) = c("Motif.ID", "TF.locus", "TF.name", "TF.family", "TF.subfamily", "src")
                                      df.motifs <- rbind(df.motifs, newrow)
                                    }    
                                    
                                    idx <- which(is.na(df.motifs$TF.locus))
                                    df.motifs <- na.omit(df.motifs)
                                    lst.pwm.motif[idx] <- NULL
                                    
                                    return(list(lst.pwm.motif=lst.pwm.motif, df.motifs=df.motifs))
                                  }
                                  
                                  
                                  
                                  
                                  ##### Conservation based ######
                                  
                                  perform_TFBS_mapping_CNS2014 <- function(v.tgs, v.th.bind_score = 0, n.cores = 15){
                                    
                                    print("compute scored CNS2014-GRN")
                                    
                                    lst.motifs <- get_cell_and_pnas_Paper_PWMs()
                                    lst.pwm.motif <- lst.motifs[[1]]
                                    df.motifs <- lst.motifs[[2]]
                                    
                                    if(b.CALC){
                                      df.CNS_2014 <- readRDS("../Datasets/df.CNS_2014.rds") 
                                    }else{
                                      df.CNS_2014 <- readRDS("Datasets/df.CNS_2014.rds") 
                                    }
                                    df.CNS_2014 <- subset(df.CNS_2014, df.CNS_2014$Target %in% v.tgs)
                                    # df.CNS_2014 <- readRDS(paste("Datasets/novelTFmotifs/df.CNS_2014_2000kb.rds", sep = ""))
                                    
                                    #df.CNS_2014 <- subset(df.CNS_2014, df.CNS_2014$end_dist_to_TSS < 1000)  
                                    vec.cns.sequences <- as.character(unique(df.CNS_2014$cns_sequence))
                                    n.cns.seq <- length(vec.cns.sequences)
                                    
                                    #system.time(for(i in 1:1){#length(lst.pwm.motif)){ 
                                    #foreach(i = 1:2) %dopar% { 
                                    cl<-makeCluster(n.cores)
                                    registerDoParallel(cl)
                                    lst.dna_binding <- foreach(i = 1:length(lst.pwm.motif), .packages=c("Biostrings")) %dopar% { 
                                      #print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
                                      df.grn.CNS2014.TFB_map <- data.frame(TF = character(), cns_sequence = character(), TF.Motif = character(), 
                                                                           bind_score = numeric(), cmp_bind_score = numeric(), BS_pos_to_TSS = numeric(), 
                                                                           stringsAsFactors = FALSE)    
                                      names(df.grn.CNS2014.TFB_map) <- c("TF", "cns_sequence", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_in_CNS")
                                      pcm <- round(100 * lst.pwm.motif[[i]])
                                      for(j in 1:n.cns.seq){
                                        hits <- matchPWM(pcm,  vec.cns.sequences[j], with.score = TRUE)
                                        nhits <- length(hits)  
                                        if(nhits >= 1){
                                          cmp_bind_score <- min(mcols(hits)$score / maxScore(pcm)) # should be >= 0.8
                                          motif.score <- mcols(hits)$score / 100
                                          if(cmp_bind_score >= 0.8){
                                            for(k in 1:nhits){
                                              if(motif.score[k] >= v.th.bind_score){
                                                newrow <- data.frame(TF =  as.character(df.motifs$TF.locus[i]), 
                                                                     cns_sequence = vec.cns.sequences[j],
                                                                     TF.Motif =  names(lst.pwm.motif)[i],
                                                                     bind_score = motif.score[k], 
                                                                     cmp_bind_score = cmp_bind_score,
                                                                     BS_pos_in_CNS = abs(end(hits)[k] + start(hits)[k])/2,
                                                                     stringsAsFactors = FALSE)
                                                df.grn.CNS2014.TFB_map <- rbind(df.grn.CNS2014.TFB_map, newrow)   
                                              }
                                            }
                                          }
                                        }
                                      }
                                      names(df.grn.CNS2014.TFB_map) <- c("TF", "cns_sequence", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_in_CNS")
                                      #saveRDS(df.grn.CNS2014.TFB_map, paste("Datasets/CNS_GRNS/tmp/df.cns2014_grn_map_",i,".rds", sep = ""))
                                      df.grn.CNS2014.TFB_map
                                    }
                                    stopCluster(cl)
                                    saveRDS(lst.dna_binding, "Datasets/novelTFmotifs/lst.dna_binding.rds")
                                    
                                    
                                    df.grn.CNS2014.TFB_map <- do.call("rbind", lst.dna_binding)
                                    
                                    # combine
                                    #   df.grn.CNS2014.TFB_map <- data.frame(TF = character(), cns_sequence = character(), TF.Motif = character(), 
                                    #                                        bind_score = numeric(), cmp_bind_score = numeric(),
                                    #                                        stringsAsFactors = FALSE)      
                                    #   for(i in 1:length(lst.pwm.motif)){ 
                                    #     print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
                                    #     df.grn.CNS2014.TFB_map <- rbind(df.grn.CNS2014.TFB_map, readRDS(paste("Datasets/CNS_GRNS/tmp/df.cns2014_grn_map_",i,".rds", sep = "")))
                                    #   }
                                    
                                    #names(df.grn.CNS2014.TFB_map) <- c("TF", "cns_sequence", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_in_CNS") 
                                    #saveRDS(df.grn.CNS2014.TFB_map, paste("Datasets/CNS_GRNS/df.grn.CNS2014.TFB_map.rds", sep = ""))
                                    #saveRDS(df.grn.CNS2014.TFB_map, paste("Datasets/CNS_GRNS/df.grn.CNS2014.TFB_map_complete.rds", sep = ""))
                                    
                                    
                                    
                                    # finalize
                                    #df.grn.CNS2014.TFB_map <- readRDS(paste("Datasets/CNS_GRNS/df.grn.CNS2014.TFB_map_complete.rds", sep = ""))
                                    df.grn.CNS2014.TFB_map <- df.grn.CNS2014.TFB_map[,c(-3,-5)]
                                    df.grn.CNS2014.TFB_map <- unique(df.grn.CNS2014.TFB_map)
                                    
                                    #df.CNS_2014["BS_pos_to_TSS"] <-  apply(df.CNS_2014[,c('start_dist_to_TSS','end_dist_to_TSS')], 1, function(x) mean(x) )
                                    #df.CNS_2014 <- df.CNS_2014[, c(-1,-2,-3,-8,-9,-10)]
                                    
                                    df.CNS_grn <- merge(df.grn.CNS2014.TFB_map, df.CNS_2014, by = "cns_sequence")
                                    
                                    
                                    #   df.CNS_grn["BS_pos_to_TSS"] <- NA
                                    #   for(i in 1:nrow(df.CNS_grn)){
                                    #     print(paste(i, "of", nrow(df.CNS_grn)))
                                    #    if(df.CNS_grn$start_dist_to_TSS[i] == df.CNS_grn$shortest_dist_to_TSS[i]){
                                    #      df.CNS_grn$BS_pos_to_TSS[i] <- df.CNS_grn$shortest_dist_to_TSS[i] + df.CNS_grn$BS_pos_in_CNS[i]
                                    #    }else{
                                    #      df.CNS_grn$BS_pos_to_TSS[i] <- df.CNS_grn$start_dist_to_TSS[i] - df.CNS_grn$BS_pos_in_CNS[i] 
                                    #    }
                                    #   }
                                    #   
                                    #saveRDS(df.CNS_grn, paste("Datasets/df.CNS_grn2014_complete.rds", sep = ""))
                                    df.CNS_grn <- unique(df.CNS_grn[,c(2,3,8,9,10)])
                                    saveRDS(df.CNS_grn, paste("Datasets/df.CNS_novel_pwm_motifs_grn.rds", sep = ""))
                                    
                                    
                                    #   vec.sequences <- unique(df.CNS_grn$cns_sequence)
                                    #   
                                    #   
                                    #   ##????
                                    #   
                                    #   foreach(i = 1:length(vec.sequences)) %dopar% { 
                                    #   #for(i in 1:length(vec.sequences)){
                                    #     df.CNS_grn.filtered <- data.frame(TF = character(), bind_score = numeric(),
                                    #                                       val = numeric(), cons_score = numeric(), Target = character(),
                                    #                                       BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE)  
                                    #     print(paste(i, "of",length(vec.sequences)))
                                    #     df.sset <- subset(df.CNS_grn, df.CNS_grn$cns_sequence == vec.sequences[i])
                                    #     if(nrow(df.sset) > 0){
                                    #       tfs <- unique(df.sset$TF)
                                    #       for(t.1 in 1:length(tfs)){
                                    #         df.sset.2 <- subset(df.sset, df.sset$TF == tfs[t.1])
                                    #         targs <- unique(df.sset.2$Target)  
                                    #         for(t.2 in 1:length(targs)){
                                    #           df.sset.3 <- subset(df.sset.2, df.sset.2$Target == targs[t.2])
                                    #           vec.pos <- unique(df.sset.3$shortest_dist_to_TSS)
                                    #           for(p in 1:length(vec.pos)){
                                    #             df.sset.4 <- subset(df.sset.3, df.sset.3$shortest_dist_to_TSS == vec.pos[p])
                                    #             df.sset.4 <- subset(df.sset.4, df.sset.4$bind_score == max(unique(df.sset.4$bind_score))) # could be more then one... filter
                                    #             df.CNS_grn.filtered <- rbind(df.CNS_grn.filtered, df.sset.4[1,])
                                    #           }
                                    #         }
                                    #       }
                                    #     }
                                    #     saveRDS(df.CNS_grn.filtered, paste("Datasets/CNS_GRNS/tmp/df.cns2014_grn_filtered_",i,".rds", sep = ""))
                                    #   }
                                    #       
                                    #   
                                    #   df.CNS_grn.filtered <- data.frame(cns_sequence = character(), TF = character(), bind_score = numeric(), BS_pos_in_CNS = numeric(),
                                    #                                      Chr = numeric(), start = numeric(), end = numeric(), val = numeric(), cons_score = numeric(), Target = character(),
                                    #                                     start_dist_to_TSS = numeric(), shortest_dist_to_TSS = numeric(), feature_position = character(), BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE) 
                                    #   for(i in 1:length(vec.sequences)){ 
                                    #     print(paste(i, "of", length(vec.sequences)))
                                    #     df.CNS_grn.filtered <- rbind(df.CNS_grn.filtered, readRDS(paste("Datasets/CNS_GRNS/tmp/df.cns2014_grn_filtered_",i,".rds", sep = "")))
                                    #   }
                                    #   names(df.CNS_grn.filtered) <- c("cns_sequence", "TF", "bind_score", "BS_pos_in_CNS", "Chr", "start", "end", "val", 
                                    #                                   "cons_score", "Target", "start_dist_to_TSS", "shortest_dist_to_TSS", "feature_position", "BS_pos_to_TSS")
                                    #   saveRDS(df.CNS_grn.filtered, paste("Datasets/CNS_GRNS/df.CNS_grn2014_filtered.rds", sep = ""))
                                    
                                  }
                                  
                                  
                                  exact_CNS2014map_seq_motifs <- function(v.tgs, v.th.seqMotif_length = 4, n.cores = 4){
                                    
                                    library(GenomicFeatures)
                                    library(biomaRt)
                                    library(Biostrings)
                                    
                                    if(b.CALC){
                                      df.seqMotifs <- read.table("../Datasets/novelTFmotifs/athamap/athmap_seqMotifs.txt", header = TRUE, stringsAsFactor = FALSE, sep = "\t")
                                    }else{
                                      df.seqMotifs <- read.table("Datasets/novelTFmotifs/athamap/athmap_seqMotifs.txt", header = TRUE, stringsAsFactor = FALSE, sep = "\t")  
                                    }
                                    
                                    df.seqMotifs <- subset(df.seqMotifs, nchar(df.seqMotifs$verified_binding_seq) >= v.th.seqMotif_length)
                                    df.seqMotifs$TF.locus <- toupper(df.seqMotifs$TF.locus)
                                    df.seqMotifs$TF.locus <- unlist(sapply(df.seqMotifs$TF.locus, function(x) gsub(" ","",x, fixed=TRUE)))
                                    #df.seqMotifs <- subset(df.seqMotifs, !df.seqMotifs$TF.locus %in% vec.regs.novel.matrix)
                                    dict_seq_motifs <- DNAStringSet(df.seqMotifs$verified_binding_seq)
                                    
                                    df.CNS_2014 <- readRDS("Datasets/df.CNS_2014.rds") 
                                    df.CNS_2014 <- subset(df.CNS_2014, df.CNS_2014$Target %in% v.tgs)
                                    # df.CNS_2014 <- readRDS(paste("Datasets/novelTFmotifs/df.CNS_2014_2000kb.rds", sep = ""))
                                    
                                    #df.CNS_2014 <- subset(df.CNS_2014, df.CNS_2014$end_dist_to_TSS < 1000)  
                                    vec.cns.sequences <- as.character(unique(df.CNS_2014$cns_sequence))
                                    n.cns.seq <- length(vec.cns.sequences)
                                    
                                    #pdict_seq_motifs <- PDict(dict_seq_motifs) 
                                    #foreach(i = 1:nrow(df.CNS_2014)) %dopar% { 
                                    for(i in 1:nrow(df.CNS_2014)){
                                      
                                      cat("Processing... ", round(i/nrow(df.CNS_2014) * 100, digits = 0) , "%", "\r") 
                                      flush.console()
                                      Target <- df.CNS_2014$Target[i]
                                      # find these motifs in promoters - get genes
                                      #sum(countPDict(dict_seq_motifs, DNAString(toString(df.CNS_2014$cns_sequence[i])),max.mismatch=0, min.mismatch=0))
                                      motifMatches <- matchPDict(dict_seq_motifs, DNAString(toString(df.CNS_2014$cns_sequence[i])),max.mismatch=0, min.mismatch=0)
                                      idxes <-  which(elementLengths(motifMatches) == 1)
                                      
                                      if(length(idxes) > 0){
                                        df.CNS_grn.seqMotifs <- data.frame(TF = character(), cns_sequence = character(), Target = character(), cons_score = numeric(),  val = numeric()) #,  BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE)
                                        
                                        #idxes <- countIndex(motifMatches)  
                                        # how many per promoter and position - 
                                        for(j in 1:length(idxes)){
                                          
                                          newrow <- data.frame(TF = df.seqMotifs$TF.locus[j], cns_sequence = df.CNS_2014$cns_sequence[i], Target = Target, 
                                                               cons_score = df.CNS_2014$cons_score[i], val = df.CNS_2014$val[i])#, BS_pos_to_TSS = BS_pos_to_TSS)
                                          df.CNS_grn.seqMotifs <- rbind(df.CNS_grn.seqMotifs, newrow)
                                          
                                          #       if(idxes[j] > 0){
                                          #         for(k in 1:idxes[j]){
                                          #           #           mean_dist_to_TSS <- round(mean(c(endIndex(motifMatches)[[j]][k], startIndex(motifMatches)[[j]][k])))
                                          #           #           if(df.CNS_2014$start_dist_to_TSS[i] == df.CNS_2014$shortest_dist_to_TSS[i]){
                                          #           #             BS_pos_to_TSS <- df.CNS_2014$shortest_dist_to_TSS[i] + mean_dist_to_TSS
                                          #           #           }else{
                                          #           #             BS_pos_to_TSS <- df.CNS_2014$start_dist_to_TSS[i] - mean_dist_to_TSS
                                          #           #           }
                                          #           newrow <- data.frame(TF = df.seqMotifs$TF.locus[j], cns_sequence = df.CNS_2014$cns_sequence[i], Target = Target, 
                                          #                                cons_score = df.CNS_2014$cons_score[i], val = df.CNS_2014$val[i])#, BS_pos_to_TSS = BS_pos_to_TSS)
                                          #           df.CNS_grn.seqMotifs <- rbind(df.CNS_grn.seqMotifs, newrow)
                                          #         }
                                          #       }
                                        }
                                        #df.CNS_grn.seqMotifs <- unique(df.CNS_grn.seqMotifs)
                                        #names(df.CNS_grn.seqMotifs) <- c("TF", "cns_sequence", "Target", "cons_score", "val", "BS_pos_to_TSS") 
                                        saveRDS(df.CNS_grn.seqMotifs, paste("Datasets/tmp/df.CNS_grn.seqMotifs_",i,".rds", sep = ""))
                                      }
                                    }
                                    
                                    df.CNS_grn.seqMotifs <- data.frame(TF = character(), Target = character(), bind_score = numeric(), BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE)       
                                    for(i in 1:nrow(df.CNS_2014)){ 
                                      if(file.exists(paste("Datasets/tmp/df.CNS_grn.seqMotifs_",i,".rds", sep = ""))){
                                        df.CNS_grn.seqMotifs <- rbind(df.CNS_grn.seqMotifs, readRDS(paste("Datasets/tmp/df.CNS_grn.seqMotifs_",i,".rds", sep = "")))
                                      }
                                    }
                                    names(df.CNS_grn.seqMotifs) <- c("TF", "cns_sequence", "Target", "cons_score", "val")#, "BS_pos_to_TSS") 
                                    #df.CNS_grn.seqMotifs <- unique(df.CNS_grn.seqMotifs)
                                    saveRDS(df.CNS_grn.seqMotifs, paste("Datasets/df.CNS_novel_seq_motifs_grn.rds", sep = ""))
                                  }
                                  
                                  
                                  
                                  combine_novel_pwm_sequence_motifs <- function(){
                                    
                                    # conservation score = 1 - only present in arabidopsis 
                                    df.cns_pwm_grn <- readRDS(paste("Datasets/df.CNS_novel_pwm_motifs_grn.rds", sep = ""))
                                    df.cns_seq_grn <- readRDS(paste("Datasets/df.CNS_novel_seq_motifs_grn.rds", sep = ""))
                                    
                                    df.cns_pwm_grn <- unique(df.cns_pwm_grn[,c(1,5,4,2)])
                                    
                                    # minimal cutoff for dna binding score is mean of distribution
                                    v.mean_bind_score <- mean(df.cns_pwm_grn$bind_score) 
                                    df.cns_pwm_grn <- subset(df.cns_pwm_grn, df.cns_pwm_grn$bind_score > v.mean_bind_score)
                                    
                                    v.seq_bind_score <- round(max(df.cns_pwm_grn$bind_score) + 1, 0)
                                    
                                    ### !!!!! COMBINE TO GET GLOBAL BINDING  ### 
                                    df.cns_seq_grn <- unique(df.cns_seq_grn[,c(1,3,4)])
                                    df.cns_seq_grn["bind_score"] <- v.seq_bind_score
                                    # identify the distance between species
                                    
                                    df.cns_grn <- rbind(df.cns_pwm_grn, df.cns_seq_grn)
                                    tfs.cns <- unique(df.cns_grn$TF)
                                    
                                    df.cns_w_counter_grn <- data.frame(TF = character(), Target = character(), cons_score = numeric(),  v.max_bindscore = numeric(), bind_counter = numeric(), v.mean_bindscore = numeric(), v.sd_bindscore = numeric()) 
                                    for(i in 1:length(tfs.cns)){
                                      cat("Processing... ", round(i/length(tfs.cns) * 100, digits = 2) , "%", "\r") 
                                      flush.console()
                                      df.cns_grn.i <- subset(df.cns_grn, df.cns_grn$TF == tfs.cns[i])
                                      tgs.cns <- unique(df.cns_grn.i$Target)
                                      for(j in 1:length(tgs.cns)){
                                        df.cns_grn.ij <- subset(df.cns_grn.i, df.cns_grn.i$TF == tfs.cns[i] & df.cns_grn.i$Target == tgs.cns[j])
                                        if(nrow(df.cns_grn.ij) > 1){
                                          newrow <- data.frame(TF = tfs.cns[i], Target = tgs.cns[j], cons_score = max(df.cns_grn.ij$cons_score), v.max_bindscore = max(df.cns_grn.ij$bind_score), bind_counter = nrow(df.cns_grn.ij), v.mean_bindscore = mean(df.cns_grn.ij$bind_score) , v.sd_bindscore = sd(df.cns_grn.ij$bind_score))
                                          df.cns_w_counter_grn <- rbind(df.cns_w_counter_grn, newrow)
                                        }
                                      } 
                                    }
                                    if(b.CALC){
                                      saveRDS(df.cns_w_counter_grn, paste("../Datasets/df.CNS_novel_motifs_w_counter.rds", sep = ""))
                                    }else{
                                      saveRDS(df.cns_w_counter_grn, paste("Datasets/df.CNS_novel_motifs_w_counter.rds", sep = ""))
                                    }
                                    
                                    
                                  }
                                  
                                  
                                  load_cns_based_grns <- function(){
                                    
                                    library(reshape2)
                                    df.cns_grn <- readRDS(paste("Datasets/df.CNS_novel_motifs_w_counter.rds", sep = ""))  
                                    
                                    mat.cns_consscore_grn <- acast(df.cns_grn, TF~Target, value.var = "cons_score")
                                    mat.cns_consscore_grn[is.na(mat.cns_consscore_grn)] <- 0
                                    class(mat.cns_consscore_grn) <- "numeric"   
                                    
                                    mat.cns_bindscore_grn <- acast(df.cns_grn, TF~Target, value.var = "v.max_bindscore")
                                    mat.cns_bindscore_grn[is.na(mat.cns_bindscore_grn)] <- 0
                                    class(mat.cns_bindscore_grn) <- "numeric"  
                                    
                                    mat.cns_counter_grn <- acast(df.cns_grn, TF~Target, value.var = "bind_counter")
                                    mat.cns_counter_grn[is.na(mat.cns_counter_grn)] <- 0
                                    class(mat.cns_counter_grn) <- "numeric"  
                                    
                                    # combine with the original cns-grn by vandepoele
                                    df.cns_grn.vandepoele <- read.table("Datasets/cns_grn_vandepoele.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
                                    
                                    df.cns_grn.vandepoele <- subset(df.cns_grn.vandepoele, !df.cns_grn.vandepoele$Target == "")
                                    mat.cns_grn.vandepoele <- acast(df.cns_grn.vandepoele, TF~Target, value.var = "cons_score")
                                    mat.cns_grn.vandepoele[is.na(mat.cns_grn.vandepoele)] <- 0
                                    class(mat.cns_grn.vandepoele) <- "numeric"   
                                    
                                    tfs.cns <- (unique(c(rownames(mat.cns_grn.vandepoele) , rownames(mat.cns_consscore_grn)))) # 356 regulators
                                    tgs.cns <- (unique(c(colnames(mat.cns_grn.vandepoele) , colnames(mat.cns_consscore_grn)))) # 12296 
                                    
                                    
                                    mat.cns_grn <- matrix(0, nrow = length(tfs.cns), ncol = length(tgs.cns), dimnames = list(tfs.cns, tgs.cns))
                                    mat.cns_grn[rownames(mat.cns_grn.vandepoele), colnames(mat.cns_grn.vandepoele)] <- mat.cns_grn.vandepoele
                                    mat.cns_grn[rownames(mat.cns_consscore_grn), colnames(mat.cns_consscore_grn)] <- mat.cns_consscore_grn
                                    
                                    return(mat.cns_grn)
                                  }
                                  
                                  
                                  cns_based_results <- function(){
                                    # conservation (can later be analyzed, because set to 1) - at least in 2 species?
                                    mat.cns_grn <- load_cns_based_grns()
                                    mat.cns_grn[mat.cns_grn > 0] <- 1 
                                    
                                    v.tfs_w_cns <- rownames(mat.cns_grn)[which(rowSums(mat.cns_grn) > 0)]
                                    
                                    hmap.tfs_w_cns <- new.env(hash=T, parent=emptyenv())
                                    for(i in 1:length(v.tfs_w_cns)){
                                      hmap.tfs_w_cns[[v.tfs_w_cns[i]]] <- i
                                    }
                                    
                                    lst.tfs_tgs_w_cns <- vector(mode = "list", length = length(v.tfs_w_cns))
                                    for(i in 1:length(v.tfs_w_cns)){
                                      cat("Processing... ", round(i/length(v.tfs_w_cns) * 100, digits = 2) , "%", "\r"); flush.console() 
                                      lst.tfs_tgs_w_cns[[i]] <- names(which(mat.cns_grn[v.tfs_w_cns[i],] == 1))
                                    }
                                    
                                    
                                    mat.cns <- matrix(0.5, nrow = length(tfs) , ncol = length(tgs), dimnames = list(tfs, tgs))
                                    mat.cns[intersect(rownames(mat.cns_grn), tfs), intersect(colnames(mat.cns_grn), tgs)] <- mat.cns_grn[intersect(rownames(mat.cns_grn), tfs), intersect(colnames(mat.cns_grn), tgs)]
                                    
                                    
                                    
                                    
                                  }
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  ############## REMAINING STUFF  ###############
                                  
                                  # compute_transitive_closure <- function(){
                                  #   
                                  #   vec.genes <- unique(c(vec.tfs, vec.tgs))
                                  #   
                                  #   #df.grn.biclust.TFB_map <- readRDS(paste("Datasets/BiClust_GRNS/df.grn.biclust.TFB_map_",e,".rds", sep = ""))
                                  #   #df.SNP_grn.filtered <- readRDS(paste("Datasets/SNP/df.grn.snp_filtered.rds", sep = ""))
                                  #   df.CNS_grn.seqMotifs <- readRDS(paste("../Datasets/CNS_GRNS/df.CNS_grn2014_SeqMotifs_filtered.rds", sep = ""))
                                  #   df.CNS_grn.filtered <- readRDS(paste("../Datasets/CNS_GRNS/df.CNS_grn2014_filtered.rds", sep = ""))
                                  #   df.grn.all_motifs_filtered <- readRDS(paste("../Datasets/novelTFmotifs/df.grn.all_motifs_filtered", sep = ""))
                                  #   df.grn.seqMotifs <- readRDS("../Datasets/novelTFmotifs/df.grn.seqMotifs.rds")
                                  #   df.grn.CNS <- readRDS(paste("../Datasets/GRMs/df.grn.CNS.rds", sep = ""))
                                  #   
                                  #   
                                  #   ### Conserved Binding 
                                  #   #   df.CNS_grn.filtered <- readRDS(paste("Datasets/CNS_GRNS/df.CNS_grn2014_filtered.rds", sep = ""))
                                  #   #   df.CNS_grn.filtered <- df.CNS_grn.filtered[,c(2,10,3)]
                                  #   #   df.CNS_grn.filtered <- subset(df.CNS_grn.filtered, df.CNS_grn.filtered$bind_score >= 6)
                                  #   #   df.CNS_grn.filtered <- df.CNS_grn.filtered[,c(1:2)]  
                                  #   #   
                                  #   #   df.CNS_grn.seqMotifs <- df.CNS_grn.seqMotifs[,c(1,3)]
                                  #   #   df.CNS_grn.seqMotifs <- unique(df.CNS_grn.seqMotifs)
                                  #   #   
                                  #   #   df.grn.CNS <- rbind(df.grn.CNS, df.CNS_grn.seqMotifs)
                                  #   #   df.grn.CNS <- rbind(df.grn.CNS, df.CNS_grn.filtered)
                                  #   #   df.grn.CNS <- unique(df.grn.CNS)
                                  #   #   
                                  #   #   
                                  #   #   df.grn.CNS <- subset(df.grn.CNS, df.grn.CNS$TF %in% rownames(mat.grn) & df.grn.CNS$Target %in% vec.genes)
                                  #   #   mat.grn.tfb.CNS <- matrix(0, nrow = length(vec.genes), ncol = length(vec.genes),
                                  #   #                             dimnames = list(vec.genes, vec.genes))
                                  #   #   
                                  #   #   for(i in 1:nrow(df.grn.CNS)){
                                  #   #     print(paste(i, "of", nrow(df.grn.CNS)))
                                  #   #     tf <- df.grn.CNS$TF[i]
                                  #   #     tg <- df.grn.CNS$Target[i]
                                  #   #     #if(mat.grn.tfb[tf, tg] < df.grn.TFB.combined$p.TFB[i]){
                                  #   #     mat.grn.tfb.CNS[tf, tg] <- 1 #df.grn.TFB.combined$p.TFB[i]  
                                  #   #     #}
                                  #   #   }
                                  #   #   saveRDS(mat.grn.tfb.CNS, "Datasets/workspace/mat.grn.tfb.CNS.all_8.rds") 
                                  #   #mat.grn.tfb.CNS <- readRDS("Datasets/workspace/mat.grn.tfb.CNS.rds") 
                                  #   #mat.grn.tfb <- readRDS("Datasets/workspace/mat.grn.tfb.probs_above_95.rds") 
                                  #   #mat.grn.tfb.CNS[(!rownames(mat.grn.tfb.CNS) %in% vec.TF_with_TFB), ] <- 0.5
                                  #   
                                  #   #### Binding probability 
                                  #   df.grn.all_motifs_filtered <- df.grn.all_motifs_filtered[,c(1,2,3)]
                                  #   
                                  #   df.CNS_grn.seqMotifs["bind_score"] <- 8
                                  #   df.CNS_grn.seqMotifs <- df.CNS_grn.seqMotifs[,c(1,3,7)]
                                  #   df.CNS_grn.seqMotifs <- unique(df.CNS_grn.seqMotifs)
                                  #   
                                  #   df.grn.seqMotifs["bind_score"] <- 7
                                  #   df.grn.seqMotifs <- df.grn.seqMotifs[,c(1,2,4)]
                                  #   df.grn.seqMotifs <- unique(df.grn.seqMotifs)
                                  #   
                                  #   df.grn.CNS["bind_score"] <- 8
                                  #   df.grn.seqMotifs <- rbind(df.grn.seqMotifs, df.grn.CNS)
                                  #   df.grn.seqMotifs <- rbind(df.grn.seqMotifs, df.CNS_grn.seqMotifs)
                                  #   
                                  #   df.grn.all_motifs_filtered <- rbind(df.grn.all_motifs_filtered, df.grn.seqMotifs)
                                  #   
                                  #   x <- df.grn.all_motifs_filtered$bind_score
                                  #   P = ecdf(x) 
                                  #   
                                  #   df.grn.TFB.combined <- df.grn.all_motifs_filtered
                                  #   df.grn.TFB.combined["p.TFB"] <- P(df.grn.all_motifs_filtered$bind_score)
                                  #   df.grn.TFB.combined <- subset(df.grn.TFB.combined, df.grn.TFB.combined$TF %in% vec.tfs & df.grn.TFB.combined$Target %in% vec.genes)
                                  #   df.grn.TFB.combined <- subset(df.grn.TFB.combined, df.grn.TFB.combined$p.TFB >= P(7))
                                  #   df.grn.TFB.combined <- unique(df.grn.TFB.combined[,c(1,2)])
                                  #   
                                  #   vec.regs.tfb <- unique(df.grn.TFB.combined$TF)
                                  #   vec.tgs.tfb <- unique(df.grn.TFB.combined$Target)
                                  #   
                                  #   mat.grn.tfb <- matrix(0, nrow = length(vec.genes), ncol = length(vec.genes),
                                  #                         dimnames = list(vec.genes, vec.genes))
                                  #   
                                  #   for(i in 1:nrow(df.grn.TFB.combined)){
                                  #     print(paste(i, "of", nrow(df.grn.TFB.combined)))
                                  #     tf <- df.grn.TFB.combined$TF[i]
                                  #     tg <- df.grn.TFB.combined$Target[i]
                                  #     mat.grn.tfb[tf, tg] <- 1 #df.grn.TFB.combined$p.TFB[i]  
                                  #   }
                                  #   saveRDS(mat.grn.tfb, "Datasets/workspace/mat.grn.tfb.probs_above_95.rds") 
                                  #   
                                  #   #mat.grn.tfb <- readRDS("Datasets/workspace/mat.grn.tfb.probs_above_95.rds") 
                                  #   #mat.grn.tfb[(!rownames(mat.grn.tfb) %in% vec.TF_with_TFB), ] <- 0.5
                                  #   
                                  #   #mat.grn.tfb <- readRDS("Datasets/workspace/mat.grn.tfb.probs_above_95.rds") 
                                  #   #mat.grn.tfb[(!rownames(mat.grn.tfb) %in% vec.TF_with_TFB), ] <- 0.5
                                  #   
                                  #   
                                  #   ###
                                  #   df.grn.CNS <- subset(df.grn.CNS, df.grn.CNS$TF %in% rownames(mat.grn) & df.grn.CNS$Target %in% vec.genes)
                                  #   mat.grn.tfb.CNS <- matrix(0, nrow = length(vec.genes), ncol = length(vec.genes),
                                  #                             dimnames = list(vec.genes, vec.genes))
                                  #   
                                  #   for(i in 1:nrow(df.grn.CNS)){
                                  #     print(paste(i, "of", nrow(df.grn.CNS)))
                                  #     tf <- df.grn.CNS$TF[i]
                                  #     tg <- df.grn.CNS$Target[i]
                                  #     #if(mat.grn.tfb[tf, tg] < df.grn.TFB.combined$p.TFB[i]){
                                  #     mat.grn.tfb.CNS[tf, tg] <- 1 #df.grn.TFB.combined$p.TFB[i]  
                                  #     #}
                                  #   }
                                  #   #saveRDS(mat.grn.tfb.CNS, "Datasets/workspace/mat.grn.tfb.CNS.rds") 
                                  #   mat.grn.tfb.CNS <- readRDS("Datasets/workspace/mat.grn.tfb.CNS.rds") 
                                  #   #mat.grn.tfb <- readRDS("Datasets/workspace/mat.grn.tfb.probs_above_95.rds") 
                                  #   mat.grn.tfb.CNS[(!rownames(mat.grn.tfb.CNS) %in% vec.TF_with_TFB), ] <- 0.5
                                  #   
                                  #   
                                  #   
                                  #   
                                  #   #df.grn.all_motifs_filtered <- rbind(df.grn.all_motifs_filtered, df.CNS_grn.seqMotifs)
                                  #   
                                  #   ####
                                  #   v.th.binding_cutoff <- 4.9
                                  #   
                                  #   
                                  #   
                                  #   #df.CNS_grn.filtered <- subset(df.CNS_grn.filtered, df.CNS_grn.filtered$bind_score >= v.th.binding_cutoff)
                                  #   #df.grn.all_motifs_filtered <- subset(df.grn.all_motifs_filtered, df.grn.all_motifs_filtered$bind_score >= v.th.binding_cutoff)
                                  #   
                                  #   df.grn.TFB.combined <- df.grn.CNS
                                  #   df.grn.TFB.combined <- rbind(df.grn.TFB.combined, df.CNS_grn.seqMotifs[,c(1,3)])
                                  #   df.grn.TFB.combined <- rbind(df.grn.TFB.combined, df.CNS_grn.filtered[,c(2,10)])
                                  #   df.grn.TFB.combined <- rbind(df.grn.TFB.combined, df.grn.all_motifs_filtered[,c(1,2)])
                                  #   df.grn.TFB.combined <- rbind(df.grn.TFB.combined, df.grn.seqMotifs[,c(1,2)])
                                  #   
                                  #   df.grn.TFB.combined <- unique(df.grn.TFB.combined)
                                  #   
                                  #   # 4.9 cutoff
                                  #   
                                  #   # + TFB (!!1/0/0.5 interconnected!!!)
                                  #   #df.grn.TFB.combined <- readRDS(paste("Datasets/GRMs/df.grn.TFB.combined.rds", sep = ""))
                                  #   
                                  #   vec.genes <- c(rownames(mat.grn), colnames(mat.grn))
                                  #   
                                  #   df.grn.TFB.combined <- subset(df.grn.TFB.combined, df.grn.TFB.combined$TF %in% rownames(mat.grn) & df.grn.TFB.combined$Target %in% vec.genes)
                                  #   vec.regs.tfb <- unique(df.grn.TFB.combined$TF)
                                  #   vec.tgs.tfb <- unique(df.grn.TFB.combined$Target)
                                  #   
                                  #   mat.grn.tfb <- matrix(0, nrow = length(vec.genes), ncol = length(vec.genes),
                                  #                         dimnames = list(vec.genes, vec.genes))
                                  #   
                                  #   for(i in 1:nrow(df.grn.TFB.combined)){
                                  #     print(paste(i, "of", nrow(df.grn.TFB.combined)))
                                  #     tf <- df.grn.TFB.combined$TF[i]
                                  #     tg <- df.grn.TFB.combined$Target[i]
                                  #     mat.grn.tfb[tf, tg] <- 1
                                  #   }
                                  #   #saveRDS(mat.grn.tfb, "Datasets/workspace/mat.grn.tfb.rds") 
                                  #   #saveRDS(mat.grn.tfb, "Datasets/workspace/mat.grn.tfb.1.rds") 
                                  #   #saveRDS(mat.grn.tfb, "Datasets/workspace/mat.grn.tfb.2.rds") 
                                  #   
                                  #   
                                  #   ###
                                  #   
                                  #   
                                  #   mat.grn.TF_TF <- mat.grn.tfb[rownames(mat.grn), rownames(mat.grn)]  
                                  #   mat.grn.TF_TF.transitive <- mat.grn.TF_TF
                                  #   for(y in 1:ncol(mat.grn.TF_TF)){
                                  #     print(y)
                                  #     tfs.y <- names(which(mat.grn.TF_TF[,y] == 1))
                                  #     if(length(tfs.y) > 0){
                                  #       for(x in 1:length(tfs.y)){
                                  #         idx.tgs <- as.numeric(which(mat.grn.TF_TF[tfs.y[x], ] == 1))
                                  #         mat.grn.TF_TF.transitive[tfs.y[x], idx.tgs] <- 1
                                  #       }  
                                  #     }
                                  #   }
                                  #   
                                  #   
                                  #   
                                  #   #mat.grn <- mat.pred.grn
                                  #   # 1 step transitive 
                                  #   mat.grn.TF_ME <- mat.grn.tfb[rownames(mat.grn), colnames(mat.grn)]  
                                  #   mat.grn.TF_ME.transitive <- mat.grn.TF_ME
                                  #   for(y in 1:ncol(mat.grn.TF_TF.transitive)){
                                  #     print(y)
                                  #     tfs.y <- names(which(mat.grn.TF_TF.transitive[,y] == 1))
                                  #     for(x in 1:length(tfs.y)){
                                  #       idx.tgs <- as.numeric(which(mat.grn.TF_ME[tfs.y[x], ] == 1))
                                  #       mat.grn.TF_ME.transitive[tfs.y[x], idx.tgs] <- 1
                                  #     }
                                  #   }
                                  #   #mat.grn.TF_ME.transitive[(!rownames(mat.grn.TF_ME.transitive) %in% vec.TF_with_TFB), ] <- 0.5
                                  #   
                                  #   
                                  #   saveRDS(mat.grn.TF_ME.transitive, "Datasets/workspace/mat.grn.TF_ME.transitive.2.rds") 
                                  #   #saveRDS(mat.grn.TF_ME.transitive, "Datasets/workspace/mat.grn.TF_ME.transitive.1.80.rds") 
                                  #   #saveRDS(mat.grn.TF_ME.transitive, "Datasets/workspace/mat.grn.TF_ME.transitive.1.95.rds") 
                                  #   
                                  #   
                                  #   
                                  #   ######
                                  #   vec.TF_with_TFB <- unique(df.grn.TFB.combined$TF)  
                                  #   mat.grn.TF_TF_with_TFB <- mat.grn.tfb[vec.TF_with_TFB, vec.TF_with_TFB]
                                  #   
                                  #   library(agop)
                                  #   strt<-Sys.time()
                                  #   mat.grn.tfb.TC <- closure_transitive(mat.grn.TF_TF_with_TFB)
                                  #   print(Sys.time()-strt)
                                  #   #saveRDS(mat.grn.tfb.TC, "Datasets/workspace/mat.grn.tfb.TC.rds") 
                                  #   #saveRDS(mat.grn.tfb.TC, "Datasets/workspace/mat.grn.tfb.TC.1.rds") 
                                  #   
                                  #   # mark a connection 
                                  #   #mat.grn.TF_TF <- mat.grn.tfb[rownames(mat.grn), rownames(mat.grn)]
                                  #   mat.grn.TF_ME <- mat.grn.tfb[rownames(mat.grn), colnames(mat.grn)]  
                                  #   mat.grn.TF_ME.transitive <- mat.grn.TF_ME
                                  #   for(y in 1:ncol(mat.grn.tfb.TC)){
                                  #     print(y)
                                  #     tfs.y <- names(which(mat.grn.tfb.TC[,y] == 1))
                                  #     for(x in 1:length(tfs.y)){
                                  #       idx.tgs <- as.numeric(which(mat.grn.TF_ME[tfs.y[x], ] == 1))
                                  #       mat.grn.TF_ME.transitive[tfs.y[x], idx.tgs] <- 1
                                  #     }
                                  #   }
                                  #   mat.grn.TF_ME.transitive[(!rownames(mat.grn.TF_ME.transitive) %in% vec.TF_with_TFB), ] <- 0.5
                                  #   
                                  #   # saveRDS(mat.grn.TF_ME.transitive, "Datasets/workspace/mat.grn.TF_ME.transitive.rds") 
                                  #   
                                  # }
                                  # 
                                  # 
                                  # 
                                  # ##### NON CONSERVATION BASED #####
                                  # 
                                  
                                  library(AnnotationDbi)
                                  library(GenomicFeatures)
                                  library(biomaRt)
                                  library(Biostrings)
                                  library(BSgenome.Athaliana.TAIR.TAIR9)
                                  genome <- BSgenome.Athaliana.TAIR.TAIR9
                                  library(TxDb.Athaliana.BioMart.plantsmart25)
                                  
                                  seqlevels(TxDb.Athaliana.BioMart.plantsmart25) <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9)
                                  transcriptCoordsByGene.GRangesList <- transcriptsBy (TxDb.Athaliana.BioMart.plantsmart25, by = "gene")
                                  isCircular(transcriptCoordsByGene.GRangesList) <- isCircular(BSgenome.Athaliana.TAIR.TAIR9)
                                  #transcriptCoordsByGene.GRangesList <- renameSeqlevels( transcriptCoordsByGene.GRangesList, c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrM", "ChrC") )
                                  transcriptCoordsByGene.GRangesList <- transcriptCoordsByGene.GRangesList[names(transcriptCoordsByGene.GRangesList) %in% "AT3G13220"]
                                  
                                  df.promSequences <- getPromoterSeq(transcriptCoordsByGene.GRangesList, genome, upstream=v.promoter_length, downstream=0)
                                  
                                  
                                  vec.genes <- c("AT5G40260", "AT1G61110")
                                  
                                  transporter AT3G13220
                                  
                                  
                                  AT3G11980
                                  MS2 ();
                                  
                                  #df.promSequences
                                  
                                  # dict_seq_motifs <- DNAStringSet("CA - AG  - TG")
                                  # 
                                  # Target <- names(df.promSequences[[i]][1])
                                  # # find these motifs in promoters - get genes
                                  # motifMatches <- matchPDict(dict_seq_motifs, DNAString(toString(df.promSequences[[1]])),max.mismatch=0, min.mismatch=0)
                                  # elementLengths(motifMatches)  
                                  
                                  
                                  #grepl("GTG",toString(df.promSequences[[2]]))
                                  sel <- gregexpr("CATATG",toString(df.promSequences[[1]]))
                                  lst.motifs <- lapply(as.list(sel[[1]]), function(m) { substr(toString(df.promSequences[[2]]), m-6, m+2)})
                                  
                                  
                                  # 
                                  # 
                                  # lst.motifs <- lapply(lst.motifs, function(m) { gregexpr("TG",m)})
                                  
                                  
                                  #attr(, 'match.length') 
                                  
                                  #attr(sel,"match.length")
                                  
                                  
                                  
                                  
                                  toString(df.promSequences[[2]])
                                  
                                  exact_map_seq_motifs <- function(vec.genes, v.promoter_length = 1000, v.th.seqMotif_length = 4){
                                    
                                    library(GenomicFeatures)
                                    library(biomaRt)
                                    library(Biostrings)
                                    library(BSgenome.Athaliana.TAIR.TAIR9)
                                    genome <- BSgenome.Athaliana.TAIR.TAIR9
                                    
                                    library(TxDb.Athaliana.BioMart.plantsmart22)
                                    transcriptCoordsByGene.GRangesList <- transcriptsBy (TxDb.Athaliana.BioMart.plantsmart22, by = "gene")
                                    transcriptCoordsByGene.GRangesList <- renameSeqlevels( transcriptCoordsByGene.GRangesList, c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrM", "ChrC") )
                                    transcriptCoordsByGene.GRangesList <- transcriptCoordsByGene.GRangesList[names(transcriptCoordsByGene.GRangesList) %in% vec.genes]
                                    df.promSequences <- getPromoterSeq(transcriptCoordsByGene.GRangesList, genome, upstream=v.promoter_length, downstream=0)
                                    
                                    
                                    
                                    df.seqMotifs <- read.table("Datasets/novelTFmotifs/athamap/athmap_seqMotifs.txt", header = TRUE, stringsAsFactor = FALSE, sep = "\t")
                                    df.seqMotifs <- subset(df.seqMotifs, nchar(df.seqMotifs$verified_binding_seq) >= v.th.seqMotif_length)
                                    df.seqMotifs$TF.locus <- toupper(df.seqMotifs$TF.locus)
                                    df.seqMotifs$TF.locus <- unlist(sapply(df.seqMotifs$TF.locus, function(x) gsub(" ","",x, fixed=TRUE)))
                                    df.seqMotifs <- subset(df.seqMotifs, !df.seqMotifs$TF.locus %in% vec.regs.novel.matrix)
                                    dict_seq_motifs <- DNAStringSet(df.seqMotifs$verified_binding_seq)
                                    #pdict_seq_motifs <- PDict(dict_seq_motifs) 
                                    foreach(i = 1:length(df.promSequences)) %dopar% { 
                                      #for(i in 1:length(df.promSequences)){
                                      df.grn.seqMotifs <- data.frame(TF = character(), seq.motif = character(), Target = character(), BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE)
                                      print(paste(i,"of",length(df.promSequences)))
                                      Target <- names(df.promSequences[[i]][1])
                                      # find these motifs in promoters - get genes
                                      motifMatches <- matchPDict(dict_seq_motifs, DNAString(toString(df.promSequences[[i]][1])),max.mismatch=0, min.mismatch=0)
                                      idxes <- countIndex(motifMatches)  
                                      # how many per promoter and position - 
                                      for(j in 1:length(idxes)){
                                        if(idxes[j] > 0){
                                          for(k in 1:idxes[j]){
                                            # SNP position in promoter
                                            mean_dist_to_TSS <- round(mean(c(endIndex(motifMatches)[[j]][k], startIndex(motifMatches)[[j]][k])))
                                            BS_pos_to_TSS <- v.promoter_length - mean_dist_to_TSS
                                            newrow <- data.frame(TF = df.seqMotifs$TF.locus[j], Target = Target, BS_pos_to_TSS = BS_pos_to_TSS)
                                            df.grn.seqMotifs <- rbind(df.grn.seqMotifs, newrow)
                                          }
                                        }
                                      }
                                      names(df.grn.seqMotifs) <- c("TF", "Target", "BS_pos_to_TSS") 
                                      saveRDS(df.grn.seqMotifs, paste("Datasets/novelTFmotifs/tmp/df.grn.seqMotifs_",i,".rds", sep = ""))
                                    }
                                    df.grn.seqMotifs <- data.frame(TF = character(), Target = character(), bind_score = numeric(), BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE)       
                                    for(i in 1:length(df.promSequences)){ 
                                      print(paste(i, "of", length(df.promSequences)))
                                      df.grn.seqMotifs <- rbind(df.grn.seqMotifs, readRDS(paste("Datasets/novelTFmotifs/tmp/df.grn.seqMotifs_",i,".rds", sep = "")))
                                    }
                                    names(df.grn.seqMotifs) <- c("TF", "Target", "BS_pos_to_TSS") 
                                    saveRDS(df.grn.seqMotifs, "Datasets/novelTFmotifs/df.grn.seqMotifs.rds")
                                  }
                                  
                                  # 
                                  # 
                                  # 
                                  # 
                                  # 
                                  # ##### mapping ###### 
                                  # 
                                  
                                  
                                  perform_TFBS_mapping_all_genes <- function(v.tgs, v.promoter_length = 1500, v.th.bind_score = 0, n.cores = 15){
                                    
                                    print("compute scored all gene GRN")
                                    
                                    library(GenomicFeatures)
                                    #source("http://bioconductor.org/biocLite.R")
                                    # biocLite("GenomicFeatures")
                                    #   biocLite("ChIPpeakAnno")
                                    #   biocLite("biomaRt")
                                    #   biocLite("Biostrings")
                                    #   install.packages("VennDiagram")
                                    #library(ChIPpeakAnno)
                                    library(biomaRt)
                                    library(Biostrings)
                                    #biocLite("BSgenome.Athaliana.TAIR.TAIR9")
                                    library(BSgenome.Athaliana.TAIR.TAIR9)
                                    genome <- BSgenome.Athaliana.TAIR.TAIR9
                                    #source("http://bioconductor.org/biocLite.R")
                                    #biocLite("TxDb.Athaliana.BioMart.plantsmart25")
                                    library(TxDb.Athaliana.BioMart.plantsmart25)
                                    
                                    seqlevels(TxDb.Athaliana.BioMart.plantsmart25) <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9)
                                    transcriptCoordsByGene.GRangesList <- transcriptsBy (TxDb.Athaliana.BioMart.plantsmart25, by = "gene")
                                    isCircular(transcriptCoordsByGene.GRangesList) <- isCircular(BSgenome.Athaliana.TAIR.TAIR9)
                                    #transcriptCoordsByGene.GRangesList <- renameSeqlevels( transcriptCoordsByGene.GRangesList, c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrM", "ChrC") )
                                    transcriptCoordsByGene.GRangesList <- transcriptCoordsByGene.GRangesList[names(transcriptCoordsByGene.GRangesList) %in% v.tgs]
                                    
                                    df.promSequences <- getPromoterSeq(transcriptCoordsByGene.GRangesList, genome, upstream=v.promoter_length, downstream=0)
                                    
                                    #install.packages("foreach")
                                    
                                    lst.motifs <- get_cell_and_pnas_Paper_PWMs()
                                    lst.pwm.motif <- lst.motifs[[1]]
                                    df.motifs <- lst.motifs[[2]]
                                    
                                    
                                    # adapt to current approach
                                    library(foreach)
                                    library(Biostrings)
                                    
                                    
                                    cl<-makeCluster(n.cores)
                                    registerDoParallel(cl)
                                    
                                    #lst.crfs <- foreach(l = 1:length(v.lambda), .packages=c("Matrix", "reshape2", "CRF", "caTools", "ROCR", "ggplot2", "igraph")) %dopar% { 
                                    
                                    lst.dna_binding <- foreach(i = 1:length(lst.pwm.motif), .packages=c("Biostrings")) %dopar% { 
                                      
                                      #for(i in 1:length(lst.pwm.motif)){
                                      print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
                                      df.grn.all_motifs <- data.frame(TF = character(), Target = character(), TF.Motif = character(), 
                                                                      bind_score = numeric(), cmp_bind_score = numeric(), BS_pos_to_TSS = numeric(),
                                                                      stringsAsFactors = FALSE)    
                                      names(df.grn.all_motifs) <- c("TF", "Target", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_to_TSS")
                                      pcm <- round(100 * lst.pwm.motif[[i]])
                                      for(j in 1:length(df.promSequences)){
                                        cat("Processing... ", round(j/length(df.promSequences) * 100, digits = 2) , "%", "\r"); flush.console() 
                                        # print(paste(j, "of", length(df.promSequences)))
                                        #hits <- matchPWM(lst.pwm.motif[[i]], DNAString(toString(df.promSequences[[j]][[1]])) , with.score = TRUE)
                                        hits <- matchPWM(pcm, DNAString(toString(df.promSequences[[j]][[1]])) , with.score = TRUE)
                                        #matchPWM(reverseComplement(pcm), DNAString(toString(df.promSequences[[j]][[1]])) , with.score = TRUE)
                                        nhits <- length(hits)  
                                        if(nhits >= 1){
                                          # cmp_bind_score <- min(mcols(hits)$score / maxScore(lst.pwm.motif[[i]])) # should be >= 0.8
                                          cmp_bind_score <- min(mcols(hits)$score / maxScore(pcm)) # should be >= 0.8
                                          motif.score <- mcols(hits)$score / 100
                                          if(cmp_bind_score >= 0.8){
                                            for(k in 1:nhits){
                                              if(motif.score[k] >= v.th.bind_score){
                                                newrow <- data.frame(TF =  as.character(df.motifs$TF.locus[i]), 
                                                                     Target = names(df.promSequences[[j]])[1],
                                                                     TF.Motif =  names(lst.pwm.motif)[i],
                                                                     bind_score = motif.score[k], 
                                                                     cmp_bind_score = cmp_bind_score,
                                                                     BS_pos_to_TSS = abs(end(hits)[k] + start(hits)[k])/2,
                                                                     stringsAsFactors = FALSE)
                                                df.grn.all_motifs <- rbind(df.grn.all_motifs, newrow)   
                                              }
                                            }
                                          }
                                        }
                                      }
                                      names(df.grn.all_motifs) <- c("TF", "Target", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_to_TSS")
                                      #saveRDS(df.grn.all_motifs, paste("Datasets/novelTFmotifs/tmp/df.grn.all_motifs_",i,".rds", sep = ""))
                                      df.grn.all_motifs
                                    }
                                    
                                    
                                    stopCluster(cl)
                                    #saveRDS(lst.dna_binding, "Datasets/novelTFmotifs/lst.dna_binding.rds")
                                    lst.dna_binding <- readRDS("Datasets/novelTFmotifs/lst.dna_binding.rds")
                                    
                                    
                                    
                                    #df.grn.dnabinding_map <- do.call("rbind", lst.dna_binding)
                                    
                                    df.grn.dnabinding_map <- data.frame(TF = character(), Target = character(), TF.Motif = character(), 
                                                                        bind_score = numeric(), cmp_bind_score = numeric(), BS_pos_to_TSS = numeric(),
                                                                        stringsAsFactors = FALSE)      
                                    for(i in 1:length(lst.pwm.motif)){ 
                                      cat("Processing... ", round(i/length(lst.pwm.motif) * 100, digits = 2) , "%", "\r") 
                                      df.grn.dnabinding_map <- rbind(df.grn.dnabinding_map, lst.dna_binding[[i]])
                                    }
                                    names(df.grn.dnabinding_map) <- c("TF", "Target", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_to_TSS") 
                                    #saveRDS(df.grn.dnabinding_map, paste("Datasets/novelTFmotifs/df.grn.dnabinding_map.rds", sep = ""))
                                    
                                    # as with cns apply a mean based filter
                                    #df.grn.dnabinding_map <- subset(df.grn.dnabinding_map, df.grn.dnabinding_map$bind_score > mean(df.grn.dnabinding_map$bind_score))
                                    df.grn.dnabinding_map <- readRDS(paste("Datasets/novelTFmotifs/df.grn.dnabinding_map.rds", sep = ""))
                                    
                                    
                                    
                                    ##
                                    
                                    l.grn <- nrow(unique(df.grn.dnabinding_map[,1:2]))
                                    df.dnabind_w_counter_grn <- data.frame(TF = character(l.grn), Target = character(l.grn), v.max_bindscore = numeric(l.grn), bind_counter = numeric(l.grn), v.mean_bindscore = numeric(l.grn), v.sd_bindscore = numeric(l.grn), stringsAsFactors = FALSE) 
                                    tfs.dna_binding <- unique(df.grn.dnabinding_map$TF)
                                    
                                    idx <- 1
                                    for(i in 1:length(tfs.dna_binding)){
                                      cat("Processing... ", round(i/length(tfs.dna_binding) * 100, digits = 2) , "%", "\r") 
                                      flush.console()
                                      df.grn.dnabinding_map.i <- subset(df.grn.dnabinding_map, df.grn.dnabinding_map$TF == tfs.dna_binding[i])
                                      tgs.dna_binding <- unique(df.grn.dnabinding_map.i$Target)
                                      for(j in 1:length(tgs.dna_binding)){
                                        df.grn.dnabinding_map.ij <- subset(df.grn.dnabinding_map.i, df.grn.dnabinding_map.i$TF == tfs.dna_binding[i] & df.grn.dnabinding_map.i$Target == tgs.dna_binding[j])
                                        if(nrow(df.grn.dnabinding_map.ij) > 1){
                                          df.dnabind_w_counter_grn$TF[idx] = tfs.dna_binding[i]
                                          df.dnabind_w_counter_grn$Target[idx] = tgs.dna_binding[j]
                                          df.dnabind_w_counter_grn$v.max_bindscore[idx] = max(df.grn.dnabinding_map.ij$bind_score)
                                          df.dnabind_w_counter_grn$bind_counter[idx] = nrow(df.grn.dnabinding_map.ij)
                                          df.dnabind_w_counter_grn$v.mean_bindscore[idx] = mean(df.grn.dnabinding_map.ij$bind_score)
                                          df.dnabind_w_counter_grn$v.sd_bindscore[idx] = sd(df.grn.dnabinding_map.ij$bind_score)
                                          idx <- idx + 1
                                        }
                                      } 
                                    }
                                    
                                    saveRDS(df.dnabind_w_counter_grn, "df.dnabind_w_counter_grn.rds")
                                    ##
                                    
                                  }
                                  
                                  
                                  load_dnabinding_based_grns <- function(){
                                    
                                    library(reshape2)
                                    
                                    df.dnabind_w_counter_grn <- readRDS("df.dnabind_w_counter_grn.rds")
                                    df.dnabind_w_counter_grn <- df.dnabind_w_counter_grn[,1:3]
                                    
                                    
                                    mat.dna_bindscore_grn <- matrix(0, nrow = length(unique(df.dnabind_w_counter_grn$TF)), ncol = length(unique(df.dnabind_w_counter_grn$Target)), dimnames = list(unique(df.dnabind_w_counter_grn$TF), unique(df.dnabind_w_counter_grn$Target)))
                                    for(i in 1:nrow(df.dnabind_w_counter_grn)){
                                      cat("Processing... ", round(i/nrow(df.dnabind_w_counter_grn) * 100, digits = 2) , "%", "\r"); flush.console() 
                                      mat.dna_bindscore_grn[df.dnabind_w_counter_grn$TF[i], df.dnabind_w_counter_grn$Target[i]] <- df.dnabind_w_counter_grn$v.max_bindscore[i]
                                    }
                                    mat.dna_bindscore_grn <- as(mat.dna_bindscore_grn, "CsparseMatrix")
                                    saveRDS(mat.dna_bindscore_grn, "Datasets/mat.dna_bindscore_grn.rds")
                                    
                                    #   
                                    #   mat.dna_bindscore_grn <- acast(df.dnabind_w_counter_grn, TF~Target, value.var = "v.max_bindscore")
                                    #   mat.dna_bindscore_grn[is.na(mat.dna_bindscore_grn)] <- 0
                                    #   class(mat.dna_bindscore_grn) <- "numeric"  
                                    
                                    
                                    
                                    #   
                                    #   mat.dna_counter_grn <- acast(df.dnabind_w_counter_grn, TF~Target, value.var = "bind_counter")
                                    #   mat.dna_counter_grn[is.na(mat.dna_counter_grn)] <- 0
                                    #   class(mat.dna_counter_grn) <- "numeric"  
                                    
                                    
                                    
                                    return(mat.cns_grn)
                                  }
                                  
                                  