# TODO: Add comment
# 
# Author: michaelbanf
###############################################################################

# load a multiset of binding motifs - last curated spring 2016
get_cell_and_pnas_Paper_PWMs <- function(){
  
  print("prepare PWM motifs and mappings")
  library(fume)
  ## create motif - gene mapping
  df.motifs <- data.frame(Motif.id = character(), TF.locus = character(), TF.name = character(), TF.family = character(), TF.subfamily = character(), src = character(), stringsAsFactors = FALSE)
  #lst.motifs <- vector(mode = "list")
  
  # Cell paper
	files <- list.files("Datasets/novelTFmotifs/Ath_CellPaper/pwms_all_motifs/")
	motif.mapping <- read.table("Datasets/novelTFmotifs/Ath_CellPaper/TF_Information_all_motifs.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
	motif.mapping$DBID <- gsub("\\-.*", "", motif.mapping$DBID)
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
  

  # JASPAR 2014
  # map tf name to .. 
  # 
  #source("http://bioconductor.org/biocLite.R")
  #biocLite("JASPAR2014")
  #biocLite("TFBSTools")
  #biocLite("RSQLite")
  #library(RSQLite)
  library(TFBSTools)
  #library(JASPAR2014)
#   opts = list()
#   opts[["tax_group"]] = "plants"
#   PFMatrixList = getMatrixSet(JASPAR2014, opts)
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


  #df.pwms.jaspar2014 <- read.table("Datasets/novelTFmotifs/Jaspar2014/pfm_plants_JASPAR_Ath_only.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
  #PFMatrixList.Ath <- PFMatrixList[vec.jaspar.motifs]
  
  # saveRDS(PFMatrixList.Ath, "Datasets/novelTFmotifs/Jaspar2014/PFMatrixList.Ath.rds")
  PFMatrixList.Ath <- readRDS("Datasets/novelTFmotifs/Jaspar2014/PFMatrixList.Ath.rds")  
  df.jaspar.motif.info <- read.table("Datasets/novelTFmotifs/Jaspar2014/jaspar2014AthMotifTable.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) 

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

	return(list(lst.pwm.motif, df.motifs))
  
}


##### mapping ###### 
perform_TFBS_mapping_all_genes <- function(vec.genes, v.promoter_length = 1000, v.th.bind_score = 4, nr.cores = 15){
  
  print("compute scored all gene GRN")
  
  library(GenomicFeatures)
  # source("http://bioconductor.org/biocLite.R")
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
  #biocLite("TxDb.Athaliana.BioMart.plantsmart21")
  library(TxDb.Athaliana.BioMart.plantsmart21)
  
  transcriptCoordsByGene.GRangesList <- transcriptsBy (TxDb.Athaliana.BioMart.plantsmart21, by = "gene")
  transcriptCoordsByGene.GRangesList <- renameSeqlevels( transcriptCoordsByGene.GRangesList, c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrM", "ChrC") )
  transcriptCoordsByGene.GRangesList <- transcriptCoordsByGene.GRangesList[names(transcriptCoordsByGene.GRangesList) %in% vec.genes]
  
  df.promSequences <- getPromoterSeq(transcriptCoordsByGene.GRangesList, genome, upstream=v.promoter_length, downstream=0)
  
  #install.packages("foreach")
  library(foreach)
  library(doMC)
  registerDoMC(nr.cores)
  
  lst.motifs <- get_cell_and_pnas_Paper_PWMs()
  lst.pwm.motif <- lst.motifs[[1]]
  df.motifs <- lst.motifs[[2]]
  
  foreach(i = 1:length(lst.pwm.motif)) %dopar% { 
  #for(i in 1:length(lst.pwm.motif)){
    #system.time(for(i in 1:2){
    print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
    df.grn.all_motifs <- data.frame(TF = character(), Target = character(), TF.Motif = character(), 
                                          bind_score = numeric(), cmp_bind_score = numeric(), BS_pos_to_TSS = numeric(),
                                          stringsAsFactors = FALSE)    
    names(df.grn.all_motifs) <- c("TF", "Target", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_to_TSS")
    pcm <- round(100 * lst.pwm.motif[[i]])
    
    for(j in 1:length(df.promSequences)){
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
    saveRDS(df.grn.all_motifs, paste("Datasets/novelTFmotifs/tmp/df.grn.all_motifs_",i,".rds", sep = ""))
  }
  # combine
  df.grn.all_motifs <- data.frame(TF = character(), Target = character(), TF.Motif = character(), 
                                        bind_score = numeric(), cmp_bind_score = numeric(), BS_pos_to_TSS = numeric(),
                                        stringsAsFactors = FALSE)      
  for(i in 1:length(lst.pwm.motif)){ 
    print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
    df.grn.all_motifs <- rbind(df.grn.all_motifs, readRDS(paste("Datasets/novelTFmotifs/tmp/df.grn.all_motifs_",i,".rds", sep = "")))
  }
  names(df.grn.all_motifs) <- c("TF", "Target", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_to_TSS") 
  saveRDS(df.grn.all_motifs, paste("Datasets/novelTFmotifs/df.grn.all_motifs.rds", sep = ""))
  
  
  # Filter
  df.grn.all_motifs <- readRDS(paste("Datasets/novelTFmotifs/df.grn.all_motifs.rds", sep = ""))
  df.grn.all_motifs <- df.grn.all_motifs[,-c(3,5)]
  names(df.grn.all_motifs) <- c("TF", "Target", "bind_score", "BS_pos_to_TSS") 
  
  vec.targs <- unique(df.grn.all_motifs$Target)
  n.targs <- length(unique(df.grn.all_motifs$Target))
  
  foreach(i = 1:n.targs) %dopar% { 
  #for(i in 1:10){#n.BSpos){ 
    #for(i in 1:length(vec.sequences)){
    df.grn.all_motifs_filtered <- data.frame(TF = character(), Target = character(), bind_score = numeric(), BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE)   
    print(paste(i, "of",n.targs))
    df.sset <- subset(df.grn.all_motifs, df.grn.all_motifs$Target == vec.targs[i])
    if(nrow(df.sset) > 0){
      tfs <- unique(df.sset$TF)
      for(t.1 in 1:length(tfs)){
        df.sset.2 <- subset(df.sset, df.sset$TF == tfs[t.1])
        targs <- unique(df.sset.2$Target)  
        for(t.2 in 1:length(targs)){
          df.sset.3 <- subset(df.sset.2, df.sset.2$Target == targs[t.2])
          vec.BSpos <- unique(df.sset.3$BS_pos_to_TSS)
          for(t.3 in 1:length(vec.BSpos)){
            df.sset.4 <- subset(df.sset.3, df.sset.3$BS_pos_to_TSS == vec.BSpos[t.3])
            df.sset.4 <- subset(df.sset.4, df.sset.4$bind_score == max(unique(df.sset.4$bind_score))) # could be more then one... filter
            df.grn.all_motifs_filtered <- rbind(df.grn.all_motifs_filtered, df.sset.4[1,])
          }
        }
      }
    }
    saveRDS(df.grn.all_motifs_filtered, paste("Datasets/novelTFmotifs/tmp/df.grn.all_motifs_filtered_",i,".rds", sep = ""))
  }
  
  df.grn.all_motifs_filtered <- data.frame(TF = character(), Target = character(), bind_score = numeric(), BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE)       
  for(i in 1:n.targs){ 
    print(paste(i, "of", n.targs))
    df.grn.all_motifs_filtered <- rbind(df.grn.all_motifs_filtered, readRDS(paste("Datasets/novelTFmotifs/tmp/df.grn.all_motifs_filtered_",i,".rds", sep = "")))
  }
  names(df.grn.all_motifs_filtered) <- c("TF", "Target", "bind_score", "BS_pos_to_TSS")
  saveRDS(df.grn.all_motifs_filtered, paste("Datasets/novelTFmotifs/df.grn.all_motifs_filtered", sep = ""))
}

exact_map_seq_motifs <- function(vec.genes, v.promoter_length = 1000, v.th.seqMotif_length = 4){
  
  library(GenomicFeatures)
  library(biomaRt)
  library(Biostrings)
  library(BSgenome.Athaliana.TAIR.TAIR9)
  genome <- BSgenome.Athaliana.TAIR.TAIR9
  
  library(TxDb.Athaliana.BioMart.plantsmart21)
  transcriptCoordsByGene.GRangesList <- transcriptsBy (TxDb.Athaliana.BioMart.plantsmart21, by = "gene")
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




##### Conservation based ######

extract_CNS2014 <- function(load_from_file = TRUE, v.promoter_length = 1000){
  
  if(!load_from_file){
    
    source("http://bioconductor.org/biocLite.R")
    #   biocLite("ChIPpeakAnno")
    #   biocLite("biomaRt")
    #   biocLite("Biostrings")
    #   install.packages("VennDiagram")
    library(ChIPpeakAnno)
    library(biomaRt)
    library(Biostrings)
    #biocLite("BSgenome.Athaliana.TAIR.TAIR9")
    library(BSgenome.Athaliana.TAIR.TAIR9)
    genome <- BSgenome.Athaliana.TAIR.TAIR9
    
    ensmart = useMart('ENSEMBL_MART_PLANT', "athaliana_eg_gene")
    df.CNS_2014 <- read.table("../Datasets/novelTFmotifs/AllFootPrintsFDR0.10_scores.bed", header = FALSE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
    names(df.CNS_2014) <- c("Chr", "start", "end", "val", "cons_score")
    g1.r <- BED2RangedData(df.CNS_2014, header=FALSE)
    annotatedData = getAnnotation(ensmart, featureType = "TSS")
    annotatedPeaks = annotatePeakInBatch(g1.r, AnnotationData= annotatedData)
    #annotatedPeaks <- annotatedPeaks[order(rownames(annotatedPeaks)),]
    annotatedPeaks.sset  = 	annotatedPeaks[!is.na(annotatedPeaks$distancetoFeature)  & 
                                             annotatedPeaks$fromOverlappingOrNearest == "NearestStart" & 
                                             annotatedPeaks$distancetoFeature < 0 &
                                             abs(annotatedPeaks$distancetoFeature) < v.promoter_length,]
    df.CNS_2014["Target"] <- NA
    df.CNS_2014["cns_sequence"] <- NA
    df.CNS_2014["start_dist_to_TSS"] <- NA
    df.CNS_2014["shortest_dist_to_TSS"] <- NA
    df.CNS_2014["feature_position"] <- NA
    
    for(i in 1:nrow(annotatedPeaks.sset)){
      print(paste(i, " of", nrow(annotatedPeaks.sset)))
      nr.row <- as.numeric(gsub("\\ .*", "",  rownames(annotatedPeaks.sset)[i]))
      chr <- df.CNS_2014$Chr[nr.row]
      start <- df.CNS_2014$start[nr.row]
      end  <- df.CNS_2014$end[nr.row]
      df.CNS_2014$Target[nr.row] <- as.character(annotatedPeaks.sset$feature[i])
      genome_by_chrom <- DNAString(genome[[chr]]) 
      df.CNS_2014$cns_sequence[nr.row] <- as.character(substring(genome_by_chrom,start,end))
      df.CNS_2014$start_dist_to_TSS[nr.row] <- abs(as.numeric(annotatedPeaks.sset$distancetoFeature[i]))
      df.CNS_2014$shortest_dist_to_TSS[nr.row]   <- abs(as.numeric(annotatedPeaks.sset$shortestDistance[i]))
      df.CNS_2014$feature_position[nr.row]   <- as.character(annotatedPeaks.sset$insideFeature[i])
    }
    df.CNS_2014 <- subset(df.CNS_2014, !is.na(df.CNS_2014$cns_sequence))
    saveRDS(df.CNS_2014, "Datasets/novelTFmotifs/df.CNS_2014.rds")
  }else{
    #df.CNS_2014 <- readRDS("Datasets/novelTFmotifs/df.CNS_2014_2000kb_sset.rds")  
    df.CNS_2014 <- readRDS("Datasets/novelTFmotifs/df.CNS_2014.rds")  
  }
}

perform_TFBS_mapping_CNS2014 <- function(vec.genes, v.th.bind_score = 4, nr.cores = 15){
  
  print("compute scored CNS2014-GRN")
  
  #install.packages("foreach")
  library(foreach)
  library(doMC)
  registerDoMC(nr.cores)
  
  lst.motifs <- get_cell_and_pnas_Paper_PWMs()
  lst.pwm.motif <- lst.motifs[[1]]
  df.motifs <- lst.motifs[[2]]
  
  df.CNS_2014 <- readRDS("Datasets/novelTFmotifs/df.CNS_2014.rds") 
  df.CNS_2014 <- subset(df.CNS_2014, df.CNS_2014$Target %in% vec.genes)
  # df.CNS_2014 <- readRDS(paste("Datasets/novelTFmotifs/df.CNS_2014_2000kb.rds", sep = ""))
  
  #df.CNS_2014 <- subset(df.CNS_2014, df.CNS_2014$end_dist_to_TSS < 1000)  
  vec.cns.sequences <- as.character(unique(df.CNS_2014$cns_sequence))
  n.cns.seq <- length(vec.cns.sequences)
  
  #system.time(for(i in 1:1){#length(lst.pwm.motif)){ 
  #foreach(i = 1:2) %dopar% { 
  foreach(i = 1:length(lst.pwm.motif)) %dopar% { 
    print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
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
    saveRDS(df.grn.CNS2014.TFB_map, paste("Datasets/CNS_GRNS/tmp/df.cns2014_grn_map_",i,".rds", sep = ""))
  }
  
  # combine
  df.grn.CNS2014.TFB_map <- data.frame(TF = character(), cns_sequence = character(), TF.Motif = character(), 
                                       bind_score = numeric(), cmp_bind_score = numeric(),
                                       stringsAsFactors = FALSE)      
  for(i in 1:length(lst.pwm.motif)){ 
    print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
    df.grn.CNS2014.TFB_map <- rbind(df.grn.CNS2014.TFB_map, readRDS(paste("Datasets/CNS_GRNS/tmp/df.cns2014_grn_map_",i,".rds", sep = "")))
  }
  names(df.grn.CNS2014.TFB_map) <- c("TF", "cns_sequence", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_in_CNS") 
  #saveRDS(df.grn.CNS2014.TFB_map, paste("Datasets/CNS_GRNS/df.grn.CNS2014.TFB_map.rds", sep = ""))
  saveRDS(df.grn.CNS2014.TFB_map, paste("Datasets/CNS_GRNS/df.grn.CNS2014.TFB_map_complete.rds", sep = ""))
  
  # finalize
  df.grn.CNS2014.TFB_map <- readRDS(paste("Datasets/CNS_GRNS/df.grn.CNS2014.TFB_map_complete.rds", sep = ""))
  df.grn.CNS2014.TFB_map <- df.grn.CNS2014.TFB_map[,c(-3,-5)]
  df.grn.CNS2014.TFB_map <- unique(df.grn.CNS2014.TFB_map)
  
  #df.CNS_2014["BS_pos_to_TSS"] <-  apply(df.CNS_2014[,c('start_dist_to_TSS','end_dist_to_TSS')], 1, function(x) mean(x) )
  #df.CNS_2014 <- df.CNS_2014[, c(-1,-2,-3,-8,-9,-10)]
  
  df.CNS_grn <- merge(df.grn.CNS2014.TFB_map, df.CNS_2014, by = "cns_sequence")
  df.CNS_grn["BS_pos_to_TSS"] <- NA
   
  for(i in 1:nrow(df.CNS_grn)){
    print(paste(i, "of", nrow(df.CNS_grn)))
   if(df.CNS_grn$start_dist_to_TSS[i] == df.CNS_grn$shortest_dist_to_TSS[i]){
     df.CNS_grn$BS_pos_to_TSS[i] <- df.CNS_grn$shortest_dist_to_TSS[i] + df.CNS_grn$BS_pos_in_CNS[i]
   }else{
     df.CNS_grn$BS_pos_to_TSS[i] <- df.CNS_grn$start_dist_to_TSS[i] - df.CNS_grn$BS_pos_in_CNS[i] 
   }
  }
  
  saveRDS(df.CNS_grn, paste("Datasets/CNS_GRNS/df.CNS_grn2014_complete.rds", sep = ""))
  
  vec.sequences <- unique(df.CNS_grn$cns_sequence)
  
  
  foreach(i = 1:length(vec.sequences)) %dopar% { 
  #for(i in 1:length(vec.sequences)){
    df.CNS_grn.filtered <- data.frame(TF = character(), bind_score = numeric(),
                                      val = numeric(), cons_score = numeric(), Target = character(),
                                      BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE)  
    print(paste(i, "of",length(vec.sequences)))
    df.sset <- subset(df.CNS_grn, df.CNS_grn$cns_sequence == vec.sequences[i])
    if(nrow(df.sset) > 0){
      tfs <- unique(df.sset$TF)
      for(t.1 in 1:length(tfs)){
        df.sset.2 <- subset(df.sset, df.sset$TF == tfs[t.1])
        targs <- unique(df.sset.2$Target)  
        for(t.2 in 1:length(targs)){
          df.sset.3 <- subset(df.sset.2, df.sset.2$Target == targs[t.2])
          vec.pos <- unique(df.sset.3$shortest_dist_to_TSS)
          for(p in 1:length(vec.pos)){
            df.sset.4 <- subset(df.sset.3, df.sset.3$shortest_dist_to_TSS == vec.pos[p])
            df.sset.4 <- subset(df.sset.4, df.sset.4$bind_score == max(unique(df.sset.4$bind_score))) # could be more then one... filter
            df.CNS_grn.filtered <- rbind(df.CNS_grn.filtered, df.sset.4[1,])
          }
        }
      }
    }
    saveRDS(df.CNS_grn.filtered, paste("Datasets/CNS_GRNS/tmp/df.cns2014_grn_filtered_",i,".rds", sep = ""))
  }
      
  
  df.CNS_grn.filtered <- data.frame(cns_sequence = character(), TF = character(), bind_score = numeric(), BS_pos_in_CNS = numeric(),
                                     Chr = numeric(), start = numeric(), end = numeric(), val = numeric(), cons_score = numeric(), Target = character(),
                                    start_dist_to_TSS = numeric(), shortest_dist_to_TSS = numeric(), feature_position = character(), BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE) 
  for(i in 1:length(vec.sequences)){ 
    print(paste(i, "of", length(vec.sequences)))
    df.CNS_grn.filtered <- rbind(df.CNS_grn.filtered, readRDS(paste("Datasets/CNS_GRNS/tmp/df.cns2014_grn_filtered_",i,".rds", sep = "")))
  }
  names(df.CNS_grn.filtered) <- c("cns_sequence", "TF", "bind_score", "BS_pos_in_CNS", "Chr", "start", "end", "val", 
                                  "cons_score", "Target", "start_dist_to_TSS", "shortest_dist_to_TSS", "feature_position", "BS_pos_to_TSS")
  saveRDS(df.CNS_grn.filtered, paste("Datasets/CNS_GRNS/df.CNS_grn2014_filtered.rds", sep = ""))
}

exact_CNS2014map_seq_motifs <- function(vec.genes, v.th.seqMotif_length = 4, nr.cores = 15){
  
  library(GenomicFeatures)
  library(biomaRt)
  library(Biostrings)
  
  library(foreach)
  library(doMC)
  registerDoMC(nr.cores)
  
  df.seqMotifs <- read.table("Datasets/novelTFmotifs/athamap/athmap_seqMotifs.txt", header = TRUE, stringsAsFactor = FALSE, sep = "\t")
  df.seqMotifs <- subset(df.seqMotifs, nchar(df.seqMotifs$verified_binding_seq) >= v.th.seqMotif_length)
  df.seqMotifs$TF.locus <- toupper(df.seqMotifs$TF.locus)
  df.seqMotifs$TF.locus <- unlist(sapply(df.seqMotifs$TF.locus, function(x) gsub(" ","",x, fixed=TRUE)))
  df.seqMotifs <- subset(df.seqMotifs, !df.seqMotifs$TF.locus %in% vec.regs.novel.matrix)
  dict_seq_motifs <- DNAStringSet(df.seqMotifs$verified_binding_seq)
  
  df.CNS_2014 <- readRDS("Datasets/novelTFmotifs/df.CNS_2014.rds") 
  df.CNS_2014 <- subset(df.CNS_2014, df.CNS_2014$Target %in% vec.genes)
  # df.CNS_2014 <- readRDS(paste("Datasets/novelTFmotifs/df.CNS_2014_2000kb.rds", sep = ""))
  
  #df.CNS_2014 <- subset(df.CNS_2014, df.CNS_2014$end_dist_to_TSS < 1000)  
  vec.cns.sequences <- as.character(unique(df.CNS_2014$cns_sequence))
  n.cns.seq <- length(vec.cns.sequences)
  
  #pdict_seq_motifs <- PDict(dict_seq_motifs) 
  #foreach(i = 1:nrow(df.CNS_2014)) %dopar% { 
  for(i in 1:nrow(df.CNS_2014)){
    df.CNS_grn.seqMotifs <- data.frame(TF = character(), cns_sequence = character(), Target = character(), cons_score = numeric(),  val = numeric(),  BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE)
    print(paste(i,"of",nrow(df.CNS_2014)))
    Target <- df.CNS_2014$Target[i]
    # find these motifs in promoters - get genes
    motifMatches <- matchPDict(dict_seq_motifs, DNAString(toString(df.CNS_2014$cns_sequence[i])),max.mismatch=0, min.mismatch=0)
    idxes <- countIndex(motifMatches)  
    # how many per promoter and position - 
    for(j in 1:length(idxes)){
      if(idxes[j] > 0){
        for(k in 1:idxes[j]){
          # SNP position in promoter
          mean_dist_to_TSS <- round(mean(c(endIndex(motifMatches)[[j]][k], startIndex(motifMatches)[[j]][k])))
          if(df.CNS_2014$start_dist_to_TSS[i] == df.CNS_2014$shortest_dist_to_TSS[i]){
            BS_pos_to_TSS <- df.CNS_2014$shortest_dist_to_TSS[i] + mean_dist_to_TSS
          }else{
            BS_pos_to_TSS <- df.CNS_2014$start_dist_to_TSS[i] - mean_dist_to_TSS
          }
          newrow <- data.frame(TF = df.seqMotifs$TF.locus[j], cns_sequence = df.CNS_2014$cns_sequence[i], Target = Target, 
                               cons_score = df.CNS_2014$cons_score[i], val = df.CNS_2014$val[i], BS_pos_to_TSS = BS_pos_to_TSS)
          df.CNS_grn.seqMotifs <- rbind(df.CNS_grn.seqMotifs, newrow)
        }
      }
    }
    names(df.CNS_grn.seqMotifs) <- c("TF", "cns_sequence", "Target", "cons_score", "val", "BS_pos_to_TSS") 
    saveRDS(df.CNS_grn.seqMotifs, paste("Datasets/CNS_GRNS/tmp/df.CNS_grn.seqMotifs_",i,".rds", sep = ""))
  }
  
  df.CNS_grn.seqMotifs <- data.frame(TF = character(), Target = character(), bind_score = numeric(), BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE)       
  for(i in 1:nrow(df.CNS_2014)){ 
    print(paste(i, "of", nrow(df.CNS_2014)))
    df.CNS_grn.seqMotifs <- rbind(df.CNS_grn.seqMotifs, readRDS(paste("Datasets/CNS_GRNS/tmp/df.CNS_grn.seqMotifs_",i,".rds", sep = "")))
  }
  names(df.CNS_grn.seqMotifs) <- c("TF", "cns_sequence", "Target", "cons_score", "val", "BS_pos_to_TSS") 
  df.CNS_grn.seqMotifs <- unique(df.CNS_grn.seqMotifs)
  saveRDS(df.CNS_grn.seqMotifs, paste("Datasets/CNS_GRNS/df.CNS_grn2014_SeqMotifs_filtered.rds", sep = ""))
}

perform_TFBS_mapping_SNP <- function(vec.genes, v.th.bind_score = 4, nr.cores = 15){
  
  print("compute scored all gene GRN")
  
  library(GenomicFeatures)
  library(biomaRt)
  library(Biostrings)
  
  library(foreach)
  library(doMC)
  registerDoMC(nr.cores)
  
  #df.snp_motifs <- read.table("Datasets/novel_snp_500_motifs.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactor = FALSE)
  df.founds <- read.table("Datasets/novel_snp_500_motifs_promoter_matches.txt",sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  vec.snp.sequences <- as.character(unique(df.founds$ext.motif))
  n.snp.seq <- length(vec.snp.sequences)

  lst.motifs <- get_cell_and_pnas_Paper_PWMs()
  lst.pwm.motif <- lst.motifs[[1]]
  df.motifs <- lst.motifs[[2]]
  
  print("perform TFB - SNP mapping")
  
  foreach(i = 1:length(lst.pwm.motif)) %dopar% { 
  #system.time(for(i in 1:1){#length(lst.pwm.motif)){
    print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
    df.grn.snp <- data.frame(TF = character(), Target = character(), TF.Motif = character(), 
                             bind_score = numeric(), cmp_bind_score = numeric(), BS_pos_in_SNP = numeric(),
                             stringsAsFactors = FALSE)    
    names(df.grn.snp) <- c("TF", "snp_sequence", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_in_SNP")
    pcm <- round(100 * lst.pwm.motif[[i]])
    ###
    for(j in 1:n.snp.seq){
      print(j)
      hits <- matchPWM(pcm, vec.snp.sequences[j], with.score = TRUE)
      nhits <- length(hits)  
      if(nhits >= 1){
        cmp_bind_score <- min(mcols(hits)$score / maxScore(pcm)) # should be >= 0.8
        motif.score <- mcols(hits)$score / 100
        if(cmp_bind_score >= 0.8){
          for(k in 1:nhits){
            if(motif.score[k] >= v.th.bind_score){
              newrow <- data.frame(TF =  as.character(df.motifs$TF.locus[i]), 
                                   snp_sequence = vec.snp.sequences[j],
                                   TF.Motif =  names(lst.pwm.motif)[i],
                                   bind_score = motif.score[k], 
                                   cmp_bind_score = cmp_bind_score,
                                   BS_pos_in_SNP = abs(end(hits)[k] + start(hits)[k])/2,
                                   stringsAsFactors = FALSE)
              df.grn.snp <- rbind(df.grn.snp, newrow)   
            }
          }
        }
      }
    }
    names(df.grn.snp) <- c("TF", "snp_sequence", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_in_SNP")
    saveRDS(df.grn.snp, paste("Datasets/SNP/tmp/df.grn.snp_",i,".rds", sep = ""))
  }
  
  # combine
  df.grn.snp <- data.frame(TF = character(), snp_sequence = character(), TF.Motif = character(), 
                           bind_score = numeric(), cmp_bind_score = numeric(), BS_pos_in_SNP = numeric(),
                           stringsAsFactors = FALSE)      
  for(i in 1:length(lst.pwm.motif)){ 
    print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
    df.grn.snp <- rbind(df.grn.snp, readRDS(paste("Datasets/SNP/tmp/df.grn.snp_",i,".rds", sep = "")))
  }
  names(df.grn.snp) <- c("TF", "snp_sequence", "TF.Motif", "bind_score", "cmp_bind_score", "BS_pos_in_SNP")

  df.grn.snp <- df.grn.snp[,c(-3,-5)]
  df.grn.snp <- unique(df.grn.snp)
  
  #df.CNS_2014["BS_pos_to_TSS"] <-  apply(df.CNS_2014[,c('start_dist_to_TSS','end_dist_to_TSS')], 1, function(x) mean(x) )
  #df.CNS_2014 <- df.CNS_2014[, c(-1,-2,-3,-8,-9,-10)]
  
  names(df.founds)[1] <- "Target"
  names(df.founds)[3] <- "snp_sequence"
  
  df.SNP_grn <- merge(df.grn.snp, df.founds, by = "snp_sequence")
  df.SNP_grn["BS_pos_to_TSS"] <- df.SNP_grn$start_dist_to_TSS - df.SNP_grn$BS_pos_in_SNP
  
  #filter  
  foreach(i = 1:length(vec.snp.sequences)) %dopar% { 
 # for(i in 1:length(vec.snp.sequences)){
    df.SNP_grn.filtered <- data.frame(snp_sequence =character(), TF = character(), bind_score = numeric(), BS_pos_in_SNP = numeric(), Target = character(),
                                      snp.motif = character(),start_dist_to_TSS = numeric(), mean_dist_to_TSS = numeric(),BS_pos_to_TSS = numeric(),
                                      stringsAsFactors = FALSE) 
    print(paste(i, "of",length(vec.snp.sequences)))
    df.sset <- subset(df.SNP_grn, df.SNP_grn$snp_sequence == vec.snp.sequences[i])
    if(nrow(df.sset) > 0){
      tfs <- unique(df.sset$TF)
      for(t.1 in 1:length(tfs)){
        df.sset.2 <- subset(df.sset, df.sset$TF == tfs[t.1])
        targs <- unique(df.sset.2$Target)  
        for(t.2 in 1:length(targs)){
          df.sset.3 <- subset(df.sset.2, df.sset.2$Target == targs[t.2])
          vec.pos <- unique(df.sset.3$start_dist_to_TSS)
          for(p in 1:length(vec.pos)){
            df.sset.4 <- subset(df.sset.3, df.sset.3$start_dist_to_TSS == vec.pos[p])
            df.sset.4 <- subset(df.sset.4, df.sset.4$bind_score == max(unique(df.sset.4$bind_score)))
            df.SNP_grn.filtered <- rbind(df.SNP_grn.filtered, df.sset.4[1,])
          }
        }
      }
    }
    saveRDS(df.SNP_grn.filtered, paste("Datasets/SNP/tmp/df.grn.snp_filtered_",i,".rds", sep = ""))
  }
  
  df.SNP_grn.filtered <- data.frame(snp_sequence =character(), TF = character(), bind_score = numeric(), BS_pos_in_SNP = numeric(), Target = character(),
                                   snp.motif = character(),start_dist_to_TSS = numeric(), mean_dist_to_TSS = numeric(),BS_pos_to_TSS = numeric(),
                                   stringsAsFactors = FALSE) 
  for(i in 1:length(vec.snp.sequences)){ 
    print(paste(i, "of", length(vec.snp.sequences)))
    df.SNP_grn.filtered <- rbind(df.SNP_grn.filtered, readRDS(paste("Datasets/SNP/tmp/df.grn.snp_filtered_",i,".rds", sep = "")))
  }
  
  names(df.SNP_grn.filtered) <- c("snp_sequence", "TF", "bind_score", "BS_pos_in_SNP", "Target", "snp.motif", "start_dist_to_TSS", "mean_dist_to_TSS", "BS_pos_to_TSS")
  saveRDS(df.SNP_grn.filtered, paste("Datasets/SNP/df.grn.snp_filtered.rds", sep = ""))
  #saveRDS(df.SNP_grn, paste("Datasets/SNP/df.grn.snp.rds", sep = ""))
}

perform_cMonkey_TFBS_mapping <- function(v.th.eVal_pwm_match = 0.05, n.top.matching.motifs = 15, nr.cores = 15, load_from_file = FALSE){
  
  # recompute the biclustering (later... )
  library(cMonkey)
  library(MotIV)
  library(foreach)
  library(doMC)
  registerDoMC(nr.cores)
  
  n.env.cMonkey <- 5
  lst.env.cMonkey <- vector(mode = "list", length = n.env.cMonkey)
  lst.env.cMonkey[[1]] <- readRDS("Datasets/BiClust_GRNS/env.cMonkey.development.rds")
  lst.env.cMonkey[[2]] <- readRDS("Datasets/BiClust_GRNS/env.cMonkey.perturbation_CoSpecificity.rds")
  lst.env.cMonkey[[3]] <- readRDS("Datasets/BiClust_GRNS/env.cMonkey.perturbation_abiotic.rds")
  lst.env.cMonkey[[4]] <- readRDS("Datasets/BiClust_GRNS/env.cMonkey.perturbation_biotic.rds")
  lst.env.cMonkey[[5]] <- readRDS("Datasets/BiClust_GRNS/env.cMonkey.perturbation_hormone.rds")
  
  lst.motifs <- get_cell_and_pnas_Paper_PWMs()
  lst.pwm.motif <- lst.motifs[[1]]
  df.motifs <- lst.motifs[[2]]
  lst.pwm.motif <- lapply(lst.pwm.motif, function(x) round(100 * x)) 
  pwm.scores <- generateDBScores(inputDB=lst.pwm.motif,cc="PCC",align="SWU",nRand=1000)
  
  
  foreach(e = 1:n.env.cMonkey) %dopar% { 
  #for(e in 1:n.env.cMonkey){
    env.cMonkey <- lst.env.cMonkey[[e]]    
    # mat.GE <- env.cMonkey$ratios$ratios
    # genes <- rownames(mat.GE)
    #conditions <- colnames(mat.GE)
    df.motifs.cMonkey <- env.cMonkey$cluster.summary()
    n.cmonkey_motifs <- 2
    for(c in 1:env.cMonkey$k.clust){
      df.grn.biclust.TFB_map <- data.frame(TF = character(), Target = character(), TF.Motif = character(), Rank.TF.Motif = numeric(), biclust_id = numeric(),
                                           nr.biclust_motif = numeric(), e.val.Motif = numeric(), p.val_per_gene = numeric(), e.val.TFB_similarity = numeric(),
                                           BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE)    
      print(paste(c,"of",env.cMonkey$k.clust))
      cluster.genes <- env.cMonkey$get.rows(c) # get genes
      cluster.conditions <- env.cMonkey$get.cols(c) # get conditions
      if(!is.null(env.cMonkey$meme.scores[[1]][[c]]$meme.out[[1]]$pssm)){
        # PWM per cluster
        lst.pwms <- vector(mode = "list", length = n.cmonkey_motifs)
        for(i in 1:length(lst.pwms)){
          lst.pwms[[i]] <- t(env.cMonkey$meme.scores[[1]][[c]]$meme.out[[i]]$pssm)
          rownames(lst.pwms[[i]]) <- c("A","C","G","T")
          colnames(lst.pwms[[i]]) <- seq(1:ncol(lst.pwms[[i]])) 
        }
        lst.pwms <- lapply(lst.pwms,trimPWMedge, threshold=1)
        names(lst.pwms) <- c("pwm1", "pwm2")
        # genes, distances to TSS, gene specific p values, motif evalue
        # top 2 matching matrices
        res <- motifMatch(inputPWM = lst.pwms, database = lst.pwm.motif, DBscores = pwm.scores, top = n.top.matching.motifs) 
        for(i in 1:length(lst.pwms)){
          n.matchin.TF_pwms <- length(res@bestMatch[[i]]@aligns) 
          if(n.matchin.TF_pwms > 0){
            e.val <- env.cMonkey$meme.scores[[1]][[c]]$meme.out[[i]]$e.value
            genes <- as.character(env.cMonkey$meme.scores[[1]][[c]]$meme.out[[i]]$posns$gene)
            #as.numeric(env.cMonkey$meme.scores[[1]][[c]]$meme.out[[i]]$posns$strand)
            dist_to_TSS_per_gene <- as.numeric(env.cMonkey$meme.scores[[1]][[c]]$meme.out[[i]]$posns$start) # distance to start
            p.val_per_gene <- as.numeric(env.cMonkey$meme.scores[[1]][[c]]$meme.out[[i]]$posns$p.value) # gene based p value
            for(k in 1:n.matchin.TF_pwms){
              motif.id <- res@bestMatch[[i]]@aligns[[k]]@TF@name
              e.val.similarity <- res@bestMatch[[i]]@aligns[[k]]@evalue
              idx <- match(motif.id, df.motifs$Motif.ID)
              for(l in 1:length(genes)){
                newrow <- data.frame(TF =  df.motifs$TF.locus[idx], 
                                     Target = genes[l],
                                     TF.Motif =  motif.id,
                                     Rank.TF.Motif = k,
                                     biclust_id = c, 
                                     nr.biclust_motif = i,
                                     e.val.TF.Motif = e.val,
                                     p.val_per_gene =  p.val_per_gene[l],
                                     e.val.TFB_similarity = e.val.similarity,
                                     BS_pos_to_TSS =  dist_to_TSS_per_gene[l],
                                     stringsAsFactors = FALSE)
                df.grn.biclust.TFB_map <- rbind(df.grn.biclust.TFB_map, newrow)  
              }
            }   
          }
        } 
      }
      names(df.grn.biclust.TFB_map) <- c("TF", "Target", "TF.Motif", "Rank.TF.Motif", "biclust_id", "nr.biclust_motif", "e.val.TF.Motif",
                                         "p.val_per_gene", "e.val.TFB_similarity", "BS_pos_to_TSS")
      
      df.grn.biclust.TFB_map.filtered <- data.frame(TF = character(), Target = character(), TF.Motif = character(), Rank.TF.Motif = numeric(), biclust_id = numeric(),
                                                    nr.biclust_motif = numeric(), e.val.Motif = numeric(), p.val_per_gene = numeric(), e.val.TFB_similarity = numeric(),
                                                    BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE) 
      
      if(nrow(df.grn.biclust.TFB_map) > 0){
        df.grn.biclust.TFB_map <- subset(df.grn.biclust.TFB_map, df.grn.biclust.TFB_map$e.val.TFB_similarity <= v.th.eVal_pwm_match)
        for(m in 1:n.cmonkey_motifs){
          df.sset <- subset(df.grn.biclust.TFB_map, df.grn.biclust.TFB_map$nr.biclust_motif == m)
          tfs <- unique(df.sset$TF)
          for(t in 1:length(tfs)){
            df.sset.2 <- subset(df.sset, df.sset$TF == tfs[t])
            targs <- unique(df.sset.2$Target)
            for(l in 1:length(targs)){
              df.sset.3 <- subset(df.sset.2, df.sset.2$Target == targs[l])
              df.sset.3 <- subset(df.sset.3, df.sset.3$e.val.TFB_similarity == min(unique(df.sset.3$e.val.TFB_similarity)))
              df.grn.biclust.TFB_map.filtered <- rbind(df.grn.biclust.TFB_map.filtered, df.sset.3[1,])
            }
          }
        }
        names(df.grn.biclust.TFB_map.filtered) <- c("TF", "Target", "TF.Motif", "Rank.TF.Motif", "biclust_id", "nr.biclust_motif", "e.val.TF.Motif",
                                                  "p.val_per_gene", "e.val.TFB_similarity", "BS_pos_to_TSS")
        saveRDS(df.grn.biclust.TFB_map.filtered, paste("Datasets/BiClust_GRNS/tmp/df.grn.biclust.TFB_map_",e,"_",c,".rds", sep = ""))
      }
    }
    # combine
    df.grn.biclust.TFB_map <- data.frame(TF = character(), Target = character(), TF.Motif = character(), Rank.TF.Motif = numeric(), biclust_id = numeric(),
                                         nr.biclust_motif = numeric(), e.val.Motif = numeric(), p.val_per_gene = numeric(), e.val.TFB_similarity = numeric(),
                                         BS_pos_to_TSS = numeric(), stringsAsFactors = FALSE)    
    for(c in 1:env.cMonkey$k.clust){
      print(paste(c,"of",env.cMonkey$k.clust))
      df.grn.biclust.TFB_map <- rbind(df.grn.biclust.TFB_map, readRDS(paste("Datasets/BiClust_GRNS/tmp/df.grn.biclust.TFB_map_",e,"_",c,".rds", sep = "")))
    }
     names(df.grn.biclust.TFB_map) <- c("TF", "Target", "TF.Motif", "Rank.TF.Motif", "biclust_id", "nr.biclust_motif", "e.val.TF.Motif",
                                        "p.val_per_gene", "e.val.TFB_similarity", "BS_pos_to_TSS")
     saveRDS(df.grn.biclust.TFB_map, paste("Datasets/BiClust_GRNS/df.grn.biclust.TFB_map_",e,".rds", sep = ""))
  }

  
#   n.env.cMonkey <- 5
#   vec.global_condition <- c("development, abiotic, biotic, hormone, tissue_specifity")
#   lst.df.grn.biclust.TFB_map <- vector(mode = "list", length = n.env.cMonkey)
#   for(e in 1:n.env.cMonkey){
#     readRDS(paste("Datasets/BiClust_GRNS/df.grn.biclust.TFB_map_",e,".rds", sep = ""))
#     
#   }
}


## REMAINING STUFF

compute_transitive_closure <- function(){
  
  vec.genes <- unique(c(vec.tfs, vec.tgs))
  
  #df.grn.biclust.TFB_map <- readRDS(paste("Datasets/BiClust_GRNS/df.grn.biclust.TFB_map_",e,".rds", sep = ""))
  #df.SNP_grn.filtered <- readRDS(paste("Datasets/SNP/df.grn.snp_filtered.rds", sep = ""))
  df.CNS_grn.seqMotifs <- readRDS(paste("../Datasets/CNS_GRNS/df.CNS_grn2014_SeqMotifs_filtered.rds", sep = ""))
  df.CNS_grn.filtered <- readRDS(paste("../Datasets/CNS_GRNS/df.CNS_grn2014_filtered.rds", sep = ""))
  df.grn.all_motifs_filtered <- readRDS(paste("../Datasets/novelTFmotifs/df.grn.all_motifs_filtered", sep = ""))
  df.grn.seqMotifs <- readRDS("../Datasets/novelTFmotifs/df.grn.seqMotifs.rds")
  df.grn.CNS <- readRDS(paste("../Datasets/GRMs/df.grn.CNS.rds", sep = ""))
  
  
  ### Conserved Binding 
  #   df.CNS_grn.filtered <- readRDS(paste("Datasets/CNS_GRNS/df.CNS_grn2014_filtered.rds", sep = ""))
  #   df.CNS_grn.filtered <- df.CNS_grn.filtered[,c(2,10,3)]
  #   df.CNS_grn.filtered <- subset(df.CNS_grn.filtered, df.CNS_grn.filtered$bind_score >= 6)
  #   df.CNS_grn.filtered <- df.CNS_grn.filtered[,c(1:2)]  
  #   
  #   df.CNS_grn.seqMotifs <- df.CNS_grn.seqMotifs[,c(1,3)]
  #   df.CNS_grn.seqMotifs <- unique(df.CNS_grn.seqMotifs)
  #   
  #   df.grn.CNS <- rbind(df.grn.CNS, df.CNS_grn.seqMotifs)
  #   df.grn.CNS <- rbind(df.grn.CNS, df.CNS_grn.filtered)
  #   df.grn.CNS <- unique(df.grn.CNS)
  #   
  #   
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
  #   saveRDS(mat.grn.tfb.CNS, "Datasets/workspace/mat.grn.tfb.CNS.all_8.rds") 
  #mat.grn.tfb.CNS <- readRDS("Datasets/workspace/mat.grn.tfb.CNS.rds") 
  #mat.grn.tfb <- readRDS("Datasets/workspace/mat.grn.tfb.probs_above_95.rds") 
  #mat.grn.tfb.CNS[(!rownames(mat.grn.tfb.CNS) %in% vec.TF_with_TFB), ] <- 0.5
  
  #### Binding probability 
  df.grn.all_motifs_filtered <- df.grn.all_motifs_filtered[,c(1,2,3)]
  
  df.CNS_grn.seqMotifs["bind_score"] <- 8
  df.CNS_grn.seqMotifs <- df.CNS_grn.seqMotifs[,c(1,3,7)]
  df.CNS_grn.seqMotifs <- unique(df.CNS_grn.seqMotifs)
  
  df.grn.seqMotifs["bind_score"] <- 7
  df.grn.seqMotifs <- df.grn.seqMotifs[,c(1,2,4)]
  df.grn.seqMotifs <- unique(df.grn.seqMotifs)
  
  df.grn.CNS["bind_score"] <- 8
  df.grn.seqMotifs <- rbind(df.grn.seqMotifs, df.grn.CNS)
  df.grn.seqMotifs <- rbind(df.grn.seqMotifs, df.CNS_grn.seqMotifs)
  
  df.grn.all_motifs_filtered <- rbind(df.grn.all_motifs_filtered, df.grn.seqMotifs)
  
  x <- df.grn.all_motifs_filtered$bind_score
  P = ecdf(x) 
  
  df.grn.TFB.combined <- df.grn.all_motifs_filtered
  df.grn.TFB.combined["p.TFB"] <- P(df.grn.all_motifs_filtered$bind_score)
  df.grn.TFB.combined <- subset(df.grn.TFB.combined, df.grn.TFB.combined$TF %in% vec.tfs & df.grn.TFB.combined$Target %in% vec.genes)
  df.grn.TFB.combined <- subset(df.grn.TFB.combined, df.grn.TFB.combined$p.TFB >= P(7))
  df.grn.TFB.combined <- unique(df.grn.TFB.combined[,c(1,2)])
  
  vec.regs.tfb <- unique(df.grn.TFB.combined$TF)
  vec.tgs.tfb <- unique(df.grn.TFB.combined$Target)
  
  mat.grn.tfb <- matrix(0, nrow = length(vec.genes), ncol = length(vec.genes),
                        dimnames = list(vec.genes, vec.genes))
  
  for(i in 1:nrow(df.grn.TFB.combined)){
    print(paste(i, "of", nrow(df.grn.TFB.combined)))
    tf <- df.grn.TFB.combined$TF[i]
    tg <- df.grn.TFB.combined$Target[i]
    mat.grn.tfb[tf, tg] <- 1 #df.grn.TFB.combined$p.TFB[i]  
  }
  saveRDS(mat.grn.tfb, "Datasets/workspace/mat.grn.tfb.probs_above_95.rds") 
  
  #mat.grn.tfb <- readRDS("Datasets/workspace/mat.grn.tfb.probs_above_95.rds") 
  #mat.grn.tfb[(!rownames(mat.grn.tfb) %in% vec.TF_with_TFB), ] <- 0.5
  
  #mat.grn.tfb <- readRDS("Datasets/workspace/mat.grn.tfb.probs_above_95.rds") 
  #mat.grn.tfb[(!rownames(mat.grn.tfb) %in% vec.TF_with_TFB), ] <- 0.5
  
  
  ###
  df.grn.CNS <- subset(df.grn.CNS, df.grn.CNS$TF %in% rownames(mat.grn) & df.grn.CNS$Target %in% vec.genes)
  mat.grn.tfb.CNS <- matrix(0, nrow = length(vec.genes), ncol = length(vec.genes),
                            dimnames = list(vec.genes, vec.genes))
  
  for(i in 1:nrow(df.grn.CNS)){
    print(paste(i, "of", nrow(df.grn.CNS)))
    tf <- df.grn.CNS$TF[i]
    tg <- df.grn.CNS$Target[i]
    #if(mat.grn.tfb[tf, tg] < df.grn.TFB.combined$p.TFB[i]){
    mat.grn.tfb.CNS[tf, tg] <- 1 #df.grn.TFB.combined$p.TFB[i]  
    #}
  }
  #saveRDS(mat.grn.tfb.CNS, "Datasets/workspace/mat.grn.tfb.CNS.rds") 
  mat.grn.tfb.CNS <- readRDS("Datasets/workspace/mat.grn.tfb.CNS.rds") 
  #mat.grn.tfb <- readRDS("Datasets/workspace/mat.grn.tfb.probs_above_95.rds") 
  mat.grn.tfb.CNS[(!rownames(mat.grn.tfb.CNS) %in% vec.TF_with_TFB), ] <- 0.5
  
  
  
  
  #df.grn.all_motifs_filtered <- rbind(df.grn.all_motifs_filtered, df.CNS_grn.seqMotifs)
  
  ####
  v.th.binding_cutoff <- 4.9
  
  
  
  #df.CNS_grn.filtered <- subset(df.CNS_grn.filtered, df.CNS_grn.filtered$bind_score >= v.th.binding_cutoff)
  #df.grn.all_motifs_filtered <- subset(df.grn.all_motifs_filtered, df.grn.all_motifs_filtered$bind_score >= v.th.binding_cutoff)
  
  df.grn.TFB.combined <- df.grn.CNS
  df.grn.TFB.combined <- rbind(df.grn.TFB.combined, df.CNS_grn.seqMotifs[,c(1,3)])
  df.grn.TFB.combined <- rbind(df.grn.TFB.combined, df.CNS_grn.filtered[,c(2,10)])
  df.grn.TFB.combined <- rbind(df.grn.TFB.combined, df.grn.all_motifs_filtered[,c(1,2)])
  df.grn.TFB.combined <- rbind(df.grn.TFB.combined, df.grn.seqMotifs[,c(1,2)])
  
  df.grn.TFB.combined <- unique(df.grn.TFB.combined)
  
  # 4.9 cutoff
  
  # + TFB (!!1/0/0.5 interconnected!!!)
  #df.grn.TFB.combined <- readRDS(paste("Datasets/GRMs/df.grn.TFB.combined.rds", sep = ""))
  
  vec.genes <- c(rownames(mat.grn), colnames(mat.grn))
  
  df.grn.TFB.combined <- subset(df.grn.TFB.combined, df.grn.TFB.combined$TF %in% rownames(mat.grn) & df.grn.TFB.combined$Target %in% vec.genes)
  vec.regs.tfb <- unique(df.grn.TFB.combined$TF)
  vec.tgs.tfb <- unique(df.grn.TFB.combined$Target)
  
  mat.grn.tfb <- matrix(0, nrow = length(vec.genes), ncol = length(vec.genes),
                        dimnames = list(vec.genes, vec.genes))
  
  for(i in 1:nrow(df.grn.TFB.combined)){
    print(paste(i, "of", nrow(df.grn.TFB.combined)))
    tf <- df.grn.TFB.combined$TF[i]
    tg <- df.grn.TFB.combined$Target[i]
    mat.grn.tfb[tf, tg] <- 1
  }
  #saveRDS(mat.grn.tfb, "Datasets/workspace/mat.grn.tfb.rds") 
  #saveRDS(mat.grn.tfb, "Datasets/workspace/mat.grn.tfb.1.rds") 
  #saveRDS(mat.grn.tfb, "Datasets/workspace/mat.grn.tfb.2.rds") 
  
  
  ###
  
  
  mat.grn.TF_TF <- mat.grn.tfb[rownames(mat.grn), rownames(mat.grn)]  
  mat.grn.TF_TF.transitive <- mat.grn.TF_TF
  for(y in 1:ncol(mat.grn.TF_TF)){
    print(y)
    tfs.y <- names(which(mat.grn.TF_TF[,y] == 1))
    if(length(tfs.y) > 0){
      for(x in 1:length(tfs.y)){
        idx.tgs <- as.numeric(which(mat.grn.TF_TF[tfs.y[x], ] == 1))
        mat.grn.TF_TF.transitive[tfs.y[x], idx.tgs] <- 1
      }  
    }
  }
  
  
  
  #mat.grn <- mat.pred.grn
  # 1 step transitive 
  mat.grn.TF_ME <- mat.grn.tfb[rownames(mat.grn), colnames(mat.grn)]  
  mat.grn.TF_ME.transitive <- mat.grn.TF_ME
  for(y in 1:ncol(mat.grn.TF_TF.transitive)){
    print(y)
    tfs.y <- names(which(mat.grn.TF_TF.transitive[,y] == 1))
    for(x in 1:length(tfs.y)){
      idx.tgs <- as.numeric(which(mat.grn.TF_ME[tfs.y[x], ] == 1))
      mat.grn.TF_ME.transitive[tfs.y[x], idx.tgs] <- 1
    }
  }
  #mat.grn.TF_ME.transitive[(!rownames(mat.grn.TF_ME.transitive) %in% vec.TF_with_TFB), ] <- 0.5
  
  
  saveRDS(mat.grn.TF_ME.transitive, "Datasets/workspace/mat.grn.TF_ME.transitive.2.rds") 
  #saveRDS(mat.grn.TF_ME.transitive, "Datasets/workspace/mat.grn.TF_ME.transitive.1.80.rds") 
  #saveRDS(mat.grn.TF_ME.transitive, "Datasets/workspace/mat.grn.TF_ME.transitive.1.95.rds") 
  
  
  
  ######
  vec.TF_with_TFB <- unique(df.grn.TFB.combined$TF)  
  mat.grn.TF_TF_with_TFB <- mat.grn.tfb[vec.TF_with_TFB, vec.TF_with_TFB]
  
  library(agop)
  strt<-Sys.time()
  mat.grn.tfb.TC <- closure_transitive(mat.grn.TF_TF_with_TFB)
  print(Sys.time()-strt)
  #saveRDS(mat.grn.tfb.TC, "Datasets/workspace/mat.grn.tfb.TC.rds") 
  #saveRDS(mat.grn.tfb.TC, "Datasets/workspace/mat.grn.tfb.TC.1.rds") 
  
  # mark a connection 
  #mat.grn.TF_TF <- mat.grn.tfb[rownames(mat.grn), rownames(mat.grn)]
  mat.grn.TF_ME <- mat.grn.tfb[rownames(mat.grn), colnames(mat.grn)]  
  mat.grn.TF_ME.transitive <- mat.grn.TF_ME
  for(y in 1:ncol(mat.grn.tfb.TC)){
    print(y)
    tfs.y <- names(which(mat.grn.tfb.TC[,y] == 1))
    for(x in 1:length(tfs.y)){
      idx.tgs <- as.numeric(which(mat.grn.TF_ME[tfs.y[x], ] == 1))
      mat.grn.TF_ME.transitive[tfs.y[x], idx.tgs] <- 1
    }
  }
  mat.grn.TF_ME.transitive[(!rownames(mat.grn.TF_ME.transitive) %in% vec.TF_with_TFB), ] <- 0.5
  
  # saveRDS(mat.grn.TF_ME.transitive, "Datasets/workspace/mat.grn.TF_ME.transitive.rds") 
  
}


# ## gene filter before proceed
# df.CNS_2014_2000kb <- readRDS(paste("Datasets/novelTFmotifs/df.CNS_2014_2000kb.rds", sep = ""))
# df.CNS_2014_2000kb_sset <- subset(df.CNS_2014_2000kb, df.CNS_2014_2000kb$gene %in% vec.genes.considered)
# saveRDS(df.CNS_2014_2000kb_sset, "Datasets/novelTFmotifs/df.CNS_2014_2000kb_sset.rds")

# 
# prepare_NCS2014_GRN <- function(load_from_file = FALSE){
# 	
#   if(!load_from_file){
#     
#     #source("http://bioconductor.org/biocLite.R")
#     #biocLite("PWMEnrich")
#     library(PWMEnrich)  
#     
#     print("compute scored CNS2014-GRN")
#     
#     df.CNS_2014 <- extract_CNS2014()
#     lst.motifs <- get_cell_and_pnas_Paper_PWMs()
#     lst.pwm.motif <- lst.motifs[[1]]
#     df.motifs <- lst.motifs[[2]]
#     
#     lst.log_pwm.motif <- vector(mode = "list", length = length(lst.pwm.motif))
#     for(i in 1:length(lst.pwm.motif)){
#       log_pwm <- new("PWM", id = names(lst.pwm.motif)[i], name = df.motifs$TF.locus[i], pfm = lst.pwm.motif[[i]], prior.params = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25), pwm = log2(lst.pwm.motif[[i]]))  
#       lst.log_pwm.motif[[i]] <-  log_pwm
#     }
#     
#     n.matches <- nrow(df.CNS_2014) * 3
#     df.grn.CNS2014 <- data.frame(TF = character(n.matches), Target = character(n.matches), TF.Motif = character(n.matches),cons_score = numeric(n.matches), align.score = numeric(n.matches),
#                                  bind_score = numeric(n.matches), p.val = numeric(n.matches), rank = numeric(n.matches), mean_dist_to_TSS = numeric(n.matches), stringsAsFactors = FALSE)    
#     names(df.grn.CNS2014) <- c("TF", "Target", "TF.Motif", "cons_score", "align_score", "bind_score", "p.val", "rank", "mean_dist_to_TSS")
#         
#     idx <- 1
#     for(j in 1:nrow(df.CNS_2014)){
#       print(paste(j, "of", nrow(df.CNS_2014)))
#       sequence <- DNAString(as.character(df.CNS_2014$cns_sequence[j]))     
#       lst.log_pwm.motif.sset <- lst.log_pwm.motif
#       for(l in 1:length(lst.log_pwm.motif.sset)){
#         if(length(sequence) < ncol(lst.log_pwm.motif.sset[[l]]$pwm)){
#           lst.log_pwm.motif.sset[[l]] <- character(0)
#         } 
#       }
#       lst.log_pwm.motif.sset <- Filter(length, lst.log_pwm.motif.sset)
#       if(length(lst.log_pwm.motif.sset) > 0){
#         res = motifEnrichment(sequence, lst.log_pwm.motif.sset, score = "affinity")
#         report = sequenceReport(res, 1)
#      
#         for(t in 1:3){
#           df.grn.CNS2014$TF[idx] <- report$target[t] 
#           df.grn.CNS2014$Target[idx] <- df.CNS_2014$gene[j]
#           df.grn.CNS2014$TF.Motif[idx] <- report$id[t] 
#           df.grn.CNS2014$cons_score[idx] <- df.CNS_2014$cons_score[j]
#           df.grn.CNS2014$align_score[idx] <- df.CNS_2014$val[j]
#           df.grn.CNS2014$bind_score[idx] <- report$raw.score[t] 
#           df.grn.CNS2014$rank[idx] <- report$rank[t]
#           df.grn.CNS2014$p.val[idx] <- report$p.value[t] 
#           df.grn.CNS2014$mean_dist_to_TSS[idx] <- mean(c(df.CNS_2014$start_dist_to_TSS[j], df.CNS_2014$end_dist_to_TSS[j]))
#           idx <- idx + 1
#         }   
#       }else{
#         idx <- idx + 3
#       }
#     }
#     df.grn.CNS2014 <- subset(df.grn.CNS2014, df.grn.CNS2014$TF != "")
#     saveRDS(df.grn.CNS2014, paste("Datasets/CNS_GRNS/df.cns2014_grn.rds", sep = ""))
#       
#    }
#     
#     
#     
#     ### --...
#     
#     print("compute scored CNS2014-GRN")
#     
#     df.CNS_2014 <- extract_CNS2014()
#    
#     system.time(for(i in 1:1){#length(lst.pwm.motif)){ 
#       print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
#         df.grn.CNS2014 <- data.frame(TF = character(n.matches), Target = character(n.matches), TF.Motif = character(n.matches), 
#                                      cons_score = numeric(n.matches), align.score = numeric(n.matches),
#                                      bind_score = numeric(n.matches), cmp_bind_score = numeric(n.matches),
#                                      mean_dist_to_TSS = numeric(n.matches), stringsAsFactors = FALSE)    
#         names(df.grn.CNS2014) <- c("TF", "Target", "TF.Motif", "cons_score", "align_score", "bind_score", "cmp_bind_score", "mean_dist_to_TSS")
#       
#       
#         for(j in 1:nrow(df.CNS_2014)){
#           print(paste(j, "of", nrow(df.CNS_2014)))
#           hits <- matchPWM(lst.pwm.motif[[i]],  as.character(df.CNS_2014$cns_sequence[j]), with.score = TRUE)
#           nhits <- length(hits)
#           
#           
#           #test.seqs <- DNAStringSet(unique(df.CNS_2014$cns_sequence)) 
#           #system.time(pwm.hits <- sapply(test.seqs, function(pseq) matchPWM(lst.pwm.motif[[i]], pseq, min.score="90%")))
#           #scores = motifScores(sequence, test, raw.scores=TRUE)
#           #head(scores[[1]])
#           #res = motifEnrichment(sequence, lst.log_pwm.motif, score = "affinity")
#           #report = sequenceReport(res, 1)
#           
#           if(nhits >= 1){
#             cmp_bind_score <- min(mcols(hits)$score / maxScore(lst.pwm.motif[[i]])) # should be >= 0.8
#             motif.score <- mcols(hits)$score
#             if(cmp_bind_score >= 0.8){
#               newrow <- data.frame(TF =  as.character(df.motifs$TF.locus[i]), Target = as.character(df.CNS_2014$gene[j]), 
#                                    TF.Motif =  lst.log_pwm.motif[[i]]$id,
#                                    cons_score = as.numeric(df.CNS_2014$cons_score[j]), 
#                                    align_score <- df.CNS_2014$val[j], 
#                                    bind_score = motif.score, 
#                                    cmp_bind_score = cmp_bind_score,
#                                    mean_dist_to_TSS = mean(c(df.CNS_2014$start_dist_to_TSS[j], df.CNS_2014$end_dist_to_TSS[j])), 
#                                    stringsAsFactors = FALSE)
#               df.grn.CNS2014 <- rbind(df.grn.CNS2014, newrow)   
#             }
#           }
#         }
#         names(df.grn.CNS2014) <- c("TF", "Target", "cons_score", "bind_score", "compared_bind_score", "start_dist", "short_dist")
#         saveRDS(df.grn.CNS2014, paste("Datasets/CNS_GRNS/tmp/df.cns2014_grn_",i,".rds", sep = ""))
#       })
#       # combine
#       df.grn.CNS2014 <- data.frame(TF = character(), Target = character(), cons_score = numeric(), bind_score = numeric(), cmp_bind_score = numeric(), 
#                                    start_dist = numeric(), short_dist = numeric(), stringsAsFactors = FALSE)    
#       for(i in 1:length(lst.pwm.motif)){ 
#         df.grn.CNS2014 <- rbind(df.grn.CNS2014, readRDS(paste("Datasets/CNS_GRNS/tmp/df.cns2014_grn_",i,".rds", sep = "")))
#       }
#       names(df.grn.CNS2014) <- c("TF", "Target", "cons_score", "bind_score", "compared_bind_score", "start_dist", "short_dist")
#       saveRDS(df.grn.CNS2014, paste("Datasets/CNS_GRNS/df.cns2014_grn.rds", sep = ""))
#   #}else{
#   #  df.grn.CNS2014 <- readRDS(paste("Datasets/CNS_GRNS/df.cns2014_grn.rds", sep = ""))
#   #}
#   return(df.grn.CNS2014)
# }
# 
# 
# 
# 
# 
# 
# extract_CNS2012 <- function(load_from_file = TRUE){
#   
#   if(!load_from_file){
#     
#     source("http://bioconductor.org/biocLite.R")
#     #   biocLite("ChIPpeakAnno")
#     #   biocLite("biomaRt")
#     #   biocLite("Biostrings")
#     #   install.packages("VennDiagram")
#     library(ChIPpeakAnno)
#     library(biomaRt)
#     library(Biostrings)
#     #biocLite("BSgenome.Athaliana.TAIR.TAIR9")
#     library(BSgenome.Athaliana.TAIR.TAIR9)
#     genome <- BSgenome.Athaliana.TAIR.TAIR9
#     
#     ensmart = useMart('ENSEMBL_MART_PLANT', "athaliana_eg_gene")
#     annotatedData = getAnnotation(ensmart, featureType = "TSS")
#     
#     df.CNS_2012 <- read.table("Datasets/novelTFmotifs/CNS_paper_2012/CNS_set.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
#     df.CNS_2012 <- subset(df.CNS_2012, df.CNS_2012$Species.Name == "Arabidopsis thaliana")
#     df.CNS_2012["start_dist_to_TSS"] <- NA
#     df.CNS_2012["end_dist_to_TSS"] <- NA
#     for(i in 1:nrow(df.CNS_2012)){
#       print(paste(i , "of", nrow(df.CNS_2012)))
#       gene <- as.character(df.CNS_2012$Gene.Identifier[i])
#       df.gene.info <- as.data.frame(annotatedData[gene,])
#       cns.region <- as.integer(unlist(strsplit(as.character(df.CNS_2012$Co.ordinates.of.CNS..Upstream.of.Gene.Identifier.TSS[i]), ":")))
#       
#       if(df.gene.info$strand == 1){
#         tss <- as.numeric(df.gene.info$start)     
#       }else{
#         tss <- as.numeric(df.gene.info$end)
#       }
#       df.CNS_2012$start_dist_to_TSS[i] <- abs(cns.region[1] - tss)
#       df.CNS_2012$end_dist_to_TSS[i] <- abs(cns.region[2] - tss)  
#     }
#     #df.CNS_2012 <- subset(df.CNS_2012, df.CNS_2012$end_dist_to_TSS < 2000)
#     saveRDS(df.CNS_2012, "Datasets/novelTFmotifs/df.CNS_2012.rds")  
#   }else{
#     df.CNS_2012 <- readRDS("Datasets/novelTFmotifs/df.CNS_2012.rds")  
#   }
# }
# 
# 
# prepare_NCS2012_GRN <- function(load_from_file = FALSE){
#   
#   df.CNS_2012 <- extract_CNS2012()
#   
#   #### --- 
#   
#   lst.motifs <- get_cell_and_pnas_Paper_PWMs()
#   lst.pwm.motif <- lst.motifs[[1]]
#   df.motifs <- lst.motifs[[2]]
#   
#   biocLite("PWMEnrich")
#   library(PWMEnrich)
#   
#   
#   library(PWMEnrich.Dmelanogaster.background)
#   data(PWMLogn.dm3.MotifDb.Dmel)
#   
#   print("compute scored CNS2012-GRN")
#   #source("http://bioconductor.org/biocLite.R")
#   # biocLite("ChIPpeakAnno")
#   # biocLite("biomaRt")
#   #biocLite("Biostrings", lib = "~/MyRLibs")
#   # install.packages("VennDiagram")
#   #library(ChIPpeakAnno)
#   #library(biomaRt)
#   library(Biostrings)
#   
#   
#   if(!load_from_file){
#     
#     system.time( # begin time measure
#       for(i in 1:length(lst.pwm.motif)){ 
#         print(paste("TF Motif ", i, "of", length(lst.pwm.motif)))
#         df.grn.CNS2012 <- data.frame(TF = character(), Target = character(), cons_score = numeric(), bind_score = numeric(), cmp_bind_score = numeric(), stringsAsFactors = FALSE)    
#         names(df.grn.CNS2012) <- c("TF", "Target", "cons_score", "bind_score", "compared_bind_score")
#         for(j in 1:nrow(df.CNS_2012)){
#           
#           
#           ids = c("bcd", "gt_FlyReg_FBgn0001150", "Kr")
#           sel.pwms = PWMLogn.dm3.MotifDb.Dmel$pwms[ids]
#           scores = motifScores(sequence, sel.pwms, raw.scores=TRUE)
#           head(scores[[1]])
#           
#           res = motifEnrichment(sequence, sel.pwms)
#           report = sequenceReport(res, 1)
#           
#           
#           hits <- matchPWM(lst.pwm.motif[[i]],  as.character(df.CNS_2012$Sequence[j]), with.score = TRUE)
#           nhits <- length(hits)
#           if(nhits >= 1){
#             cmp_bind_score <- min(mcols(hits)$score / maxScore(lst.pwm.motif[[i]])) # should be >= 0.8
#             motif.score <- mcols(hits)$score
#             if(cmp_bind_score >= 0.8){
#               newrow <- data.frame(TF =  as.character(df.motifs$TF.locus[i]), Target = as.character(df.CNS_2012$Gene.Identifier[j]), 
#                                    cons_score = as.numeric(df.CNS_2012$Conservation.Score[j]), bind_score = motif.score, 
#                                    cmp_bind_score = cmp_bind_score, stringsAsFactors = FALSE)
#               df.grn.CNS2012 <- rbind(df.grn.CNS2012, newrow)
#             }
#           }
#         }
#         names(df.grn.CNS2012) <- c("TF", "Target", "cons_score", "bind_score", "compared_bind_score")
#         saveRDS(df.grn.CNS2012, paste("Datasets/CNS_GRNS/tmp/df.cns2012_grn_",i,".rds", sep = ""))
#       }
#     ) # end time measure
#     
#     # combine
#     df.grn.CNS2012 <- data.frame(TF = character(), Target = character(), cons_score = numeric(), bind_score = numeric(), cmp_bind_score = numeric(), stringsAsFactors = FALSE)    
#     for(i in 1:length(lst.pwm.motif)){ 
#       df.grn.CNS2012 <- rbind(df.grn.CNS2012, readRDS(paste("Datasets/CNS_GRNS/tmp/df.cns2012_grn_",i,".rds", sep = "")))
#     }
#     names(df.grn.CNS2012) <- c("TF", "Target", "cons_score", "bind_score", "compared_bind_score")
#   }
#   saveRDS(df.grn.CNS2012, paste("Datasets/CNS_GRNS/df.cns2012_grn.rds", sep = ""))
# }
# 
# 
# 
# 
# 
# compute_TFBinding_Matrix_complete <- function( note ="", load_from_file = TRUE){
# 	
# 	tfbs <- read.csv("Datasets/TFBS/TFBS_conservation.csv", header = TRUE, sep = ",")
# 	TFs <- unique(as.character(tfbs$Transcription.factor))
# 	targets <- unique(as.character(tfbs$Target.gene))
# 	
# 	if(!load_from_file){
# 		tfbs <- read.csv("Datasets/TFBS/TFBS_conservation.csv", header = TRUE, sep = ",")
# 		
# 		mat.tfb <- matrix(0, nrow = length(TFs), ncol = length(targets))
# 		rownames(mat.tfb) <- TFs
# 		colnames(mat.tfb)  <- targets
# 		
# 		mat.tfb.cons <- matrix(0, nrow = length(TFs), ncol = length(targets))
# 		rownames(mat.tfb.cons) <- TFs
# 		colnames(mat.tfb.cons)  <- targets
# 		
# 		n <- nrow(tfbs)
# 		print("compute transcription factor binding matrix - conservation")
# 		for(i in 1:n){		
# 			print(paste(i,"of",n))
# 			reg <- as.character(tfbs$Transcription.factor[i])
# 			targ <- as.character(tfbs$Target.gene[i])
# 			w.cons <- as.numeric(tfbs$Species.conservation[i])
# 			
# 			mat.tfb.cons[reg, targ] <- w.cons
# 			mat.tfb[reg, targ] <- 1
# 				
# 		}
# 		saveRDS(mat.tfb, file = paste("Datasets/r_objects/mat.tfb.",note,"1.rds", sep =""))
# 		saveRDS(mat.tfb.cons, file = paste("Datasets/r_objects/mat.tfb.",note,"2.rds", sep =""))
# 	}else{
# 		mat.tfb <- readRDS(file = paste("Datasets/r_objects/mat.tfb.",note,"1.rds", sep =""))
# 		mat.tfb.cons <- readRDS(file = paste("Datasets/r_objects/mat.tfb.",note,"2.rds", sep =""))
# 	}
# 	
# 	return(list(mat.tfb, mat.tfb.cons))	
# }
# 
# 
# 
# #if(source == "pcc"){
# #	tfbs <- read.csv("Datasets/TFBS/TFBS_PCC.csv", header = TRUE, sep = ",")
# #}
# 
# #' Compute transcription factor binding adjacency matrix (based on motif conservation evidence)
# #'
# #' This function 
# #' 
# #' @param writeResults - write results to file (default = TRUE)
# #' @param ... expressions evaluated in the context of \code{df} and 
# #'   then fed to \code{\link{order}}
# #' @keywords manip
# #' @export
# #' @examples
# #' lst.mat.tfb <- compute_TFB_matrix_conservation(vec.genes)
# #' mat.tfb       <- lst.mat.tfb[[1]]
# #' mat.tfb.cons  <- lst.mat.tfb[[2]]
# compute_TFB_matrix_conservation <- function(vec.rel.genes, note ="", load_from_file = TRUE){
# 
# 	if(!load_from_file){
# 		tfbs <- read.csv("Datasets/TFBS/TFBS_conservation.csv", header = TRUE, sep = ",")
# 		
# 		mat.tfb <- matrix(0, nrow = length(vec.rel.genes), ncol = length(vec.rel.genes))
# 		rownames(mat.tfb) <- vec.rel.genes 
# 		colnames(mat.tfb) <- vec.rel.genes 
# 		
# 		mat.tfb.cons <- matrix(0, nrow = length(vec.rel.genes), ncol = length(vec.rel.genes))
# 		rownames(mat.tfb.cons) <- vec.rel.genes 
# 		colnames(mat.tfb.cons) <- vec.rel.genes 
# 		
# 		n <- nrow(tfbs)
# 		print("compute transcription factor binding matrix - conservation")
# 		for(i in 1:n){		
# 				print(paste(i,"of",n))
# 				reg <- as.character(tfbs$Transcription.factor[i])
# 				targ <- as.character(tfbs$Target.gene[i])
# 				w.cons <- as.numeric(tfbs$Species.conservation[i])
# 				if(is.element(reg, vec.rel.genes) && is.element(targ, vec.rel.genes)){
# 					mat.tfb.cons[reg, targ] <- w.cons
# 					mat.tfb[reg, targ] <- 1
# 				}	
# 			}
# 			saveRDS(mat.tfb, file = paste("Datasets/r_objects/mat.tfb.",note,"1.rds", sep =""))
# 			saveRDS(mat.tfb.cons, file = paste("Datasets/r_objects/mat.tfb.",note,"2.rds", sep =""))
# 		}else{
# 			mat.tfb <- readRDS(file = paste("Datasets/r_objects/mat.tfb.",note,"1.rds", sep =""))
# 			mat.tfb.cons <- readRDS(file = paste("Datasets/r_objects/mat.tfb.",note,"2.rds", sep =""))
# 		}
# 		
# 	return(list(mat.tfb, mat.tfb.cons))	
# }
# 
# 
# #' Compute transcription factor binding adjacency matrix (based on AGRIS evidence)
# #'
# #' This function 
# #' 
# #' @param writeResults - write results to file (default = TRUE)
# #' @param ... expressions evaluated in the context of \code{df} and 
# #'   then fed to \code{\link{order}}
# #' @keywords manip
# #' @export
# #' @examples
# #' mat.tfb <- compute_TFB_matrix_agris(vec.genes)
# compute_TFB_matrix_agris <- function(vec.rel.genes, df.geneloci_TF, writeResults = TRUE,  reset = FALSE){
# 	if(reset){
# 		bindingSites <- loadTFBS("Datasets/BindingSite.tbl")
# 	
# 		bindingSites$TF.Family[bindingSites$TF.Family == "MYB-RELATED"] <- "MYB_related"
# 		bindingSites$TF.Family[bindingSites$TF.Family == "BHLH"] <- "bHLH"
# 		bindingSites$TF.Family[bindingSites$TF.Family == "E2F-DP"] <- "E2F/DP"
# 		bindingSites$TF.Family[bindingSites$TF.Family == "HB"] <- "HB-PHD"
# 		bindingSites$TF.Family[bindingSites$TF.Family == "AP2-EREBP"] <- "AP2"
# 		bindingSites$TF.Family[bindingSites$TF.Family == "BZIP"] <- "bZIP"
# 		write.csv(bindingSites,"Datasets/TFBS/BindingSite_renamed.csv", row.names=FALSE)
# 	}else{
# 		bindingSites <- read.csv("Datasets/TFBS/BindingSite_renamed.csv")
# 	}
# 	
# 	mat.TFBS <- matrix(0, nrow = length(vec.rel.genes), ncol = length(vec.rel.genes))
# 	rownames(mat.TFBS) <- vec.rel.genes
# 	colnames(mat.TFBS) <- vec.rel.genes
# 	print("compute transcription factor binding matrix - AGRIS")
# 	for(i in 1:length(unique(bindingSites$TF.Family))){	
# 		tf.fam <- as.character(unique(bindingSites$TF.Family)[i])
# 		subset.TFBS <- subset(bindingSites, TF.Family == tf.fam)		
# 		TFs.binding <- subset(df.geneloci_TF, df.geneloci_TF$Family == tf.fam)
# 		idx.i <- match(TFs.binding$Locus, vec.rel.genes)
# 		idx.j <- match(subset.TFBS$Promoter.Locus, vec.rel.genes)
# 		idx.i  <- idx.i[!is.na(idx.i)]
# 		idx.j  <- idx.j[!is.na(idx.j)]
# 		for(m in 1:length(idx.i)){
# 			for(n in 1:length(idx.j)){
# 				mat.TFBS[idx.i[m], idx.j[n]] <- 1		
# 			}
# 		}	
# 	}
# 	if(writeResults){
# 		saveRDS(mat.TFBS, file = "Datasets/r_objects/mat.tfb.agris.rds")
# 	}
# 	return(mat.TFBS)
# }
# 
# 
# #' Extract metabolic enzyme encoding genes by metabolic pathway
# #'
# #' This function identifies all metabolic enzyme encoding genes belonging to each pathway in AraCyc
# #' returns lists of pwys and their metabolic enzyme encoding genes, a list of pw.names and pw.ids
# #' @param writeResults - write results to file (default = TRUE)
# #' @param ... expressions evaluated in the context of \code{df} and 
# #'   then fed to \code{\link{order}}
# #' @keywords manip
# #' @export
# #' @examples
# #' lst.pw.MEs <- extract_ME_by_PWY()
# #' pw.enzymes <- lst.pw.MEs[[1]]
# #' pw.ids     <- lst.pw.MEs[[2]]
# #' pw.names   <- lst.pw.MEs[[3]]
# compute_tfbs_between_modules <- function(lst.geneModules, vec.rel.genes, mat.TFBS, idx.dset = 1, idx.partition = 6, note ="", load_from_file = TRUE){ # idx_per_cluster
# 	
# 	if(!load_from_file){
# 		n.modules <- length(lst.geneModules)
# 		idx <- vector(mode = "list", n.modules)
# 		for(i in 1:n.modules){
# 		  id <- as.numeric(match(lst.geneModules[[i]], vec.rel.genes))
# 		  idx[[i]] <- id[!is.na(id)] 
# 		}
# 		mat.TFBS.weights.norm <- matrix(0, nrow = n.modules, ncol = n.modules)
# 		rownames(mat.TFBS.weights.norm) <- seq(1:n.modules)
# 		colnames(mat.TFBS.weights.norm) <- seq(1:n.modules)
# 	  
# 		mat.TFBS.weights <- matrix(0, nrow = n.modules, ncol = n.modules)
# 		rownames(mat.TFBS.weights) <- seq(1:n.modules)
# 		colnames(mat.TFBS.weights) <- seq(1:n.modules)		
# 		print("compute Transcription Factor Binding between Modules")
# 		for(i in 1:n.modules){
# 			for(j in 1:n.modules){	
# 			  	idx_per_cluster <- expand.grid(idx[[i]], idx[[j]])
# 				v.TFBS <- numeric()	
# 				v.TFBS <- apply( idx_per_cluster , 1 , function(x) mat.TFBS[x[1] ,x[2]])
# 				n.genes <- (length(lst.geneModules[[i]]) * length(lst.geneModules[[j]]))
# 			  	mat.TFBS.weights.norm[i,j] <- (sum(v.TFBS) / n.genes)
# 				mat.TFBS.weights[i,j] <- sum(v.TFBS)
# 			}
# 		}
# 		saveRDS(mat.TFBS.weights.norm, file = paste("Datasets/r_objects/mat.tfb.modules.norm.",note,idx.dset,"_",idx.partition,".rds", sep =""))
# 		saveRDS(mat.TFBS.weights, file = paste("Datasets/r_objects/mat.tfb.modules.",note, idx.dset,"_",idx.partition,".rds", sep =""))
# 		return(list(mat.TFBS.weights.norm, mat.TFBS.weights))	
# 	}else{
# 		mat.TFBS.weights.norm <- readRDS(file = paste("Datasets/r_objects/mat.tfb.modules.norm.",note, idx.dset,"_",idx.partition,".rds", sep =""))
# 		mat.TFBS.weights <- readRDS(file = paste("Datasets/r_objects/mat.tfb.modules.",note, idx.dset,"_",idx.partition,".rds", sep =""))
# 		return(list(mat.TFBS.weights.norm, mat.TFBS.weights))
# 	}
# }
# 
# 
# 
# #### FINAL EVALUATION METHOD?
# compute_TFB_PWY_enrichment <- function(vec.rel.genes, lst.pw.enzymes, vec.pw.ids, vec.pw.names, lst.pw.reactions.id, lst.pw.reactions.enzymes, writeResults = TRUE){
#   
#   pw.names <- readRDS(file = "Datasets/r_objects/pw_names.rds")
#   pw.ids <- readRDS(file = "Datasets/r_objects/pw_ids.rds")
#   pw.enzymes <- readRDS(file = "Datasets/r_objects/pw_enzymes.rds")
#   pw.reactions.id  <- readRDS(file = "Datasets/r_objects/pw_reactions.id.rds")
#   pw.reactions.enzymes <- readRDS(file = "Datasets/r_objects/pw_reactions.enzymes.rds")
#   
#   tfbs <- read.csv("Datasets/TFBS/TFBS_conservation.csv", header = TRUE, sep = ",")
#   tfs <- as.character(unique(tfbs$Transcription.factor))
#   
#   tf_pwy_regs <- data.frame(TF = character(), PWY.ID = character(), REG.RATIO = numeric(),  NR.REGS = numeric(), PWY.SIZE = numeric(), PWY.NAME = character())
#   
#   n.genome <- 27000
# 
#   mat.pValues <- matrix(1, nrow=length(tfs), ncol= length(pw.ids ))
#   colnames(mat.pValues) <-pw.ids 
#   rownames(mat.pValues) <- tfs
#   
#   for(i in 1:length(tfs)){
# 	  
# 	    print(paste(i,"of", length(tfs)))
# 	    tfbs.sset <- subset(tfbs, tfbs$Transcription.factor == tfs[i])
# 	    
# 		targ.genes <- as.character(tfbs.sset$Target.gene)
# 		n.targ.genes <- length(targ.genes)
# 	
# 	    for(j in 1:length(pw.reactions.id)){
# 	    	
# 			n.tg.enzymes <- numeric() 
# 			n.pwy.knowns <- 0
# 			
# 			for(r in 1:length(pw.reactions.enzymes[[j]])){
# 				
# 				if(pw.reactions.enzymes[[j]][[r]][1] != "unknown"){
# 					n.pwy.knowns <- n.pwy.knowns + 1
# 					n.reg.reaction.genes <- length(intersect(pw.reactions.enzymes[[j]][[r]], targ.genes))
# 					n.tg.enzymes <- c(n.tg.enzymes, n.reg.reaction.genes)	
# 					
# 					# for every individual reaction per pathways
# 					#n.reaction.genes <- length(pw.reactions.enzymes[[j]][[r]])
# 					#n.non_reg.pwy.genes <- n.reaction.genes  - n.reg.reaction.genes
# 					#n.non_targ.genes <- n.genome - n.targ.genes
# 					#counts = (matrix(data = c(n.reg.reaction.genes, n.non_reg.pwy.genes, n.targ.genes, n.non_targ.genes), nrow = 2))
# 					#p.ind.reaction <- fisher.test(counts)	
# 				}
# 			}
# 	
# 			n.pwy.reactions <- n.pwy.knowns
# 			n.reg.pwy.reactions <- length(which(n.tg.enzymes > 0))
# 			n.non_reg.pwy.reactions <- n.pwy.reactions - n.reg.pwy.reactions
# 			n.non_targ.genes <- n.genome - n.targ.genes
# 			
# 			counts = (matrix(data = c(n.reg.pwy.reactions, n.non_reg.pwy.reactions, n.targ.genes, n.non_targ.genes), nrow = 2))
# 			mat.pValues[i,j] <- fisher.test(counts)$p.value
# 			}
# 		}
# 	
# 		vec.tmp <- as.numeric(mat.pValues)
# 		p.values.adjusted <- p.adjust(vec.tmp, method= "BH")
# 		mat.pValues<- matrix(p.values.adjusted, nrow=length(tfs), byrow = T)
# 		colnames(mat.pValues) <-pw.ids 
# 		rownames(mat.pValues) <- tfs
# 		
# 		
# 		# COMPLETE 			pwy genes 		non pwy genes
# 		# reg. genes	
# 		# non reg genes
# 		#n.pwy.genes <- 6 #length(which(pw.enzymes[[i]] != "unknown"))
# 		#n.reg.pwy.genes <- 4 # length(intersect(pw.enzymes[[i]], targ.genes))
# 		#n.non_reg.pwy.genes <- n.pwy.genes - n.reg.pwy.genes
# 		#n.non_targ.genes <- n.genome - n.targ.genes
# 		#counts = (matrix(data = c(n.reg.pwy.genes, n.non_reg.pwy.genes, n.targ.genes, n.non_targ.genes), nrow = 2))
# 		#p.complete <- fisher.test(counts)
# 		
# 
# 		pwy.targets <- intersect(as.character(tfbs.sset$Target.gene), as.character(pw.enzymes[[j]]))		
# 		size.pwy <- length(as.character(pw.enzymes[[j]]))
# 		nr.reg <- length(pwy.targets)
# 		ratio.reg <- nr.reg / size.pwy
# 		
# 		if(nr.reg >= 1){
# 			tf_pwy_regs <- rbind(tf_pwy_regs, data.frame(TF = as.character(tfs[i]), PWY.ID = as.character(pw.ids[[j]]), REG.RATIO = ratio.reg, NR.REGS = nr.reg , PWY.SIZE = size.pwy, PWY.NAME = as.character(pw.names[[j]])))
# 		}
# 	
#   return(tf_pwy_regs)
#   write.table(tf_pwy_regs, "/Users/michaelbanf/Documents/postdoctoral_work/programs/pwy_enrichment.xls", row.names = FALSE, sep = "\