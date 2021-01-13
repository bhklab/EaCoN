######################################
#### Load and Clean the ASCAT RDS ####
loadBestFitRDS <- function(sample, gamma, segmenter){
  RDS.file <- file.path(sample, toupper(segmenter), 'ASCN', 
                        paste0('gamma', format(gamma, nsmall=2)),
                        paste(sample, 'ASCN', toupper(segmenter), 'RDS', sep="."))
  rds <- readRDS(RDS.file)
  return(rds)
}

loadL2R <- function(sample, segmenter){
  RDS.file <- file.path(sample, toupper(segmenter), 'L2R', 
                        paste(sample, 'SEG', toupper(segmenter), 'RDS', sep="."))
  rds <- readRDS(RDS.file)
  l2r <- cleanL2Rdat(rds)
  return(l2r)
}

cleanL2Rdat <- function(data){
  l2r.segments <- data$cbs$nocut
  l2r.segments$Value <- l2r.segments$Log2Ratio
  
  g.cut <- data$meta$eacon[["L2R-segments-gain-cutoff"]]
  l.cut <- data$meta$eacon[["L2R-segments-loss-cutoff"]]
  
  gain.idx <- which(l2r.segments$Value > g.cut)
  loss.idx <- which(l2r.segments$Value < l.cut)
  normal.idx <- which(l2r.segments$Value >= l.cut & l2r.segments$Value <= g.cut)
  
  l2r.seg.obj <- list(pos = l2r.segments, 
                      idx = list(gain = gain.idx, 
                                 loss = loss.idx, 
                                 normal = normal.idx), 
                      cutval = c(l.cut, g.cut))
  seg.col <- list(gain = "blue", outscale.gain = "midnightblue", loss = "red", outscale.loss = "darkred", normal = "black")
  
  genome.pkg = data$meta$basic$genome.pkg
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)
  cs <- EaCoN:::chromobjector(BSg.obj)
  
  l2r.chr <- as.integer(unname(unlist(cs$chrom2chr[as.character(data$data$SNPpos$chrs)])))
  
  l2r.value <- data.frame(Chr = l2r.chr,
                          Start = as.integer(data$data$SNPpos$pos),
                          End = as.integer(data$data$SNPpos$pos),
                          Value = data$data$Tumor_LogR_wins[,1],
                          stringsAsFactors = FALSE)
  samplename=data$meta$basic$samplename
  list("l2r.value"=l2r.value,
       "l2r.seg.obj"=l2r.seg.obj,
       "seg.col"=seg.col,
       "genome.pkg"=genome.pkg,
       "sample.name"=samplename)
}

cleanGR <- function(gr0){
  for(i in c(1:ncol(elementMetadata(gr0)))){
    icol <- elementMetadata(gr0)[,i]
    if(class(icol)=='numeric'){
      elementMetadata(gr0)[,i] <- round(icol, 3)
    }
  }
  gr0$TCN <- rowSums(as.matrix(elementMetadata(gr0)[,c('nMajor', 'nMinor')]))
  gr0$seg.mean <- round(log2(rowSums(as.matrix(elementMetadata(gr0)[,c('nAraw', 'nBraw')])) / 2),3)
  gr0$seg.mean[gr0$seg.mean < log2(1/50)] <- round(log2(1/50), 2)
  gr0
}

getGenes <- function(genome.build="hg19", make.into.gr=FALSE){
  switch(genome.build,
         hg18={
           suppressPackageStartupMessages(require(TxDb.Hsapiens.UCSC.hg18.knownGene))
           if(!exists("TxDb.Hsapiens.UCSC.hg18.knownGene")) stop("Requires TxDb.Hsapiens.UCSC.hg18.knownGene")
           package <- TxDb.Hsapiens.UCSC.hg18.knownGene
         },
         hg19={
           suppressPackageStartupMessages(require(TxDb.Hsapiens.UCSC.hg19.knownGene))
           if(!exists("TxDb.Hsapiens.UCSC.hg19.knownGene")) stop("Requires TxDb.Hsapiens.UCSC.hg19.knownGene")
           package <- TxDb.Hsapiens.UCSC.hg19.knownGene
         },
         hg38={
           suppressPackageStartupMessages(require(TxDb.Hsapiens.UCSC.hg38.knownGene))
           if(!exists("TxDb.Hsapiens.UCSC.hg38.knownGene")) stop("Requires TxDb.Hsapiens.UCSC.hg38.knownGene")
           package <- TxDb.Hsapiens.UCSC.hg38.knownGene
         },
         stop("genome must be 'hg19' or 'hg38'"))
  
  if(make.into.gr){
    genes0 <- genes(package)
    idx <- rep(seq_along(genes0), elementNROWS(genes0$gene_id))
    genes <- granges(genes0)[idx]
    genes$gene_id = unlist(genes0$gene_id)
    genes
  } else {
    list(txdb=package,
         txdb.genes=genes(package))
  }
}

getMapping <- function(in.col='ENTREZID', 
                       out.cols=c("SYMBOL", "ENSEMBL")){
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  gene.map <- as.data.frame(sapply(out.cols, function(oc){
    mapIds(org.Hs.eg.db, keys=keys(org.Hs.eg.db, in.col),
           keytype="ENTREZID", column=oc, multiVals = 'first')
  }))
  gene.map$ENTREZID <- keys(org.Hs.eg.db, in.col)
  # gene.map <- select(org.Hs.eg.db, keys=keys(org.Hs.eg.db, in.col), 
  #                    keytype="ENTREZID", columns=out.cols)
  
  gene.map[,1] <- as.character(gene.map[,1])
  gene.map[,2] <- as.character(gene.map[,2])
  gene.map[,3] <- as.character(gene.map[,3])
  gene.map
}

fitOptimalGamma <- function(fit.val, sample=NULL, gamma.meta=NULL, pancan.ploidy=NULL){
  # #### Test Data
  # fit.val <- all.fits[[2]]$fit
  # sample <- all.fits[[2]]$sample
  # meta <- ccle.meta[,c('SNP arrays', 'tcga_code')]
  # colnames(meta) <- c('sample', 'TCGA_code')
  # ####
  
  if(is.null(pancan.ploidy)){
    ## Load the pancan.ploidy dataset
    pancan.dir <- "/mnt/work1/users/pughlab/references/TCGA/TCGA_Pancan_ploidyseg/cleaned"
    pancan.data <- readRDS(file=file.path(pancan.dir, 'pancanPloidy.RDS'))
    pancan.ploidy <- pancan.data$breaks$ploidy
  }
  
  ### Validation checks
  if(is.null(gamma.meta)){
    warning("A meta file is needed to estimate Gamma, otherwise average pancan ploidy is used")
    meta.col.check <- FALSE
  } else {
    meta.col.check <- all(c("sample", "TCGA_code") %in% colnames(gamma.meta))
    tcga.code <- 'AVG'
  }
  
  if(meta.col.check){
    sample.row.idx <- match(sample, gamma.meta$sample)
    if(is.na(sample.row.idx)){
      warning(paste0(sample, " was not found in the metadata"))
      tcga.code <- 'AVG'
    } else {
      tcga.code <- gamma.meta[sample.row.idx,]$TCGA_code
      if(is.na(tcga.code)){
        warning(paste0(sample, " does not have an annotated TCGA onco-code"))
        tcga.code <- 'AVG'
      } 
    }
  } else {
    warning("Could not find columns 'sample' and 'TCGA_code' in metafile, using average pancan ploidy.")
    tcga.code <- 'AVG'
  }
  
  if(tcga.code == 'COAD/READ'){
    tcga.code <- 'COAD'
  }
  if(!any(grepl(tcga.code, colnames(pancan.ploidy)))){
    warning(paste0("Could not find any TCGA code to match: ", tcga.code))
    tcga.code <- 'AVG'
  }
  
  ### Calculate best fit:
  ploidy.prior <- pancan.ploidy[,c('breaks', tcga.code)]
  ## Smooth the ploidy curve using a loess fit
  ploidy.prior$loess <- suppressWarnings(
    loess(formula = paste(tcga.code, "breaks", sep = "~"),
          data = ploidy.prior,
          degree = 1, ncmax= 200,
          span = 0.05)$fitted
  )
  ploidy.prior$breaks <- round(ploidy.prior$breaks, 1)
  # plot(ploidy.prior[,c('breaks', tcga.code)], type='p')
  # lines(ploidy.prior[,c('breaks', 'loess')], col="blue")
  fit.val$ploidy.round <- round(fit.val$psi,1)
  fit.val.priors <- merge(fit.val, ploidy.prior, by.x='ploidy.round', by.y='breaks', all.x=T)
  
  ## Create a score: (Goodnes_of_Fit * ploidy_likelihood)
  fit.val.priors$score <- with(fit.val.priors, GoF * loess)
  fit.val.priors <- fit.val.priors[order(fit.val.priors$gamma),]
  
  ## Select the score with the largest gamma
  max.score.idx <- which.max(fit.val.priors$score)
  return(max.score.idx[length(max.score.idx)])
}

ASCAT.selectBestFit <- function(fit.val, gamma.method='GoF', ...){
  idx <- switch(gamma.method,
                GoF=which.max(fit.val$GoF),
                psi=which.min(fit.val$psi),
                score=fitOptimalGamma(fit.val, ...),
                dev=temp())
  gamma <- fit.val$gamma[idx]
  tmsg(paste0("Gamma: ", gamma, " [method=", gamma.method, "]"))
  
  return(gamma)
}

genWindowedBed <- function(bin.size=1000000, seq.style="UCSC"){
  suppressPackageStartupMessages(require(BSgenome.Hsapiens.UCSC.hg19))
  chrs <- seqlengths(Hsapiens)[paste0("chr", c(1:22,"X", "Y"))]
  
  ## Construct intervals across the genome of a certain bin size
  start.points <- seq(1, 500000000, by=bin.size)
  grl <- lapply(names(chrs), function(chr.id){
    chr <- chrs[chr.id]
    ir <- IRanges(start=start.points[start.points < chr], width=bin.size)
    end(ir[length(ir),]) <- chr
    gr <- GRanges(seqnames = chr.id, ir)
    gr
  })
  
  ## Assemble all GRanges and set seq level style
  grl <- as(grl, "GRangesList")
  suppressWarnings(seqlevelsStyle(grl) <- seq.style)
  gr <- unlist(grl)
  gr
}

flagMultiBins <- function(olaps){
  dup.idx <- which(duplicated(subjectHits(olaps), fromLast=TRUE))
  dup.idx <- c(dup.idx, which(duplicated(subjectHits(olaps), fromLast=FALSE)))
  dup.idx <- sort(dup.idx)
  return(dup.idx)
}

reduceMultiBins <- function(cnv, dup.idx, olaps, reduce='median'){
  dup.df <- as.data.frame(olaps[dup.idx,])
  dup.spl <- split(dup.df, dup.df$subjectHits)
  
  dup.em <- lapply(dup.spl, function(i,reduce='median') {
    em.mat <- as.matrix(mcols(cnv[i$queryHits,]))
    storage.mode(em.mat) <- "numeric"
    switch(reduce,
           mean=colMeans(em.mat),
           min=do.call(pmin, lapply(1:nrow(em.mat), function(j) em.mat[j,])),
           max=do.call(pmax, lapply(1:nrow(em.mat), function(j) em.mat[j,])),
           median=apply(em.mat, 2, median))
    
  }, reduce=reduce)
  dup.em <- as.data.frame(do.call(rbind, dup.em))
  dup.em$sample <- unique(cnv$sample)
  return(dup.em)
}

populateNewMcols <- function(ref, cnv, dup.em, olaps, dup.idx){
  em  <- matrix(nrow=length(ref), 
                ncol=(ncol(mcols(cnv)) + ncol(mcols(ref))), 
                dimnames = list(NULL,
                                c(colnames(mcols(ref)),
                                  colnames(mcols(cnv)))))
  dedup.olaps <- olaps[-dup.idx,]
  em <- as.data.frame(em)
  em[subjectHits(dedup.olaps),] <- as.data.frame(cbind(mcols(ref)[subjectHits(dedup.olaps),],
                                                       mcols(cnv)[queryHits(dedup.olaps),]))
  em[unique(subjectHits(olaps)[dup.idx]),] <- dup.em
  return(em)
}

segmentCNVs <- function(cnv, bed, reduce='mean', feature.id='bins', l2r.dat=NULL){
  olaps = findOverlaps(cnv, bed)
  # Flag BED bins that map to multiple CNVs
  dup.idx <- EaCoN:::flagMultiBins(olaps)
  # Use a summary metric (Default=mean) to reduce the CNV information that
  # spans multiple bed windows
  dup.em <- EaCoN:::reduceMultiBins(cnv, dup.idx, olaps, reduce=reduce)
  # Initialize a metadata matrix and populate it for the BED GRanges object
  em <- EaCoN:::populateNewMcols(bed, cnv, dup.em, olaps, dup.idx)
  
  if(!is.null(l2r.dat)){
    ## Developmental: Include the block running median
    require(dplyr)
    l2rraw.gr <- makeGRangesFromDataFrame(l2r.dat$l2r.value, keep.extra.columns = TRUE)
    l2r.gr <- makeGRangesFromDataFrame(l2r.dat$l2r.seg.obj$pos, keep.extra.columns = TRUE)[,'Log2Ratio']
    seqlevels(l2rraw.gr)  <- seqlevels(l2r.gr) <- c(1:22, "X", "Y")
    seqlevelsStyle(l2rraw.gr) <- seqlevelsStyle(l2r.gr) <- 'UCSC'

    ref <- bed
    mcols(ref) <- em
    
    ## Append Log2Ratio segmented values
    olaps = GenomicRanges::findOverlaps(l2r.gr, ref)
    dup.idx <- EaCoN:::flagMultiBins(olaps)
    dup.em <- EaCoN:::reduceMultiBins(l2r.gr, dup.idx, olaps, reduce=reduce)
    em <- EaCoN:::populateNewMcols(ref, l2r.gr, dup.em, olaps, dup.idx)
    em$Log2Ratio <- round(em$Log2Ratio, 3)
    
    ## Append Log2Ratio unsegmented values
    olaps <- findOverlaps(l2rraw.gr, bed)
    lrr.vals <- mcols(l2rraw.gr)[,1]
    lrr.meds <- sapply(split(lrr.vals, subjectHits(olaps)), median, na.rm=TRUE)
    em$L2Rraw <- NA
    em[unique(subjectHits(olaps)),]$L2Rraw <- round(lrr.meds,3)
  }
  
  # Append metadata and return
  bed$ID <- em$ID <- paste0(feature.id, "_", c(1:nrow(em)))
  
  return(list(seg=bed, genes=em))
}

.assignEntrezToSegment <- function(cnv0, anno){
  olaps = findOverlaps(cnv0, anno)
  mcols(olaps)$gene_id = anno$gene_id[subjectHits(olaps)]  # Fixed the code here
  cnv_factor = factor(queryHits(olaps), levels=seq_len(queryLength(olaps)))
  gene.id <- splitAsList(mcols(olaps)$gene_id, cnv_factor)
  return(gene.id)
}

.splitSegmentByGene <- function(cnv0, cols){
  seg.entrez <- apply(as.data.frame(mcols(cnv0)), 1, function(i){
    ids <- unlist(strsplit(x = as.character(unlist(i[['gene_id']])), split=","))
    segs <- do.call(rbind, replicate(length(ids), round(unlist(i[cols]),3), simplify = FALSE))
    
    as.data.frame(cbind(segs, 'ENTREZ'=ids))
  })
  seg.entrez <- do.call(rbind, seg.entrez)
  if(any(duplicated(seg.entrez$ENTREZ))) seg.entrez <- seg.entrez[-which(duplicated(seg.entrez$ENTREZ)),]
  return(seg.entrez)
}

annotateCNVs <- function(cnv, txdb, anno=NULL,
                         cols=c("seg.mean", "nA", "nB"),
                         l2r.dat=NULL){
  stopifnot(is(cnv, "GRanges"), is(txdb, "TxDb"))
  
  ## Assign EntrezID to each segment 
  if(is.null(anno)) anno = genes(txdb)
  cnv$gene_id = EaCoN:::.assignEntrezToSegment(cnv, anno)
  if(!is.null(l2r.dat)) l2r.dat$gene_id <-  EaCoN:::.assignEntrezToSegment(l2r.dat, anno)
  
  ## Split the segments by genes
  seg.entrez <- EaCoN:::.splitSegmentByGene(cnv, cols)
  if(!is.null(l2r.dat)) {
    l2r.entrez <- EaCoN:::.splitSegmentByGene(l2r.dat, cols='Log2Ratio')
    seg.l2r.entrez <- merge(seg.entrez, l2r.entrez, by='ENTREZ', all.x=TRUE)
    seg.entrez <- seg.l2r.entrez
  }
  
  ## Map ensembl and HUGO IDs to the ENTREZ ids
  seg.anno <- merge(seg.entrez, EaCoN:::getMapping(),
                    by.x="ENTREZ", by.y="ENTREZID", all.x=TRUE)
  if(any(duplicated(seg.anno$ENTREZ))) seg.anno <- seg.anno[-which(duplicated(seg.anno$ENTREZ)),]
  cols <- colnames(seg.entrez)[grep("ENTREZ", colnames(seg.entrez), invert = TRUE)]
  for(each.col in cols){
    seg.anno[,each.col] <- as.numeric(as.character(seg.anno[,each.col]))
  }
  
  list("seg"=cnv, "genes"=seg.anno)  
}

#' Annotate CNV
#' @description This function annotates the output of the ASCN ASCAT file. It scans
#' the output directory for the *.gammaEval.txt to select the best "Goodness of Fit"
#' score, then load in the corresponding RDS file.  Using this RDS file, it extracts
#' and calculates the absolute allele-specific copy-number (modalA, modalB, modalAB),
#' as well as the L2R copy-ratios associated with them (nAraw, nBraw, seg.mean). It 
#' then annotates each segment and attributes a value to each gene in the knownGenes 
#' data structure
#'
#' @examples
#'  annotateRDS(fit.val.df, 'SampleX', 'ASCAT', build='hg19', bin.size=50000)
#' 
#' @param fit.val A dataframe of the *.gammaEval.txt file
#' @param sample Sample name 
#' @param segmenter Segmenter used (only works with ASCAT)
#' @param build Genome build (only works with hg19) [Default=hg19]
#' @param bin.size Bin-size for feature creation [Default=50000]
#' @param ... 
#'
#' @return Returns a list containing 3 elements: 'seg', 'genes', and 'bins'
#' for annotated datastructure fo a seg file, genes per row, and bins per row
#' @export
annotateRDS <- function(fit.val, sample, segmenter, build='hg19', 
                        bin.size=50000, feature.set=c('bins'), load.l2r=TRUE, ...){
  print(paste0("Bin size: ", bin.size))
  ## Assemble the ASCAT Seg file into a CNV GRanges object
  # load.l2r=TRUE
  # bin.size=5000
  # sample <- all.fits[r['start']:r['end']][[1]]$sample
  # fit.val <- all.fits[r['start']:r['end']][[1]]$fit
  # gamma <- EaCoN:::ASCAT.selectBestFit(fit.val, sample=sample, gamma.method='score',
  #                              gamma.meta=meta, pancan.ploidy=pancan.ploidy)
  # segmenter <- 'ASCAT'
  # build <- 'hg19'
  # tmsg=EaCoN:::tmsg
  # EaCoN:::
  gamma <- EaCoN:::ASCAT.selectBestFit(fit.val, sample=sample, ...)
  my.data <- EaCoN:::loadBestFitRDS(sample, gamma, segmenter)
  l2r.data <- if(load.l2r) EaCoN:::loadL2R(sample, segmenter) else NULL
  genes <- EaCoN:::getGenes(build)
  # EaCoN:::EaCoN.l2rplot.geno(l2r = l2r.data$l2r.value,
  #                            seg = l2r.data$l2r.seg.obj, 
  #                            seg.col = l2r.data$seg.col,
  #                            seg.type = "block", seg.normal = TRUE, 
  #                            genome.pkg = l2r.data$genome.pkg,
  #                            title = paste0(l2r.data$samplename, " L2R"),
  #                            ylim = c(-1.5,1.5))

  
  tmsg(paste0("Annotating sample: ", sample, " [gamma:", gamma, "]..."))
  cnv <- makeGRangesFromDataFrame(my.data$segments_raw, keep.extra.columns=TRUE, 
                                  start.field='startpos', end.field='endpos')
  cnv <- EaCoN:::cleanGR(cnv)

  ## Annotate the CNVs based on:
  cl.anno <- list()
  
  # Genes
  cols <- c('nMajor', 'nMinor', 'nAraw', 'nBraw', 'TCN', 'seg.mean')
  if(!is.null(l2r.data)){
    l2r.gr <- GenomicRanges::makeGRangesFromDataFrame(l2r.data$l2r.seg.obj$pos, keep.extra.columns = TRUE)
    seqlevels(l2r.gr) <- c(1:22, "X", "Y")
    seqlevelsStyle(l2r.gr) <- 'UCSC'
  }
  cl.anno[['genes']] <- suppressMessages(EaCoN:::annotateCNVs(cnv, genes$txdb, 
                                                      anno=genes$txdb.genes, cols=cols,
                                                      l2r.dat=l2r.gr))
  cl.anno$genes$genes <- cl.anno$genes$genes[-which(is.na(cl.anno$genes$genes$SYMBOL)),]
  
  # Raw seg
  cl.anno[['seg']] <- as.data.frame(cnv)

  feature.anno <- lapply(feature.set, function(fset){
    switch(fset,
           bins={
             tmsg("Assembling 'Bin' features...")
             windowed.bed <- EaCoN:::genWindowedBed(bin.size=bin.size)
             EaCoN:::segmentCNVs(cnv, bed=windowed.bed, reduce='median', feature.id = fset, l2r.dat=l2r.data)
           },
           tads={
             tmsg("Assembling 'TAD' features...")
             data(consensusTAD)
             segmentCNVs(cnv, tad.gr, reduce='median', feature.id = fset)
           },
           cres={
             tmsg("Assembling 'CRE' features...")
             data(geneCRE)
             segmentCNVs(cnv, cre.gr, reduce='median', feature.id = fset)
           })
  })
  names(feature.anno) <- feature.set
  
  cl.anno <- append(cl.anno, feature.anno)
  
  return(cl.anno)
}

#' Batch wrapper for annotateRDS()
#'
#' @param all.fits Sample named-list of all samples containing two elements:
#'  fit[data.frame]: data.frame of the *.gammaEval.txt file
#'  sample[character]: sample name
#' @param segmenter Segmenter used (E.g. ASCAT) (Only ASCAT works currently)
#' @param nthread Max number of threads to use [Default=1]
#' @param cluster.type Cluster type [Default=PSOCK]
#' @param ... 
#'
#' @examples
#'     gr.cnv <- annotateRDS.Batch(all.fits, toupper(segmenter), nthread=3)
#' 
#' @return Returns annotateRDS() objects
#' @export
annotateRDS.Batch <- function(all.fits, segmenter, nthread = 1, 
                              cluster.type = "PSOCK", bin.size=50000, ...){
  if (length(all.fits) < nthread) nthread <- length(all.fits)
  `%dopar%` <- foreach::"%dopar%"
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  grcnv.batch <- foreach::foreach(r = all.fits, 
                                  .inorder = TRUE, 
                                  .errorhandling = "stop",
                                  .export = ls(globalenv())) %dopar% {
                                    annotateRDS(r$fit, r$sample, segmenter, build='hg19', bin.size=bin.size, ...)
                                  }
  parallel::stopCluster(cl)
  names(grcnv.batch) <- sapply(all.fits, function(x) x$sample)
  return(grcnv.batch)
}

############################################
#### Building cBioportal Output Objects ####
.appendToCbioSeg <- function(cbio.path, cbio.file, seg, raw=FALSE, overwrite=NULL){
  # if(raw){
  #   cbio.file['file'] <- paste0("RAW", cbio.file['file'])
  # }
  exist.seg <- read.table(file.path(cbio.path, cbio.file['file']), sep="\t", header=T,
                          stringsAsFactors = F, check.names = F, fill=F)
  exist.spl <- split(exist.seg, f=exist.seg$ID)
  if(!is.null(overwrite)) {
    ov.idx <- sapply(overwrite, function(id) grep(id, names(exist.spl)))
    tmsg(paste0("Overwriting samples: ", paste(names(exist.spl)[ov.idx], collapse=",")))
    exist.spl[ov.idx] <- NULL
  }
  
  new.spl <- split(seg, f=seg$ID)
  new.ids <- which(!names(new.spl) %in% names(exist.spl))
  
  if(length(new.ids) > 0){
    tmsg(paste0("New samples being added to cBio Seg file : ", 
                paste(names(new.spl)[new.ids], collapse=",")))
    seg <- do.call(rbind, append(new.spl[new.ids], exist.spl))
  } else {
    tmsg("No new samples to add.  If you want to replace an existing sample, please specify 
         using overwrite=c('SampleA', 'SampleB')")
  }
  return(seg)
}

.appendToCbioMat <- function(cbio.path, cbio.file, cnv.mat, overwrite=NULL){
  exist.cna <- read.table(file.path(cbio.path, cbio.file['file']), sep="\t", header=T,
                          stringsAsFactors = F, check.names = F, fill=F)
  if(!is.null(overwrite)) {
    ov.idx <- sapply(overwrite, function(id) grep(id, colnames(exist.cna)))
    tmsg(paste0("Overwriting samples: ", paste(colnames(exist.cna)[ov.idx], collapse=",")))
    exist.cna <- exist.cna[,-ov.idx]
  }
  new.cols <- which(!colnames(cnv.mat) %in% colnames(exist.cna))
  
  if(length(new.cols) > 0){
    tmsg(paste0("New samples being added to cBio CNA matrices: ", 
                paste(colnames(cnv.mat)[new.cols], collapse=",")))
    cnv.mat <- suppressWarnings(merge(exist.cna, cnv.mat[,c(1,2, new.cols)], 
                                      by=c('Hugo_Symbol', 'Entrez_Gene_Id'), all.x=TRUE))
  } else {
    tmsg("No new samples to add.  If you want to replace an existing sample, please specify 
         using overwrite=c('SampleA', 'SampleB')")
  }
  return(cnv.mat)
}

#' cBio-formatted table builder 
#' @description Generates cBioportal-style formatted tsv's and seg files as per
#' details from https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats
#' 
#' @param gr.cnv List output from annotateRDS.Batch() function
#' @param cbio.path Relative path to the cBio output path
#' @param pattern Regex pattern for existing linear_CNA and _CNA files [Default="_CNA"]
#' @param cbio.cna.file Names of existing data_CNA.txt files [Defaul=NULL]
#' @param cbio.linear.file Names of existing data_linear_CNA.txt files [Defaul=NULL]
#' @param cbio.seg.file Names of existing data_cna_hg19.seg files [Defaul=NULL]
#' @param amp.thresh Absolute CN to start calling AMP from GAINs [Default=5]
#' @param add.on.to.existing Builds on existing cBio data found [Default=TRUE]
#' @param overwrite Sample ID to overwrite in cBio output if already existing [Default=NULL]
#' @param ... 
#'
#' @return NULL
#' @export
#'
#' @examples
#'     buildCbioOut(gr.cnv, cbio.path="./out/cBio", overwrite=c('YT_4941', '5637_3858'))
buildCbioOut <- function(gr.cnv, cbio.path="./out/cBio", pattern="_CNA", 
                         cbio.cna.file=NULL, cbio.linear.file=NULL, cbio.seg.file=NULL, cbio.RAWseg.file=NULL,
                         amp.thresh=5, add.on.to.existing=TRUE, ...){
  .checkFile <- function(cbio.path, file.id, pat){
    idx <- grep(list.files(cbio.path), pattern=pat, perl=TRUE)[1]
    if(!is.na(idx)){
      cbio.file <- list.files(cbio.path)[idx]
      exists.stat <- TRUE
    } else {
      cbio.file <- file.id
      exists.stat <- FALSE
    }
    c('file'=cbio.file, 'exists'=exists.stat)
  }
  
  .adjustCnaMat <- function(cnv.mat, mat.type, amp.thresh=NULL, ord=NULL){
    tcn.mat <- cnv.mat[,-c(1,2),drop=FALSE]
    if(mat.type=='CNA'){
      ## Convert Total CN to the -2, -1, 0, 1, 2, standards
      tcn.mat[tcn.mat >= amp.thresh] <- amp.thresh
      tcn.mat <- tcn.mat - 2
      for(i in c(1:(amp.thresh-3))){ tcn.mat[tcn.mat == i] <- 1 }
      tcn.mat[tcn.mat == (amp.thresh-2)] <- 2
    } else if(mat.type=='linear') {
      tcn.mat <- round(tcn.mat, 2)
    }
    
    ## Recombine the CNV mat
    colnames(tcn.mat) <- names(gr.cnv)
    colnames(cnv.mat)[c(1,2)] <- c('Hugo_Symbol', 'Entrez_Gene_Id')
    cnv.mat <- cbind(cnv.mat[,c(1:2)], tcn.mat)
    
    ## Set the order
    if(is.null(ord)){
      cnv.mat <- cnv.mat[order(cnv.mat$Hugo_Symbol),]
    } else {
      cnv.mat <- cnv.mat[match(ord, cnv.mat$Hugo_Symbol),]
      na.idx <- apply(cnv.mat, 1, function(x) all(is.na(x)))
      if(any(na.idx)) cnv.mat <- cnv.mat[-which(na.idx),]
    }
    return(cnv.mat)
  }
  
  tmsg(paste0("Building a cBioportal Object..."))
  ## Locating existing cBioportal Objects
  suppressWarnings(dir.create(cbio.path, recursive = TRUE))
  if(is.null(cbio.cna.file)){
    cbio.cna.file <- .checkFile(cbio.path, 'data_CNA.txt',   pat=paste0("(?<!linear)", pattern))
  }
  if(is.null(cbio.linear.file)){
    cbio.linear.file <- .checkFile(cbio.path, 'data_linear_CNA.txt',   pat=paste0("linear", pattern))
  }
  if(is.null(cbio.seg.file)){
    cbio.seg.file <- .checkFile(cbio.path, 'data_cna_hg19.seg',   pat='^data.*seg$')
  }
  if(is.null(cbio.RAWseg.file)){
    cbio.RAWseg.file <- .checkFile(cbio.path, 'RAWdata_cna_hg19.seg',   pat='^RAW.*seg$')
  }
  
  ## Create the data_CNA.txt file: https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#discrete-copy-number-data
  cna.mat <- suppressWarnings(Reduce(function(x,y) merge(x,y, by=c('SYMBOL', 'ENTREZ'), 
                                                         all.x=TRUE, all.y=TRUE),
                                     lapply(gr.cnv, function(cnv){ cnv$genes$genes[,c('SYMBOL', 'ENTREZ', 'TCN')]})))
  cna.mat <- .adjustCnaMat(cna.mat, mat.type = 'CNA', amp.thresh = amp.thresh, ord = NULL)
  
  ## Create the data_linear_CNA.txt file: https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#continuous-copy-number-data
  linear.mat <- suppressWarnings(Reduce(function(x,y) merge(x,y, by=c('SYMBOL', 'ENTREZ'), all.x=TRUE),
                                        lapply(gr.cnv, function(cnv){ cnv$genes$genes[,c('SYMBOL', 'ENTREZ', 'seg.mean')]})))
  linear.mat <- .adjustCnaMat(linear.mat, mat.type = 'linear', ord = cna.mat$Hugo_Symbol)
  
  ## Create the data_cna_hg19.seg data file: https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#segmented-data
  segs <- do.call("rbind", lapply(gr.cnv, function(x) x$seg))
  segs$ID <- gsub(".[0-9]*$", "", rownames(segs))
  raw.segs <- segs
  segs <- segs[,c('ID', 'seqnames', 'start', 'end', 'width', 'seg.mean')]
  colnames(segs) <- c('ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean')
  
  ## If existing cBio objects exist, append to the existing data structure
  if(add.on.to.existing){
    if(as.logical(cbio.cna.file['exists'])){
      cna.mat <- .appendToCbioMat(cbio.path, cbio.cna.file, cna.mat, ...)
    }
    if(as.logical(cbio.linear.file['exists'])){
      linear.mat <- .appendToCbioMat(cbio.path, cbio.linear.file, linear.mat, ...)
    }
    if(as.logical(cbio.seg.file['exists'])){
      segs <- .appendToCbioSeg(cbio.path, cbio.seg.file, segs, ...)
    }
    if(as.logical(cbio.RAWseg.file['exists'])){
      raw.segs <- .appendToCbioSeg(cbio.path, cbio.RAWseg.file, raw.segs, raw=TRUE, ...)
    }
  }
  
  ## Write cBioportal Matrices
  .write <- function(...){
    write.table(..., sep="\t", col.names=TRUE, row.names=FALSE, quote=F)
  }
  .write(x=raw.segs, file=file.path(cbio.path, cbio.RAWseg.file['file']))
  .write(x=segs, file=file.path(cbio.path, cbio.seg.file['file']))
  .write(x=cna.mat, file=file.path(cbio.path, cbio.cna.file['file']))
  .write(x=linear.mat, file=file.path(cbio.path, cbio.linear.file['file']))
}

#########################################
#### Building Expression Sets (ESet) ####
.overlapMetaWithExprs <- function(exprs, meta=NULL, pre.regex="^", post.regex="(.cel)?$"){
  ## Identify which columns in the metadata contains the IDs that match the exprs matrix IDs
  if(!is.null(meta)){
    #meta <- gdsc.meta
    meta.idx <- lapply(colnames(meta), function(col.id){
      match.idxs <- sapply(colnames(exprs), function(c.id){
        grep(pattern = paste0(pre.regex, c.id, post.regex), x = meta[,col.id], ignore.case = T)
      })
      return(unlist(match.idxs))
    })
    idlen <- sapply(meta.idx, length)
    col.id.idx <- which.max(idlen)
  } else {
    idlen <- 0
    col.id.idx <- 1
  }
  
  
  if(idlen[col.id.idx] > 0){
    # If a column containing Sample IDs is found...
    tmsg(paste0("Building phenodata on ", idlen[col.id.idx], "/", ncol(exprs), " samples"))
    if(idlen[col.id.idx] > ncol(exprs)){
      multi.idx <- which(table(names(meta.idx[[col.id.idx]])) > 1)
      stop(paste0("One or more sample(s) (", 
                  paste(names(meta.idx[[col.id.idx]])[multi.idx], collapse=","),
                  ") was/were mapped to multiple values in the metadata"))
    }
    meta$phenodata.id <- NA
    meta$phenodata.id[meta.idx[[col.id.idx]]] <- colnames(exprs)
    
    tmsg(paste0("Meta column containing SampleIDs: ", colnames(meta)[col.id.idx]))
    
    ## Identify the remaining metadata features
    other.col.idx <- c(1:ncol(meta))[-col.id.idx]
    other.col.id <- colnames(meta)[other.col.idx]
    tmsg(paste0("Misc meta columns: ", paste(other.col.id, collapse=",")))
    
    ## Order and compose the meta data structure
    m.meta <- merge(data.frame(phenodata.id=colnames(exprs)), meta, 
                    by='phenodata.id', all.x=TRUE)[,-1]
    rownames(m.meta) <- colnames(exprs)
  } else {
    # If no meta data is present or no columns contained sample IDs that matched
    tmsg("WARNING: No appropriate meta file was provided. Please make sure sample IDs are stored in one column")
    
    m.meta <- data.frame(phenodata.id=colnames(exprs))
    rownames(m.meta) <- colnames(exprs)
  }
  
  return(m.meta)
}

.createEsetEnv <- function(mats, exprs.id='seg.mean'){
  eset.env <- new.env()
  exprs.idx <- grep(exprs.id, names(mats))
  nonexprs.idx <- c(1:length(mats))[-exprs.idx]
  
  assign("exprs", mats[[exprs.idx]], envir=eset.env)
  for(idx in nonexprs.idx){
    assign(names(mats)[idx], mats[[idx]], envir=eset.env)
  }
  
  return(eset.env)
}

reduceEsetMats <- function(gene.lrr, cols, features='SYMBOL', ord=FALSE,
                           keys=c("ENTREZ", "SYMBOL", "ENSEMBL")){
  mt <- lapply(cols, function(each.col, features){
    m <- suppressWarnings(Reduce(f=function(x,y) merge(x,y,by=keys),
                                 lapply(gene.lrr, function(i) i[['genes']][,c(keys, each.col)])))
    if(ord) m <- m[match(gene.lrr[[1]][['genes']][,keys], m[,keys]),]
    if(any(duplicated(m[,features]))) m <- m[-which(duplicated(m[,features])),]
    if(any(is.na(m[,features]))) m <- m[-which(is.na(m[,features])),]
    rownames(m) <- m[,features]
    m <- m[,-c(1:length(keys)),drop=FALSE]
    colnames(m) <- names(gene.lrr) 
    as.matrix(m)
  }, features=features)
  mt
}

#' PSet builder for PharmacoGX
#' @description Builds the PSet that is used in PharmacoGX.  It is defined to interface
#' with the output from annotateRDS.Batch() function
#'
#' @param gr.cnv List output from annotateRDS.Batch() function
#' @param anno.name The ID of the group you are analyzing (e.g. 'GDSC')
#' @param pset.path Relative path to the PSet output path
#' @param cols The columns from gr.cnv to build into assayData() objects
#' @param meta A dataframe where one column contains the sample IDs found in gr.cnv
#' @param ... 
#'
#' @return NULL
#' @export
#'
#' @examples 
#'     buildPSetOut(gr.cnv, "CGP", pset.path, meta=cell.line.anno)
buildPSetOut <- function(gr.cnv, anno.name, pset.path, 
                         cols=c('seg.mean', 'nAraw', 'nBraw', 'nMinor', 'nMajor', 'TCN'), 
                         verbose=T, seg.id='seg.mean', out.idx=NULL, ...){
  dir.create(pset.path, recursive = T, showWarnings = F)
  
  #### Assemble assayData environment ####
  if(verbose) print("Building assayData [Genes]...")
  genes.cols <- cols[which(cols %in% colnames(gr.cnv[[1]]$genes$genes))]
  gene.mats <- EaCoN:::reduceEsetMats(lapply(gr.cnv, function(x) x$genes), 
                                      genes.cols, keys='SYMBOL', features='SYMBOL', ord=TRUE)
  names(gene.mats) <- genes.cols
  gene.env <- EaCoN:::.createEsetEnv(gene.mats, seg.id)
    
  #### Assemble featureData #### 
  if(verbose) print("Assembling featureData [Genes]...")
  gene.fdata <- AnnotatedDataFrame(data=as.data.frame(matrix(nrow=nrow(gene.env$exprs),ncol=0)),
                                   varMetadata=data.frame(labelDescription=c()))
  rownames(gene.fdata) <- rownames(gene.env$exprs)
  
  #### Assemble PhenoData ####
  if(verbose) print("Assembling phenoData...")
  # Use exprs(assayData)
  meta <- EaCoN:::.overlapMetaWithExprs(exprs=gene.env$exprs, ...)
  cl.phenoData <- new("AnnotatedDataFrame", data=meta)
  
  #### Assemble the eset #### 
  gene.eset <- ExpressionSet(assayData=gene.env,
                             phenoData=cl.phenoData,
                             annotation=anno.name,
                             featureData=gene.fdata)
  
  
  
  feature.sets <- names(gr.cnv[[1]])
  feature.sets <- feature.sets[-grep("genes|seg", feature.sets)]
  f.esets <- lapply(feature.sets, function(f, pheno, anno){
    if(verbose) print(paste0("Assembling PSet for ", f))
    ## Build assayData
    assay.mats <- EaCoN:::reduceEsetMats(lapply(gr.cnv, function(x) x[[f]]), 
                                 cols, features='ID', keys='ID', ord=TRUE)
    names(assay.mats) <- cols
    assay.env <- EaCoN:::.createEsetEnv(assay.mats, seg.id)
    
    ## Build featureData
    fdata <- AnnotatedDataFrame(data=as.data.frame(gr.cnv[[1]][[f]])[,1:6],
                                varMetadata=data.frame(labelDescription=c(
                                  "Chromosome", "start", "end", 
                                  "width", "Strand", "feature_id")))
    rownames(fdata) <- rownames(assay.env$exprs)
    
    ExpressionSet(assayData=assay.env,
                  phenoData=pheno,
                  annotation=anno,
                  featureData=fdata)
  }, pheno=cl.phenoData, anno=anno.name)
  names(f.esets) <- feature.sets
  
  
  
  #### output segmented data ####
  if(length(out.idx) > 1) out.idx <- paste(out.idx, collapse="-")
  save(gene.eset, file=file.path(pset.path, paste0(anno.name, "_gene_ESet.", out.idx, ".RData")))
  sapply(feature.sets, function(f){
    out.eset <- f.esets[[f]]
    save(out.eset, file=file.path(pset.path, paste0(anno.name, "_", f, 
                                                    "_ESet.", out.idx, ".RData")))
  })
}




.junk <- function(){
  stop("This function should never be called")
  ### Temp
  fp <- file.path(sample, "ASCAT", "ASCN", paste0(sample, ".gammaEval.txt"))
  fit.val <- read.table(fp, sep="\t", header=TRUE, stringsAsFactors = F,
                        check.names = F, fill=T)
  ###
  
  all.fits <- list()
  all.fits[[sample]] <- list("fit"=fit.val,
                             "sample"=sample)
  
  gr.cnv <- lapply(all.fits, function(x) annotateRDS(x$fit, x$sample, segmenter, 
                                                     build='hg19', bin.size=50000))
  cbio.path=file.path("out", "cBio")
  buildCbioOut(gr.cnv, cbio.path="./out/cBio", overwrite=sample)
  
  pset.path=file.path("out", "PSet")
  buildPSetOut(gr.cnv, "CGP", pset.path, meta=cell.line.anno)
}
