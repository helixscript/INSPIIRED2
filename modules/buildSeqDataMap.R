#!/usr/bin/env -S Rscript --vanilla
for (p in c('argparse', 'tidyverse', 'ShortRead', 'parallel', 'data.table')) suppressPackageStartupMessages(library(p, character.only = TRUE))

parser <- ArgumentParser()
parser$add_argument("--outputDir",        type = "character", required = TRUE,       help = "Directory for output files.")
parser$add_argument("--inputData",        type = "character", required = TRUE,       help = "Path to sequencing data.")
parser$add_argument("--softwareRoot",     type = "character", required = TRUE,       help = "Path to INSPIIRED2 installation.")
parser$add_argument("--threads",          type = "integer",   default  = 50,         help = "Number of threads to use.")
parser$add_argument("--fileTag",          type = "character", default  = "testHMMs", help = "String appended to output files in the output directory.")
parser$add_argument("--ramDiskPath",      type = "character", default  = "/dev/shm", help = "Path to system ramdisk file system. Will default to output directory if ramdisk file system is not supported.")
parser$add_argument("--mapWidth",         type = "integer",   default  = 100,        help = "Width of map in NT.")
parser$add_argument("--startPosBinWidth", type = "integer",   default  = 3,          help = "Bin size of alignment start positions before overflow bin.")
parser$add_argument("--nBinRows",         type = "integer",   default  = 5000,       help = "Number of heat map rows.")
parser$add_argument("--outputImgHeight",  type = "double",    default  = 5,          help = "Height, in inches, of outpout image.")

runModule <- function() {
  startModule()
  yaml::write_yaml(args, file.path(args$outputDir, paste0(args$fileTag, '.yml')))
  
  on.exit({
    unlink(args$tmpDir, recursive = TRUE, force = TRUE)
    unlink(args$logDir, recursive = TRUE, force = TRUE)
    unlink(args$ramDisk, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  makeSeqHeatmap <- function(x,nBins=5000,sortReads=TRUE,sortWidth=20,
                             genomicStart=50,genomicWidth=11,
                             consensusHigh=.75,grayFloor=.25,
                             grayDark=.25,grayLight=.90){
    nReads <- as.numeric(length(x)); nPos <- unique(width(x)); bases <- c("A","C","G","T")
    if(length(nPos)!=1) stop("All sequences must have the same width.")
    if(nBins>nReads) stop("nBins cannot exceed number of reads.")
    
    ord <- seq_along(x)
    if(sortReads){
      k1 <- subseq(x,start=1,width=min(sortWidth,nPos))
      if(!is.null(genomicStart) && genomicStart<=nPos){
        k2 <- subseq(x,start=genomicStart,width=min(genomicWidth,nPos-genomicStart+1))
        ord <- order(xscat(k1,k2))
      } else ord <- order(k1)
      x <- x[ord]
    }
    
    br <- floor(seq(0,nReads,length.out=nBins+1))
    consensus <- character(nBins*nPos); agreement <- numeric(nBins*nPos)
    
    for(b in seq_len(nBins)){
      idx <- (br[b]+1):br[b+1]
      cm <- consensusMatrix(x[idx],baseOnly=FALSE)
      cnt <- cm[bases,seq_len(nPos),drop=FALSE]
      win <- max.col(t(cnt),ties.method="first")
      ii <- ((b-1)*nPos+1):(b*nPos)
      consensus[ii] <- bases[win]
      agreement[ii] <- cnt[cbind(win,seq_len(nPos))]/length(idx)
    }
    
    d <- data.table(bin=rep(seq_len(nBins),each=nPos),
                    position=rep(seq_len(nPos),nBins),
                    consensus=consensus,
                    agreement=agreement)
    
    lo <- round(grayFloor*100); hi <- round(consensusHigh*100)-1
    grayLevels <- paste0("g",lo:hi)
    grayCols <- grDevices::gray(seq(grayDark,grayLight,length.out=length(grayLevels)))
    names(grayCols) <- grayLevels
    
    d[,grayPct:=pmax(lo,pmin(hi,round(agreement*100)))]
    d[,tile:=ifelse(agreement>=consensusHigh,consensus,paste0("g",grayPct))]
    
    cols <- c(A="green3", C="gold", T="dodgerblue3", G="red3", grayCols)
    mid <- round((lo+hi)/2)
    legendBreaks <- c("A","C","G","T",paste0("g",hi),paste0("g",mid),paste0("g",lo))
    legendLabels <- c("A","C","G","T",  
                      paste0("<",round(consensusHigh*100),"% to ", paste0("\u2265",mid,"%")),
                      paste0('<',mid,"% to ",paste0("\u2265",lo,"%")),
                      paste0("<",lo,"%"))
    
    p <- ggplot(d,aes(position,bin,fill=tile))+
      geom_raster()+
      scale_fill_manual(values=cols,breaks=legendBreaks,labels=legendLabels,name="Base / agreement")+
      scale_y_continuous(expand=c(0,0),labels=scales::comma)+
      scale_x_continuous(expand=c(0,0),breaks=seq(0,nPos,10))+
      labs(x="Nucleotide position",y="Read bin")+
      theme_minimal(base_size=10)+
      theme(panel.grid=element_blank())
    
    list(data=d,plot=p,order=ord)
  }
   
   o <- subseq(readFastq(args$inputData)@sread, 1, args$mapWidth)
   p <- makeSeqHeatmap(o, nBins = args$nBinRows)
   
   suppressMessages(ggsave(file.path(args$outputDir, paste0(args$fileTag, '.png')), p$plot, units = 'in', height = args$outputImgHeight))
}

args <- parser$parse_args()
source(file.path(args$softwareRoot, 'lib', 'common.R'))

tryCatch({
  runModule()
}, error = function(e) {
  cat("ERROR: ", conditionMessage(e), "\n", sep = "", file = stderr())
  flush(stderr())
  quit(save = "no", status = 1, runLast = FALSE)
})