
# packages ----------------------------------------------------------------

library(Rsamtools)
library(data.table)

# set variables depending on input data -----------------------------------

num_cores <- 16 
setDTthreads(num_cores)

#barcode expected list, example for Smart-seq3xpress data
xpress_annot <- fread("barcode_annotation.txt")
bcs <- xpress_annot[TSO=="M462_TSO" & cell_source=="HEK"]$XC_DNBPE

taglist <- c("BC", "UB","GE","GI") #for current zUMIs versions
taglist <- c("XC", "XM","XT") #for bam files from old zUMIs versions

bam_path <- "SmartseqLV.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam"

# loading function --------------------------------------------------------

load_reads <- function(featfile, taglist, bcs, cores){
  idxstats <- Rsamtools::idxstatsBam(featfile)
  if("*" %in% idxstats$seqnames){
    idxstats <- idxstats[idxstats$seqnames != "*", ]
    idxstats$seqnames <- as.character(idxstats$seqnames)
  }
  
  if(Rsamtools::testPairedEndBam(featfile)){
    is_PE <- TRUE
  }else{
    is_PE <- FALSE
  }
  
  rsamtools_reads <- parallel::mclapply(1:nrow(idxstats), function(x) {
    if(is_PE){
      parms <- ScanBamParam(tag=taglist,
                            what="pos",
                            flag = scanBamFlag(isFirstMateRead = TRUE),
                            tagFilter = list(BC = bcs),
                            which = GRanges(seqnames = idxstats[x,"seqnames"], ranges = IRanges(1,idxstats[x,"seqlength"])))
      dat <- scanBam(file = featfile, param = parms)
    }else{
      parms <- ScanBamParam(tag=taglist,
                            what="pos",
                            tagFilter = list(BC = bcs),
                            which = GRanges(seqnames = idxstats[x,"seqnames"], ranges = IRanges(1,idxstats[x,"seqlength"])))
      dat <- scanBam(file = featfile, param = parms)
    }
    
    if("GI" %in% taglist){
      dt <- data.table(RG = dat[[1]]$tag$BC, UB = dat[[1]]$tag$UB, GE = dat[[1]]$tag$GE, GEin = dat[[1]]$tag$GI)
    }else{
      dt <- data.table(RG = dat[[1]]$tag$XC, UB = dat[[1]]$tag$XM, GE = dat[[1]]$tag$XT)
    }
    return(dt)
  }, mc.cores = cores)
  rsamtools_reads <- rbindlist(rsamtools_reads, fill = TRUE, use.names = TRUE)
  
  if("GI" %in% taglist){
    rsamtools_reads[ , ftype :="NA"][
      is.na(GEin)==F, ftype :="intron"][
        is.na(GE)==F  , ftype:="exon"][
          is.na(GE)     , GE:=GEin][
            ,GEin:=NULL ]
  }
  rsamtools_reads <- rsamtools_reads[GE!="NA"]
  return(rsamtools_reads)
}

# xpress added on 19072022 ------------------------------------------------

xpress_dat_hek <- load_reads(bam_path, taglist, bcs, num_cores) #load reads
xpress_dat_hek <- xpress_dat_hek[UB != "",.N, by = c("RG","GE","UB")] # count up reads per UMI
xpress_dat_hek <- xpress_dat_hek[grep("ENSG",GE)] #make sure only typical gene IDs are present (eg. no spike ins)
fwrite(xpress_dat_hek, file = "smartseq3xpress_hek293_data.hd1.txt.gz", quote = F, sep ="\t")
