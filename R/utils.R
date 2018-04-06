# adjust chromosome names, such that each has the "chr" prefix and the mitochondrial
# chromosome is named chrM
adjChrom <- function(chrs){
	chrs <- as.character(chrs)
	prependChr <- c(1:99, "X", "Y", "M", "MT")
	do.prepend <- chrs %in% prependChr
	if (any(do.prepend)){
		chrs[do.prepend] <- paste0("chr", chrs[do.prepend])
	}

	#MT --> M
	do.mt <- chrs == "chrMT"
	if (any(do.mt)){
		chrs[do.mt] <- "chrM"
	}
	return(chrs)
}

# convert a bigbed to a bed file using system call
convert.bb.to.bed <- function(inFn,outFn){
	system(paste("bigBedToBed",inFn,outFn))
}


# wrapper around import.gff (rtracklayer) to import GFF file to data frame
import.gff.df <- function(fn, ...){
	require(rtracklayer)
	gr <- import.gff(fn, ...)
	gr.df <- data.frame(gr)
	colnames(gr.df)[1:5] <- c("chrom", "start", "end", "gr.width", "gr.strand")
	gr.df[,"chrom"] <- as.character(gr.df[,"chrom"])
	return(gr.df)
}

# wrapper around import.gff (rtracklayer) to import GFF file to data frame
import.bb.df <- function(fn, ...){
	require(rtracklayer)
	tmpPrefix <- file.path(tempdir(), gsub("\\..*$", "", basename(fn)))
	bbFn <- paste0(tmpPrefix, ".bb")
	bedFn <- paste0(tmpPrefix, ".bed")
	download.file(fn, bbFn)
	convert.bb.to.bed(bbFn, bedFn)
	gr <- import(bedFn, "BED")
	gr.df <- data.frame(gr)
	colnames(gr.df)[1:5] <- c("chrom", "start", "end", "gr.width", "gr.strand")
	gr.df[,"chrom"] <- as.character(gr.df[,"chrom"])
	return(gr.df)
}
