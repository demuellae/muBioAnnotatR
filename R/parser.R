
getParseResult <- function(regionSets, metadata){
	res <- list(
		regionSets = regionSets,
		metadata = metadata
	)
	class(res) <- "ParsedAnnotation"
	return(res)
}

parse.bed3 <- function(fn, assembly, metadata){
	df <- read.table(fn, header=FALSE, sep="\t")
	gr <- df2granges(df, ids=NULL, strand.col=NULL, coord.format="B0RE", assembly=assembly, doSort=TRUE, adjNumChromNames=TRUE)
	rs <- list(gr)
	md <- metadata
	return(getParseResult(rs, md))
}

parse.ensembl.regbuild <- function(fn, assembly, metadata){
	df <- import.gff.df(fn)
	assembly4gr <- assembly
	df[, "chrom"] <- adjChrom(df[, "chrom"])
	if (assembly=="hg38"){
		assembly4gr <- "hg38_chr"
	}

	# df[,"featureType_norm"] <- normalize.str(df[,"feature_type"])
	gr <- df2granges(df, ids=df$ID, strand.col="gr.strand", coord.format="B1RI", assembly=assembly4gr, doSort=TRUE)
	featureTypes <- sort(unique(elementMetadata(gr)[, "feature_type"]))
	if (length(featureTypes) > 1){
		rs <- lapply(featureTypes, FUN=function(ft){gr[elementMetadata(gr)[, "feature_type"]==ft]})
		names(rs) <- featureTypes
		md <- do.call("rbind", rep(list(metadata), length(featureTypes)))
		md[,"name"] <- paste0(featureTypes)
		md[,"description"] <- paste0(md[,"description"], " - ", featureTypes)
	} else {
		rs <- list(gr)
		md <- metadata
	}

	return(getParseResult(rs, md))	
}

parse.ensembl.regbuild.bp <- function(fn, assembly, metadata){
	df <- import.bb.df(fn)
	assembly4gr <- assembly
	df[, "chrom"] <- adjChrom(df[, "chrom"])
	if (assembly=="hg38"){
		assembly4gr <- "hg38_chr"
	}

	# df[,"featureType_norm"] <- normalize.str(df[,"feature_type"])
	gr <- df2granges(df, ids=df$name, strand.col="gr.strand", coord.format="B1RI", assembly=assembly4gr, doSort=TRUE)
	elementMetadata(gr)[,"feature_type"] <- gsub("^(.*)_.*$", "\\1",elementMetadata(gr)[,"name"])
	featureTypes <- sort(unique(elementMetadata(gr)[, "feature_type"]))
	if (length(featureTypes) > 1){
		rs <- lapply(featureTypes, FUN=function(ft){gr[elementMetadata(gr)[, "feature_type"]==ft]})
		names(rs) <- featureTypes
		md <- do.call("rbind", rep(list(metadata), length(featureTypes)))
		md[,"name"] <- paste0(featureTypes)
		md[,"description"] <- paste0(md[,"description"], " - ", featureTypes)
	} else {
		rs <- list(gr)
		md <- metadata
	}

	return(getParseResult(rs, md))	
}

parse.gencode <- function(fn, assembly, metadata, addPromoters=TRUE, promoters.up=1500, promoters.down=500){
	# legalGeneTypes <- NULL
	legalGeneTypes <- c("protein_coding", "antisense", "lincRNA", "pseudogene", "processed_transcript", "miRNA", "rRNA", "snRNA", "snoRNA", "misc_RNA")
	df <- import.gff.df(fn)
	assembly4gr <- assembly
	if (assembly=="hg38"){
		df[, "chrom"] <- adjChrom(df[, "chrom"])
		assembly4gr <- "hg38_chr"
	}
	df.genes <- df[df[,"type"]=="gene",]
	# df.genes <- df[df[,"type"]=="gene" & df[,"gene_status"]=="KNOWN",]

	gr <- df2granges(df.genes, ids=df.genes$ID, strand.col="gr.strand", coord.format="B1RI", assembly=assembly4gr, doSort=TRUE)
	geneTypes <- sort(unique(elementMetadata(gr)[, "gene_type"]))
	if (!is.null(legalGeneTypes)){
		geneTypes <- intersect(geneTypes, legalGeneTypes)
	}
	if (length(geneTypes) > 1){
		rs <- lapply(geneTypes, FUN=function(ft){gr[elementMetadata(gr)[, "gene_type"]==ft]})
		names(rs) <- geneTypes
		md <- do.call("rbind", rep(list(metadata), length(geneTypes)))
		md[,"name"] <- paste0("genes_", geneTypes)
		md[,"description"] <- paste0(md[,"description"], " - ", geneTypes)
	} else {
		rs <- list(gr)
		md <- metadata
	}

	if (addPromoters){
		rs.prom <- lapply(rs, FUN=function(x){promoters(x, upstream=promoters.up, downstream=promoters.down)})
		md.prom <- md
		md.prom[,"name"] <- gsub("^genes_", "promoters_", md.prom[,"name"])
		md.prom[,"description"] <- paste0(md.prom[,"description"], " - Promoters [-", promoters.up, "bp:+", promoters.down, "bp]")
		names(rs.prom) <- md.prom[,"name"]
		rs <- c(rs, rs.prom)
		md <- rbind(md, md.prom)
	}
	return(getParseResult(rs, md))	
}

parse.ucsc.cgis <- function(fn, assembly, metadata){
	assembly4gr <- assembly
	if (assembly=="hg38"){
		assembly4gr <- "hg38_chr"
	}

	require(rtracklayer)
	bs <- browserSession()
	genome(bs) <- assembly
	tt <- getTable(ucscTableQuery(bs, track="CpG Islands", table="cpgIslandExt"))
	gr <- df2granges(tt, ids=tt[,"name"], chrom.col="chrom", start.col="chromStart", end.col="chromEnd", coord.format="B0RE", assembly=assembly4gr, adjNumChromNames=TRUE, doSort=TRUE)
	
	rs <- list(gr)
	md <- metadata
	return(getParseResult(rs, md))
}

parse.ucsc.repeats <- function(fn, assembly, metadata){
	# legalRepeatTypes <- NULL
	legalRepeatTypes <- c("DNA", "LINE", "LTR", "SINE", "Satellite", "Simple_repeat", "RNA", "rRNA", "scRNA", "snRNA", "tRNA", "Other", "Unknown", "Low_complexity")

	assembly4gr <- assembly
	if (assembly=="hg38"){
		assembly4gr <- "hg38_chr"
	}

	require(rtracklayer)
	bs <- browserSession()
	genome(bs) <- assembly
	tt <- getTable(ucscTableQuery(bs, track="RepeatMasker", table="rmsk"))
	gr <- df2granges(tt, ids=tt[,"id"], chrom.col="genoName", start.col="genoStart", end.col="genoEnd", strand.col="strand", coord.format="B0RE", assembly=assembly4gr, adjNumChromNames=TRUE, doSort=TRUE)
	
	repeatTypes <- sort(unique(elementMetadata(gr)[, "repClass"]))
	if (!is.null(legalRepeatTypes)){
		repeatTypes <- intersect(repeatTypes, legalRepeatTypes)
	}
	if (length(repeatTypes) > 1){
		rs <- lapply(repeatTypes, FUN=function(ft){gr[elementMetadata(gr)[, "repClass"]==ft]})
		names(rs) <- repeatTypes
		md <- do.call("rbind", rep(list(metadata), length(repeatTypes)))
		md[,"name"] <- paste0(repeatTypes)
		md[,"description"] <- paste0(md[,"description"], " - ", repeatTypes)
	} else {
		rs <- list(gr)
		md <- metadata
	}

	return(getParseResult(rs, md))
}

parse.motifmatchr <- function(fn, assembly, metadata){
	require(motifmatchr)
	require(muRtools)
	require(ChrAccR) # prepareMotifmatchr

	assembly4gr <- assembly
	if (assembly=="hg38"){
		assembly4gr <- "hg38_chr"
	}
	seqlengths <- getSeqlengths4assembly(assembly4gr)
	validChroms <- names(seqlengths)

	motifs <- strsplit(fn, ";")[[1]]

	mmArgs <- prepareMotifmatchr(assembly4gr, motifs)
	
	# for each chromosome, get the motif matches
	gr <- do.call("c", lapply(names(mmArgs$genome), FUN=function(chromName){
		# logger.status(c("chr:", chromName))#TODO: remove
		mmRes <- matchMotifs(mmArgs$motifs, mmArgs$genome[[chromName]], out="positions")
		# convert the resulting IRangesList to GRanges
		return(do.call("c", lapply(names(mmRes), FUN=function(motifName){
			x <- unlist(mmRes[[motifName]])
			rr <- GRanges()
			if (length(x) > 0){
				elementMetadata(x)[,"motifName"] <- motifName
				strands <- elementMetadata(x)[,"strand"]
				elementMetadata(x)[,"strand"] <- NULL
				rr <- GRanges(seqnames=chromName,ranges=x,strand=strands)
			}
			seqlevels(rr) <- validChroms
			seqlengths(rr) <- seqlengths
			genome(rr) <- assembly4gr
			return(rr)
		})))
	}))
	elementMetadata(gr)[,"motifName"] <- factor(elementMetadata(gr)[,"motifName"])

	# sort
	oo <- order(as.integer(seqnames(gr)), start(gr), end(gr), as.integer(strand(gr)))
	gr <- gr[oo]

	rs <- list(gr)
	md <- metadata
	return(getParseResult(rs, md))
}
