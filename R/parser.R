
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

	motifNames <- sort(as.character(unique(elementMetadata(gr)[, "motifName"])))
	if (length(motifNames) > 1){
	# if (FALSE){
		rs <- lapply(motifNames, FUN=function(mn){gr[elementMetadata(gr)[, "motifName"]==mn]})
		names(rs) <- motifNames
		md <- do.call("rbind", rep(list(metadata), length(motifNames)))
		md[,"name"] <- paste0(motifNames)
		md[,"description"] <- paste0(md[,"description"], " - ", motifNames)
	} else {
		rs <- list(gr)
		md <- metadata
	}
	return(getParseResult(rs, md))
}

parse.mumbach2017.hichip.supptab <- function(fn, assembly, metadata){
	require(openxlsx)
	require(muRtools) #grLiftOver
	fUrl <- "https://media.nature.com/original/nature-assets/ng/journal/v49/n11/extref/ng.3963-S4.xlsx"
	cellTab <- data.frame(
		cellType=c("mESC", "Tnaive", "Th17", "Treg", "GM12878", "K562", "MyLa", "HCASMC"),
		assembly=c("mm9", "hg19", "hg19", "hg19", "hg19", "hg19", "hg19", "hg19"),
		antibody=c("H3K27ac", "H3K27ac", "H3K27ac", "H3K27ac", "H3K27ac", "H3K27ac", "H3K27ac", "H3K27ac"),
		sheetName=c("mES H3K27ac Loops", "Naive H3K27ac Loops", "Th17 H3K27ac Loops", "Treg H3K27ac Loops", "GM12878 H3K27ac Loops", "K562 H3K27ac Loops", "My-La H3K27ac Loops", "HCASMC H3K27ac Loops"),
		stringsAsFactors=FALSE
	)
	compatAssemblies <- list(mm9=c("mm9", "mm10"), mm10=c("mm9", "mm10"), hg19=c("hg19", "hg38"), hg38=c("hg19", "hg38"))


	fn <- tempfile(fileext=".xlsx")
	download.file(fUrl, fn)

	cellTab <- cellTab[cellTab$assembly %in% compatAssemblies[[assembly]],,drop=FALSE]
	Ncells <- nrow(cellTab)
	if (Ncells < 1) logger.error("Found no compatible assembly")
	grl <- list()
	for (i in 1:Ncells){
		dfj <- read.xlsx(fn, cellTab[i, "sheetName"])
		dfx <- dfj[,1:3]
		colnames(dfx) <- c("chrom", "chromStart", "chromEnd")
		dfx[,"loop"] <- paste0("loop", 1:nrow(dfj))
		dfy <- dfj[,4:6]
		colnames(dfy) <- c("chrom", "chromStart", "chromEnd")
		dfy[,"loop"] <- paste0("loop", 1:nrow(dfj))
		df <- rbind(dfx, dfy)
		df[,"loop"] <- factor(df[,"loop"], levels=paste0("loop", 1:nrow(dfj)))

		df[, "chrom"] <- adjChrom(df[, "chrom"])
		assembly4gr <- as.character(cellTab[i, "assembly"])
		gr <- df2granges(df, ids=df$loop, coord.format="B1RI", assembly=assembly4gr, doSort=TRUE)

		# liftover if assembly does not match
		if (assembly != assembly4gr){
			gr <- grLiftOver(gr, assembly)
			nOcc <- table(names(gr))
			gr <- gr[nOcc[names(gr)] == 2] # discard regions where one matching partner could not be mapped by liftOver
		}
		grl <- c(grl, list(gr))
	}
	names(grl) <- cellTab[,"cellType"]
	
	md <- do.call("rbind", rep(list(metadata), Ncells))
	md[,"name"] <- paste0(cellTab[,"cellType"])
	md[,"description"] <- paste0(md[,"description"], " - ", cellTab[,"cellType"], " - ", cellTab[,"antibody"])

	return(getParseResult(grl, md))	
}

parse.hnisz2013.se.supp <- function(fn, assembly, metadata){
	require(muRtools) #grLiftOver
	fUrl <- fn

	fn <- tempfile(fileext=".zip")
	download.file(fUrl, fn)
	exDir <- tempfile(fileext="unzip")
	unzip(fn, exdir=exDir)
	bedFns <- list.files(exDir)
	cellTypes <- gsub("\\.bed$", "", bedFns)

	grl <- list()
	for (i in seq_along(bedFns)){
		fn <- bedFns[i]
		rl <- readLines(file.path(exDir, fn))
		iStart <- grep("^track.+color=255,0,0$", rl)
		tt <- read.table(file.path(exDir, fn), skip=iStart, quote="", sep="\t", stringsAsFactors=FALSE)
		gr <- df2granges(tt, ids=NULL, chrom.col=1L, start.col=2L, end.col=3L, strand.col=NULL, coord.format="B1RI", assembly="hg19")
		colnames(elementMetadata(gr)) <- c("name", "count")
		elementMetadata(gr)[,"cellType"] <- cellTypes[i]
		if (assembly != "hg19"){
			gr <- grLiftOver(gr, assembly)
		}
		grl <- c(grl, list(gr))
	}
	
	md <- do.call("rbind", rep(list(metadata), length(cellTypes)))
	md[,"name"] <- cellTypes
	md[,"description"] <- paste0(md[,"description"], " - ", cellTypes)

	return(getParseResult(grl, md))	
}

parse.gwas.catalog.snps <- function(fn, assembly, metadata){
	require(muRtools)

	if (is.element(assembly, c("hg19", "GRCh37"))){
		require(SNPlocs.Hsapiens.dbSNP144.GRCh37)
		snpDb <- SNPlocs.Hsapiens.dbSNP144.GRCh37
	} else if (is.element(assembly, c("hg38", "GRCh38"))){
		require(SNPlocs.Hsapiens.dbSNP144.GRCh38)
		snpDb <- SNPlocs.Hsapiens.dbSNP144.GRCh38
	} else {
		logger.error(c("GWAS Catalog parsing is not supported for assembly", assembly))
	}
	# genomeObj <- getGenomeObject(assembly, adjChrNames=TRUE)
	# fn <- "http://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
	fUrl <- fn
	fn <- tempfile(fileext=".tsv")
	download.file(fUrl, fn)

	tt <- readTab(fn)
	tt <- tt[,!(colnames(tt) %in% c("CHR_ID", "CHR_POS"))]
	tt <- tt[grepl("^rs", tt$SNPS),]

	#problem: some rows in tt have multiple SNPs --> expand tt to contain multiple rows for these cases
	snpSplit <- strsplit(tt$SNPS, ";")
	tt <- do.call("rbind", lapply(seq_along(snpSplit), FUN=function(i){
		res <- do.call("rbind", rep(list(tt[i, , drop=FALSE]), length(snpSplit[[i]])))
		res[,"SNPS"] <- normalize.str(snpSplit[[i]])
		return(res)
	}))
	tt <- tt[grepl("^rs[0-9]+$", tt$SNPS),]


	snpLocs <- GRanges(snpsById(snpDb, tt$SNPS, ifnotfound="warn"))
	snpLocs <- setGenomeProps(snpLocs, assembly)

	idx <-  match(tt$SNPS, elementMetadata(snpLocs)[,"RefSNP_id"])
	isUnmapped <- is.na(idx)
	if (any(isUnmapped)){
		logger.warning(c("The following SNP IDs could not be matched to coordinates:", paste(unique(tt$SNPS[is.na(idx)]), collapse=",")))
		tt <- tt[!isUnmapped,]
		idx <- idx[!isUnmapped]
	}

	gr <- snpLocs[idx]
	elementMetadata(gr) <- cbind(tt, elementMetadata(gr))

	return(getParseResult(list(gr), metadata))	
}

parse.pics.snps.farh2015.supp <- function(fn, assembly, metadata){
	require(readxl)
	require(muRtools)

	if (is.element(assembly, c("hg19", "GRCh37"))){
		require(SNPlocs.Hsapiens.dbSNP144.GRCh37)
		snpDb <- SNPlocs.Hsapiens.dbSNP144.GRCh37
	} else if (is.element(assembly, c("hg38", "GRCh38"))){
		require(SNPlocs.Hsapiens.dbSNP144.GRCh38)
		snpDb <- SNPlocs.Hsapiens.dbSNP144.GRCh38
	} else {
		logger.error(c("PICS parsing is not supported for assembly", assembly))
	}
	# genomeObj <- getGenomeObject(assembly, adjChrNames=TRUE)
	# fn <- "https://media.nature.com/original/nature-assets/nature/journal/v518/n7539/extref/nature13835-s1.xls"
	fUrl <- fn

	fn <- tempfile(fileext=".xls")
	download.file(fUrl, fn)

	tt <- read_excel(fn)
	tt <- tt[,!(colnames(tt) %in% c("chr", "pos"))]

	snpLocs <- GRanges(snpsById(snpDb, tt$SNP, ifnotfound="warn"))
	snpLocs <- setGenomeProps(snpLocs, assembly)

	idx <-  match(tt$SNP, elementMetadata(snpLocs)[,"RefSNP_id"])
	isUnmapped <- is.na(idx)
	if (any(isUnmapped)){
		logger.warning(c("The following SNP IDs could not be matched to coordinates:", paste(unique(tt$SNP[is.na(idx)]), collapse=",")))
		tt <- tt[!isUnmapped,]
		idx <- idx[!isUnmapped]
	}

	gr <- snpLocs[idx]
	elementMetadata(gr) <- cbind(tt, elementMetadata(gr))

	grl <- split(gr, elementMetadata(gr)[,"Disease"])
	diseases <- names(grl)
	
	md <- do.call("rbind", rep(list(metadata), length(diseases)))
	md[,"name"] <- diseases
	md[,"description"] <- paste0(md[,"description"], " - Disease:", diseases)

	return(getParseResult(grl, md))	
}

parse.gtex.eqtl <- function(fn, assembly, metadata){
	require(muRtools)

	if (assembly=="hg38"){
		# dbSnpPkg <- "SNPlocs.Hsapiens.dbSNP149.GRCh38"
		dbSnpPkg <- "SNPlocs.Hsapiens.dbSNP151.GRCh38"
		require(dbSnpPkg, character.only=TRUE)
		snpObj <- eval(parse(text=dbSnpPkg))
	} else {
		stop(c("GTEx SNP annotation not supported for assembly:", assembly))
	}

	# fn <- ""
	if (fn == "v7"){
		annotUrl <- "https://storage.googleapis.com/gtex_analysis_v7/reference/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz"
		eqtlUrl <- "https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz"
	} else {
		stop(c("Unknown GTEx version:", fn))
	}

	downFn <- tempfile(fileext=".txt.gz")
	download.file(annotUrl, downFn)
	varAnnot <- readTab(downFn)
	rownames(varAnnot) <- varAnnot[,"variant_id"]

	downFn <- tempfile(fileext=".tar.gz")
	download.file(eqtlUrl, downFn)
	destDir <- tempfile() #destDir <- "/scratch/users/muellerf/temp/blubb"
	dir.create(destDir)
	utils::untar(fn, exdir=destDir)
	varDir <- file.path(destDir, gsub("\\.tar\\.gz", "", basename(eqtlUrl)))
	linkFns <- list.files(varDir, pattern="signif_variant_gene_pairs.txt.gz")
	tissues <- gsub("\\.v[0-9]\\.signif_variant_gene_pairs\\.txt\\.gz", "", linkFns)

	qtlTab <- do.call("rbind", lapply(1:length(linkFns), FUN=function(i){
		data.frame(
			tissue=tissues[i],
			readTab(file.path(varDir, linkFns[i]))
		)
	}))

	qtlTab[,"gene_id_short"] <- gsub("\\.[0-9]+$", "", qtlTab[,"gene_id"])

	qtlTab[,"rsId"] <- varAnnot[qtlTab[,"variant_id"], "rs_id_dbSNP147_GRCh37p13"]
	qtlTab[qtlTab[,"rsId"] %in% ".","rsId"] <- NA

	qtlTab <- qtlTab[!is.na(qtlTab[,"rsId"]),]

	gr <- snpsById(snpObj, qtlTab[,"rsId"], ifnotfound="warning")

	if (genome(gr)[1]!=assembly){
		# adjust chromosome names
		prependChr <- !grepl("chr", seqlevels(gr))
		if (any(prependChr)){
			seqlevels(gr)[prependChr] <- paste0("chr", seqlevels(gr)[prependChr])
		}
		seqlevels(gr)[seqlevels(gr)=="chrMT"] <- "chrM"
		gr <- muRtools::setGenomeProps(gr, assembly, onlyMainChrs=TRUE)
	}
	snpIds2match <- elementMetadata(gr)[,"RefSNP_id"]

	idx <- qtlTab[,"rsId"] %in% snpIds2match
	elementMetadata(gr) <- cbind(elementMetadata(gr), qtlTab[idx,])
	gr <- muRtools::sortGr(gr)

	grl <- split(gr, elementMetadata(gr)[,"tissue"])
	tissues <- names(grl)
	# saveRDS(grl, file.path("/oak/stanford/groups/wjg/muellerf/data/GTEx", paste0("GTEx_eQTL_signif_variant_gene_pairs", "_v7", "_hg38", ".rds")))
	
	md <- do.call("rbind", rep(list(metadata), length(tissues)))
	md[,"name"] <- tissues
	md[,"description"] <- paste0(md[,"description"], " - tissue:", tissues)

	return(getParseResult(grl, md))	
}

parse.vista.enhancers <- function(fn, assembly, metadata){
	require(muRtools) #grLiftOver
	fUrl <- "https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?page_size=100;show=1;search.result=yes;page=1;form=search;search.form=no;action=search;search.sequence=1"

	parse.one <- function(string,result){
		m <- do.call(rbind,lapply(seq_along(string),function(i){
		st <- attr(result,"capture.start")[i,]
		substring(string[i],st,st+attr(result,"capture.length")[i,]-1)
		}))
		colnames(m) <- attr(result,"capture.names")
		m
	}

	fn_unformatted <- tempfile(fileext=".fax")
	fn <- tempfile(fileext=".fa")
	download.file(fUrl, fn_unformatted)
	
	# replace <pre> tag
	system2("sed", c("'1 s/^<pre>//'", fn_unformatted, ">", fn))

	lls <- readLines(fn)
	lls <- gsub("^>", "", lls[grepl("^>", lls)])
	re <- paste0("^",
		"(?P<organism>(Human|Mouse))\\|",
		"(?P<chr>chr.+):(?P<chrStart>[0-9]+)-(?P<chrEnd>[0-9]+) \\| ",
		"(?P<element>element.+) \\| ",
		"(?P<status>(positive|negative)) \\|?",
		"(?P<organ> .+)?",
		"$"
	)
	isMatch <- grepl(re, lls, perl=TRUE)
	if (sum(!isMatch) > 0){
		logger.error(c("The following sample ids do not match the regular expression for annotation:", paste(lls[!isMatch], collapse=", ")))
	}
	rem <- regexpr(re, lls, perl=TRUE, ignore.case=FALSE)
	parsed <- parse.one(lls, rem)
	parsed[parsed==""] <- NA
	sannot <- data.frame(
		organism = factor(parsed[,"organism"]),
		chr = parsed[,"chr"],
		chrStart = as.integer(parsed[,"chrStart"]),
		chrEnd = as.integer(parsed[,"chrEnd"]),
		elementName = gsub(" ", "", parsed[,"element"]),
		status = factor(parsed[,"status"]),
		organ=gsub("^ ?\\| ?", "", parsed[,"organ"]),
		stringsAsFactors=FALSE
	)
	organList <- lapply(sannot[,"organ"], FUN=function(s){
		if (is.na(s)) return(NULL)
		ss <- strsplit(s, " | ", fixed=TRUE)[[1]]
		res <- paste(gsub("^(.+)\\[([0-9]+)\\/([0-9]+)\\]$", "\\2", ss), "of", gsub("^(.+)\\[([0-9]+)\\/([0-9]+)\\]$", "\\3", ss))
		names(res) <- gsub("^(.+)\\[([0-9]+)\\/([0-9]+)\\]$", "\\1", ss)
		return(res)
	})
	allOrgs <- sort(unique(do.call("c", lapply(organList, FUN=function(x){
		if (is.null(x)) return(NULL)
		return(names(x))
	}))))
	occTab <- t(sapply(organList, FUN=function(x){
		if (is.null(x)) return(rep(FALSE, length(allOrgs)))
		allOrgs %in% names(x)
	}))
	colnames(occTab) <- allOrgs

	if (grepl("^hg", assembly)){
		# human
		idx <- which(sannot[,"organism"]=="Human")
		refAss <- "hg19"
	} else if (grepl("^mm", assembly)){
		# mouse
		idx <- which(sannot[,"organism"]=="Mouse")
		refAss <- "mm9"
	} else {
		logger.error(c("Unsupported assembly:", assembly))
	}
	sannot <- sannot[idx, ]
	gr <- df2granges(sannot, chrom.col="chr", start.col="chrStart", end.col="chrEnd", strand.col=NULL, coord.format="B1RI", assembly=refAss)
	names(gr) <- NULL
	if (assembly != refAss){
		gr <- grLiftOver(gr, assembly)
		# discard elements that could not be liftOvered
		idx <- idx[sannot[,"elementName"] %in% elementMetadata(gr)[,"elementName"]]
	}
	rm(sannot)

	occTab <- occTab[idx, ]
	occTab <- occTab[, colSums(occTab) > 0]

	featureTypes <- colnames(occTab)
	orgGrl <- lapply(featureTypes, FUN=function(on){
		gr[occTab[,on]]
	})
	names(orgGrl) <- featureTypes

	md <- do.call("rbind", rep(list(metadata), length(featureTypes)))
	md[,"name"] <- paste0(featureTypes)
	md[,"description"] <- paste0(md[,"description"], " - ", featureTypes)
	return(getParseResult(orgGrl, md))
}

