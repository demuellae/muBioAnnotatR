#' getTfMotifClusters_altius
#'
#' Retrieve TF motif model clustering from the database provided 
#'
#' @param assembly genome assembly. Currently there is only a database for hg38
#' @param datasetDir directory where the data has been downloaded to.
#'                 In future releases, this can be \code{NULL} and the data will be downloaded from the website
#' @return An S3 datastructure containing genomic occurrences of motif clusters and cluster annotations
#'
#' @details
#' Data from https://resources.altius.org/~jvierstra/projects/motif-clustering/
#' and https://github.com/jvierstra/motif-clustering
#' @author Fabian Mueller
#' @noRd
getTfMotifClusters_altius <- function(
		assembly="hg38",
		datasetDir="/oak/stanford/groups/wjg/muellerf/resources/TFmotifs/AltiusMotifModelClusters/v1.0"
	){
	if (!is.element(assembly, c("hg38"))){
		logger.error(c("Invalid genome assembly:", assembly))
	}
	
	if (!dir.exists(datasetDir)){
		logger.error(c("datasetDir does not exist:", datasetDir))
	}
	fns <- c(
		clusterBed=list.files(datasetDir, pattern=paste0(assembly, "\\.archetype_motifs.+\\.bed(\\.gz)?$"))[1],
		clusterAnnot=list.files(datasetDir, pattern="archetype_clusters.tsv")[1],
		motifAnnot=list.files(datasetDir, pattern="motifs.tsv")[1]
	)

	cAnnot <- readTab(file.path(datasetDir, fns["clusterAnnot"]))
	cAnnot[,"cid"] <- make.unique(normalize.str(gsub("/", "_", cAnnot[,"Name"])), sep="_")
	rownames(cAnnot) <- cAnnot[,"cid"]
	mAnnot <- readTab(file.path(datasetDir, fns["motifAnnot"]))
	mAnnot[,"cid"] <- cAnnot[match(mAnnot[,"Cluster"], cAnnot[,"Cluster_ID"]),"cid"]
	rownames(mAnnot) <- mAnnot[,"Motif"]

	motifGeneStr <- toupper(sapply(strsplit(mAnnot[,"Motif"], split="[_.]"), FUN=function(x){x[1]}))
	mAnnot[,"geneSymbols"] <- motifGeneStr
	# split by + (Jaspar fusion proteins?)

	gr <- df2granges(
		readTab(file.path(datasetDir, fns["clusterBed"]), header=FALSE),
		ids=NULL,
		chrom.col=1L, start.col=2L, end.col=3L, strand.col=6L,
		coord.format="B1RI",
		assembly=assembly,
		doSort=TRUE,
		adjNumChromNames=TRUE
	)
	colnames(elementMetadata(gr)) <- c("clusterName", "score", "motifName", "number")
	elementMetadata(gr)[,"cid"] <- mAnnot[elementMetadata(gr)[,"motifName"],"cid"]
	grl <- split(gr, elementMetadata(gr)[,"cid"])
	grl <- grl[rownames(cAnnot)]

	res <- list(
		genome=assembly,
		clusterAnnot = cAnnot,
		motifAnnot = mAnnot,
		clusterOcc = grl
	)
	class(res) <- "TfMotifClusters"
	# saveRDS(res, file.path(datasetDir, paste0("tfMotifClusters_", assembly, ".rds")))
	return(res)
}
