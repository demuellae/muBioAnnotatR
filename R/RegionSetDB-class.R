#' RegionSetDB
#'
#' A class for storing a database of genomic region annotations
#'
#' @section Slots:
#' \describe{
#'   \item{\code{parserSummary}}{
#'		A summary table containing information for parsing region sets
#'   }
#'   \item{\code{sampleAnnot}}{
#'		Sample annotation Table
#'   }
#'   \item{\code{genomes}}{
#'		Character vector of genome assemblies contained in the database
#'   }
#' }
#'
#' @section Methods:
#' \describe{
#'    \item{\code{\link{blubb,RegionSetDB-method}}}{
#'      some function
#'    }
#' }
#'
#' @name RegionSetDB-class
#' @rdname RegionSetDB-class
#' @author Fabian Mueller
#' @exportClass RegionSetDB
setClass("RegionSetDB",
	slots = list(
		parserSummary = "data.frame",
		regionSets    = "list",
		regionSetMd   = "data.frame",
		dbDir         = "characterOrNULL",
		loadingStatus = "character"
	),
	package = "muBioAnnotatR"
)
setMethod("initialize","RegionSetDB",
	function(
		.Object
	) {
		emptyMdTab <- data.frame(
			pid                    = character(0),
			name                   = character(0),
			description            = character(0),
			assembly               = character(0),
			collection             = character(0),
			collection_description = character(0),
			source                 = character(0),
			parser                 = character(0),
			file                   = character(0),
			stringsAsFactors = FALSE
		)
		.Object@parserSummary <- emptyMdTab
		.Object@regionSets    <- list()
		.Object@regionSetMd   <- emptyMdTab
		.Object@dbDir         <- NULL
		.Object@loadingStatus <- character(0)
		.Object
	}
)

#' @param parserAnnot A data frame containing annotations for the region sets to be parsed.
#'                    requires the following columns: \code{"pid", "name", "description", "assembly", "collection", "collection_description", "source", "parser", "file"}
#' @name RegionSetDB
#' @rdname RegionSetDB-class
#' @aliases intialize,RegionSetDB-method
#' @export
RegionSetDB <- function(parserAnnot=NULL){
	obj <- new("RegionSetDB")
	if (!is.null(parserAnnot)){
		obj <- addRegionSetFromDataFrame(obj, parserAnnot)
	}
	return(obj)
}

################################################################################
# Getters and Setters
################################################################################
if (!isGeneric("genomes")) {
	setGeneric(
		"genomes",
		function(.object) standardGeneric("genomes"),
		signature=c(".object")
	)
}
#' genomes-methods
#'
#' Return genome assemblies contained in the dataset
#'
#' @param .object \code{\linkS4class{RegionSetDB}} object
#' @return Character vector of genome assemblies contained in the dataset
#'
#' @rdname genomes-RegionSetDB-method
#' @docType methods
#' @aliases genomes
#' @aliases genomes,RegionSetDB-method
#' @author Fabian Mueller
#' @export
setMethod("genomes",
	signature(
		.object="RegionSetDB"
	),
	function(
		.object
	) {
		return(sort(unique(.object@regionSetMd[,"assembly"])))
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("collections")) {
	setGeneric(
		"collections",
		function(.object, ...) standardGeneric("collections"),
		signature=c(".object")
	)
}
#' collections-methods
#'
#' Return collectionss contained in the dataset
#'
#' @param .object \code{\linkS4class{RegionSetDB}} object
#' @param genomeNames optional vector of genome assemblies for which the collection should be retrieved
#' @return Character vector of collection names contained in the dataset
#'
#' @rdname collections-RegionSetDB-method
#' @docType methods
#' @aliases collections
#' @aliases collections,RegionSetDB-method
#' @author Fabian Mueller
#' @export
setMethod("collections",
	signature(
		.object="RegionSetDB"
	),
	function(
		.object,
		genomeNames=genomes(.object)
	) {
		tt <- .object@regionSetMd[.object@regionSetMd[,"assembly"] %in% genomeNames,]
		return(sort(unique(tt[,"collection"])))
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("regionSetNames")) {
	setGeneric(
		"regionSetNames",
		function(.object, ...) standardGeneric("regionSetNames"),
		signature=c(".object")
	)
}
#' regionSetNames-methods
#'
#' Return region set names contained in the dataset
#'
#' @param .object \code{\linkS4class{RegionSetDB}} object
#' @param genomeNames optional vector of genome assemblies for which the names should be retrieved
#' @param collectionNames optional vector of collections for which the names should be retrieved
#' @return Character vector of region set names contained in the dataset
#'
#' @rdname regionSetNames-RegionSetDB-method
#' @docType methods
#' @aliases regionSetNames
#' @aliases regionSetNames,RegionSetDB-method
#' @author Fabian Mueller
#' @export
setMethod("regionSetNames",
	signature(
		.object="RegionSetDB"
	),
	function(
		.object,
		genomeNames=genomes(.object),
		collectionNames=collections(.object)
	) {
		tt <- .object@regionSetMd[(.object@regionSetMd[,"assembly"] %in% genomeNames) & (.object@regionSetMd[,"collection"] %in% collectionNames), ]
		return(sort(unique(tt[,"name"])))
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("regionSetGr")) {
	setGeneric(
		"regionSetGr",
		function(.object, ...) standardGeneric("regionSetGr"),
		signature=c(".object")
	)
}
#' regionSetGr-methods
#'
#' Return region set \code{GRanges} object for a given name, collection and genome
#'
#' @param .object \code{\linkS4class{regionSetGrDB}} object
#' @param name    name of the region set to be retrieved (character)
#' @param col     name of the collection the region set is in (character)
#' @param gen     genome assembly for the region set (character)
#' @return \code{GRanges} object containing the region set or \code{NULL} if it does not exist
#'
#' @rdname regionSetGr-regionSetGrDB-method
#' @docType methods
#' @aliases regionSetGr
#' @aliases regionSetGr,regionSetGrDB-method
#' @author Fabian Mueller
#' @export
setMethod("regionSetGr",
	signature(
		.object="RegionSetDB"
	),
	function(
		.object,
		name,
		col,
		gen
	) {
		i <- which(.object@regionSetMd[,"assembly"] == gen & .object@regionSetMd[,"collection"] == col & .object@regionSetMd[,"name"] == name)
		if (length(i) < 1){
			logger.warning(c("There is no region set", name, "for collection", col, "and genome", gen, "--> returning NULL"))
			return(NULL)
		}
		if (length(i) > 1){
			logger.error(c("Invalid RegionSetDB objet: there are multiple instances of region set", name, "for collection", col, "and genome", gen, "--> returning NULL"))
			return(NULL)
		}
		if (.object@loadingStatus[i] == "unloaded"){
			logger.status(c("Loading region set from disk"))
			.object <- loadRegionSets(.object, indx=i) #TODO: is there some way, we can use concurrent programming to actually modify the input .object and make subsequent loading of the same region set faster?
		}
		res <- .object@regionSets[[i]]
		#non-standard genome names ("_chr" suffix for added "chr" prefix to standard genomes)
		if (genome(res)[1]==paste0(gen, "_chr")) {
			logger.info('Removed "_chr" suffix from genome name')
			genome(res) <- gen
		}
		return(res)
	}
)
#-------------------------------------------------------------------------------
#' length-methods
#'
#' Return the number of region sets in the DB
#'
#' @param x \code{\linkS4class{RegionSetDB}} object
#' @return integer indicating the number of region sets in the DB
#'
#' @rdname length-RegionSetDB-method
#' @docType methods
#' @aliases length,RegionSetDB-method
#' @author Fabian Mueller
#' @export
setMethod("length",
	signature(
		x="RegionSetDB"
	),
	function(
		x
	) {
		return(length(x@regionSets))
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("dbDir")) {
	setGeneric(
		"dbDir",
		function(.object) standardGeneric("dbDir"),
		signature=c(".object")
	)
}
if (!isGeneric("dbDir<-")) {
	setGeneric(
		"dbDir<-",
		function(.object, value) standardGeneric("dbDir<-"),
		signature=c(".object", "value")
	)
}
#' dbDir-methods
#'
#' Retrieve or set the directories were the object and region sets are stored
#'
#' @param .object \code{\linkS4class{RegionSetDB}} object
#'
#' @rdname dbDir-RegionSetDB-method
#' @docType methods
#' @aliases dbDir
#' @aliases dbDir,RegionSetDB-method
#' @author Fabian Mueller
#' @export
setMethod("dbDir", signature(.object="RegionSetDB"),
	function(.object){
		.object@dbDir
	}
)
#' @aliases dbDir
#' @aliases dbDir,RegionSetDB-method
#' @export
setReplaceMethod("dbDir", signature(.object="RegionSetDB"),
	function(.object, value){
		if (!is.null(value) && !is.character(value)){
			logger.error("Invalid setting for dbDir")
		}
		.object@dbDir <- value
		invisible(.object)
	}
)
################################################################################
# Adding and removing region sets
################################################################################
if (!isGeneric("addRegionSet")) {
	setGeneric(
		"addRegionSet",
		function(.object, ...) standardGeneric("addRegionSet"),
		signature=c(".object")
	)
}
#' addRegionSet-methods
#'
#' parses a region set and adds it to the \code{\linkS4class{RegionSetDB}} object
#'
#' @param .object \code{\linkS4class{RegionSetDB}} object
#' @param pid         a unique region set id used for the parser (character)
#' @param name        a name for the region set (character)
#' @param description a brief description of the region set (character)
#' @param assembly    the genome assembly of the region set (character)
#' @param collection  a collection or group this region set belongs to (character)
#' @param collection_description a brief description of the collection/group (character)
#' @param source      a reference or source from where the region set was obtained (character)
#' @param parser      name of the parser function used to import the region set (character)
#' @param file        the filename or input argument to parser function (character)
#' 
#' @return a modified \code{\linkS4class{RegionSetDB}} object with the region set added
#'
#' @rdname addRegionSet-RegionSetDB-method
#' @docType methods
#' @aliases addRegionSet
#' @aliases addRegionSet,RegionSetDB-method
#' @author Fabian Mueller
#' @export
setMethod("addRegionSet",
	signature(
		.object="RegionSetDB"
	),
	function(
		.object,
		pid,
		name,
		description,
		assembly,
		collection,
		collection_description,
		source,
		parser,
		file
	) {
		md <- data.frame(
			pid                    = pid,
			name                   = name,
			description            = description,
			assembly               = assembly,
			collection             = collection,
			collection_description = collection_description,
			source                 = source,
			parser                 = parser,
			file                   = file,
			stringsAsFactors = FALSE
		)

		argList <- list(fn=file, assembly=assembly, metadata=md)
		logger.start(c("Parsing annotation:", pid))
		parsedAnnot <- do.call(parser, argList)
		logger.completed()

		md <- rbind(.object@regionSetMd, parsedAnnot$metadata)

		#check if the 'name' in each genome-collection combination is unique
		chkList <-  lapply(genomes(.object), FUN=function(gg){
			lapply(collections(.object, genomeNames=gg), FUN=function(cc){
				rr <- regionSetNames(.object, genomeNames=gg, collectionNames=cc)
				if (any(duplicated(rr))){
					logger.error(c("Duplicate region set names in collection", cc, "for genome", gg, ":", paste(rr[duplicated(rr)], collapse=", ")))
				}
				return(rr)
			})
		})

		.object@parserSummary <- rbind(.object@parserSummary, md)
		.object@regionSets    <- c(.object@regionSets, parsedAnnot$regionSets)
		.object@regionSetMd   <- md
		.object@loadingStatus <- c(.object@loadingStatus, rep("loaded", nrow(md)))

		return(.object)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("addRegionSetFromDataFrame")) {
	setGeneric(
		"addRegionSetFromDataFrame",
		function(.object, ...) standardGeneric("addRegionSetFromDataFrame"),
		signature=c(".object")
	)
}
#' addRegionSetFromDataFrameFromDataFrame-methods
#'
#' parses multiple region sets and adds them to the \code{\linkS4class{RegionSetDB}} object. The parsing parameters are read from a data.frame
#'
#' @param .object \code{\linkS4class{RegionSetDB}} object
#' @param df      a data.frame object containing the required columns (see \code{\link{addRegionSet,RegionSetDB-metho}} for a brief description of the required columns)
#' 
#' @return a modified \code{\linkS4class{RegionSetDB}} object with the region sets added
#'
#' @rdname addRegionSetFromDataFrame-RegionSetDB-method
#' @docType methods
#' @aliases addRegionSetFromDataFrame
#' @aliases addRegionSetFromDataFrame,RegionSetDB-method
#' @author Fabian Mueller
#' @export
setMethod("addRegionSetFromDataFrame",
	signature(
		.object="RegionSetDB"
	),
	function(
		.object,
		df
	) {
		reqCols <- colnames(.object@parserSummary)

		missingCols <- setdiff(reqCols, colnames(df))
		if (length(missingCols)){
			logger.error(c("Missing columns in annotation summary table:", paste(missingCols, collapse=",")))
		}
		if (any(duplicated(df[,"pid"]))){
			logger.error(c("Duplicate ids in parser annotation"))
		}

		df <- df[, reqCols]

		for (i in 1:nrow(df)){
			.object <- addRegionSet(.object, pid=df[i, "pid"], name=df[i, "name"], description=df[i, "description"], assembly=df[i, "assembly"], collection=df[i, "collection"], collection_description=df[i, "collection_description"], source=df[i, "source"], parser=df[i, "parser"], file=df[i, "file"])
		}

		return(.object)
	}
)
################################################################################
# Saving and loading RegionSetDB objects
################################################################################
if (!isGeneric("unloadRegionSets")) {
	setGeneric(
		"unloadRegionSets",
		function(.object, ...) standardGeneric("unloadRegionSets"),
		signature=c(".object")
	)
}
#' unloadRegionSets-methods
#'
#' unload region sets from the \code{\linkS4class{RegionSetDB}} object
#'
#' @param .object \code{\linkS4class{RegionSetDB}} object
#' @param indx    indices of the region sets to be unloaded
#' @return a modified \code{\linkS4class{RegionSetDB}} object with the region sets unloaded
#'
#' @rdname unloadRegionSets-RegionSetDB-method
#' @docType methods
#' @aliases unloadRegionSets
#' @aliases unloadRegionSets,RegionSetDB-method
#' @author Fabian Mueller
#' @noRd
setMethod("unloadRegionSets",
	signature(
		.object="RegionSetDB"
	),
	function(
		.object,
		indx=NULL
	) {
		if (is.null(indx)) indx <- seq_along(.object)
		if (is.logical(indx)) indx <- which(indx)

		.object@loadingStatus[indx] <- "unloaded"
		.object@regionSets[indx] <- rep(list(NULL), length(indx))
		return(.object)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("loadRegionSets")) {
	setGeneric(
		"loadRegionSets",
		function(.object, ...) standardGeneric("loadRegionSets"),
		signature=c(".object")
	)
}
#' loadRegionSets-methods
#'
#' load region sets for the \code{\linkS4class{RegionSetDB}} object from disk
#'
#' @param .object \code{\linkS4class{RegionSetDB}} object
#' @param indx    indices of the region sets to be unloaded
#' @param genome  only load region sets for a specified genome assembly
#' @return a modified \code{\linkS4class{RegionSetDB}} object with the region sets loaded
#'
#' @rdname loadRegionSets-RegionSetDB-method
#' @docType methods
#' @aliases loadRegionSets
#' @aliases loadRegionSets,RegionSetDB-method
#' @author Fabian Mueller
#' @noRd
setMethod("loadRegionSets",
	signature(
		.object="RegionSetDB"
	),
	function(
		.object,
		indx=NULL,
		genome=NULL
	) {
		if (is.null(indx)) indx <- seq_along(.object)
		if (is.logical(indx)) indx <- which(indx)
		if (!is.null(genome)){
			indx.genome <- which(.object@regionSetMd[,"assembly"]==genome)
			indx <- intersect(indx, indx.genome)
		}

		#load region sets if they have not all been loaded
		if (any(.object@loadingStatus[indx] == "unloaded")){
			rsDir <- file.path(dbDir(.object), "region_sets")
			if (!dir.exists(rsDir)){
				logger.error(c("Cannot load stored region sets from", rsDir, "- directory does not exist"))
			}
		
			.object@loadingStatus[indx] <- "loaded"
			for (i in indx){
				.object@regionSets[[i]] <- readRDS(file.path(rsDir, paste0("rs", i, ".rds")))
			}
		}
		return(.object)
	}
)
#-------------------------------------------------------------------------------
#' saveRegionSetDB
#' 
#' Save a RegionSetDB dataset to disk for later loading
#' @param .object \code{\linkS4class{RegionSetDB}} object
#' @param path    destination to save the object to
#' @return nothing of particular interest
#' @author Fabian Mueller
#' @export
saveRegionSetDB <- function(.object, path=dbDir(.object)){
	if (dir.exists(path)){
		logger.error("Could not save object. Path already exists")
	}
	if (is.null(path)){
		logger.error("Could not save object. Specify a non-NULL path")
	}
	dir.create(path, recursive=FALSE)
	dbDir(.object) <- path

	rsDir <- file.path(path, "region_sets")
	dir.create(rsDir)
	for (i in seq_along(.object@regionSets)){
		saveRDS(.object@regionSets[[i]], file.path(rsDir, paste0("rs", i, ".rds")))
	}

	obj.red <- unloadRegionSets(.object)
	dsFn <- file.path(path, "rsdb.rds")
	saveRDS(obj.red, dsFn)
	invisible(NULL)
}
#-------------------------------------------------------------------------------
#' loadRegionSetDB
#' 
#' Load a RegionSetDB dataset from disk
#' @param path    Location of saved \code{\linkS4class{RegionSetDB}} object
#' @param loadRegionSets should the individual region sets be loaded into the object. If not this will be done as they are required in a lazy-load fashion
#' @return \code{\linkS4class{RegionSetDB}} object
#' @author Fabian Mueller
#' @export
loadRegionSetDB <- function(path, loadRegionSets=FALSE){
	if (!dir.exists(path)){
		logger.error(c("Could not load object. Path does not exist:", path))
	}
	dsFn <- file.path(path, "rsdb.rds")
	.object <- readRDS(dsFn)
	if (loadRegionSets){
		.object <- loadRegionSets(.object)
	}
	return(.object)
}
