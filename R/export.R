
if (!isGeneric("export2Lola")) {
	setGeneric(
		"export2Lola",
		function(.object, ...) standardGeneric("export2Lola"),
		signature=c(".object")
	)
}
#' export2Lola-methods
#'
#' Export a \code{\linkS4class{RegionSetDB}} object as LOLA database
#'
#' @param .object \code{\linkS4class{RegionSetDB}} object
#' @param destDir destination directory
#' @param doIndex index the LOLA database for faster re-loading
#' @return nothing of particular interest
#'
#' @rdname export2Lola-RegionSetDB-method
#' @docType methods
#' @aliases export2Lola
#' @aliases export2Lola,RegionSetDB-method
#' @author Fabian Mueller
#' @noRd
setMethod("export2Lola",
	signature(
		.object="RegionSetDB"
	),
	function(
		.object,
		destDir,
		doIndex=TRUE
	) {
		require(LOLA)
		# map from summary cols to LOLA cols
		LOLA.INDEX.MAP.COLS <- c(
			name = "description",
			description = "description_long",
			source = "dataSource",
			cellType = "cellType",
			tissue = "tissue",
			antibody = "antibody",
			treatment = "treatment"
		)
		#reverse map for colelction info
		LOLA.COLS.TO.SUMMARY.COLS.COLLECTION <- c(
			description = "collection_description",
			source = "source"
		)

		.object <- loadRegionSets(.object)
		annotSummary <- .object@regionSetMd

		if (!dir.exists(destDir)){
			dir.create(destDir)
		} else {
			logger.error(c("LOLA directory already exists:", destDir))
		}

		assemblies <- unique(annotSummary[,"assembly"])

		logger.status(c("Preparing directory structure"))
		collections <- unique(annotSummary[,"collection"])
		for (cc in collections){
			curAss <- unique(annotSummary[annotSummary[,"collection"]==cc,"assembly"])
			for (assembly in curAss){
				assDir <- file.path(destDir, assembly)
				if (!dir.exists(assDir)){
					dir.create(assDir)
				}
				colDir <- file.path(assDir, cc)
				dir.create(colDir)

				x <- annotSummary[annotSummary[,"collection"]==cc & annotSummary[,"assembly"]==assembly, LOLA.COLS.TO.SUMMARY.COLS.COLLECTION[c("source", "description")]]
				if (nrow(x) > 1){
					if (sum(duplicated(x))!=(nrow(x)-1)){
						logger.error(c("Collection annotation", cc, "for assembly", assembly, "not uniquely defined"))
					}
					x <- x[1,,drop=FALSE]
				}
				colInfo <- data.frame(
					collector = "createAnnotationDB.R",
					date = format(Sys.time(), "%Y-%m-%d"),
					source = x[[LOLA.COLS.TO.SUMMARY.COLS.COLLECTION["source"]]],
					description = x[[LOLA.COLS.TO.SUMMARY.COLS.COLLECTION["description"]]],
					stringsAsFactors = FALSE
				)
				fn <- file.path(colDir, "collection.txt")
				write.table(colInfo, fn, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
			}
		}

		#map column names to LOLA index column names
		cns.mapped <- LOLA.INDEX.MAP.COLS[colnames(annotSummary)]
		cns.mapped.notna <- !is.na(cns.mapped)
		annotSummary.lola <- annotSummary
		annotSummary.lola[["filename"]] <- NA
		#map column names to LOLA index column names			
		if (any(cns.mapped.notna)){
			colnames(annotSummary.lola)[cns.mapped.notna] <- cns.mapped[cns.mapped.notna]
		}

		cns.lola <- intersect(c("filename", LOLA.INDEX.MAP.COLS), colnames(annotSummary.lola))
		cns.other <- setdiff(colnames(annotSummary.lola), LOLA.INDEX.MAP.COLS)
		#reorder columns, s.t. LOLA relevant information comes first
		annotSummary.lola <- annotSummary.lola[,c(cns.lola, cns.other)]

		for (i in seq_along(.object)){
			regionSetsInfo <- annotSummary.lola[i,,drop=FALSE]
			Nrs <- nrow(regionSetsInfo)

			regionSetsInfo[["filename"]] <- paste0("regionSet_", 1:Nrs, ".bed")			

			curAss <- regionSetsInfo[1,"assembly"]
			curCol <- regionSetsInfo[1,"collection"]
			colDir <- file.path(destDir, curAss, curCol)
			regDir <- file.path(colDir, "regions")

			fn <- file.path(colDir, "index.txt")
			# check if there is already an index file. If so, multiple annotations are mapped to the same collection
			# get continuing running numbers for file names
			doAdd <- file.exists(fn)
			if (doAdd){
				oldIndexTab <- read.table(fn, sep="\t", header=TRUE)
				maxFnInd <- max(as.integer(sub("^regionSet_([0-9]+).bed", "\\1", oldIndexTab[,"filename"])))
				indsOffset <- (maxFnInd+1):(maxFnInd+Nrs)
				regionSetsInfo[["filename"]] <- paste0("regionSet_", indsOffset, ".bed")
				write.table(regionSetsInfo, fn, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
			} else {
				write.table(regionSetsInfo, fn, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
				dir.create(regDir)
			}
			for (j in 1:Nrs){
				granges2bed(.object@regionSets[[i]], file.path(regDir, regionSetsInfo[j, "filename"]), score=NULL, addAnnotCols=TRUE, colNames=FALSE, doSort=TRUE)
			}
		}

		regDbs <- NULL
		if (doIndex){
			require(LOLA)
			regDbs <- list()
			for (aa in assemblies){
				assDir <- file.path(destDir, assembly)
				regDb <- loadRegionDB(assDir)
				regDbs <- c(regDbs, list(regDb))
			}
			names(regDbs) <- assemblies
		}
		invisible(regDbs)
	}
)
