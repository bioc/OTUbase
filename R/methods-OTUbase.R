
#### Access to assignmentData annotated data frame
### meant to be generic for both OTU and Classification
setMethod("assignmentData",
		signature=".OTUbase",
		function(object) slot(object, "assignmentData"))


setReplaceMethod("assignmentData",
                 signature=signature(
                   object=".OTUbase", value="AnnotatedDataFrame"),
                 function(object, value){
			object@assignmentData <- value
			object
                 })

setMethod("aData", ".OTUbase",
		function(object) pData(assignmentData(object)))


setReplaceMethod("aData",
                 signature=signature(
                   object=".OTUbase", value="data.frame"),
                 function(object, value){
 			pd <- assignmentData(object)
			pData(pd) <- value
			assignmentData(object) <- pd
			object
                 })



setMethod("assignmentLabels",
		signature=".OTUbase",
		function(object) varLabels(assignmentData(object)))

setReplaceMethod("assignmentLabels",
                 signature=signature(
                   object=".OTUbase", value="character"),
                 function(object, value) {
			pd <- assignmentData(object)
			varLabels(pd) <- value
			object@assignmentData <- pd
			object
                 })


############# End accessors for assignment data

#### Access to sampleData annotated data frame
setMethod(sampleData, ".OTUbase",
		function(object) slot(object, "sampleData"))


setReplaceMethod("sampleData",
                 signature=signature(
                   object=".OTUbase", value="AnnotatedDataFrame"),
                 function(object, value) {
			object@sampleData <- value
			object
                 })


setMethod(sData, ".OTUbase",
		function(object) pData(sampleData(object)))


setReplaceMethod("sData",
                 signature=signature(
                   object=".OTUbase", value="data.frame"),
                 function(object, value) {
			pd <- sampleData(object)
			pData(pd) <- value
			sampleData(object) <- pd
			object
                 })


setMethod("sampleLabels",
		signature=".OTUbase",
		function(object) varLabels(sampleData(object)))


setReplaceMethod("sampleLabels",
                 signature=signature(
                   object=".OTUbase", value="character"),
                 function(object, value) {
			pd <- sampleData(object)
			varLabels(pd) <- value
			object@sampleData <- pd
			object
                 })

############# END accessors for sampleData

setMethod(sampleID, ".OTUbase",
		function(object, ...) slot(object, "sampleID"))

setMethod(sampleNames, ".OTUbase",
		function(object) rownames(sData(object)))

setMethod(nsamples, ".OTUbase",
		function(object, ...) length(sampleNames(object)))

setMethod(assignmentNames, ".OTUbase",
		function(object, ...) rownames(aData(object)))

### allows for more information on the id line, pulls just the id
setMethod(seqnames, ".OTUbase",
		function(object){
			idline<-as.character(id(object))
			id<-1:length(idline)
			seqnames<-sapply(id,function(i){unlist(strsplit(idline[[i]],split=" "))[1]})	
			seqnames
		})
