## put generic functions here

setGeneric("otuID",function(object, ...) standardGeneric("otuID"))
setGeneric("sampleID",function(object, ...) standardGeneric("sampleID"))

setGeneric("nsamples",function(object, ...) standardGeneric("nsamples"))
setGeneric("notus",function(object, ...) standardGeneric("notus"))
setGeneric("seqnames", function(object, ...) standardGeneric("seqnames"))

setGeneric("aData", function(object, ...) standardGeneric("aData"))
setGeneric("sData", function(object, ...) standardGeneric("sData"))
setGeneric("aData<-", function(object, value) standardGeneric("aData<-"))
setGeneric("sData<-", function(object, value) standardGeneric("sData<-"))


setGeneric("sampleData", function(object, ...) standardGeneric("sampleData"))
setGeneric("sampleData<-", function(object, value) standardGeneric("sampleData<-"))

setGeneric("assignmentData", function(object, ...) standardGeneric("assignmentData"))
setGeneric("assignmentData<-", function(object, value) standardGeneric("assignmentData<-"))

setGeneric("abundance", function(object, ...) standardGeneric("abundance"))

setGeneric("tax",function(object, ...) standardGeneric("tax"))
setGeneric("tax<-", function(object, value) standardGeneric("tax<-"))

setGeneric("assignmentLabels", function(object, ...) standardGeneric("assignmentLabels"))
setGeneric("assignmentLabels<-", function(object, value) standardGeneric("assignmentLabels<-"))

setGeneric("sampleLabels", function(object, ...) standardGeneric("sampleLabels"))
setGeneric("sampleLabels<-", function(object, value) standardGeneric("sampleLabels<-"))

setGeneric("assignmentNames", function(object, ...) standardGeneric("assignmentNames"))


setGeneric("clusterSamples", function(object, ...) standardGeneric("clusterSamples"))
