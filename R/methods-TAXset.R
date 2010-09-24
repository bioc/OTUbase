### RDP
setMethod(tax, ".TAXset",
		function(object) slot(object, "tax"))


setReplaceMethod("tax",
                 signature=signature(
                   object=".TAXset", value="data.frame"),
                 function(object, value) {
			object@tax <- value
			object
                 })

setMethod(id,"TAXsetB",
		function(object) slot(object, "id"))


## show
setMethod(show, ".TAXset", 
		function(object) {
			cat("Class:", class(object),"\n")
			if ((class(object)=="TAXsetF")||(class(object)=="TAXsetQ")){
				cat("Number of Sequences:", length(sread(object)),"reads\n")
				cat("Sequence Width:", paste(range(width(sread(object))),collapse=".."),"cycles\n")}
			
			cat("Number of Samples:",nsamples(object),"\n")
			if (length(sampleNames(sampleData(object)))>0) {
				cat("sampleData: T")
				cat("\tncol:", ncol(sData(object)),"\n")
			}else cat("sampleData: F\n")
			if (length(sampleNames(assignmentData(object)))>0){ 
				cat("assignmentData: T")
				cat("\tncol:", ncol(aData(object)),"\n")
			}
			else cat("assignmentData: F\n")
			
		})


setMethod(abundance, ".TAXset",
		function(object,taxCol, weighted=F, collab, sampleCol, ...){
			if (missing(taxCol)) stop("Must provide tax column for abundance.")

			## make list of samples or sample characteristics
			if (!missing(sampleCol)){
			    a<-sData(object)[,sampleCol]
			    t<-match(sampleID(object), row.names(sData(object)))
			    s<-a[t]
			}else s<-sampleID(object)


			o<-as.character(tax(object)[,taxCol])
			
			abund<-table(o,s)
			
			if (weighted){abund<-apply(abund,2,function(j) j/sum(j))}
			
			if (!missing(collab)) colnames(abund)<-sData(object)[,collab][match(colnames(abund),row.names(sData(object)))]
			
			return(abund)
		}
)


setMethod(clusterSamples, ".TAXset",
		function(object, taxCol, collab, distmethod='bray', clustermethod='complete', assignmentCol,...){
			if (missing(taxCol)) stop("Missing taxCol")
			if (!missing(collab) && !missing(assignmentCol)){
				## eventually change to allow assignmentCol. Currently does nothing new.
				a<-t(abundance(object, collab=collab, taxCol=taxCol, weighted=T))
			}else if (!missing(collab)){
				a<-t(abundance(object, collab=collab, taxCol=taxCol, weighted=T))
			}else if (!missing(assignmentCol)){
				a<-t(abundance(object, taxCol=taxCol, weighted=T))
			}else a<-t(abundance(object, taxCol=taxCol, weighted=T))
			
			d<-vegdist(a, method=distmethod)
			rn<-row.names(as.matrix(d))
			clust<-hclust(d,method=clustermethod)
			plot(clust,labels=rn, sub=NA, xlab=NA,...)
			return(clust)
		}
)



