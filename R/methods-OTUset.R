## put OTUset methods here
# Access slots
setMethod(id,"OTUsetB",
		function(object) slot(object, "id"))

setMethod(otuID,".OTUset",
		function(object) slot(object, "otuID"))


setMethod(notus, ".OTUset",
		function(object) length(unique(otuID(object))))


# show
setMethod(show, ".OTUset", 
    function(object) {
	cat("Class:", class(object),"\n")
	if ((class(object)=="OTUsetF")||(class(object)=="OTUsetQ")){
		cat("Number of Sequences:", length(sread(object)),"reads\n")
		cat("Sequence Width:", paste(range(width(sread(object))),collapse=".."),"cycles\n")
	}

        cat("Number of OTUs:",notus(object),"\n")
	cat("Number of Samples:",nsamples(object),"\n")
	if (length(sampleNames((object)))>0) {
		cat("sampleData: T")
		cat("\tncol:", ncol(sData(object)),"\n")
	}else cat("sampleData: F\n")
	if (length(assignmentNames(object))>0){ 
		cat("assignmentData: T")
		cat("\tncol:", ncol(aData(object)),"\n")
	}else cat("assignmentData: F\n")

})


## summary information
setMethod(abundance, ".OTUset",
    function(object,assignmentCol, sampleCol, weighted=F,collab,...){
        s<-sampleID(object)


	## make list of OTUs or assignments
        if (!missing(assignmentCol)){
            g<-aData(object)[,assignmentCol]
            c<-match(otuID(object),row.names(aData(object)))
            o<-g[c]
        }else o<-otuID(object)

	## make list of samples or sample characteristics
	if (!missing(sampleCol)){
	    a<-sData(object)[,sampleCol]
	    t<-match(sampleID(object), row.names(sData(object)))
	    s<-a[t]
	}else s<-sampleID(object)

        abund<-table(o,s)

        if (weighted){abund<-apply(abund,2,function(j) j/sum(j))}

        if (!missing(collab)) colnames(abund)<-sData(object)[,collab][match(colnames(abund),row.names(sData(object)))]

        return(abund)
    }
)


setMethod(clusterSamples, ".OTUset",
		function(object, assignmentCol, collab, distmethod='bray', clustermethod='complete',...){
			if (!missing(collab) && !missing(assignmentCol)){
				a<-t(abundance(object, collab=collab, weighted=T,assignmentCol=assignmentCol))
			}else if (!missing(collab)){
				a<-t(abundance(object, collab=collab, weighted=T))
			}else if (!missing(assignmentCol)){
				a<-t(abundance(object, weighted=T, assignmentCol=assignmentCol))
			}else a<-t(abundance(object, weighted=T))
			
			d<-vegdist(a, method=distmethod)
			rn<-row.names(as.matrix(d))
			clust<-hclust(d,method=clustermethod)
			plot(clust,labels=rn, sub=NA, xlab=NA,...)
			return(clust)
		}
)

