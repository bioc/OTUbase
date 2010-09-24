# plots and summary functions

otusize<-function(object){
    y <- table(na.exclude(otuID(object)))
}

sharedotus<-function(object){
    y <- tapply(sampleID(object),otuID(object),function(x)length(unique(x)))
}

otuspersample<-function(object){
    y<-tapply(otuID(object),sampleID(object),function(x) length(unique(x)))
}

seqspersample<-function(object){
    y <- table(na.exclude(sampleID(object)))
}

# Calculate Richness
# Use functions from package "vegan"
o_estimateR<-function(object,...){
	r<-apply(abundance(object, weighted=F,...),2,function(i){estimateR(i,...)})
	r
}
o_diversity<-function(object, ...){
	r<-diversity(t(abundance(object, weighted=F,...)),...)
	r
}
##############################
# plots samples by number of sequences and estimated richness
rseqplot<-function(object,...){
    x<-seqspersample(object)
    y<-o_estimateR(object)[2,]
    c<-match(names(x),names(y))
    x<-as.numeric(x)
    y<-as.numeric(y[c])
    plot(x,y,ylab="Estimated Richness",xlab="Number Sequences",main="Samples plotted by number of sequences and estimated richness",...)
}




# plots samples by number of sequences and number of otus
otuseqplot<-function(object,...){
    x<-seqspersample(object)
    y<-otuspersample(object)
    c<-match(names(x),names(y))
    x<-as.numeric(x)
    y<-as.numeric(y[c])
    plot(x,y,ylab="Number OTUs",xlab="Number Sequences",main="Samples plotted by number of sequences and number of OTUs",...)
    
}


