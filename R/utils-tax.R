## read in taxonomic data and create object of class. Currently only reads in rdp fixed formated data
readTAXset<-function(dirPath=".", taxfile, namefile, fastafile, qualfile, samplefile, sampleADF, assignmentADF, sADF.names, aADF.names, type='rdpfix',...){
    
    if(!missing(sampleADF)& missing(samplefile)) stop("Sample file is required to read in sampleADF")
    
    ## create object, read in fasta and qual if available

    if (!missing(qualfile)& !missing(fastafile)){
        fasta<-readFasta(file.path(dirPath,fastafile))
        qualfile<-readQual(dirPath,sread(fasta),qualfile)
        object<-new("TAXsetQ",fasta,quality=qualfile)
    }else if (!missing(fastafile)){
        fasta<-readFasta(file.path(dirPath,fastafile))
        object<-new("TAXsetF",fasta)
    }else if (missing(qualfile) & missing(fastafile)){
        object<-new("TAXsetB")
    }else {
        stop("Error: must pass either Fasta, Fasta and Qual, or neither")
    }

    ## read in rdpfile if available
    ## Can add functions to read in other file types here.
    
    if (type=="rdpfix") preTax<-.readrdpfix(dirPath=dirPath, rdpfile=taxfile)
    else preTax<-.readunknowntax(dirPath=dirPath, taxfile=taxfile)
    
    ## expand preTax (if namesfile is given)

    if (!missing(namefile)){
        names<-read.table(file.path(dirPath,namefile), head=F, stringsAsFactors=F, sep="\t")
        completeTax<-.addUnique(tax=preTax, names=names)
    }else { 
        message("Namefile not given. All sequences assumed to be in taxfile.")
        completeTax<-preTax
    }


    ## fill tax slot
    ######## BUG ####### If fastafile has fewer sequences then are in the taxfile, the extra taxfile names will be left out.
    if (missing(fastafile)){ 
        object@id<-BStringSet(row.names(completeTax))
        row.names(completeTax)<-NULL
        object@tax<-completeTax
    }else {
        message("Warning: If fastafile has fewer sequences then are in the taxfile, the extra taxfile names will be left out.")
        completeTax<-completeTax[match(as.character(id(object)),row.names(completeTax)),]
        row.names(completeTax)<-NULL
        object@tax<-completeTax
    }

    ## read in sampleID

    if (!missing(samplefile)){
        object@sampleID <-.readsample(id(object),file.path(dirPath,samplefile))
    }

    ## read in sampleData

    if(!missing(sampleADF)){
        object@sampleData<-.readADF(ADF=sampleADF,dirPath=dirPath,ADFnames=sADF.names,ID=sampleID(object))
    }

    ### figure out how to read in taxData. What should be required? Nothing. Reading in as simple annotated Data Frame
    if(!missing(assignmentADF)){
        object@assignmentData<-read.AnnotatedDataFrame(assignmentADF,dirPath,row.names=aADF.names)
    }
    return(object)

}


## function to add unique sequences back into RDP output. Requires fixed RDP ouput processed by .readRDP() function and names file in mothur format.

.addUnique<-function(tax, names){
    ## get only names present in taxonomy
    names<-names[match(row.names(tax), names$V1),]

    parsed <- strsplit(names[,2],",",fixed=T)
    newdata <- cbind(rep(names[,1],times=sapply(parsed,length)),unlist(parsed))
    fulltax <- cbind(newdata[,2],tax[match(newdata[,1],row.names(tax)),])
    row.names(fulltax)<-fulltax[,1]
    fulltax<-fulltax[,-1]
    return(fulltax)
}



## function to read in rdp fixed taxonomy files. Returns dataframe with rownames as sequence ids

.readrdpfix<-function(dirPath, rdpfile){
    ids<-read.table(file.path(dirPath, rdpfile), head=F, stringsAsFactors=F, sep="\t")$V1
    preTax<-.readRDP(id=ids, dirPath=dirPath, ADF=rdpfile)
    return(preTax)
}

.readunknowntax<-function(dirPath, taxfile,...){
        message("taxfile is assumed to be readable with read.table(). First column must contain sequence names.")
        message("Use carefully. An odd file type may lead to warningless errors.")
        preTax<-as.data.frame(read.table(taxfile, stringsAsFactors=F,...))
        row.names(preTax)<-preTax[,1]
        preTax<-cbind(preTax, "")
        preTax<-preTax[,-1]
	return(preTax)
}


## can only subset by samples. Need to add subsetting by assignment
subTAXset<-function(object,samples){
	main_index<-!is.na(match(sampleID(object), samples))

	new_sampleID<-sampleID(object)[main_index]
	new_tax<-tax(object)[main_index,]
	new_id<-id(object)[main_index]

	padf_c<-!is.na(match(sampleNames(sampleData(object)),unique(as.character(new_sampleID))))
	new_sampleData<-sampleData(object)[padf_c,]

	## not changeing assignmentData object for TAXset subsets.
	new_assignmentData<-assignmentData(object)

	if (class(object)=="TAXsetB") new_object<-new("TAXsetB", sampleID=new_sampleID, tax=new_tax, id=new_id, sampleData=new_sampleData, assignmentData=new_assignmentData)
	if (class(object)=="TAXsetF") {
		new_sread<-sread(object)[main_index]
		new_object<-new("TAXsetF", sampleID=new_sampleID, tax=new_tax, id=new_id, sread=new_sread, sampleData=new_sampleData, assignmentData=new_assignmentData)
	}
	if (class(object)=="TAXsetQ") {
		new_sread<-sread(object)[main_index]
		new_qual<-quality(object)[main_index]
		new_object<-new("TAXsetQ", sampleID=new_sampleID, tax=new_tax, id=new_id, sread=new_sread, quality=new_qual, sampleData=new_sampleData, assignmentData=new_assignmentData)
	}
	return(new_object)
}	



