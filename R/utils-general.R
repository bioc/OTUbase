
# Matt's readotu function (for Mothur otu files)
.readmothurotu <- function(ids, filename,level){
    if (!is.character(level)) stop("Level must be written as a character string")
    if (length(levels) > 1 ) stop("must specify only one level")
    otu_data <- strsplit(readLines(filename),split="\t")
    levels <- sapply(otu_data,function(x) x[1])
    otu_data <- otu_data[[match(level,levels)]]
    if (length(otu_data) == 0) stop("something must have gone wrong")

    nclust <- otu_data[2]
    otu_data<- otu_data[-c(1,2)]
    split_otu <- strsplit(otu_data,split=",")

    otu_lens <- sapply(split_otu,length)
    otu_id <- rep(1:nclust,times=otu_lens)
    otu_id <- paste("otu",otu_id,sep="")
    otu_obj <- unlist(split_otu)
    
    seqnames<-.getseqname(ids)
    ord <- match(seqnames,otu_obj)
    otulist<-otu_id[ord]
    otulist<-as.vector(sapply(otulist,function(d){if (is.na(d)){d<-"NA"}else{d<-d}}))
    return(as.character(otulist))    
}

# Function to read in cd-hit files
.readcdhitotu <- function(ids, filename){
    # read in data
    filedata<-scan(file=filename,strip.white=T, what=character(), sep="\n")

    # manipulate data into workable format
    d<-sapply(unlist(strsplit(paste(filedata,collapse="bbbaaarrr"), split="Cluster")), function (i) strsplit(i, split="bbbaaarrr"))
    d<-d[-1]
    names(d)<-1:length(d)
    lapply(d, function(i) i[-1])->d
    d[[length(d)]][length(d[[length(d)]])+1]<-">"
    d<-lapply(d, function(i) i[-length(i)])

    split_otu<-lapply(d, function(i){
        i<-sapply(i, function(j) unlist(strsplit(j, split=">|\\.\\.\\."))[2])
        names(i)<-character(length(i))
        i
    })
    nclust<-length(split_otu)
    # convert list of clusters into vector of cluster names
    otu_lens <- sapply(split_otu,length)
    otu_id <- rep(1:nclust,times=otu_lens)
    otu_id <- paste("otu",otu_id,sep="")
    otu_obj <- unlist(split_otu)
    seqnames<-.getseqname(ids)

    ord <- match(seqnames,otu_obj)
    otulist<-otu_id[ord]
    otulist<-as.vector(sapply(otulist,function(d){if (is.na(d)){d<-"NA"}else{d<-d}}))
    return(as.character(otulist))    
}

####################
# Function to extract subset of OTUset object given sample or otu names
subOTUset<-function(object,samples, otus){
	if (!missing(samples) & !missing(otus)){
		temp_object<-subOTUset(object,samples=samples)
		new_object<-subOTUset(temp_object,otus=otus)
		return(new_object)
	}else if (!missing(samples)){
		main_index<-!is.na(match(sampleID(object), samples))
	}else if (!missing(otus)){
		main_index<-!is.na(match(otuID(object), otus))
	}

	new_sampleID<-sampleID(object)[main_index]
	new_otuID<-otuID(object)[main_index]
	new_id<-id(object)[main_index]

	padf_c<-!is.na(match(sampleNames(sampleData(object)),unique(as.character(new_sampleID))))
	new_sampleData<-sampleData(object)[padf_c,]
	fadf_c<-!is.na(match(sampleNames(assignmentData(object)),unique(as.character(new_otuID))))
	new_assignmentData<-assignmentData(object)[fadf_c,]

	if (class(object)=="OTUsetB") new_object<-new("OTUsetB", sampleID=new_sampleID, otuID=new_otuID, id=new_id, sampleData=new_sampleData, assignmentData=new_assignmentData)
	if (class(object)=="OTUsetF") {
		new_sread<-sread(object)[main_index]
		new_object<-new("OTUsetF", sampleID=new_sampleID, otuID=new_otuID, id=new_id, sread=new_sread, sampleData=new_sampleData, assignmentData=new_assignmentData)
	}
	if (class(object)=="OTUsetQ") {
		new_sread<-sread(object)[main_index]
		new_qual<-quality(object)[main_index]
		new_object<-new("OTUsetQ", sampleID=new_sampleID, otuID=new_otuID, id=new_id, sread=new_sread, quality=new_qual, sampleData=new_sampleData, assignmentData=new_assignmentData)
	}
	return(new_object)
}	


#################################
#method to read in group file (mothur format) to fill slot group
.readsample<-function(ids,filename){
    groupfile<-readLines(filename)        #reads in group list
    grouplist<-unlist(strsplit(groupfile,"\t"))            #breaks it into seqname and group
    seqnames<-.getseqname(ids)
    sample<-grouplist[match(seqnames,grouplist)+1]    #puts group names in correct order in slot
    sample<-as.vector(sapply(sample,function(d){if (is.na(d)){d<-"NA"}else{d<-d}}))
    return(as.character((sample)))
}

################################
#function to check ADF for correct names (to link to other slots)
#currently if error the object is not created. Modify this to just not read in ADF.    
.readADF<-function(ADF,dirPath,ADFnames,ID){
    #reads in AnnotatedDataFrames
    if (class(ADF)=="AnnotatedDataFrame") adf<-ADF
    else if (is.character(ADF)){
        if (!missing(ADFnames)){
            adf<-read.AnnotatedDataFrame(ADF,dirPath,row.names=ADFnames)
        }else {    #look for column of data frame with sample names
            message("Looking for column with names")
            adf_table<-read.table(file.path(dirPath,ADF),sep="\t",head=T)
            # Match each column to sampleID to find matches
            tftable<-apply(adf_table,2,function(i){!is.na(match(as.character(i),as.character(ID)))})
            #number of matches in each column
            mcol<-apply(tftable,2,function(i){length(na.exclude(match(i,TRUE)))})
            if (length(mcol[mcol>0])==0) stop("Error: Could not find names in AnnotatedDataFrame")
            if (length(mcol[mcol>0])>1) stop("Error: Multiple possible name columns in AnnotatedDataFrame, please specify column number")
            if (length(mcol[mcol>0])==1) {
                ADFnames<-match(TRUE,mcol>0)
                adf<-read.AnnotatedDataFrame(ADF,dirPath,row.names=ADFnames)
            }
        }
    }
    if(exists("adf")) {
        #checks for matching of sample names in ADF and object
        if (!is.na(match(TRUE,(!is.na(match(sampleNames(adf),as.character(ID))))))) radf<-adf
        else stop("Sample names in AnnotatedDataFrame do not match ID")
    }
    return(radf)
}

	
#function to get sample names given column of sampleData and value
.adf_extract<-function(adf,colnum,value,exact){
	if (!exact){ #use grep to find all occurances
		c<-sapply(value,function(i){grep(pData(adf)[,colnum],pattern=i,ignore.case=T)})
	}else if (exact){ #use match to find exact matches
		c<-sapply(value, function(i){
			m<-match(pData(adf)[,colnum],i)
			(1:length(m))[!is.na(m)]})
	}
	i<-unique(as.numeric(unlist(c)))
	return(sampleNames(adf)[i])
}


###### functions to get sample names that have specific value in sampleData dataframe
getSamples<-function(object,colnum=1,value=".",exact=F){
	adf<-sampleData(object)
	return(.adf_extract(adf, colnum=colnum,value=value,exact=exact))
}

###### functions to get otu names that have specific value in featureData dataframe
getOTUs<-function(object,colnum=1,value=".",exact=F){
	adf<-assignmentData(object)
	return(.adf_extract(adf, colnum=colnum,value=value,exact=exact))
}


#combining functions into method that will accept directory and read in data
readOTUset<-function(dirPath = '.',otufile,level="0.03",fastafile,qualfile,samplefile,sampleADF,assignmentADF,sADF.names,aADF.names, rdp=F, otufiletype="mothur"){

    if(!missing(sampleADF)& missing(samplefile)) stop("Sample file is required to read in sampleADF")

    if (!missing(qualfile)& !missing(fastafile)){
        fasta<-readFasta(file.path(dirPath,fastafile))
        qualfile<-readQual(dirPath,sread(fasta),qualfile)
        object<-new("OTUsetQ",fasta,quality=qualfile)
    }else if (!missing(fastafile)){
        fasta<-readFasta(file.path(dirPath,fastafile))
        object<-new("OTUsetF",fasta)
    }else if (missing(qualfile) & missing(fastafile)){
        object<-new("OTUsetB")
        if (otufiletype=="mothur"){
            object@id<-.idfrommothurotu(file.path(dirPath,otufile),level)
        }
        else if (otufiletype=="cdhit"){
            object@id<-.idfromcdhitotu(file.path(dirPath, otufile))
        }
        else {
            stop("Error: OTU file type not recognized")
        }
    }else {
        stop("Error: must pass either Fasta, Fasta and Qual, or neither")
    }

    if (!missing(otufile)){
        if (otufiletype=="mothur"){
        	object@otuID <- .readmothurotu(id(object),file.path(dirPath,otufile),level)
	}
	else if (otufiletype=="cdhit"){
		object@otuID <- .readcdhitotu(id(object), file.path(dirPath,otufile))
	}
	else {
		stop("Error: OTU file type not recognized")
	}
        #object<-readotu(object,file.path(dirPath,otufile),level)
        if (!missing(samplefile)){
            object@sampleID <-.readsample(id(object),file.path(dirPath,samplefile))
        }
    }
    

    if(!missing(sampleADF)){
        object@sampleData<-.readADF(ADF=sampleADF,dirPath=dirPath,ADFnames=sADF.names,ID=sampleID(object))
    }
    
    # Add option to read in RDP classifier (fixed) data from OTUrep fasta directly

    if (!missing(assignmentADF)&(rdp==F)){
        object@assignmentData<-as(.readADF(ADF=assignmentADF,dirPath=dirPath,ADFnames=aADF.names,ID=otuID(object)), "AnnotatedDataFrame")
    }else if (!missing(assignmentADF)&(rdp==T)){
	 object@assignmentData<-as(.readRDP(ADF=assignmentADF, dirPath=dirPath, id=as.character(id(object)), otuid=otuID(object)), "AnnotatedDataFrame")
    }
    object
}

.readRDP<-function(ADF, dirPath, id, otuid){
    message("Note: Classification file must be RDP's fixed taxonomy.")
    adf<-read.table(file.path(dirPath,ADF), head=F, sep="\t", stringsAsFactors=F)

    ## get rid of NA columns in RDP output. I'm not sure why they are there.
    adf$V3<-"Root"
    adf$V4<-"norank"
    adf$V5<-"1.0"

    # strips off extra numbers on the sequence ID added by Mothur in versions before 1.13.0
    adf$V1<-sapply(strsplit(as.character(adf$V1),split="\\|"), function(i) i[1])

    ## get rid of useless "-" column
    adf[,-2]->adf

    ## label column names
    colnames(adf)<-c("otuID", "root", " ", "root_score", "domain", " ", "domain_score", "phylum", " ", "phylum_score", "class", " ", "class_score", "order", " ", "order_score", "family", " ", "family_score", "genus", " ", "genus_score")
    
    ## get rid of crap columns
    adf[,-c(3,6,9,12,15,18,21)]->adf

    ## match rdp seqnames with seqnames in object, then replace seq name with otuid
    if (!missing(otuid)) adf$otuID<-otuid[match(adf$otuID, id)]

    adf[,1]->row.names(adf)
    adf[,-1]->adf
    return(adf)
}
    
    

# extracts sequence id from mothur otufile (used when no fastafile is present)
.idfrommothurotu<-function(filename,level){
    if (!is.character(level)) stop("Level must be written as a character string")
    if (length(levels) > 1 ) stop("must specify only one level")
    otu_data <- strsplit(readLines(filename),split="\t")
    levels <- sapply(otu_data,function(x) x[1])
    otu_data <- otu_data[[match(level,levels)]]
    if (length(data) == 0) stop("something  must have gone wrong")
    
    otu_data<- otu_data[-c(1,2)]
    split_otu <- unique(unlist(strsplit(otu_data,split=",")))

    BStringSet(split_otu,use.names=F)
}

# extracts sequence id from cdhit otufile
.idfromcdhitotu<-function(filename){
    otu_data<-scan(file=filename,strip.white=T, what=character(), sep="\n")
    d<-sapply(unlist(strsplit(paste(otu_data,collapse="bbbaaarrr"), split="Cluster")), function (i) strsplit(i, split="bbbaaarrr"))
    d<-d[-1]
    names(d)<-1:length(d)
    d<-lapply(d, function(i) i[-1])
    d[[length(d)]][length(d[[length(d)]])+1]<-">"
    d<-lapply(d, function(i) i[-length(i)])
    otu_data<-lapply(d, function(i){
        i<-sapply(i, function(j) unlist(strsplit(j, split=">|\\.\\.\\."))[2])
        names(i)<-character(length(i))
        i
    })
    otu_data<-as.vector(unlist(otu_data))
    BStringSet(otu_data, use.names=F)
}

# gets seqname from id line
.getseqname<-function(ids){
    ids<-as.character(ids)
    seqnames<-as.vector(sapply(ids,function(i){unlist(strsplit(i,split=" "))[1]}))
    seqnames
}



.getBar <- function(width=getOption("width"))
  paste(rep("=", width), collapse="")


