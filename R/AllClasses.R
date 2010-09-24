## put classes here

## class not exported by ShortRead
setClass(".ShortReadBase")



## base class for both OTUset and RDPset
setClass(".OTUbase",
		representation(
		    sampleID="character", ## number of reads long
		    assignmentData="AnnotatedDataFrame", ## number of rows = number of assignments (ie OTUs)
		    sampleData="AnnotatedDataFrame") ## number of rows = number of samples
)
#base OTU class
setClass(".OTUset", 
		contains = c(".OTUbase"),
		representation(
            otuID="character") ## number of reads long
)

#OTU class including reads
setClass("OTUsetF",
    contains = c(".OTUset", "ShortRead")
)

#OTU class including quality scores
setClass("OTUsetQ",
    contains = c("ShortReadQ",".OTUset")	
)
#class without read or quality information
setClass("OTUsetB",
    contains = c(".OTUset"),
    representation(
        id="BStringSet")
)


## Base class for RDP type data
setClass(".TAXset",
	contains = c(".OTUbase"),
    representation(
        tax = "data.frame")
)

# RDP class without read or quality information
setClass("TAXsetB",
    contains = c(".TAXset"),
    representation(
        id="BStringSet")
)

# RDP class including reads
setClass("TAXsetF",
    contains = c(".TAXset", "ShortRead")
)

# RDP class including quality scores
setClass("TAXsetQ",
    contains = c(".TAXset", "ShortReadQ")
)
