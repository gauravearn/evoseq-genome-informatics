# Here's a good place to put your top-level package documentation

.onAttach <- function (lib, pkgname="evoseq") {
    ## Put stuff here you want to run when your package is loaded
    genecoding <- function(genomealigned, ids, strand) {
    library(dplyr)
    library(stringr)
    # Author Gaurav Sablok
    # Universitat Potsdam
    # Date: 2024-5-1
    # a gene virualizer for the genome annotation coming from the proteome
    # alignments. You can visualize the specific genes, strands or you can visualize
    # all the annotations at the specific level. will prepare the file for the all 
    # the coding regions for the visualization. 
    system("cat *.gff | grep -v ^# > normalizedfilemRNA.txt")
    readfile <- paste(getwd(), "normalizedfilemRNA.txt", sep = "/")
    fileinput <- read.csv(readfile, sep = "\t")
    sectionedfile <- fileinput[c(1,2,3,4,5,6,7)]
    filteredfile <- dplyr::filter(sectionedfile, sectionedfile[3] == "CDS")
    colnames(filteredfile) <- c("col1", "col2", "col3", "col4", "col5", "col6", "col7")
    inputids <- ids
    if (strand == "positive" || ids) {
        posfilterids <- dplyr::filter(filteredfile, filteredfile["col1"] == inputids)
        posfilterstrand <- dplyr::filter(posfilterids, posfilterids["col7"] == "+")
        idvector <- vector(length = length(posfilterids$col1))
        idtype <- vector(length=length(posfilterids$col1))
        idstart <- vector(length=length(posfilterids$col1))
        idend <- vector(length=length(posfilterids$col1))
        idorientation <- vector(length=length(posfilterids$col1))
        newids <- rep(inputids,length(posfilterids$col1))
        for (i in seq_along(posfilterstrand[1])){
                 idvector <- posfilterstrand[1][i] }
        for (i in seq_along(posfilterstrand[1])){
                   idtype <- posfilterstrand[2][i] }
        for (i in seq_along(posfilterstrand[4])){
                   idstart <- posfilterstrand[4][i]}
        for (i in seq_along(posfilterstrand[5])){
                   idend <- postfilterstrand[5][i]}
        for (i in seq_along(postfilterstrand[7])){
             idorientation <- postfilterstrand[7][i]}
       idtype <- str_replace_all(idtype, "miniprot", "mRNA")    
       data <- cbind(newids, idtype, idstart, idend, idorientation)
       colnames(data)
       colnames(data) <- c("newids", "idtype", "idstart", "idend", "strand")
       data["orientation"] <- data["strand"]
       updatestrand <- as.list(ifelse(data[,"strand"] == "+", "forward"))
       updateorientation <- as.list(ifelse(data[,"strand"] == "+", "1"))
       finaldata <- cbind(data, updatestrand, updateorientation)
       colnames(finaldata) <- c("molecule", "gene", "start", "end", "strand", "orientation", "strand", "orientation")
       return(finaldata[c(1,2,3,4,7,8)]) 
    } 
    if (strand == "negative" || ids) {
        negfilterids <- dplyr::filter(filteredfile, filteredfile["col1"] == inputids)
        negfilterstrand <- dplyr::filter(filterids, negfilterids["col7"] == "-")
    idvector <- vector(length = length(posfilterids$col1))
        idtype <- vector(length=length(posfilterids$col1))
        idstart <- vector(length=length(posfilterids$col1))
        idend <- vector(length=length(posfilterids$col1))
        idorientation <- vector(length=length(posfilterids$col1))
        newids <- rep(inputids,length(posfilterids$col1))
        for (i in seq_along(posfilterstrand[1])){
                 idvector <- posfilterstrand[1][i] }
        for (i in seq_along(posfilterstrand[1])){
                   idtype <- posfilterstrand[2][i] }
        for (i in seq_along(posfilterstrand[4])){
                   idstart <- posfilterstrand[4][i]}
        for (i in seq_along(posfilterstrand[5])){
                   idend <- postfilterstrand[5][i]}
        for (i in seq_along(postfilterstrand[7])){
             idorientation <- postfilterstrand[7][i]}
       idtype <- str_replace_all(idtype, "miniprot", "mRNA")    
       data <- cbind(newids, idtype, idstart, idend, idorientation)
       colnames(data)
       colnames(data) <- c("newids", "idtype", "idstart", "idend", "strand")
       data["orientation"] <- data["strand"]
       updatestrand <- as.list(ifelse(data[,"strand"] == "+", "forward"))
       updateorientation <- as.list(ifelse(data[,"strand"] == "+", "1"))
       finaldata <- cbind(data, updatestrand, updateorientation)
       colnames(finaldata) <- c("molecule", "gene", "start", "end", "strand", "orientation", "strand", "orientation") 
       return(finaldata[c(1,2,3,4,7,8)])
    }  
}   
}
