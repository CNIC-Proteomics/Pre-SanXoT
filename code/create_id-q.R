library("optparse")

# get input parameters
option_list = list(
  make_option(c("-i", "--indir"), type="character", default=NULL,
              help="dataset directory", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="Fr_*", 
              help="Pattern of folders that contains the MSFs [default= %default]", metavar="character"),
  make_option(c("-s", "--search"), type="integer", default=2,
              help="Search engine [default= %default]", metavar="character"),
  make_option(c("-d", "--daemon"), type="logical", default=TRUE,
              help="Daemon used [default= %default]", metavar="character"),
  make_option(c("-c", "--channels"), type="integer", default=NULL,
              help="Channel IDs used in the experiment. Eg. '2,5:8,10,12:15'", metavar="character"),
  make_option(c("-l", "--label"), type="character", default="TMT", 
              help="Type of label used [default= %default]", metavar="character"),
  make_option(c("-t", "--tags"), type="character", default="ALL", 
              help="Tags Used in the Experiment [default= %default]", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="output directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

opt$indir <- "D:/projects/Pre-SanXoT/pratio_FractionsSeparated/Pre-SanXoT_2"
opt$pattern <- "Fr_*"
opt$search <- 2
opt$daemon <- TRUE
opt$label <- "TMT"
opt$tags <- "ALL"
opt$channels <- source(textConnection( paste("c(",gsub("\\-", ":", "1:10"),")")) )$value 
opt$outdir <- "D:/projects/Pre-SanXoT/pratio_FractionsSeparated/Pre-SanXoT_2"

if ( is.null(opt$indir) || is.null(opt$pattern) || is.null(opt$search) || is.null(opt$daemon) || is.null(opt$channels) || is.null(opt$outdir) ){
  print_help(opt_parser)
  stop("All arguments must be supplied.n", call.=FALSE)
}

# create workspace
dir.create(opt$outdir, showWarnings = FALSE)

print( "input parameters: ")
print( opt )


files <- list.files(path = opt$indir,pattern="_results", full.names = TRUE)
ID_all <- do.call("rbind", lapply(files, read.table, header = TRUE))
write.table(ID_all, file = paste(opt$outdir,"/ID-all.txt",sep=""), sep="\t", row.names = FALSE)

# Pattern of folders that contains the MSFs
Patern=c("Fr_*")

# Channels used in the Experiments
ChannelID=c(1:10)

# Type of label used
Typeoflabel=c("TMT")

# Tags Used in the Experiment (All is "ALL")
TagsUsed=c("ALL")

# Control Tag
ControlTag=c("131")

# Mean Tag Calculation
MeanCalculation=c("TRUE")

# Mean Tags
MeanTags=c("130_C","131","129_C","129_N")

# First Tag
FirstTag=c("126")


# Number of comparatives within the Experiment
Comparatives=c("10")

# To Absolute Quantification (TRUE = Absolute Quantification, FALSE = Relative Quantification or BOTH = Both)
Absolute=c("BOTH")

# if (length(Expto)<2) {
  
  k<-q_all
  x<-ID_all
  x$Raw_FirstScan<-do.call(paste, c(x[c("RAWFile","FirstScan")], sep = ""))
  k$Raw_FirstScan<-do.call(paste, c(k[c("FileName","FirstScan")], sep = ""))
  x$Raw_FirstScan<-as.character(x$Raw_FirstScan)
  k$Raw_FirstScan<-as.character(k$Raw_FirstScan)
  
  if (Typeoflabel=="TMT"){
    if (TagsUsed=="ALL"){
      q<-k[,c("Raw_FirstScan","X126","X127_N","X127_C","X128_N","X128_C","X129_N","X129_C","X130_N","X130_C","X131")]
    } else {
      q<-k[,colnames_TMT]
    }
  }else{
    
    if (TagsUsed=="ALL"){
      q<-k[,c("Raw_FirstScan","X113","X114","X115","X116","X117","X118","X119","X121")]
    } else {
      q<-k[,colnames_iTRAQ]
    }
  }
  
  if (MeanTags=="ALL"){
    if (Typeoflabel == "TMT") {
      MeanTags<-c("X126","X127_N","X127_C","X128_N","X128_C","X129_N","X129_C","X130_N","X130_C","X131")
    } else {
      MeanTags<-c("X113","X114","X115","X116","X117","X118","X119","X121")
    }
  }

  all<-merge(x,q)
  FirstTagIndex=as.numeric(grep(paste0("X",FirstTag), colnames(all)))
  CalcIndex=trunc(seq(FirstTagIndex, by=(length(ChannelID)/as.numeric(Comparatives)), len = as.numeric(Comparatives)),1)
  
  for (i in CalcIndex){
    
    ControlIndex=as.numeric(grep(paste0("X",ControlTag), colnames(all)))
    
    if (MeanCalculation == "TRUE") {
      all$Mean <- rowMeans(all[,paste0("X",MeanTags)])
      MeanIndex=as.numeric(grep("Mean", colnames(all)))
      all$newcolumn <- log2(all[,i]/all$Mean)
      l <- substring(colnames(all)[i],2)
      colnames(all)[ncol(all)] <- paste0("Xs_",l,"_Mean")
    } else {
      all$newcolumn <- log2(all[,i]/all[,ControlIndex])
      l <- substring(colnames(all)[i],2)
      colnames(all)[ncol(all)] <- paste0("Xs_",l,"_",ControlTag)
    }
    
    if (Absolute == "TRUE") {
      all$newcolumn <- all[,c(i)]
      colnames(all)[ncol(all)] <- paste0("Vs_",l,"_ABS")
    }
    
    if (Absolute == "FALSE") {
      if (MeanCalculation == "TRUE") {
        all$newcolumn <- apply(all[,c(i,MeanIndex)], 1, max)
        colnames(all)[ncol(all)] <- paste0("Vs_",l,"_Mean")
      } else {
        all$newcolumn <- apply(all[,c(i,ControlIndex)], 1, max)
        colnames(all)[ncol(all)] <- paste0("Vs_",l,"_",ControlTag)
      }
    }
    
    if (Absolute == "BOTH") {
      all$newcolumn <- all[,c(i)]
      colnames(all)[ncol(all)] <- paste0("Vs_",l,"_ABS")
      if (MeanCalculation == "TRUE") {
        all$newcolumn <- apply(all[,c(i,MeanIndex)], 1, max)
        colnames(all)[ncol(all)] <- paste0("Vs_",l,"_Mean")
      } else {
        all$newcolumn <- apply(all[,c(i,ControlIndex)], 1, max)
        colnames(all)[ncol(all)] <- paste0("Vs_",l,"_",ControlTag)
      }
    }
  }
  write.table(all, file = paste(WD,"/",j,"/Pre-SanXoT/ID-q.txt",sep=""), sep="\t", row.names = FALSE)
  
# } else {
#   
#   for (j in Expto){
#     
#     files <- list.files(path = paste(WD,"/",j,"/Pre-SanXoT",sep=""),pattern="Q-all.xls", full.names = TRUE)
#     
#     k<-read.table(files, header=TRUE, sep=",")
#     
#     files <- list.files(path = paste(WD,"/",j,"/Pre-SanXoT",sep=""),pattern="ID-all.txt", full.names = TRUE)
#     
#     x<-read.table(files, header=TRUE, sep="\t")
#     
#     x$Raw_FirstScan<-do.call(paste, c(x[c("RAWFile","FirstScan")], sep = ""))
#     k$Raw_FirstScan<-do.call(paste, c(k[c("FileName","FirstScan")], sep = ""))
#     x$Raw_FirstScan<-as.character(x$Raw_FirstScan)
#     k$Raw_FirstScan<-as.character(k$Raw_FirstScan)
#     
#     if (Typeoflabel=="TMT"){
#       
#       if (TagsUsed=="ALL"){
#         
#         q<-k[,c("Raw_FirstScan","X126","X127_N","X127_C","X128_N","X128_C","X129_N","X129_C","X130_N","X130_C","X131")]
#         
#       } else {
#         
#         q<-k[,colnames_TMT]}
#       
#     }else{
#       
#       if (TagsUsed=="ALL"){
#         
#         q<-k[,c("Raw_FirstScan","X113","X114","X115","X116","X117","X118","X119","X121")]
#         
#       } else {
#         
#         q<-k[,colnames_iTRAQ]}}
#     
#     all<-merge(x,q)
#     
#     FirstTagIndex=as.numeric(grep(paste0("X",FirstTag), colnames(all)))
#     
#     CalcIndex=trunc(seq(FirstTagIndex, by=(length(ChannelID)/as.numeric(Comparatives)), len = as.numeric(Comparatives)),1)
#     
#     if (MeanTags=="ALL"){
#       
#       if (Typeoflabel == "TMT"){
#         
#         MeanTags<-c("X126","X127_N","X127_C","X128_N","X128_C","X129_N","X129_C","X130_N","X130_C","X131")
#         
#       } else {
#         
#         MeanTags<-c("X113","X114","X115","X116","X117","X118","X119","X121")}}
#     
#     for (i in CalcIndex){
#       
#       ControlIndex=as.numeric(grep(paste0("X",ControlTag), colnames(all)))
#       
#       if (MeanCalculation == "TRUE") {
#         
#         all$Mean <- rowMeans(all[,paste0("X",MeanTags)])
#         
#         MeanIndex=as.numeric(grep("Mean", colnames(all)))
#         
#         all$newcolumn <- log2(all[,i]/all$Mean)
#         
#         l <- substring(colnames(all)[i],2)
#         
#         colnames(all)[ncol(all)] <- paste0("Xs_",l,"_Mean")
#         
#       } else {
#         
#         all$newcolumn <- log2(all[,i]/all[,ControlIndex])
#         
#         l <- substring(colnames(all)[i],2)
#         
#         colnames(all)[ncol(all)] <- paste0("Xs_",l,"_",ControlTag)}
#       
#       if (Absolute == "TRUE"){
#         
#         all$newcolumn <- all[,c(i)]
#         
#         colnames(all)[ncol(all)] <- paste0("Vs_",l,"_ABS")}
#       
#       if (Absolute == "FALSE"){
#         
#         if (MeanCalculation == "TRUE"){
#           
#           all$newcolumn <- apply(all[,c(i,MeanIndex)], 1, max)
#           
#           colnames(all)[ncol(all)] <- paste0("Vs_",l,"_Mean")
#           
#         } else {
#           
#           all$newcolumn <- apply(all[,c(i,ControlIndex)], 1, max)
#           
#           colnames(all)[ncol(all)] <- paste0("Vs_",l,"_",ControlTag)}}
#       
#       if (Absolute == "BOTH"){
#         
#         all$newcolumn <- all[,c(i)]
#         
#         colnames(all)[ncol(all)] <- paste0("Vs_",l,"_ABS")
#         
#         if (MeanCalculation == "TRUE"){
#           
#           all$newcolumn <- apply(all[,c(i,MeanIndex)], 1, max)
#           
#           colnames(all)[ncol(all)] <- paste0("Vs_",l,"_Mean")
#           
#         } else {
#           
#           all$newcolumn <- apply(all[,c(i,ControlIndex)], 1, max)
#           
#           colnames(all)[ncol(all)] <- paste0("Vs_",l,"_",ControlTag)
#         }
#       }
#     }}
# }
# write.table(all, file = paste(WD,"/",j,"/Pre-SanXoT/ID-q.txt",sep=""), sep="\t", row.names = FALSE)