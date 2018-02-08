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


# create the Q-all
files <- list.files(path = opt$indir, pattern="*.msf.csv", full.names = TRUE)
all_q <- do.call("rbind", lapply(files, read.csv, header = TRUE))
if (opt$daemon == TRUE | opt$search == 2) {
  all_q$FileName<-paste(all_q$FileName,".raw",sep="")
} else {
  all_q$FileName<-substring(all_q$FileName,1,(nchar(as.character(all_q$FileName))-4))
  all_q$FileName<-paste(all_q$FileName,".raw",sep="")
}
write.table(all_q, file = paste(opt$outdir,"/Q-all.txt",sep=""), sep="\t", row.names = FALSE)


# calculation
y<-all_q
q_all<-data.frame()

for (i in opt$channels) {
  TMT<-y[,"QuanChannelID",drop=FALSE]==i
  z<-y[TMT,][,,drop=FALSE]
  TMTgood<-complete.cases(z)	#posicion de NaN
  a<-z[TMTgood,][,,drop=FALSE]
  c<-a[,c("FirstScan","Height1","FileName")]
  colnames(c)=c("FirstScan",i,"FileName")
  if ( nrow(q_all)==0 ) {
    q_all <- c
  }
  else {
    q_all<-merge(q_all,c)
  }
}

source("global.R")

list[] <- extract_tags(opt$label, opt$tags)

# # extract tags  
# if ( opt$label == "TMT" ) {
#   if ( opt$tags == "ALL" ) {
#     colnames(q_all)=c("FirstScan","FileName","X126","X127_N","X127_C","X128_N","X128_C","X129_N","X129_C","X130_N","X130_C","X131")
#   }
#   else {
#     colnames_TMT=c("FirstScan","FileName")
#     TagsUsed=paste0("X",opt$tags)
#     colnames_TMT=append(colnames_TMT, TagsUsed)
#     colnames(q_all)=colnames_TMT
#     colnames_TMT=c("Raw_FirstScan")
#     colnames_TMT=append(colnames_TMT, TagsUsed)
#   }
# } else {
#   if (TagsUsed == "ALL"){
#       colnames(q_all)=c("FirstScan","FileName","X113","X114","X115","X116","X117","X118","X119","X121")
#   }
#   else {
#     colnames_iTRAQ=c("FirstScan","FileName")
#     TagsUsed=paste0("X",opt$tags)
#     colnames_iTRAQ=append(colnames_iTRAQ, TagsUsed)
#     colnames(q_all)=colnames_iTRAQ
#     colnames_iTRAQ=c("Raw_FirstScan")
#     colnames_iTRAQ=append(colnames_iTRAQ, TagsUsed)
#   }
# }
write.table(q_all, file = paste(opt$outdir,"/Q-all.xls",sep=""), sep=",", row.names = FALSE)

