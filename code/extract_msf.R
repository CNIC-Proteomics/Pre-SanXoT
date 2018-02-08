library("RSQLite")
library("optparse")

# get input parameters
option_list = list(
  make_option(c("-i", "--indir"), type="character", default=NULL,
              help="dataset directory", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="Fr_*", 
              help="Pattern of folders that contains the MSFs [default= %default]", metavar="character"),
  make_option(c("-s", "--search"), type="integer", default=2,
              help="Search engine [default= %default]", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="output directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# BEGIN: IMPORTANT!!!! HARD-CORE inputs!!!
opt$indir <- "D:/projects/Pre-SanXoT/pratio_FractionsSeparated/MSF"
opt$pattern <- "Fr_*"
opt$search <- 2
opt$outdir <- "D:/projects/Pre-SanXoT/pratio_FractionsSeparated/Pre-SanXoT_2"
# END: IMPORTANT!!!! HARD-CORE inputs!!!

if ( is.null(opt$indir) || is.null(opt$outdir) ){
  print_help(opt_parser)
  stop("All arguments must be supplied.n", call.=FALSE)
}

print( "input parameters: ")
print( opt )

# create workspace if not exit
dir.create(opt$outdir, showWarnings = FALSE)

# recreate the list.dirs function!!
list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE, full.names=FALSE, ignore.case=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs, full.names=TRUE, recursive=FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}

# extract the MSF files and the result files
MSFfolders <- list.dirs(path = opt$indir, pattern=opt$pattern)
for (k in MSFfolders) {
  idir = paste(opt$indir,"/",k,sep="")
  files <- list.files(path = idir,pattern="*.msf")
  
  for (i in files) {
    dbfile = paste(idir,"/",i,sep="")
    print( paste0("extract dbfile: ", dbfile) )
    
    db=dbConnect(SQLite(), dbname= dbfile)
    if(opt$search == 2) {
      data=dbGetQuery(conn = db,
                      "SELECT [SpectrumHeaders].[FirstScan],
                      [ReporterIonQuanResults].[Mass] AS [Mass2],
                      [ReporterIonQuanResults].[Height] AS [Height1],
                      [SpectrumHeaders].[RetentionTime],
                      [ReporterIonQuanResults].[QuanChannelID],
                      [MassPeaks].[MassPeakID],
                      [Workflows].[WorkflowName] AS [FileName]
                      FROM [ReporterIonQuanResults]
                      INNER JOIN [SpectrumHeaders] ON [ReporterIonQuanResults].[SpectrumID] =
                      [SpectrumHeaders].[SpectrumID]
                      INNER JOIN [MassPeaks] ON [MassPeaks].[MassPeakID] =
                      [SpectrumHeaders].[MassPeakID]
                      INNER JOIN [WorkflowInputFiles] ON [MassPeaks].[FileID] =
                      [WorkflowInputFiles].[FileID]
                      INNER JOIN [Workflows] ON [WorkflowInputFiles].[WorkflowID] =
                      [Workflows].[WorkflowID]
                      WHERE [ReporterIonQuanResults].[Mass] > 0")
      } else {
        data=dbGetQuery(conn = db,
                        "SELECT [SpectrumHeaders].[FirstScan],
                        [ReporterIonQuanResults].[Mass] AS [Mass2],
                        [ReporterIonQuanResults].[Height] AS [Height1],
                        [SpectrumHeaders].[RetentionTime],
                        [ReporterIonQuanResults].[QuanChannelID],
                        [MassPeaks].[MassPeakID],
                        [WorkflowInfo].[WorkflowName] AS [FileName]
                        FROM [ReporterIonQuanResults]
                        INNER JOIN [SpectrumHeaders] ON [ReporterIonQuanResults].[SpectrumID] =
                        [SpectrumHeaders].[SpectrumID]
                        INNER JOIN [MassPeaks] ON [MassPeaks].[MassPeakID] =
                        [SpectrumHeaders].[MassPeakID]
                        INNER JOIN [FileInfos] ON [MassPeaks].[FileID] = [FileInfos].[FileID],
                        [WorkflowInfo]
                        WHERE [ReporterIonQuanResults].[Mass] > 0")
      }
      i <- substr(i, 1, nchar(i) - 4)
      write.csv(data, file=paste(opt$outdir,"/",i,".msf.csv",sep=""),row.names=FALSE)
  }

  # copy MSF results
  files <- list.files(path = idir, pattern="_results", full.names = TRUE)
  print( paste0("extract result file: ", files) )  
  if (length(files) > 0) {
    ID_all<- read.table(files, sep="\t",comment.char = "!",quote = "¿", header = TRUE)
    files <- list.files(path = idir, pattern="_results")
    write.table(ID_all, file = paste(opt$outdir,"/", k, files,sep=""), sep="\t", row.names = FALSE)
  }

} # end MSFfolder

dbDisconnect(db)

