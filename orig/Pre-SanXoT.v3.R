(WD <- getwd())
if (!is.null(WD)) setwd(WD)
WD <- "D:/projects/Pre-SanXoT/orig/"
source(paste0(WD,"/Config.txt"))

################################################################################################################

# Pre-SanXoT

################################################################################################################

list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
    # use full.names=TRUE to pass to file.info
    all <- list.files(path, pattern, all.dirs,
                      full.names=TRUE, recursive=FALSE, ignore.case)
    dirs <- all[file.info(all)$isdir]
    # determine whether to return full names or just dir names
    if(isTRUE(full.names))
        return(dirs)
    else
        return(basename(dirs))
}

MSFfolders <- list.dirs(path = paste0(WD,"/",Expto,"/MSF"), pattern=Patern)

library("RSQLite")

for (j in Expto){
    
  for (k in MSFfolders){

    files <- list.files(path = paste(WD,"/",j,"/MSF/",k,sep=""),pattern="*.msf")

    for (i in files) {
        
        db=dbConnect(SQLite(), dbname=paste(WD,"/",j,"/MSF/",k,"/",i,sep=""))
        
        if(SearchEngine=="2"){
        
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
        # BEGIN: jmrc
        dir.create(paste(WD,"/",j,"/Pre-SanXoT/",sep=""), showWarnings = FALSE)
        # END: jmrc
        write.csv(data, file=paste(WD,"/",j,"/Pre-SanXoT/",i,".csv",sep=""),row.names=FALSE)
    }
  }
}

##############################
# CREATE ID-Q -------------- #
##############################

for (j in Expto){

    files <- list.files(path = paste(WD,"/",j,"/Pre-SanXoT",sep=""),pattern="*.csv", full.names = TRUE)
    
    all_q <- do.call("rbind", lapply(files, read.csv, header = TRUE))
    
    if (Daemon == "TRUE" | SearchEngine=="2"){
        
        all_q$FileName<-paste(all_q$FileName,".raw",sep="")
        
    } else {
    
    all_q$FileName<-substring(all_q$FileName,1,(nchar(as.character(all_q$FileName))-4))
    
    all_q$FileName<-paste(all_q$FileName,".raw",sep="")}
    
    write.table(all_q, file = paste(WD,"/",j,"/Pre-SanXoT/Q-all.txt",sep=""), sep="\t", row.names = FALSE)}

if (length(Expto)<2) {
  
  y<-all_q
  q_all<-data.frame()
  for (i in ChannelID){
    TMT<-y[,"QuanChannelID",drop=FALSE]==i
    z<-y[TMT,][,,drop=FALSE]
    TMTgood<-complete.cases(z)	#posicion de NaN
    a<-z[TMTgood,][,,drop=FALSE]
    c<-a[,c("FirstScan","Height1","FileName")]
    colnames(c)=c("FirstScan",i,"FileName")
    
    if (nrow(q_all)==0){
      q_all<-c
    } else {q_all<-merge(q_all,c)}
  }
  
  if (Typeoflabel=="TMT"){
    if (TagsUsed=="ALL"){
      colnames(q_all)=c("FirstScan","FileName","X126","X127_N","X127_C","X128_N","X128_C","X129_N","X129_C","X130_N","X130_C","X131")
    } else {
      colnames_TMT=c("FirstScan","FileName")
      TagsUsed=paste0("X",TagsUsed)
      colnames_TMT=append(colnames_TMT, TagsUsed)
      colnames(q_all)=colnames_TMT
      colnames_TMT=c("Raw_FirstScan")
      colnames_TMT=append(colnames_TMT, TagsUsed)}
  }else{
    if (TagsUsed=="ALL"){
      colnames(q_all)=c("FirstScan","FileName","X113","X114","X115","X116","X117","X118","X119","X121")
    } else {
      colnames_iTRAQ=c("FirstScan","FileName")
      TagsUsed=paste0("X",TagsUsed)
      colnames_iTRAQ=append(colnames_iTRAQ, TagsUsed)
      colnames(q_all)=colnames_iTRAQ
      colnames_iTRAQ=c("Raw_FirstScan")
      colnames_iTRAQ=append(colnames_iTRAQ, TagsUsed)
    }
  }
  
  write.table(q_all, file = paste(WD,"/",j,"/Pre-SanXoT/Q-all.xls",sep=""), sep=",", row.names = FALSE)
  
} else {
  
  for (j in Expto){
    
    files <- list.files(path = paste(WD,"/",j,"/Pre-SanXoT",sep=""),pattern="Q-all.txt", full.names = TRUE)
    y<-read.table(files, header=TRUE, sep="\t")
    q_all<-data.frame()
    
    for (i in ChannelID){
        TMT<-y[,"QuanChannelID",drop=FALSE]==i
        z<-y[TMT,][,,drop=FALSE]
        TMTgood<-complete.cases(z)	#posicion de NaN
        a<-z[TMTgood,][,,drop=FALSE]
        c<-a[,c("FirstScan","Height1","FileName")]
        colnames(c)=c("FirstScan",i,"FileName")
        
        if (nrow(q_all)==0){
            
            q_all<-c
            
        } else {q_all<-merge(q_all,c)}
    }
    
    if (Typeoflabel=="TMT"){
      if (TagsUsed=="ALL"){
          colnames(q_all)=c("FirstScan","FileName","X126","X127_N","X127_C","X128_N","X128_C","X129_N","X129_C","X130_N","X130_C","X131")
      } else {
          colnames_TMT=c("FirstScan","FileName")
          TagsUsed=paste0("X",TagsUsed)
          colnames_TMT=append(colnames_TMT, TagsUsed)
          colnames(q_all)=colnames_TMT
          colnames_TMT=c("Raw_FirstScan")
          colnames_TMT=append(colnames_TMT, TagsUsed)
      }
    }else{
      if (TagsUsed=="ALL") {
        colnames(q_all)=c("FirstScan","FileName","X113","X114","X115","X116","X117","X118","X119","X121")
      } else {
        colnames_iTRAQ=c("FirstScan","FileName")
        TagsUsed=paste0("X",TagsUsed)
        colnames_iTRAQ=append(colnames_iTRAQ, TagsUsed)
        colnames(q_all)=colnames_iTRAQ
        colnames_iTRAQ=c("Raw_FirstScan")
        colnames_iTRAQ=append(colnames_iTRAQ, TagsUsed)
      }
    }
    write.table(q_all, file = paste(WD,"/",j,"/Pre-SanXoT/Q-all.xls",sep=""), sep=",", row.names = FALSE)
  }
}

################################
# CREATE ID-ALL -------------- #
################################
for (j in Expto){
    
  for (k in MSFfolders) {
    files <- list.files(path = paste(WD,"/",j,"/MSF/",k,sep=""),pattern="_results", full.names = TRUE)
    if (length(files) > 0) {
      ID_all<- read.table(files, sep="\t",comment.char = "¡",quote = "¿", header = TRUE)
      files <- list.files(path = paste(WD,"/",j,"/MSF/",k,sep=""),pattern="_results")
      write.table(ID_all, file = paste(WD,"/",j,"/Pre-SanXoT/",k,files,sep=""), sep="\t", row.names = FALSE)
    }
  }
  
  files <- list.files(path = paste(WD,"/",j,"/Pre-SanXoT",sep=""),pattern="_results", full.names = TRUE)
  ID_all <- do.call("rbind", lapply(files, read.table, header = TRUE))
  write.table(ID_all, file = paste(WD,"/",j,"/Pre-SanXoT/ID-all.txt",sep=""), sep="\t", row.names = FALSE)
  print( "_resutls" )
  print( files )
  # file.remove(files)
  files <- list.files(path = paste(WD,"/",j,"/Pre-SanXoT",sep=""),pattern="*.csv", full.names = TRUE)
  print( "csv" )
  print( files )
  # file.remove(files)
}

##############################
# CREATE ID-q -------------- #
##############################
if (length(Expto)<2) {
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
    if (Typeoflabel == "TMT"){
      MeanTags<-c("X126","X127_N","X127_C","X128_N","X128_C","X129_N","X129_C","X130_N","X130_C","X131")
    } else {
      MeanTags<-c("X113","X114","X115","X116","X117","X118","X119","X121")
    }
  }
  
  all<-merge(x,q)
  FirstTagIndex=as.numeric(grep(paste0("X",FirstTag), colnames(all)))
  CalcIndex=trunc(seq(FirstTagIndex, by=(length(ChannelID)/as.numeric(Comparatives)), len = as.numeric(Comparatives)),1)
  
  # BEGIN: jmrc TEMPORAL!!
  i = 24
  all$Mean <- rowMeans(all[,paste0("X",MeanTags)])
  MeanIndex=as.numeric(grep("Mean", colnames(all)))
  # END 
  for (i in CalcIndex){
    
    ControlIndex=as.numeric(grep(paste0("X",ControlTag), colnames(all)))
    
    if (MeanCalculation == "TRUE") {
      # all$Mean <- rowMeans(all[,paste0("X",MeanTags)])
      # MeanIndex=as.numeric(grep("Mean", colnames(all)))
      all$newcolumn <- log2(all[,i]/all$Mean)
      l <- substring(colnames(all)[i],2)
      colnames(all)[ncol(all)] <- paste0("Xs_",l,"_Mean")
    }
    else {
      all$newcolumn <- log2(all[,i]/all[,ControlIndex])
      l <- substring(colnames(all)[i],2)
      colnames(all)[ncol(all)] <- paste0("Xs_",l,"_",ControlTag)
    }
    
    if (Absolute == "TRUE"){
      all$newcolumn <- all[,c(i)]
      colnames(all)[ncol(all)] <- paste0("Vs_",l,"_ABS")
    }
    if (Absolute == "FALSE"){
      if (MeanCalculation == "TRUE"){
        all$newcolumn <- apply(all[,c(i,MeanIndex)], 1, max)
        colnames(all)[ncol(all)] <- paste0("Vs_",l,"_Mean")
      }
      else {
        all$newcolumn <- apply(all[,c(i,ControlIndex)], 1, max)
        colnames(all)[ncol(all)] <- paste0("Vs_",l,"_",ControlTag)
      }
    }
    if (Absolute == "BOTH"){
      all$newcolumn <- all[,c(i)]
      colnames(all)[ncol(all)] <- paste0("Vs_",l,"_ABS")
      if (MeanCalculation == "TRUE"){
        all$newcolumn <- apply(all[,c(i,MeanIndex)], 1, max)
        colnames(all)[ncol(all)] <- paste0("Vs_",l,"_Mean")
      }
      else {
        all$newcolumn <- apply(all[,c(i,ControlIndex)], 1, max)
        colnames(all)[ncol(all)] <- paste0("Vs_",l,"_",ControlTag)
      }
    }
  } # end for calcIndex
  
  print( paste(WD,"/",j,"/Pre-SanXoT/ID-q.txt",sep="") )
  write.table(all, file = paste(WD,"/",j,"/Pre-SanXoT/ID-q.txt",sep=""), sep="\t", row.names = FALSE)
  
}
# else {
  
  # for (j in Expto){
  #   
  #   files <- list.files(path = paste(WD,"/",j,"/Pre-SanXoT",sep=""),pattern="Q-all.xls", full.names = TRUE)
  #   
  #   k<-read.table(files, header=TRUE, sep=",")
  #   
  #   files <- list.files(path = paste(WD,"/",j,"/Pre-SanXoT",sep=""),pattern="ID-all.txt", full.names = TRUE)
  #   
  #   x<-read.table(files, header=TRUE, sep="\t")
  #   
  # x$Raw_FirstScan<-do.call(paste, c(x[c("RAWFile","FirstScan")], sep = ""))
  # k$Raw_FirstScan<-do.call(paste, c(k[c("FileName","FirstScan")], sep = ""))
  # x$Raw_FirstScan<-as.character(x$Raw_FirstScan)
  # k$Raw_FirstScan<-as.character(k$Raw_FirstScan)
  # 
  # if (Typeoflabel=="TMT"){
  #     
  #     if (TagsUsed=="ALL"){
  #         
  #         q<-k[,c("Raw_FirstScan","X126","X127_N","X127_C","X128_N","X128_C","X129_N","X129_C","X130_N","X130_C","X131")]
  #         
  #     } else {
  #         
  #         q<-k[,colnames_TMT]}
  #     
  # }else{
  #     
  #     if (TagsUsed=="ALL"){
  #         
  #         q<-k[,c("Raw_FirstScan","X113","X114","X115","X116","X117","X118","X119","X121")]
  #         
  #     } else {
  #         
  #         q<-k[,colnames_iTRAQ]}}
  #     
  # all<-merge(x,q)
  #     
  # FirstTagIndex=as.numeric(grep(paste0("X",FirstTag), colnames(all)))
  #     
  # CalcIndex=trunc(seq(FirstTagIndex, by=(length(ChannelID)/as.numeric(Comparatives)), len = as.numeric(Comparatives)),1)
  # 
  # if (MeanTags=="ALL"){
  #     
  #     if (Typeoflabel == "TMT"){
  #         
  #         MeanTags<-c("X126","X127_N","X127_C","X128_N","X128_C","X129_N","X129_C","X130_N","X130_C","X131")
  #         
  #     } else {
  #         
  #         MeanTags<-c("X113","X114","X115","X116","X117","X118","X119","X121")}}
  # 
  # for (i in CalcIndex){
  # 
  #     ControlIndex=as.numeric(grep(paste0("X",ControlTag), colnames(all)))
  #     
  #     if (MeanCalculation == "TRUE") {
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
  #     } else {
  #     
  #     all$newcolumn <- log2(all[,i]/all[,ControlIndex])
  #     
  #     l <- substring(colnames(all)[i],2)
  #     
  #     colnames(all)[ncol(all)] <- paste0("Xs_",l,"_",ControlTag)}
  #     
  #     if (Absolute == "TRUE"){
  #         
  #         all$newcolumn <- all[,c(i)]
  #         
  #         colnames(all)[ncol(all)] <- paste0("Vs_",l,"_ABS")}
  #         
  #     if (Absolute == "FALSE"){
  #         
  #         if (MeanCalculation == "TRUE"){
  #             
  #             all$newcolumn <- apply(all[,c(i,MeanIndex)], 1, max)
  #             
  #             colnames(all)[ncol(all)] <- paste0("Vs_",l,"_Mean")
  #             
  #         } else {
  #         
  #         all$newcolumn <- apply(all[,c(i,ControlIndex)], 1, max)
  #         
  #         colnames(all)[ncol(all)] <- paste0("Vs_",l,"_",ControlTag)}}
  #     
  #     if (Absolute == "BOTH"){
  #         
  #         all$newcolumn <- all[,c(i)]
  #         
  #         colnames(all)[ncol(all)] <- paste0("Vs_",l,"_ABS")
  #         
  #         if (MeanCalculation == "TRUE"){
  #             
  #             all$newcolumn <- apply(all[,c(i,MeanIndex)], 1, max)
  #             
  #             colnames(all)[ncol(all)] <- paste0("Vs_",l,"_Mean")
  #             
  #         } else {
  #             
  #         all$newcolumn <- apply(all[,c(i,ControlIndex)], 1, max)
  #             
  #         colnames(all)[ncol(all)] <- paste0("Vs_",l,"_",ControlTag)}}}
  #   }
  # }
  
  # BEGIN: jmrc
  # if (Random == "YES"){
  #   
  #   for (i in CalcIndex){
  #     
  #     for (m in CalcIndex){
  #       
  #       all$newcolumn <- log2(all[,i]/all[,m])
  #       
  #       l <- substring(colnames(all)[i],2)
  #       
  #       o <- substring(colnames(all)[m],2)
  #       
  #       colnames(all)[ncol(all)] <- paste0("Xs_",l,"_",o)
  #       
  #       all$newcolumn <- apply(all[,c(i,m)], 1, max)
  #       
  #       colnames(all)[ncol(all)] <- paste0("Vs_",l,"_",o)}}}
  # END: jmrc  
    
# write.table(all, file = paste(WD,"/",j,"/Pre-SanXoT/ID-q.txt",sep=""), sep="\t", row.names = FALSE)

################################################################################################################

# Tag file Maker
# BEGIN: jmrc
# NOTE!! I don't know if it is used!!!
# END: jmrc

################################################################################################################
# 
# if (TagsUsed=="ALL"){
#   
#   if (Typeoflabel == "TMT"){
#     
#     TagsUsed<-c("126","127_N","127_C","128_N","128_C","129_N","129_C","130_N","130_C","131")
#     
#   } else {
#     
#     TagsUsed<-c("113","114","115","116","117","118","119","121")}}
# 
# Tag<-c()
# 
# for (i in TagsUsed){
#   for (j in Expto){
#     
#     tags_temp<-paste(j,i,sep="_")
#     
#     if (NROW(Tag)==0){
#       
#       Tag<-tags_temp
#       
#     } else {Tag<-rbind(Tag,tags_temp)}}}
# 
# Path<-c()
# 
# for (i in TagsUsed){
#   for (j in Expto){
#     
#     path_temp<-paste(WD,"/",j,"/SanXoT/",i,"/data/Q2A_lowerNormV.xls",sep="")
#     
#     if (NROW(Tag)==0){
#       
#       Path<-path_temp
#       
#     } else {Path<-rbind(Path,path_temp)}}}
# 
# Tag<-as.data.frame(Tag)
# row.names(Tag) <- NULL
# colnames(Tag) <- "Tag"
# 
# Path<-as.data.frame(Path)
# row.names(Path) <- NULL
# colnames(Path) <- "Path"
# 
# Tag_file_temp<-cbind(Tag,Path)
# 
# # BEGIN: jmrc
# IntegrationDir = paste(WD,"/",Expto,"/Integration",sep="")
# dir.create(IntegrationDir, showWarnings = FALSE)
# # END: jmrc
# # BEGIN: jmrc
# #write.table(Tag_file_temp, file=paste0(WD,"/Integration/Tag_file_temp.txt"),sep="\t", row.names = FALSE)
# write.table(Tag_file_temp, file=paste0(IntegrationDir,"/Tag_file_temp.txt"),sep="\t", row.names = FALSE)
# # END: jmrc
# 
# if (NOI == 1){
#   
#   Tag_file<-Tag_file_temp[Tag_file_temp$Tag %in% Integration,]
#   
#   # BEGIN: jmrc
#   #write.table(Tag_file, file=paste0(WD,"/Integration/Tag_file.txt"),sep="\t", row.names = FALSE)
#   write.table(Tag_file, file=paste0(IntegrationDir,"/Tag_file.txt"),sep="\t", row.names = FALSE)
#   # END: jmrc
#   
# } else {
#   
#   for (j in Integrations) {
#     
#     tag<-paste(j,get(j),sep="_")
#     
#     for (i in tag){
#       
#       tag<-substring(tag,(nchar(j)+2),nchar(tag))
#       
#       tag<-paste(Expto,tag,sep="_")
#       
#       Tag_file<-Tag_file_temp[Tag_file_temp$Tag %in% tag, ]
#       
#       if (nrow(Tag_file)>1){
# 
#         # BEGIN: jmrc
#         #write.table(Tag_file, file=paste0(WD,"/Integration/",j,"_Tag_file.txt"),sep="\t", row.names = FALSE)
#         write.table(Tag_file, file=paste0(IntegrationDir,"/",j,"_Tag_file.txt"),sep="\t", row.names = FALSE)
#         # END: jmrc
#         
#       }
#     }
#   }
# }
# 
