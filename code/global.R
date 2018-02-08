# extract tags  

extract_tags <- function(label, tags) {
  if ( label == "TMT" ) {
    if ( tags == "ALL" ) {
      colnames(q_all)=c("FirstScan","FileName","X126","X127_N","X127_C","X128_N","X128_C","X129_N","X129_C","X130_N","X130_C","X131")
    }
    else {
      colnames_TMT=c("FirstScan","FileName")
      TagsUsed=paste0("X",tags)
      colnames_TMT=append(colnames_TMT, TagsUsed)
      colnames(q_all)=colnames_TMT
      colnames_TMT=c("Raw_FirstScan")
      colnames_TMT=append(colnames_TMT, TagsUsed)
    }
  } else {
    if (TagsUsed == "ALL"){
      colnames(q_all)=c("FirstScan","FileName","X113","X114","X115","X116","X117","X118","X119","X121")
    }
    else {
      colnames_iTRAQ=c("FirstScan","FileName")
      TagsUsed=paste0("X",tags)
      colnames_iTRAQ=append(colnames_iTRAQ, TagsUsed)
      colnames(q_all)=colnames_iTRAQ
      colnames_iTRAQ=c("Raw_FirstScan")
      colnames_iTRAQ=append(colnames_iTRAQ, TagsUsed)
    }
  }
}

