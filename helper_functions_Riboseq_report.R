library(glue)
library(types)
library(stringr)
library(dplyr)
library(UpSetR)


calculate_background_threshold <- function(unique_region_type = ? character){
  #get directory list of subdirectories of the current working directory
  directories<-list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
  
  #get random intersect counts of the region type of interest
  randomfiles=list()
  for(i in directories){
    if(!(identical(list.files(i,pattern=paste0("*", unique_region_type, "_random_intersect_counts_relative_sorted.bed")), character(0)))){
      randomfiles=c(randomfiles,paste0(i,"/",list.files(i,pattern=paste0("*", unique_region_type, "_random_intersect_counts_relative_sorted.bed"))))
    }
  }
  print(randomfiles)
  
  #read in intersections with random alignments and calculate a threshold for each
  randomdataframes <- lapply(randomfiles, read.csv, header = FALSE, sep="\t")
  background=list()
  for(randomframe in randomdataframes){
    colnames(randomframe)=c("ID","start","stop","read_count", "relative_count")
    CI_upper_lim=t.test(randomframe$relative_count)[[4]][[2]]
    #for the whole relative counts of the intersections with the random regions
    #take upper tail of the 95% CI: Threshold to be significantly different 
    #from the null distribution
    background = c(background,CI_upper_lim)
  }
  return(background)
}


print_thresholds <- function(unique_region_type = ? character, background){
  #get directory list of subdirectories of the current working directory
  directories<-list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
  #obtain names of the datasets for which the thresholds are calculated
  randomSetnames=list()
  for(i in directories){
    if(!(identical(list.files(i,pattern=paste0("*", unique_region_type, "_random_intersect_counts_relative_sorted.bed")), character(0)))){
      randomSetnames=c(randomSetnames,list.files(i,pattern=paste0("*", unique_region_type, "_random_intersect_counts_relative_sorted.bed")))
    }
  }
  randomnames=stringr::str_replace(randomSetnames, pattern = "_random_intersect_counts_relative_sorted.bed", replacement = "")
  
  printbackground=data.frame(x=background)
  printbackground=as.data.frame(t(printbackground))
  rownames(printbackground)=randomnames
  colnames(printbackground)="Threshold"
  print(kable(printbackground, caption="Thresholds for each dataset"))
}


get_intersect_files <- function(unique_region_type = ? character){
  directories=list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
  files=list()
  setnames=list()
  for(i in directories){
    if(!(identical(list.files(i,pattern=paste0("*", unique_region_type, "_intersect_counts_relative_sorted.bed")), character(0)))){
      files=c(files,paste0(i,"/",list.files(i,pattern=paste0("*", unique_region_type, "_intersect_counts_relative_sorted.bed"))))
      setnames=c(setnames,list.files(i,pattern=paste0("*", unique_region_type, "_intersect_counts_relative_sorted.bed")))
    }
  }
  dataframes <- lapply(unlist(files), read.csv, header = FALSE, sep="\t")
  
  datalist=list()
  for(a in dataframes){
    datalist=c(datalist,list(a))
  }
  dataframes=datalist
  names(dataframes) <- stringr::str_replace(setnames, pattern = "_intersect_counts_relative_sorted.bed", replacement = "")
  
  return(dataframes)
}


count_unique_regions_above_threhsold <- function(dataframes, threshold){
  
  relevantregionscount=list()
  relevantregions=list()
  i=1
  for(frame in dataframes){
    colnames(frame)=c("ID","start","stop","ORF","read_count",  "relative_count")
    temp=0
    c=1
    for(count in frame$relative_count){
      if(count>=threshold[[i]]){
        if(frame$read_count[[c]] > 1){
          temp = temp + 1
        }
      }
      c = c + 1
    }
    i = i + 1
    relevantregionscount = c(relevantregionscount,temp)
  }
  return(relevantregionscount)
}


print_unique_region_above_t <- function(relevantregionscount, dataframes){
  printframe=data.frame(x=relevantregionscount)
  printframe=as.data.frame(t(printframe))
  rownames(printframe)=names(dataframes)
  colnames(printframe)=paste0("Number of unique regions with relative count >= threshold")
  print(kable(printframe, caption="Regions above the threshold"))
  print(printframe)
}


get_top_5_unique_regions <- function(dataframes_list, threshold){
  upsetlist=list()
  j=1
  for(f in dataframes_list){
    bed=f
    bed = filter(bed, V6 > threshold[[j]])
    list<- rep(NA, length(bed[,1]))
    for(i in 1:length(bed[,1])){
      list[i]=paste0(bed[i,1],":",bed[i,2],":",bed[i,3])
    }
    upsetlist=c(upsetlist,list(list))
    colnames(bed)=c("ID","start","stop","ORF","read_count",  "relative_count")
    print(kable(bed[1:5,], caption=names(dataframes_list)[j], row.names = FALSE))
    j = j + 1
  }
  return(upsetlist)
}
