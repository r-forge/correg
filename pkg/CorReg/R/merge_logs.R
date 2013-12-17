#' merge Rstudio Cran log files in a dataset
#' @param logdirectory the path where files are placed
#' @param packages vector of the names of the packages to keep in the global file(NULL will keep every packages)
#' @param destdir directory for export
#' @export
merge_logs<-function(logdirectory="CRANlogs",packages="CorReg",destdir=NULL){
   #later version would allow to choose a period
   #start=as.Date('2012-10-01'),today=format(Sys.Date(), "%Y-%m-%d")
   file_list <- list.files(logdirectory, full.names=TRUE,pattern=".csv.gz")
   logs <- list()
   for (file in file_list){
      print(paste("Reading", file, "..."))
      logs[[file]] <- read.table(file, header = TRUE, sep = ",", quote = "\"",
                                 dec = ".", fill = TRUE, comment.char = "", as.is=TRUE)
      logs[[file]]=logs[[file]][logs[[file]]$package%in%packages,]
   }
   
   # rbind together all files
   dat <- rbindlist(logs)
   
   # add some keys and define variable types
   dat[, date:=as.Date(date)]
   dat[, package:=factor(package)]
   dat[, country:=factor(country)]
   dat[, weekday:=weekdays(date)]
   dat[, week:=strftime(as.POSIXlt(date),format="%Y-%W")]
   
   setkey(dat, package, date, week, country)
   if(is.null(destdir)){
      save(data.frame(dat), file=paste(logdirectory,"/CRANlogs",paste(packages,collapse=""),".RData",sep=""))
      print(paste("load(",paste(logdirectory,"/CRANlogs",packages,".RData",sep=""),")"))
      
   }else{
      filename=paste(destdir,"/CRANlogs",paste(packages,collapse=""),".csv",sep="")
      write.csv(dat,file=filename)
   }
   return(dat)
}

