#'managing CRAN logs. 
#' @description code adapted from www.nicebread.de finally-tracking-cran-packages-downloads (original code is under the FreeBSD license, by Felix Schonbrodt)
#' @param start starting date
#' @param today last date
#' @param logdirectory the path where files will be placed
#' @param global boolean to define wether a global file will be created or not
#' @param packages vector of the names of the packages to keep in the global file(NULL will keep every packages)
#' @export
#logdirectory="L:/packages/CRANlogs"
log_manager<-function(start=as.Date('2012-10-01'),today=format(Sys.Date()-1, "%Y-%m-%d"),logdirectory="CRANlogs", global=FALSE,packages=NULL){   
   today=as.Date(today)
   start=as.Date(start)
   print(start)
   print(today)
   all_days <- seq(start, today, by = 'day')
   
   year <- as.POSIXlt(all_days)$year + 1900
   urls <- paste0('http://cran-logs.rstudio.com/', year, '/', all_days, '.csv.gz')
   
   # only download the files you don't have:
   missing_days <- setdiff(as.character(all_days), tools::file_path_sans_ext(dir(logdirectory), TRUE))
   if(length(missing_days)>0){
      dir.create(logdirectory)
      for (i in 1:length(missing_days)) {
         print(paste0(i, "/", length(missing_days)))
         download.file(urls[i], paste0(logdirectory,'/', missing_days[i], '.csv.gz'))
         
      }
      if(global){#we reduce the file
         
         file_list <- list.files(logdirectory, full.names=TRUE)
         
         logs <- list()
         for (file in file_list) {
            print(paste("Reading", file, "..."))
            logs[[file]] <- read.table(file, header = TRUE, sep = ",", quote = "\"",
                                       dec = ".", fill = TRUE, comment.char = "", as.is=TRUE)
            logs[[file]]=logs[[file]][logs[[file]]$package%in%packages,]
         }
         
         # rbind together all files
#          require(data.table)
         dat <- rbindlist(logs)
         
         # add some keys and define variable types
         dat[, date:=as.Date(date)]
         dat[, package:=factor(package)]
         dat[, country:=factor(country)]
         dat[, weekday:=weekdays(date)]
         dat[, week:=strftime(as.POSIXlt(date),format="%Y-%W")]
         
         setkey(dat, package, date, week, country)
         
         save(dat, file=paste(logdirectory,"/CRANlogs",packages,".RData",sep=""))
      }
   }else{
      print("nothing new")
   }
}