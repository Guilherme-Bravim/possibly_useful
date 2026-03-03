################################################################################
##                      FUNCAO PARA DOSAGEM ALELICA
################################################################################


allelic.dosage<-function(data,missing="NA",callrate=0.90,maf=0.05){
library(progress)

    pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = 3, clear = FALSE, width=60
  )
  
  results<-list()
  start.time<-Sys.time() 
  #-----------------------  
  # Alleic Dosage  
  #-----------------------
  pb <- progress_bar$new(
    format = "Progress of the allelic dosage process: [:bar] :percent eta: :eta",
    total = ncol(data), clear = FALSE, width= 100)
  
  
  transform<-function(data){
    
    
    data.1<- data[data != missing]
    
    if(length(data.1)==0){
      data.transform<- data # aqui vai remover todos marcadores que nao possuem dados
      pb$tick()
      return(data.transform)
    }else {
      
      data.transform<- data
      resume <- data.frame(table(data.1))
      
      alle.type <- function(x) {
        alleles <- unique(unlist(strsplit(as.character(x), "")))
        paste0(alleles, collapse = "")
      }
      
      resume$alle.type <- sapply(resume$data.1, alle.type)
      freq.alle <- data.frame(table(unlist(strsplit(data.1 , "")))/sum(table(unlist(strsplit(data.1 , "")))))
      colnames(freq.alle)<-c("alle.type","Freq.alle")
      resume<-merge(resume,freq.alle,by="alle.type",all=T)
      
      resume$dosage <- NA
      resume[which(is.na(resume$Freq.alle)), "dosage"] <- 1
      if (max(resume$Freq.alle, na.rm = T) != min(resume$Freq.alle, na.rm =T)) {
        resume[which.max(resume$Freq.alle), "dosage"] <- 0 # possibilidade de inserir na funcao para o usuario indicar a forma
        resume[which.min(resume$Freq.alle), "dosage"] <- 2
      } else{
        resume[1, "dosage"] <- 0 # possibilidade de inserir na funcao para o usuario indicar a forma
        resume[3, "dosage"] <- 2
      }
      
      dom<-as.character(resume[which(resume$dosage==0),2])
      het<-as.character(resume[which(resume$dosage==1),2])
      rec<-as.character(resume[which(resume$dosage==2),2])
      
      data.transform[data.transform == dom] <- 0
      data.transform[data.transform == het] <- 1
      data.transform[data.transform == rec] <- 2
      
      pb$tick() 
      return(data.transform)
    }
  }
  
  
  encode <- as.data.frame(apply(data,2, transform))
  results[["E"]]<-encode
  
  #-----------------------  
  # Quality control
  #----------------------- 
  pb2 <- progress_bar$new(
    format = "Progress of the quality control process: [:bar] :percent eta: :eta",
    total = 4, clear = FALSE, width= 100)
  
  data.control<-results[["E"]]  
  pb2$tick() 
  # CALL RATE   
  CR<-lapply(data.control,function(x){
    cr<-1-(length(which(x==missing))/sum(table(x)))
    return(cr)
  })
  
  inf<-data.frame(which(unlist(CR) < callrate))
  
  if(exists("inf") & nrow(inf)>=1){
    colnames(inf)<-"markers"
    M.cr<-data.control[,-c(inf[,"markers"])]  
  }else{M.cr<- data.control}
  
  results[["M.CallRate"]]<-M.cr
  results[["Markers.removed"]]<-ncol(data.control)-ncol(M.cr)
  pb2$tick()  
  # MAF 
  MAF<-lapply(M.cr,function(x){
    data.filter<-x[x != missing]
    maf<-((2*length(which(data.filter==0)))+length(which(data.filter==1)))/(2*sum(table(data.filter)))
    return(maf)
  })
  
  inf.maf<-data.frame(which(unlist(MAF) > 1-maf |  unlist(MAF) < maf))
  
  if(exists("inf.maf") & nrow(inf.maf)>=1){
    colnames(inf.maf)<-"markers"
    M.maf<-M.cr[,-c(inf.maf[,"markers"])]  
  }else{M.maf<- M.cr}
  
  results[["M.MAF"]]<-M.maf
  results[["Markers.removed.maf"]]<-ncol(M.cr)-ncol(M.maf)
  pb2$tick() 
  
  # data imputation
  # Criterio de imputacao: 2p <= 0.5 : 0 / 0.5 < 2p <= 1.5 : 1 / 1.5 < 2p <= 2 = 2  
  
  data.i<-lapply(as.data.frame(M.maf),function(x){
    
    data.imp<-x
      if(is.na(missing)){
    data.filter<-x[!is.na(missing)]   
      }else{
    data.filter<-x[x != missing] }
        
    p<-((2*length(which(data.filter==0)))+length(which(data.filter==1)))/(2*sum(table(data.filter)))
    p2<-2*p
    
    if (p2 <= 0.5) {
      w <- 0
    } else if (p2 > 0.5 & p2 <= 1.5) {
      w <- 1
    } else {
      w <- 2
    }
    
    if(missing=="NA"){
      data.imp[is.na(data.imp)]<-w
      data.imp<-as.numeric(data.imp)
      return(data.imp)  

    }else{
      data.imp[which(data.imp==missing)]<-w
      data.imp<-as.numeric(data.imp)
      return(data.imp)    
    }
  
  })
  
  data.imput<-do.call(cbind,data.i)
  rownames(data.imput)<- rownames(results[["M.MAF"]])
  results[["M"]]<-data.imput
  pb2$tick() 
  ##
  
  end.time<-Sys.time()
  diff.time <- difftime(end.time,start.time)
  
  cat("################################################################################\n\n")
  cat("                     Analysis completed                 \n\n")
  cat(paste("Analysis performed on", Sys.Date(), "at",format(start.time, "%H:%M:%S"),"\n\n"))
  cat(paste("Number of initial markers:", ncol(data.control),"\n\n"))
  cat(paste("Call Rate:", callrate, "-- Number of markers removed:",results[["Markers.removed"]],"\n\n"))
  cat(paste("MAF:", maf, "-- Number of markers removed:",results[["Markers.removed.maf"]],"\n\n"))
  cat(paste("Total number of markers removed by Call Rate and MAF:",results[["Markers.removed"]]+results[["Markers.removed.maf"]],"\n\n"))
  cat(paste("Time for analysis processing:",round(diff.time[[1]],2),units(diff.time),"\n\n"))
  cat("################################################################################")
  return(results)
}






























