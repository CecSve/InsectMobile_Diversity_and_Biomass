rm(list=ls())

Global_ID <- "DE"

Insect_env <- function(custom_path_basis = "script"){
  
  if(custom_path_basis == "script"){
    
    path_basis <- dirname(rstudioapi::getSourceEditorContext()$path) #gives u the position of this script (in windows)
    
  }else{
    path_basis <- custom_path_basis
  }
  
  path_basis <<- path_basis
  


  if(.Platform$OS.type == "unix") {
      dir.create(file.path(path_basis, "R_Library_unix"), showWarnings = FALSE)
      library_path <<- paste0(path_basis,"/R_Library_unix")
      
  }else{
      dir.create(file.path(path_basis, "R_Library"), showWarnings = FALSE)
      library_path <<- paste0(path_basis,"/R_Library")
      
  }
  .libPaths( c( library_path, .libPaths()) )
  
  
  WD_BASIS <<- paste0(path_basis,"/Proj_",Global_ID)
  
  dir.create(file.path(WD_BASIS), showWarnings = FALSE)
  

  dir.create(file.path(WD_BASIS, ".temp"), showWarnings = FALSE)
  dir.create(file.path(WD_BASIS, "InsektMobile_data"), showWarnings = FALSE)
  
  dir.create(file.path(WD_BASIS, "WD"), showWarnings = FALSE)
  dir.create(file.path(WD_BASIS, "Plot"), showWarnings = FALSE)
  dir.create(file.path(paste0(WD_BASIS, "/Plot"),"Distribution"), showWarnings = FALSE)
  
  setwd(paste0(WD_BASIS,"/WD"))
  
  return(WD_BASIS)
  
}

Insect_env()

VH_excel <- function(filename, path){
  
  library(readxl)
  
  t <- paste0(path,"/",filename)
  
  names_dta <- excel_sheets(t)
  
  dta <- lapply(excel_sheets(t), read_excel, path =t)
  
  for(i in 1:length(dta)){
    dta[[i]] <- as.data.frame(dta[[i]])
  }
  
  names(dta) <- names_dta
  
  c <- i <- 1
  
  for(i in 1:length(dta)){
    t <- as.data.frame(matrix(NA, nrow = ncol(dta[[i]]),ncol = 5))
    t$V1 <- colnames(dta[[i]]); colnames(t) <- c("colname","mean","min","max", "levels")
    
    for(c in 1:ncol(dta[[i]])){
      t[t$colname == colnames(dta[[i]][c]),]$levels <- nrow(unique(dta[[i]][c]))
      
      if(is.numeric(dta[[i]][c][1])){
        t[t$colname == colnames(dta[[i]][c]),]$mean <- sum(dta[[i]][c])/nrow(dta[[i]][c])
        t[t$colname == colnames(dta[[i]][c]),]$min <- min(dta[[i]][c], na.rm = TRUE)
        t[t$colname == colnames(dta[[i]][c]),]$max <- max(dta[[i]][c], na.rm = TRUE)
      }
    }
  }
  return(dta)
}

# 2018 Import#####

library(data.table)

dta <- VH_excel(filename = "Dry weight_Insektenmobil_Julianas lab book_16_05_2019-she (1).xlsx",path = paste0(WD_BASIS,"/InsektMobile_data"))[[1]]

dta <- dta[,colnames(dta) %in% c("TOP","Sample ID","Size","Dry mass (mg)")]

dta <- dta[!is.na(dta$`Sample ID`),]

dta_bit <- VH_excel(filename = "DE_2018samples_03032020.xlsx",path = paste0(WD_BASIS,"/InsektMobile_data"))[[2]]

dta <- merge(dta, dta_bit, by.x = "TOP", by.y = "PCRID", all.x = TRUE)


dta <- dta[,colnames(dta) %in% c("TOP","Sample ID","Size","Dry mass (mg)") | data.table::like(colnames(dta),"Qubit")]

t <- data.frame(do.call('rbind', strsplit(dta$`Sample ID`,'_',fixed=TRUE)))


trans_env <- data.frame("Index" = unique(t$X1),
                        "Code" = c("A","Bl","BLF","U","F","G","W"))


dta$new_sample_ID <-NA
dta$new_PCR_ID <- NA
count<-0
count_S<-0

r <- 1
for(r in 1:nrow(dta)){
  
  new_sample_ID_r <- c()
  if(t$X1[r] %in% c("Blank","Blank filter")){
    
    new_sample_ID_r <- t$X1[r] 
    
    # new_sample_ID_r <- paste0(new_sample_ID_r,"_",dta$Size[r])
    count <- count + 1
    dta$new_PCR_ID[r] <- paste0("Blank18_",count) 
  }else{
    
    new_sample_ID_r <- trans_env[trans_env$Index == t$X1[r],]$Code
    
    trans_no <- data.frame("new_nr" = 1:9,
                           "old_nr" = c("01","02","03","04","05","06","07","08","09"))
    
    if(t$X2[r] %in% trans_no$old_nr){
      t$X2[r] <- trans_no[trans_no$old_nr == t$X2[r],]$new_nr
    }
    
    new_sample_ID_r <- paste0(new_sample_ID_r,"_",t$X2[r])
    
    if(t$X3[r] == "12-15"){
      
      new_sample_ID_r <- paste0(new_sample_ID_r,"_","A")
    }else{
      #if "17-20"
      new_sample_ID_r <- paste0(new_sample_ID_r,"_","B")
      
    }#end(if else timepoint)
    
    
    if(dta$Size[r] == "-10mm"){
      new_sample_ID_r <- paste0(new_sample_ID_r,"_","S")
      
    }else{
      if(dta$Size[r] == "+10mm"){
        new_sample_ID_r <- paste0(new_sample_ID_r,"_","L")
      }
    }
    
    count_S <- count_S + 1
    
    dta$new_PCR_ID[r] <- paste0("S18_",count_S) 
    
  }#end if else Blank noblank
  
  new_sample_ID_r <- paste0(new_sample_ID_r,"_","2018")
  
  dta$new_sample_ID[r] <- new_sample_ID_r
  
  
}#end r

dta <- dta[,colnames(dta) %in% c("new_PCR_ID","new_sample_ID","Dry mass (mg)") | data.table::like(colnames(dta),"Qubit")]

t <- dta[,data.table::like(colnames(dta),"Qubit")]


dta <- dta[,c("new_PCR_ID","new_sample_ID","Dry mass (mg)")]

dta$Qubit <- t




setwd(paste0(WD_BASIS,"/WD")); fwrite(dta,paste0("DE_2018samples_03032020.csv"))









# 2019 Import#####

library(data.table)
library(stringr)

dta <- VH_excel(filename = "SampleList_Germany_2019.xlsx",path = paste0(WD_BASIS,"/InsektMobile_data"))[[5]]

dta <- dta[,colnames(dta) %in% c("no.","Sample ID","Size","Dry mass (mg)")]

dta_bit <- VH_excel(filename = "3y4DNAdilutions.xlsx",path = paste0(WD_BASIS,"/InsektMobile_data"))[[1]]

dta$Qubit <- dta_bit[,data.table::like(colnames(dta_bit),"Qubit")][1:nrow(dta)]

t <- data.frame(do.call('rbind', strsplit(dta$`Sample ID`,'_',fixed=TRUE)))

unique(t$X1)


dta$new_sample_ID <- NA
dta$new_PCR_ID <- NA
#manuel correction
t[t$X2 == "02:",]$X2 <- "02"
t[t$X2 == " 16",]$X2 <- "16"
t[t$X4 == "17-20<10mm",]$X4 <- "17-20"
t[t$X4 == "12",]$X4 <- "12-15"


r <- 20
for(r in 1:nrow(dta)){
  
  new_sample_ID_r <- c()
  if(t$X1[r] %in% c("Blank","blank")){
    
    new_sample_ID_r <- "Blank"
    
    #error left
    #tr <- gsub(".", "", as.character(dta$Size[r])); tr
    #str_replace(dta$Size[1], "mm", "_")
    
    #new_sample_ID_r <- paste0(new_sample_ID_r,"_",dta$Size[r])
    
    #new_sample_ID_r <-  str_replace(new_sample_ID_r, ".0", "0")
    #new_sample_ID_r <-  str_replace(new_sample_ID_r, ".2", "2")
    
    
    
  }else{
    
    if(t$X1[r] %in% c("Feucht","feucht","feuch")){
      new_sample_ID_r <- ("F")
      
    }else{
      if(t$X1[r] %in% c("Wald","wald")){
        new_sample_ID_r <- ("W")
        
      }else{
        if(t$X1[r] %in% c("agrar")){
          new_sample_ID_r <- "A"
          
        }else{
          if(data.table::like(t$X1[r],"gr")){
            new_sample_ID_r <- "G"
            
          }else{
            if(t$X1[r] %in% c("urban","Urban")){
              
              new_sample_ID_r <- "U"
              
              
            }
          }
        }
      }
    }
    
    
    
    
    
    trans <- data.frame("new_nr" = 1:9,
                        "old_nr" = c("01","02","03","04","05","06","07","08","09"))
    
    if(t$X2[r] %in% trans$old_nr){
      t$X2[r] <- trans[trans$old_nr == t$X2[r],]$new_nr
    }
    
    new_sample_ID_r <- paste0(new_sample_ID_r,"_",t$X2[r])
    
    
    
    
    if(t$X4[r] == "12-15"){
      
      new_sample_ID_r <- paste0(new_sample_ID_r,"_","A")
      
    }else{
      #if "17-20"
      new_sample_ID_r <- paste0(new_sample_ID_r,"_","B")
      
    }#end(if else timepoint)
    
    
    if(dta$Size[r] == "<10mm"){
      new_sample_ID_r <- paste0(new_sample_ID_r,"_","S")
      
    }else{
      if(dta$Size[r] == ">10mm"){
        new_sample_ID_r <- paste0(new_sample_ID_r,"_","L")
      }
    }
    #dta$new_PCR_ID[r] <- paste0("S19_",r) 
    
  }#end if ele Blank noblank
  
  new_sample_ID_r <- paste0(new_sample_ID_r,"_","2019")
  
  dta$new_sample_ID[r] <- new_sample_ID_r
  
  dta$new_PCR_ID[r] <- paste0("S19_",r) 
  
}#end r

dta <- dta[,colnames(dta) %in% c("new_PCR_ID","new_sample_ID","Dry mass (mg)","Qubit")]

dta <- dta[,c("new_PCR_ID","new_sample_ID","Dry mass (mg)","Qubit")]

setwd(paste0(WD_BASIS,"/WD")); fwrite(dta,paste0("SampleList_Germany_2019.csv"))









#Combine 2018 & 2019#####

library(data.table)

setwd(paste0(WD_BASIS,"/WD")); dta_18 <- fread(paste0("DE_2018samples_03032020.csv"))

setwd(paste0(WD_BASIS,"/WD"));dta_19 <- fread(paste0("SampleList_Germany_2019.csv"))

dta <- rbind(dta_18,dta_19)

setwd(paste0(WD_BASIS,"/WD")); fwrite(dta,paste0("03_lab_data_cleaned_DE.csv"))
