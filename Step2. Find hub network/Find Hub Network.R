
######positive sort posi lambda2
rm(list=ls())
library(MASS)
library(glmnet) 

file_name=paste("file path", sep="")
BNC_SDPC_case1_4<- read.table(file_name, header=TRUE, sep=",")

BNC_SDPC_case1 <- matrix(NA,nrow=1,ncol=3)
index<-1

max_no <- dim(BNC_SDPC_case1_4)[1]

for( i in 1:max_no )
{
  
  Y_gene<-as.character(BNC_SDPC_case1_4[i,1])
  X_gene<-as.character(BNC_SDPC_case1_4[i,2])
  value_gene<-as.character(BNC_SDPC_case1_4[i,3])
  
  ############################
  
  BNC_SDPC_case1 <- rbind( BNC_SDPC_case1 , c(Y_gene,X_gene,value_gene) ) 
}

############################

BNC_SDPC_case1<-BNC_SDPC_case1[-1,]
file_name=paste("file path", sep="")
write.table(BNC_SDPC_case1,row.names=FALSE,col.names=FALSE, file=file_name,sep=",")

file_name=paste("file path", sep="")
BNC_SDPC_case1<- as.matrix(read.table(file_name, header=FALSE, sep=","))

exist_genes <- unique( c(BNC_SDPC_case1[,1],BNC_SDPC_case1[,2]) )

genes_neighbor_num <- matrix(NA,ncol=2)
genes_neighbor_list <- matrix(NA,ncol=100)

max_num <- length(exist_genes)

for( index in exist_genes )
{
  temp_TF <- index == BNC_SDPC_case1[,1]
  neighbor_1 <- BNC_SDPC_case1[temp_TF,2]
  
  temp_TF <- index == BNC_SDPC_case1[,2]
  neighbor_2 <- BNC_SDPC_case1[temp_TF,1]
  
  neighbors <- unique( c( neighbor_1,neighbor_2 ) )
  
  genes_neighbor_num <- rbind( genes_neighbor_num , c( index , length(neighbors) ) )
  
  index_neighbors <- c( index , neighbors )
  temp_genes_neighbor_list <- matrix("",ncol=100)
  temp_genes_neighbor_list[1:length(index_neighbors)] <- index_neighbors
  
  genes_neighbor_list <- rbind( genes_neighbor_list , temp_genes_neighbor_list )
  
}

genes_neighbor_list <- genes_neighbor_list[-1,]
genes_neighbor_num <- genes_neighbor_num[-1,]

genes_neighbor_total <- cbind( genes_neighbor_num , genes_neighbor_list )

col_names <- c("hub_gene","degree","hub_gene",paste("X",c(1:99),sep="") )


file_name=paste("file path", sep="")
write.table( genes_neighbor_total ,row.names=FALSE,col.names=col_names, file=file_name,sep=",")


file_name=paste("file path", sep="")
write.table( genes_neighbor_total ,row.names=FALSE,col.names=col_names, file=file_name,sep=",")

######positive sort posi lambda4
rm(list=ls())
library(MASS)
library(glmnet) 

file_name=paste("file path", sep="")
BNC_SDPC_case1_4<- read.table(file_name, header=TRUE, sep=",")

BNC_SDPC_case1 <- matrix(NA,nrow=1,ncol=3)
index<-1

max_no <- dim(BNC_SDPC_case1_4)[1]

for( i in 1:max_no )
{
  
  Y_gene<-as.character(BNC_SDPC_case1_4[i,1])
  X_gene<-as.character(BNC_SDPC_case1_4[i,2])
  value_gene<-as.character(BNC_SDPC_case1_4[i,3])
  
  ############################
  
  BNC_SDPC_case1 <- rbind( BNC_SDPC_case1 , c(Y_gene,X_gene,value_gene) ) 
}

############################

BNC_SDPC_case1<-BNC_SDPC_case1[-1,]
file_name=paste("file path", sep="")
write.table(BNC_SDPC_case1,row.names=FALSE,col.names=FALSE, file=file_name,sep=",")

file_name=paste("file path", sep="")
BNC_SDPC_case1<- as.matrix(read.table(file_name, header=FALSE, sep=","))

exist_genes <- unique( c(BNC_SDPC_case1[,1],BNC_SDPC_case1[,2]) )

genes_neighbor_num <- matrix(NA,ncol=2)
genes_neighbor_list <- matrix(NA,ncol=100)

max_num <- length(exist_genes)

for( index in exist_genes )
{
  temp_TF <- index == BNC_SDPC_case1[,1]
  neighbor_1 <- BNC_SDPC_case1[temp_TF,2]
  
  temp_TF <- index == BNC_SDPC_case1[,2]
  neighbor_2 <- BNC_SDPC_case1[temp_TF,1]
  
  neighbors <- unique( c( neighbor_1,neighbor_2 ) )
  
  genes_neighbor_num <- rbind( genes_neighbor_num , c( index , length(neighbors) ) )
  
  index_neighbors <- c( index , neighbors )
  temp_genes_neighbor_list <- matrix("",ncol=100)
  temp_genes_neighbor_list[1:length(index_neighbors)] <- index_neighbors
  
  genes_neighbor_list <- rbind( genes_neighbor_list , temp_genes_neighbor_list )
  
}

genes_neighbor_list <- genes_neighbor_list[-1,]
genes_neighbor_num <- genes_neighbor_num[-1,]

genes_neighbor_total <- cbind( genes_neighbor_num , genes_neighbor_list )

col_names <- c("hub_gene","degree","hub_gene",paste("X",c(1:99),sep="") )


file_name=paste("file path", sep="")
write.table( genes_neighbor_total ,row.names=FALSE,col.names=col_names, file=file_name,sep=",")

######negative sort posi lambda2
rm(list=ls())
library(MASS)
library(glmnet) 

file_name=paste("file path", sep="")
BNC_SDPC_case1_4<- read.table(file_name, header=TRUE, sep=",")

BNC_SDPC_case1 <- matrix(NA,nrow=1,ncol=3)
index<-1

max_no <- dim(BNC_SDPC_case1_4)[1]

for( i in 1:max_no )
{
  
  Y_gene<-as.character(BNC_SDPC_case1_4[i,1])
  X_gene<-as.character(BNC_SDPC_case1_4[i,2])
  value_gene<-as.character(BNC_SDPC_case1_4[i,3])
  
  ############################
  
  BNC_SDPC_case1 <- rbind( BNC_SDPC_case1 , c(Y_gene,X_gene,value_gene) ) 
}

############################

BNC_SDPC_case1<-BNC_SDPC_case1[-1,]
file_name=paste("file path", sep="")
write.table(BNC_SDPC_case1,row.names=FALSE,col.names=FALSE, file=file_name,sep=",")

file_name=paste("file path", sep="")
BNC_SDPC_case1<- as.matrix(read.table(file_name, header=FALSE, sep=","))

exist_genes <- unique( c(BNC_SDPC_case1[,1],BNC_SDPC_case1[,2]) )

genes_neighbor_num <- matrix(NA,ncol=2)
genes_neighbor_list <- matrix(NA,ncol=100)

max_num <- length(exist_genes)

for( index in exist_genes )
{
  temp_TF <- index == BNC_SDPC_case1[,1]
  neighbor_1 <- BNC_SDPC_case1[temp_TF,2]
  
  temp_TF <- index == BNC_SDPC_case1[,2]
  neighbor_2 <- BNC_SDPC_case1[temp_TF,1]
  
  neighbors <- unique( c( neighbor_1,neighbor_2 ) )
  
  genes_neighbor_num <- rbind( genes_neighbor_num , c( index , length(neighbors) ) )
  
  index_neighbors <- c( index , neighbors )
  temp_genes_neighbor_list <- matrix("",ncol=100)
  temp_genes_neighbor_list[1:length(index_neighbors)] <- index_neighbors
  
  genes_neighbor_list <- rbind( genes_neighbor_list , temp_genes_neighbor_list )
  
}

genes_neighbor_list <- genes_neighbor_list[-1,]
genes_neighbor_num <- genes_neighbor_num[-1,]

genes_neighbor_total <- cbind( genes_neighbor_num , genes_neighbor_list )

col_names <- c("hub_gene","degree","hub_gene",paste("X",c(1:99),sep="") )


file_name=paste("file path", sep="")
write.table( genes_neighbor_total ,row.names=FALSE,col.names=col_names, file=file_name,sep=",")

######negative sort posi lambda4
rm(list=ls())
library(MASS)
library(glmnet) 

file_name=paste("file path", sep="")
BNC_SDPC_case1_4<- read.table(file_name, header=TRUE, sep=",")

BNC_SDPC_case1 <- matrix(NA,nrow=1,ncol=3)
index<-1

max_no <- dim(BNC_SDPC_case1_4)[1]

for( i in 1:max_no )
{
  
  Y_gene<-as.character(BNC_SDPC_case1_4[i,1])
  X_gene<-as.character(BNC_SDPC_case1_4[i,2])
  value_gene<-as.character(BNC_SDPC_case1_4[i,3])
  
  ############################
  
  BNC_SDPC_case1 <- rbind( BNC_SDPC_case1 , c(Y_gene,X_gene,value_gene) ) 
}

############################

BNC_SDPC_case1<-BNC_SDPC_case1[-1,]
file_name=paste("file path", sep="")
write.table(BNC_SDPC_case1,row.names=FALSE,col.names=FALSE, file=file_name,sep=",")

file_name=paste("file path", sep="")
BNC_SDPC_case1<- as.matrix(read.table(file_name, header=FALSE, sep=","))

exist_genes <- unique( c(BNC_SDPC_case1[,1],BNC_SDPC_case1[,2]) )

genes_neighbor_num <- matrix(NA,ncol=2)
genes_neighbor_list <- matrix(NA,ncol=100)

max_num <- length(exist_genes)

for( index in exist_genes )
{
  temp_TF <- index == BNC_SDPC_case1[,1]
  neighbor_1 <- BNC_SDPC_case1[temp_TF,2]
  
  temp_TF <- index == BNC_SDPC_case1[,2]
  neighbor_2 <- BNC_SDPC_case1[temp_TF,1]
  
  neighbors <- unique( c( neighbor_1,neighbor_2 ) )
  
  genes_neighbor_num <- rbind( genes_neighbor_num , c( index , length(neighbors) ) )
  
  index_neighbors <- c( index , neighbors )
  temp_genes_neighbor_list <- matrix("",ncol=100)
  temp_genes_neighbor_list[1:length(index_neighbors)] <- index_neighbors
  
  genes_neighbor_list <- rbind( genes_neighbor_list , temp_genes_neighbor_list )
  
}

genes_neighbor_list <- genes_neighbor_list[-1,]
genes_neighbor_num <- genes_neighbor_num[-1,]

genes_neighbor_total <- cbind( genes_neighbor_num , genes_neighbor_list )

col_names <- c("hub_gene","degree","hub_gene",paste("X",c(1:99),sep="") )


file_name=paste("file path", sep="")
write.table( genes_neighbor_total ,row.names=FALSE,col.names=col_names, file=file_name,sep=",")


