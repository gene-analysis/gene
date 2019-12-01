rm(list=ls())

library(glmnet)
library(igraph)
library(NetSwan)

high_grade_colon_mRNA <- read.csv("file path", row.names = 1)
low_grade_colon_mRNA <- read.csv("file path", row.names = 1)

## transformation
high_grade_colon_mRNA_1 <- log2(high_grade_colon_mRNA+1)
low_grade_colon_mRNA_1 <- log2(low_grade_colon_mRNA+1)

Samples_mean<-apply(high_grade_colon_mRNA_1,2,mean)   #Standardize the data
Samples_sd<-apply(high_grade_colon_mRNA_1,2,sd)
for( i in 1:length(Samples_mean) )
{
  high_grade_colon_mRNA_1[ ,i]<-(high_grade_colon_mRNA_1[ ,i]-Samples_mean[i])/Samples_sd[i]
}

Samples_mean<-apply(low_grade_colon_mRNA_1,2,mean)   #Standardize the data
Samples_sd<-apply(low_grade_colon_mRNA_1,2,sd)
for( i in 1:length(Samples_mean) )
{
  low_grade_colon_mRNA_1[ ,i]<-(low_grade_colon_mRNA_1[ ,i]-Samples_mean[i])/Samples_sd[i]
}

low_grade_colon_mRNA_standardize <- low_grade_colon_mRNA_1
high_grade_colon_mRNA_standardize <- high_grade_colon_mRNA_1


high_list <- as.data.frame(read.csv("file path"))
high_list_name <- colnames(high_list)
high_grade_colon_mRNA_standardize <- subset(high_grade_colon_mRNA_standardize, select=high_list_name)
low_grade_colon_mRNA_standardize <- subset(low_grade_colon_mRNA_standardize,select=high_list_name)


########### network estimation 1 ##########
Sample_size <- 0.5*(dim(low_grade_colon_mRNA_standardize)[1]+dim(high_grade_colon_mRNA_standardize)[1])
nV <- dim(low_grade_colon_mRNA_standardize)[2]

ALPHA<-0.1
penalty_lambda<-round(2/Sample_size^0.5 *qnorm(ALPHA/2/nV^2, mean = 0, sd = 1, lower.tail = FALSE),5)
penalty_lambda_2<-penalty_lambda/2

##low grade
low_grade_colon_mRNA_standardize <- as.matrix(low_grade_colon_mRNA_standardize)

result_table <- matrix(NA,nrow=1,ncol=3)
colnames(result_table) <- c("Y_gene","X_gene","coef")

for(index in c(1:nV))
{  
  yy<-low_grade_colon_mRNA_standardize[,index]                   # Decide yy (child) and xx (parents) in the LASSO type penalized LR
  xx<-low_grade_colon_mRNA_standardize[,-index]
  
  y_node<-colnames(low_grade_colon_mRNA_standardize)[index] 
  x_node<-colnames(low_grade_colon_mRNA_standardize)[-index]   
  
  glmnet_fit<-glmnet(xx, yy,family=c("gaussian"),lambda=penalty_lambda_2)
  inter_coef_list<-coef(glmnet_fit)
  coef_list<-inter_coef_list[-1]
  
  temp_TF <- coef_list!=0
  selected_coef_list <- coef_list[temp_TF]
  selected_x_node <- x_node[temp_TF]
  
  if( length(selected_x_node)!=0 )
  { 
    temp_result_table <- cbind(y_node,selected_x_node,selected_coef_list) 
    result_table <- rbind(result_table,temp_result_table)
  }
}   

result_table_low_grade <- result_table[-1,]

 
 # high grade
high_grade_colon_mRNA_standardize <- as.matrix(high_grade_colon_mRNA_standardize)
 
result_table <- matrix(NA,nrow=1,ncol=3)
colnames(result_table) <- c("Y_gene","X_gene","coef")
 
for(index in c(1:nV))
{
  yy<-high_grade_colon_mRNA_standardize[,index]                   # Decide yy (child) and xx (parents) in the LASSO type penalized LR
  xx<-high_grade_colon_mRNA_standardize[,-index]
  
  y_node<-colnames(high_grade_colon_mRNA_standardize)[index]
  x_node<-colnames(high_grade_colon_mRNA_standardize)[-index]
  
  glmnet_fit<-glmnet(xx, yy,family=c("gaussian"),lambda=penalty_lambda_2)
  inter_coef_list<-coef(glmnet_fit)
  coef_list<-inter_coef_list[-1]
  
  temp_TF <- coef_list!=0
  selected_coef_list <- coef_list[temp_TF]
  selected_x_node <- x_node[temp_TF]
  
  if( length(selected_x_node)!=0 )
  {
    temp_result_table <- cbind(y_node,selected_x_node,selected_coef_list)
    result_table <- rbind(result_table,temp_result_table)
  }
  
}
 
result_table_high_grade <- result_table[-1,]
 
write.csv(result_table_high_grade, "file path")
write.csv(result_table_low_grade, "file path")


########### network estimation 2 : lambda*0.5 ##########
Sample_size <- 0.5*(dim(low_grade_colon_mRNA_standardize)[1]+dim(high_grade_colon_mRNA_standardize)[1])
nV <- dim(low_grade_colon_mRNA_standardize)[2]

ALPHA<-0.1
penalty_lambda<-round(2/Sample_size^0.5 *qnorm(ALPHA/2/nV^2, mean = 0, sd = 1, lower.tail = FALSE),5)
penalty_lambda_2<-penalty_lambda/2
penalty_lambda_3<-penalty_lambda_2/2

# low grade
low_grade_colon_mRNA_standardize <- as.matrix(low_grade_colon_mRNA_standardize)

result_table <- matrix(NA,nrow=1,ncol=3)
colnames(result_table) <- c("Y_gene","X_gene","coef")

for(index in c(1:nV))
{  
  yy<-low_grade_colon_mRNA_standardize[,index]                   # Decide yy (child) and xx (parents) in the LASSO type penalized LR
  xx<-low_grade_colon_mRNA_standardize[,-index]
  
  y_node<-colnames(low_grade_colon_mRNA_standardize)[index] 
  x_node<-colnames(low_grade_colon_mRNA_standardize)[-index]   
  
  glmnet_fit<-glmnet(xx, yy,family=c("gaussian"),lambda=penalty_lambda_2)
  inter_coef_list<-coef(glmnet_fit)
  coef_list<-inter_coef_list[-1]
  
  temp_TF <- coef_list!=0
  selected_coef_list <- coef_list[temp_TF]
  selected_x_node <- x_node[temp_TF]
  
  if( length(selected_x_node)!=0 )
  { 
    temp_result_table <- cbind(y_node,selected_x_node,selected_coef_list) 
    result_table <- rbind(result_table,temp_result_table)
  }
}   

result_table_low_grade_lambda3 <- result_table[-1,]


# high grade
high_grade_colon_mRNA_standardize <- as.matrix(high_grade_colon_mRNA_standardize)

result_table <- matrix(NA,nrow=1,ncol=3)
colnames(result_table) <- c("Y_gene","X_gene","coef")

for(index in c(1:nV))
{  
  yy<-high_grade_colon_mRNA_standardize[,index]                   # Decide yy (child) and xx (parents) in the LASSO type penalized LR
  xx<-high_grade_colon_mRNA_standardize[,-index]
  
  y_node<-colnames(high_grade_colon_mRNA_standardize)[index] 
  x_node<-colnames(high_grade_colon_mRNA_standardize)[-index]   
  
  glmnet_fit<-glmnet(xx, yy,family=c("gaussian"),lambda=penalty_lambda_3)
  inter_coef_list<-coef(glmnet_fit)
  coef_list<-inter_coef_list[-1]
  
  temp_TF <- coef_list!=0
  selected_coef_list <- coef_list[temp_TF]
  selected_x_node <- x_node[temp_TF]
  
  if( length(selected_x_node)!=0 )
  { 
    temp_result_table <- cbind(y_node,selected_x_node,selected_coef_list) 
    result_table <- rbind(result_table,temp_result_table)
  }
  
}   

result_table_high_grade_lambda2 <- result_table[-1,]

write.csv(result_table_high_grade_lambda2, "file path")
write.csv(result_table_low_grade_lambda3, "file path")


########### network estimation 3 : lambda*0.25 ##########
Sample_size <- 0.5*(dim(low_grade_colon_mRNA_standardize)[1]+dim(high_grade_colon_mRNA_standardize)[1])
nV <- dim(low_grade_colon_mRNA_standardize)[2]

ALPHA<-0.1
penalty_lambda<-round(2/Sample_size^0.5 *qnorm(ALPHA/2/nV^2, mean = 0, sd = 1, lower.tail = FALSE),5)
penalty_lambda_2<-penalty_lambda/2
penalty_lambda_3<-penalty_lambda_2/2
penalty_lambda_4<-penalty_lambda_3/2

# low grade
low_grade_colon_mRNA_standardize <- as.matrix(low_grade_colon_mRNA_standardize)

result_table <- matrix(NA,nrow=1,ncol=3)
colnames(result_table) <- c("Y_gene","X_gene","coef")

for(index in c(1:nV))
{  
  yy<-low_grade_colon_mRNA_standardize[,index]                   # Decide yy (child) and xx (parents) in the LASSO type penalized LR
  xx<-low_grade_colon_mRNA_standardize[,-index]
  
  y_node<-colnames(low_grade_colon_mRNA_standardize)[index] 
  x_node<-colnames(low_grade_colon_mRNA_standardize)[-index]   
  
  glmnet_fit<-glmnet(xx, yy,family=c("gaussian"),lambda=penalty_lambda_4)
  inter_coef_list<-coef(glmnet_fit)
  coef_list<-inter_coef_list[-1]
  
  temp_TF <- coef_list!=0
  selected_coef_list <- coef_list[temp_TF]
  selected_x_node <- x_node[temp_TF]
  
  if( length(selected_x_node)!=0 )
  { 
    temp_result_table <- cbind(y_node,selected_x_node,selected_coef_list) 
    result_table <- rbind(result_table,temp_result_table)
  }
}   

result_table_low_grade_lambda4 <- result_table[-1,]


# high grade
high_grade_colon_mRNA_standardize <- as.matrix(high_grade_colon_mRNA_standardize)

result_table <- matrix(NA,nrow=1,ncol=3)
colnames(result_table) <- c("Y_gene","X_gene","coef")

for(index in c(1:nV))
{  
  yy<-high_grade_colon_mRNA_standardize[,index]                   # Decide yy (child) and xx (parents) in the LASSO type penalized LR
  xx<-high_grade_colon_mRNA_standardize[,-index]
  
  y_node<-colnames(high_grade_colon_mRNA_standardize)[index] 
  x_node<-colnames(high_grade_colon_mRNA_standardize)[-index]   
  
  glmnet_fit<-glmnet(xx, yy,family=c("gaussian"),lambda=penalty_lambda_4)
  inter_coef_list<-coef(glmnet_fit)
  coef_list<-inter_coef_list[-1]
  
  temp_TF <- coef_list!=0
  selected_coef_list <- coef_list[temp_TF]
  selected_x_node <- x_node[temp_TF]
  
  if( length(selected_x_node)!=0 )
  { 
    temp_result_table <- cbind(y_node,selected_x_node,selected_coef_list) 
    result_table <- rbind(result_table,temp_result_table)
  }
  
}   

result_table_high_grade_lambda4 <- result_table[-1,]

write.csv(result_table_high_grade_lambda4, "file path")
write.csv(result_table_low_grade_lambda4, "file path")

