
################ network estimation lambda  ###########################

Sample_size <- 0.5*(dim(normal_colon_mRNA_standardize)[1]+dim(tumor_colon_mRNA_standardize)[1])
nV <- dim(normal_colon_mRNA_standardize)[2]

ALPHA<-0.1
penalty_lambda<-round(2/Sample_size^0.5 *qnorm(ALPHA/2/nV^2, mean = 0, sd = 1, lower.tail = FALSE),5)
penalty_lambda_2<-penalty_lambda/2
penalty_lambda_4<-penalty_lambda_2/2

##normal
normal_colon_mRNA_standardize <- as.matrix(normal_colon_mRNA_standardize)

result_table <- matrix(NA,nrow=1,ncol=3)
colnames(result_table) <- c("Y_gene","X_gene","coef")

for(index in c(1:nV))
{  
  yy<-normal_colon_mRNA_standardize[,index]                   # Decide yy (child) and xx (parents) in the LASSO type penalized LR
  xx<-normal_colon_mRNA_standardize[,-index]
  
  y_node<-colnames(normal_colon_mRNA_standardize)[index] 
  x_node<-colnames(normal_colon_mRNA_standardize)[-index]   
  
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

result_table_normal <- result_table[-1,]


##tumor
tumor_colon_mRNA_standardize <- as.matrix(tumor_colon_mRNA_standardize)

result_table <- matrix(NA,nrow=1,ncol=3)
colnames(result_table) <- c("Y_gene","X_gene","coef")
nV1 <- dim(tumor_colon_mRNA_standardize)[2]
for(index in c(1:nV1))
{  
  yy<-tumor_colon_mRNA_standardize[,index]                   # Decide yy (child) and xx (parents) in the LASSO type penalized LR
  xx<-tumor_colon_mRNA_standardize[,-index]
  
  y_node<-colnames(tumor_colon_mRNA_standardize)[index] 
  x_node<-colnames(tumor_colon_mRNA_standardize)[-index]   
  
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

result_table_tumor <- result_table[-1,]
