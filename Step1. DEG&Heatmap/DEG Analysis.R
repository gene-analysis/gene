##comparison Degree centrality
rm(list=ls())
library(sqldf)
file_name_tumor=paste("file path", sep="")
tumor_Data<-read.table(file_name_tumor, header=TRUE, sep=",")
colnames(tumor_Data)

DEG_list <- as.data.frame(read.csv("file path"))

group<-tumor_Data[,"gene"]

Group_X_list <- c("high")
Group_Y_list <- c("low")

GROUP_X_PART<-group==Group_X_list
GROUP_Y_PART<-group==Group_Y_list

#ratio_GROUP_X_PART <- sum(GROUP_X_PART)/(sum(GROUP_X_PART)+sum(GROUP_Y_PART))
#ratio_GROUP_Y_PART <- sum(GROUP_Y_PART)/(sum(GROUP_X_PART)+sum(GROUP_Y_PART))

summary_table<-matrix("",nrow=1,ncol=5)

summary_group_1 <-c("DEG","median(log)","","","")
summary_group_2 <-c("","Lymph node-positive","Lymph node-negative","p-value","Relative gene expression level")
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "i"

for (i in DEG_list[c(1:1000),1]){
node_info <- as.data.frame(tumor_Data[,1])
colnames(node_info) <- c("name")
node_info <- sqldf('SELECT a.name, b.type
                     FROM node_info a
                     LEFT JOIN DEG_list b USING(name)')

new_i<- log2(tumor_Data[,i]+1)
gene_name<-sprintf("%s",i)

temp_summary_1<-cbind(median(new_i[GROUP_X_PART],na.rm=TRUE),median(new_i[GROUP_Y_PART],na.rm=TRUE))
level<- ifelse((median(new_i[GROUP_X_PART],na.rm=TRUE)-median(new_i[GROUP_Y_PART],na.rm=TRUE))>0,"Up","Down")

summary_group_1 <-cbind(gene_name,temp_summary_1,wilcox.test(new_i[GROUP_X_PART],new_i[GROUP_Y_PART])$p.value,level)
summary_table<-rbind(summary_table,summary_group_1)
}

file_name=paste("file path", sep="")
write.table(summary_table,row.names=FALSE,col.names=FALSE, file=file_name,sep=",")


###################################################

### "TEAD3"

TEAD3 <- log(tumor_Data[,"TEAD3"]+1)

temp_summary<-cbind(mean(TEAD3[GROUP_X_PART],na.rm=TRUE),-sd(TEAD3[GROUP_X_PART],na.rm=TRUE),mean(TEAD3[GROUP_Y_PART],na.rm=TRUE),-sd(TEAD3[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(TEAD3[GROUP_X_PART],na.rm=TRUE),"",median(TEAD3[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(TEAD3,na.rm=TRUE),-sd(TEAD3,na.rm=TRUE) )
temp_total_2<-cbind( median(TEAD3,na.rm=TRUE),"" )

summary_group_1 <-cbind("TEAD3","mean (SD)",temp_summary,t.test(TEAD3[GROUP_X_PART],TEAD3[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("TEAD3","median",temp_summary_2,wilcox.test(TEAD3[GROUP_X_PART],TEAD3[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "INTS10"

INTS10 <- log(tumor_Data[,"INTS10"]+1)

temp_summary<-cbind(mean(INTS10[GROUP_X_PART],na.rm=TRUE),-sd(INTS10[GROUP_X_PART],na.rm=TRUE),mean(INTS10[GROUP_Y_PART],na.rm=TRUE),-sd(INTS10[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(INTS10[GROUP_X_PART],na.rm=TRUE),"",median(INTS10[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(INTS10,na.rm=TRUE),-sd(INTS10,na.rm=TRUE) )
temp_total_2<-cbind( median(INTS10,na.rm=TRUE),"" )

summary_group_1 <-cbind("INTS10","mean (SD)",temp_summary,t.test(INTS10[GROUP_X_PART],INTS10[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("INTS10","median",temp_summary_2,wilcox.test(INTS10[GROUP_X_PART],INTS10[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "AGPAT5"

AGPAT5 <- log(tumor_Data[,"AGPAT5"]+1)

temp_summary<-cbind(mean(AGPAT5[GROUP_X_PART],na.rm=TRUE),-sd(AGPAT5[GROUP_X_PART],na.rm=TRUE),mean(AGPAT5[GROUP_Y_PART],na.rm=TRUE),-sd(AGPAT5[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(AGPAT5[GROUP_X_PART],na.rm=TRUE),"",median(AGPAT5[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(AGPAT5,na.rm=TRUE),-sd(AGPAT5,na.rm=TRUE) )
temp_total_2<-cbind( median(AGPAT5,na.rm=TRUE),"" )

summary_group_1 <-cbind("AGPAT5","mean (SD)",temp_summary,t.test(AGPAT5[GROUP_X_PART],AGPAT5[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("AGPAT5","median",temp_summary_2,wilcox.test(AGPAT5[GROUP_X_PART],AGPAT5[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "NAT1"

NAT1 <- log(tumor_Data[,"NAT1"]+1)

temp_summary<-cbind(mean(NAT1[GROUP_X_PART],na.rm=TRUE),-sd(NAT1[GROUP_X_PART],na.rm=TRUE),mean(NAT1[GROUP_Y_PART],na.rm=TRUE),-sd(NAT1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(NAT1[GROUP_X_PART],na.rm=TRUE),"",median(NAT1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(NAT1,na.rm=TRUE),-sd(NAT1,na.rm=TRUE) )
temp_total_2<-cbind( median(NAT1,na.rm=TRUE),"" )

summary_group_1 <-cbind("NAT1","mean (SD)",temp_summary,t.test(NAT1[GROUP_X_PART],NAT1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("NAT1","median",temp_summary_2,wilcox.test(NAT1[GROUP_X_PART],NAT1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "MINPP1"

MINPP1 <- log(tumor_Data[,"MINPP1"]+1)

temp_summary<-cbind(mean(MINPP1[GROUP_X_PART],na.rm=TRUE),-sd(MINPP1[GROUP_X_PART],na.rm=TRUE),mean(MINPP1[GROUP_Y_PART],na.rm=TRUE),-sd(MINPP1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(MINPP1[GROUP_X_PART],na.rm=TRUE),"",median(MINPP1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(MINPP1,na.rm=TRUE),-sd(MINPP1,na.rm=TRUE) )
temp_total_2<-cbind( median(MINPP1,na.rm=TRUE),"" )

summary_group_1 <-cbind("MINPP1","mean (SD)",temp_summary,t.test(MINPP1[GROUP_X_PART],MINPP1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("MINPP1","median",temp_summary_2,wilcox.test(MINPP1[GROUP_X_PART],MINPP1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "ERI1"

ERI1 <- log(tumor_Data[,"ERI1"]+1)

temp_summary<-cbind(mean(ERI1[GROUP_X_PART],na.rm=TRUE),-sd(ERI1[GROUP_X_PART],na.rm=TRUE),mean(ERI1[GROUP_Y_PART],na.rm=TRUE),-sd(ERI1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(ERI1[GROUP_X_PART],na.rm=TRUE),"",median(ERI1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(ERI1,na.rm=TRUE),-sd(ERI1,na.rm=TRUE) )
temp_total_2<-cbind( median(ERI1,na.rm=TRUE),"" )

summary_group_1 <-cbind("ERI1","mean (SD)",temp_summary,t.test(ERI1[GROUP_X_PART],ERI1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("ERI1","median",temp_summary_2,wilcox.test(ERI1[GROUP_X_PART],ERI1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "PBK"

PBK <- log(tumor_Data[,"PBK"]+1)

temp_summary<-cbind(mean(PBK[GROUP_X_PART],na.rm=TRUE),-sd(PBK[GROUP_X_PART],na.rm=TRUE),mean(PBK[GROUP_Y_PART],na.rm=TRUE),-sd(PBK[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(PBK[GROUP_X_PART],na.rm=TRUE),"",median(PBK[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(PBK,na.rm=TRUE),-sd(PBK,na.rm=TRUE) )
temp_total_2<-cbind( median(PBK,na.rm=TRUE),"" )

summary_group_1 <-cbind("PBK","mean (SD)",temp_summary,t.test(PBK[GROUP_X_PART],PBK[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("PBK","median",temp_summary_2,wilcox.test(PBK[GROUP_X_PART],PBK[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "GSR"

GSR <- log(tumor_Data[,"GSR"]+1)

temp_summary<-cbind(mean(GSR[GROUP_X_PART],na.rm=TRUE),-sd(GSR[GROUP_X_PART],na.rm=TRUE),mean(GSR[GROUP_Y_PART],na.rm=TRUE),-sd(GSR[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(GSR[GROUP_X_PART],na.rm=TRUE),"",median(GSR[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(GSR,na.rm=TRUE),-sd(GSR,na.rm=TRUE) )
temp_total_2<-cbind( median(GSR,na.rm=TRUE),"" )

summary_group_1 <-cbind("GSR","mean (SD)",temp_summary,t.test(GSR[GROUP_X_PART],GSR[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("GSR","median",temp_summary_2,wilcox.test(GSR[GROUP_X_PART],GSR[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "RGL2"

RGL2 <- log(tumor_Data[,"RGL2"]+1)

temp_summary<-cbind(mean(RGL2[GROUP_X_PART],na.rm=TRUE),-sd(RGL2[GROUP_X_PART],na.rm=TRUE),mean(RGL2[GROUP_Y_PART],na.rm=TRUE),-sd(RGL2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(RGL2[GROUP_X_PART],na.rm=TRUE),"",median(RGL2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(RGL2,na.rm=TRUE),-sd(RGL2,na.rm=TRUE) )
temp_total_2<-cbind( median(RGL2,na.rm=TRUE),"" )

summary_group_1 <-cbind("RGL2","mean (SD)",temp_summary,t.test(RGL2[GROUP_X_PART],RGL2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("RGL2","median",temp_summary_2,wilcox.test(RGL2[GROUP_X_PART],RGL2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "PDE12"

PDE12 <- log(tumor_Data[,"PDE12"]+1)

temp_summary<-cbind(mean(PDE12[GROUP_X_PART],na.rm=TRUE),-sd(PDE12[GROUP_X_PART],na.rm=TRUE),mean(PDE12[GROUP_Y_PART],na.rm=TRUE),-sd(PDE12[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(PDE12[GROUP_X_PART],na.rm=TRUE),"",median(PDE12[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(PDE12,na.rm=TRUE),-sd(PDE12,na.rm=TRUE) )
temp_total_2<-cbind( median(RGL2,na.rm=TRUE),"" )

summary_group_1 <-cbind("PDE12","mean (SD)",temp_summary,t.test(PDE12[GROUP_X_PART],PDE12[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("PDE12","median",temp_summary_2,wilcox.test(PDE12[GROUP_X_PART],PDE12[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "CDCA2"

CDCA2 <- log(tumor_Data[,"CDCA2"]+1)

temp_summary<-cbind(mean(CDCA2[GROUP_X_PART],na.rm=TRUE),-sd(CDCA2[GROUP_X_PART],na.rm=TRUE),mean(CDCA2[GROUP_Y_PART],na.rm=TRUE),-sd(CDCA2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(CDCA2[GROUP_X_PART],na.rm=TRUE),"",median(CDCA2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(CDCA2,na.rm=TRUE),-sd(CDCA2,na.rm=TRUE) )
temp_total_2<-cbind( median(CDCA2,na.rm=TRUE),"" )

summary_group_1 <-cbind("CDCA2","mean (SD)",temp_summary,t.test(CDCA2[GROUP_X_PART],CDCA2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("CDCA2","median",temp_summary_2,wilcox.test(CDCA2[GROUP_X_PART],CDCA2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "COQ2"

COQ2 <- log(tumor_Data[,"COQ2"]+1)

temp_summary<-cbind(mean(COQ2[GROUP_X_PART],na.rm=TRUE),-sd(COQ2[GROUP_X_PART],na.rm=TRUE),mean(COQ2[GROUP_Y_PART],na.rm=TRUE),-sd(COQ2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(COQ2[GROUP_X_PART],na.rm=TRUE),"",median(COQ2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(COQ2,na.rm=TRUE),-sd(COQ2,na.rm=TRUE) )
temp_total_2<-cbind( median(COQ2,na.rm=TRUE),"" )

summary_group_1 <-cbind("COQ2","mean (SD)",temp_summary,t.test(COQ2[GROUP_X_PART],COQ2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("COQ2","median",temp_summary_2,wilcox.test(COQ2[GROUP_X_PART],COQ2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "ESCO2"

ESCO2 <- log(tumor_Data[,"ESCO2"]+1)

temp_summary<-cbind(mean(ESCO2[GROUP_X_PART],na.rm=TRUE),-sd(ESCO2[GROUP_X_PART],na.rm=TRUE),mean(ESCO2[GROUP_Y_PART],na.rm=TRUE),-sd(ESCO2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(ESCO2[GROUP_X_PART],na.rm=TRUE),"",median(ESCO2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(ESCO2,na.rm=TRUE),-sd(ESCO2,na.rm=TRUE) )
temp_total_2<-cbind( median(ESCO2,na.rm=TRUE),"" )

summary_group_1 <-cbind("ESCO2","mean (SD)",temp_summary,t.test(ESCO2[GROUP_X_PART],ESCO2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("ESCO2","median",temp_summary_2,wilcox.test(ESCO2[GROUP_X_PART],ESCO2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "LYSMD2"

LYSMD2 <- log(tumor_Data[,"LYSMD2"]+1)

temp_summary<-cbind(mean(LYSMD2[GROUP_X_PART],na.rm=TRUE),-sd(LYSMD2[GROUP_X_PART],na.rm=TRUE),mean(LYSMD2[GROUP_Y_PART],na.rm=TRUE),-sd(LYSMD2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(LYSMD2[GROUP_X_PART],na.rm=TRUE),"",median(LYSMD2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(LYSMD2,na.rm=TRUE),-sd(LYSMD2,na.rm=TRUE) )
temp_total_2<-cbind( median(LYSMD2,na.rm=TRUE),"" )

summary_group_1 <-cbind("LYSMD2","mean (SD)",temp_summary,t.test(LYSMD2[GROUP_X_PART],LYSMD2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("LYSMD2","median",temp_summary_2,wilcox.test(LYSMD2[GROUP_X_PART],LYSMD2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "NUDT6"

NUDT6 <- log(tumor_Data[,"NUDT6"]+1)

temp_summary<-cbind(mean(NUDT6[GROUP_X_PART],na.rm=TRUE),-sd(NUDT6[GROUP_X_PART],na.rm=TRUE),mean(NUDT6[GROUP_Y_PART],na.rm=TRUE),-sd(NUDT6[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(NUDT6[GROUP_X_PART],na.rm=TRUE),"",median(NUDT6[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(NUDT6,na.rm=TRUE),-sd(NUDT6,na.rm=TRUE) )
temp_total_2<-cbind( median(NUDT6,na.rm=TRUE),"" )

summary_group_1 <-cbind("NUDT6","mean (SD)",temp_summary,t.test(NUDT6[GROUP_X_PART],NUDT6[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("NUDT6","median",temp_summary_2,wilcox.test(NUDT6[GROUP_X_PART],NUDT6[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "ITFG3"

ITFG3 <- log(tumor_Data[,"ITFG3"]+1)

temp_summary<-cbind(mean(ITFG3[GROUP_X_PART],na.rm=TRUE),-sd(ITFG3[GROUP_X_PART],na.rm=TRUE),mean(ITFG3[GROUP_Y_PART],na.rm=TRUE),-sd(ITFG3[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(ITFG3[GROUP_X_PART],na.rm=TRUE),"",median(ITFG3[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(ITFG3,na.rm=TRUE),-sd(ITFG3,na.rm=TRUE) )
temp_total_2<-cbind( median(ITFG3,na.rm=TRUE),"" )

summary_group_1 <-cbind("ITFG3","mean (SD)",temp_summary,t.test(ITFG3[GROUP_X_PART],ITFG3[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("ITFG3","median",temp_summary_2,wilcox.test(ITFG3[GROUP_X_PART],ITFG3[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "NEIL2"

NEIL2<- log(tumor_Data[,"NEIL2"]+1)

temp_summary<-cbind(mean(NEIL2[GROUP_X_PART],na.rm=TRUE),-sd(NEIL2[GROUP_X_PART],na.rm=TRUE),mean(NEIL2[GROUP_Y_PART],na.rm=TRUE),-sd(NEIL2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(NEIL2[GROUP_X_PART],na.rm=TRUE),"",median(NEIL2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(NEIL2,na.rm=TRUE),-sd(NEIL2,na.rm=TRUE) )
temp_total_2<-cbind( median(NEIL2,na.rm=TRUE),"" )

summary_group_1 <-cbind("NEIL2","mean (SD)",temp_summary,t.test(NEIL2[GROUP_X_PART],NEIL2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("NEIL2","median",temp_summary_2,wilcox.test(NEIL2[GROUP_X_PART],NEIL2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)




### "ATP6V1B2"

ATP6V1B2 <- log(tumor_Data[,"ATP6V1B2"]+1)

temp_summary<-cbind(mean(ATP6V1B2[GROUP_X_PART],na.rm=TRUE),-sd(ATP6V1B2[GROUP_X_PART],na.rm=TRUE),mean(ATP6V1B2[GROUP_Y_PART],na.rm=TRUE),-sd(ATP6V1B2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(ATP6V1B2[GROUP_X_PART],na.rm=TRUE),"",median(ATP6V1B2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(ATP6V1B2,na.rm=TRUE),-sd(ATP6V1B2,na.rm=TRUE) )
temp_total_2<-cbind( median(ATP6V1B2,na.rm=TRUE),"" )

summary_group_1 <-cbind("ATP6V1B2","mean (SD)",temp_summary,t.test(ATP6V1B2[GROUP_X_PART],ATP6V1B2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("ATP6V1B2","median",temp_summary_2,wilcox.test(ATP6V1B2[GROUP_X_PART],ATP6V1B2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "PINX1"

PINX1 <- log(tumor_Data[,"PINX1"]+1)

temp_summary<-cbind(mean(PINX1[GROUP_X_PART],na.rm=TRUE),-sd(PINX1[GROUP_X_PART],na.rm=TRUE),mean(PINX1[GROUP_Y_PART],na.rm=TRUE),-sd(PINX1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(PINX1[GROUP_X_PART],na.rm=TRUE),"",median(PINX1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(PINX1,na.rm=TRUE),-sd(PINX1,na.rm=TRUE) )
temp_total_2<-cbind( median(PINX1,na.rm=TRUE),"" )

summary_group_1 <-cbind("PINX1","mean (SD)",temp_summary,t.test(PINX1[GROUP_X_PART],PINX1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("PINX1","median",temp_summary_2,wilcox.test(PINX1[GROUP_X_PART],PINX1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "MDM2"

MDM2 <- log(tumor_Data[,"MDM2"]+1)

temp_summary<-cbind(mean(MDM2[GROUP_X_PART],na.rm=TRUE),-sd(MDM2[GROUP_X_PART],na.rm=TRUE),mean(MDM2[GROUP_Y_PART],na.rm=TRUE),-sd(MDM2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(MDM2[GROUP_X_PART],na.rm=TRUE),"",median(MDM2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(MDM2,na.rm=TRUE),-sd(MDM2,na.rm=TRUE) )
temp_total_2<-cbind( median(MDM2,na.rm=TRUE),"" )

summary_group_1 <-cbind("MDM2","mean (SD)",temp_summary,t.test(MDM2[GROUP_X_PART],MDM2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("MDM2","median",temp_summary_2,wilcox.test(MDM2[GROUP_X_PART],MDM2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "FASLG"

FASLG <- log(tumor_Data[,"FASLG"]+1)

temp_summary<-cbind(mean(FASLG[GROUP_X_PART],na.rm=TRUE),-sd(FASLG[GROUP_X_PART],na.rm=TRUE),mean(FASLG[GROUP_Y_PART],na.rm=TRUE),-sd(FASLG[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(FASLG[GROUP_X_PART],na.rm=TRUE),"",median(FASLG[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(FASLG,na.rm=TRUE),-sd(FASLG,na.rm=TRUE) )
temp_total_2<-cbind( median(FASLG,na.rm=TRUE),"" )

summary_group_1 <-cbind("FASLG","mean (SD)",temp_summary,t.test(FASLG[GROUP_X_PART],FASLG[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("FASLG","median",temp_summary_2,wilcox.test(FASLG[GROUP_X_PART],FASLG[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "C6orf15"

C6orf15 <- log(tumor_Data[,"C6orf15"]+1)

temp_summary<-cbind(mean(C6orf15[GROUP_X_PART],na.rm=TRUE),-sd(C6orf15[GROUP_X_PART],na.rm=TRUE),mean(C6orf15[GROUP_Y_PART],na.rm=TRUE),-sd(C6orf15[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(C6orf15[GROUP_X_PART],na.rm=TRUE),"",median(C6orf15[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(C6orf15,na.rm=TRUE),-sd(C6orf15,na.rm=TRUE) )
temp_total_2<-cbind( median(C6orf15,na.rm=TRUE),"" )

summary_group_1 <-cbind("C6orf15","mean (SD)",temp_summary,t.test(C6orf15[GROUP_X_PART],C6orf15[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("C6orf15","median",temp_summary_2,wilcox.test(C6orf15[GROUP_X_PART],C6orf15[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "MCPH1"

MCPH1 <- log(tumor_Data[,"MCPH1"]+1)

temp_summary<-cbind(mean(MCPH1[GROUP_X_PART],na.rm=TRUE),-sd(MCPH1[GROUP_X_PART],na.rm=TRUE),mean(MCPH1[GROUP_Y_PART],na.rm=TRUE),-sd(MCPH1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(MCPH1[GROUP_X_PART],na.rm=TRUE),"",median(MCPH1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(MCPH1,na.rm=TRUE),-sd(MCPH1,na.rm=TRUE) )
temp_total_2<-cbind( median(MCPH1,na.rm=TRUE),"" )

summary_group_1 <-cbind("MCPH1","mean (SD)",temp_summary,t.test(MCPH1[GROUP_X_PART],MCPH1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("MCPH1","median",temp_summary_2,wilcox.test(MCPH1[GROUP_X_PART],MCPH1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "BAT3"

BAT3 <- log(tumor_Data[,"BAT3"]+1)

temp_summary<-cbind(mean(BAT3[GROUP_X_PART],na.rm=TRUE),-sd(BAT3[GROUP_X_PART],na.rm=TRUE),mean(BAT3[GROUP_Y_PART],na.rm=TRUE),-sd(BAT3[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(BAT3[GROUP_X_PART],na.rm=TRUE),"",median(BAT3[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(BAT3,na.rm=TRUE),-sd(BAT3,na.rm=TRUE) )
temp_total_2<-cbind( median(BAT3,na.rm=TRUE),"" )

summary_group_1 <-cbind("BAT3","mean (SD)",temp_summary,t.test(BAT3[GROUP_X_PART],BAT3[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("BAT3","median",temp_summary_2,wilcox.test(BAT3[GROUP_X_PART],BAT3[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "ATF6B"

ATF6B <- log(tumor_Data[,"ATF6B"]+1)

temp_summary<-cbind(mean(ATF6B[GROUP_X_PART],na.rm=TRUE),-sd(ATF6B[GROUP_X_PART],na.rm=TRUE),mean(ATF6B[GROUP_Y_PART],na.rm=TRUE),-sd(ATF6B[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(ATF6B[GROUP_X_PART],na.rm=TRUE),"",median(ATF6B[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(ATF6B,na.rm=TRUE),-sd(ATF6B,na.rm=TRUE) )
temp_total_2<-cbind( median(ATF6B,na.rm=TRUE),"" )

summary_group_1 <-cbind("ATF6B","mean (SD)",temp_summary,t.test(ATF6B[GROUP_X_PART],ATF6B[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("ATF6B","median",temp_summary_2,wilcox.test(ATF6B[GROUP_X_PART],ATF6B[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "FAS"

FAS <- log(tumor_Data[,"FAS"]+1)

temp_summary<-cbind(mean(FAS[GROUP_X_PART],na.rm=TRUE),-sd(FAS[GROUP_X_PART],na.rm=TRUE),mean(FAS[GROUP_Y_PART],na.rm=TRUE),-sd(FAS[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(FAS[GROUP_X_PART],na.rm=TRUE),"",median(FAS[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(FAS,na.rm=TRUE),-sd(FAS,na.rm=TRUE) )
temp_total_2<-cbind( median(FAS,na.rm=TRUE),"" )

summary_group_1 <-cbind("FAS","mean (SD)",temp_summary,t.test(FAS[GROUP_X_PART],FAS[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("FAS","median",temp_summary_2,wilcox.test(FAS[GROUP_X_PART],FAS[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "ME2"

ME2 <- log(tumor_Data[,"ME2"]+1)

temp_summary<-cbind(mean(ME2[GROUP_X_PART],na.rm=TRUE),-sd(ME2[GROUP_X_PART],na.rm=TRUE),mean(ME2[GROUP_Y_PART],na.rm=TRUE),-sd(ME2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(ME2[GROUP_X_PART],na.rm=TRUE),"",median(ME2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(ME2,na.rm=TRUE),-sd(ME2,na.rm=TRUE) )
temp_total_2<-cbind( median(ME2,na.rm=TRUE),"" )

summary_group_1 <-cbind("ME2","mean (SD)",temp_summary,t.test(ME2[GROUP_X_PART],ME2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("ME2","median",temp_summary_2,wilcox.test(ME2[GROUP_X_PART],ME2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "CNOT7"

CNOT7 <- log(tumor_Data[,"CNOT7"]+1)

temp_summary<-cbind(mean(CNOT7[GROUP_X_PART],na.rm=TRUE),-sd(CNOT7[GROUP_X_PART],na.rm=TRUE),mean(CNOT7[GROUP_Y_PART],na.rm=TRUE),-sd(CNOT7[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(CNOT7[GROUP_X_PART],na.rm=TRUE),"",median(CNOT7[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(CNOT7,na.rm=TRUE),-sd(CNOT7,na.rm=TRUE) )
temp_total_2<-cbind( median(CNOT7,na.rm=TRUE),"" )

summary_group_1 <-cbind("CNOT7","mean (SD)",temp_summary,t.test(CNOT7[GROUP_X_PART],CNOT7[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("CNOT7","median",temp_summary_2,wilcox.test(CNOT7[GROUP_X_PART],CNOT7[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "ASAH1"

ASAH1 <- log(tumor_Data[,"ASAH1"]+1)

temp_summary<-cbind(mean(ASAH1[GROUP_X_PART],na.rm=TRUE),-sd(ASAH1[GROUP_X_PART],na.rm=TRUE),mean(ASAH1[GROUP_Y_PART],na.rm=TRUE),-sd(ASAH1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(ASAH1[GROUP_X_PART],na.rm=TRUE),"",median(ASAH1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(ASAH1,na.rm=TRUE),-sd(ASAH1,na.rm=TRUE) )
temp_total_2<-cbind( median(ASAH1,na.rm=TRUE),"" )

summary_group_1 <-cbind("ASAH1","mean (SD)",temp_summary,t.test(ASAH1[GROUP_X_PART],ASAH1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("ASAH1","median",temp_summary_2,wilcox.test(ASAH1[GROUP_X_PART],ASAH1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "CASP1"

CASP1 <- log(tumor_Data[,"CASP1"]+1)

temp_summary<-cbind(mean(CASP1[GROUP_X_PART],na.rm=TRUE),-sd(CASP1[GROUP_X_PART],na.rm=TRUE),mean(CASP1[GROUP_Y_PART],na.rm=TRUE),-sd(CASP1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(CASP1[GROUP_X_PART],na.rm=TRUE),"",median(CASP1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(CASP1,na.rm=TRUE),-sd(CASP1,na.rm=TRUE) )
temp_total_2<-cbind( median(CASP1,na.rm=TRUE),"" )

summary_group_1 <-cbind("CASP1","mean (SD)",temp_summary,t.test(CASP1[GROUP_X_PART],CASP1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("CASP1","median",temp_summary_2,wilcox.test(CASP1[GROUP_X_PART],CASP1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "RARA"

RARA <- log(tumor_Data[,"RARA"]+1)

temp_summary<-cbind(mean(RARA[GROUP_X_PART],na.rm=TRUE),-sd(RARA[GROUP_X_PART],na.rm=TRUE),mean(RARA[GROUP_Y_PART],na.rm=TRUE),-sd(RARA[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(RARA[GROUP_X_PART],na.rm=TRUE),"",median(RARA[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(RARA,na.rm=TRUE),-sd(RARA,na.rm=TRUE) )
temp_total_2<-cbind( median(RARA,na.rm=TRUE),"" )

summary_group_1 <-cbind("RARA","mean (SD)",temp_summary,t.test(RARA[GROUP_X_PART],RARA[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("RARA","median",temp_summary_2,wilcox.test(RARA[GROUP_X_PART],RARA[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "LEMD2"

LEMD2 <- log(tumor_Data[,"LEMD2"]+1)

temp_summary<-cbind(mean(LEMD2[GROUP_X_PART],na.rm=TRUE),-sd(LEMD2[GROUP_X_PART],na.rm=TRUE),mean(LEMD2[GROUP_Y_PART],na.rm=TRUE),-sd(LEMD2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(LEMD2[GROUP_X_PART],na.rm=TRUE),"",median(LEMD2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(LEMD2,na.rm=TRUE),-sd(LEMD2,na.rm=TRUE) )
temp_total_2<-cbind( median(LEMD2,na.rm=TRUE),"" )

summary_group_1 <-cbind("LEMD2","mean (SD)",temp_summary,t.test(LEMD2[GROUP_X_PART],LEMD2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("LEMD2","median",temp_summary_2,wilcox.test(LEMD2[GROUP_X_PART],LEMD2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "TCHH"

TCHH <- log(tumor_Data[,"TCHH"]+1)

temp_summary<-cbind(mean(TCHH[GROUP_X_PART],na.rm=TRUE),-sd(TCHH[GROUP_X_PART],na.rm=TRUE),mean(TCHH[GROUP_Y_PART],na.rm=TRUE),-sd(TCHH[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(TCHH[GROUP_X_PART],na.rm=TRUE),"",median(TCHH[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(TCHH,na.rm=TRUE),-sd(TCHH,na.rm=TRUE) )
temp_total_2<-cbind( median(TCHH,na.rm=TRUE),"" )

summary_group_1 <-cbind("TCHH","mean (SD)",temp_summary,t.test(TCHH[GROUP_X_PART],TCHH[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("TCHH","median",temp_summary_2,wilcox.test(TCHH[GROUP_X_PART],TCHH[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "XPO7"

XPO7 <- log(tumor_Data[,"XPO7"]+1)

temp_summary<-cbind(mean(XPO7[GROUP_X_PART],na.rm=TRUE),-sd(XPO7[GROUP_X_PART],na.rm=TRUE),mean(XPO7[GROUP_Y_PART],na.rm=TRUE),-sd(XPO7[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(XPO7[GROUP_X_PART],na.rm=TRUE),"",median(XPO7[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(XPO7,na.rm=TRUE),-sd(XPO7,na.rm=TRUE) )
temp_total_2<-cbind( median(XPO7,na.rm=TRUE),"" )

summary_group_1 <-cbind("XPO7","mean (SD)",temp_summary,t.test(XPO7[GROUP_X_PART],XPO7[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("XPO7","median",temp_summary_2,wilcox.test(XPO7[GROUP_X_PART],XPO7[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "CXCL11"

CXCL11 <- log(tumor_Data[,"CXCL11"]+1)

temp_summary<-cbind(mean(CXCL11[GROUP_X_PART],na.rm=TRUE),-sd(CXCL11[GROUP_X_PART],na.rm=TRUE),mean(CXCL11[GROUP_Y_PART],na.rm=TRUE),-sd(CXCL11[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(CXCL11[GROUP_X_PART],na.rm=TRUE),"",median(CXCL11[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(CXCL11,na.rm=TRUE),-sd(CXCL11,na.rm=TRUE) )
temp_total_2<-cbind( median(CXCL11,na.rm=TRUE),"" )

summary_group_1 <-cbind("CXCL11","mean (SD)",temp_summary,t.test(CXCL11[GROUP_X_PART],CXCL11[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("CXCL11","median",temp_summary_2,wilcox.test(CXCL11[GROUP_X_PART],CXCL11[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "PHLDB3"

PHLDB3 <- log(tumor_Data[,"PHLDB3"]+1)

temp_summary<-cbind(mean(PHLDB3[GROUP_X_PART],na.rm=TRUE),-sd(PHLDB3[GROUP_X_PART],na.rm=TRUE),mean(PHLDB3[GROUP_Y_PART],na.rm=TRUE),-sd(PHLDB3[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(PHLDB3[GROUP_X_PART],na.rm=TRUE),"",median(PHLDB3[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(PHLDB3,na.rm=TRUE),-sd(PHLDB3,na.rm=TRUE) )
temp_total_2<-cbind( median(PHLDB3,na.rm=TRUE),"" )

summary_group_1 <-cbind("PHLDB3","mean (SD)",temp_summary,t.test(PHLDB3[GROUP_X_PART],PHLDB3[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("PHLDB3","median",temp_summary_2,wilcox.test(PHLDB3[GROUP_X_PART],PHLDB3[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "CARKD"

CARKD <- log(tumor_Data[,"CARKD"]+1)

temp_summary<-cbind(mean(CARKD[GROUP_X_PART],na.rm=TRUE),-sd(CARKD[GROUP_X_PART],na.rm=TRUE),mean(CARKD[GROUP_Y_PART],na.rm=TRUE),-sd(CARKD[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(CARKD[GROUP_X_PART],na.rm=TRUE),"",median(CARKD[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(CARKD,na.rm=TRUE),-sd(CARKD,na.rm=TRUE) )
temp_total_2<-cbind( median(CARKD,na.rm=TRUE),"" )

summary_group_1 <-cbind("CARKD","mean (SD)",temp_summary,t.test(CARKD[GROUP_X_PART],CARKD[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("CARKD","median",temp_summary_2,wilcox.test(CARKD[GROUP_X_PART],CARKD[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "CAMK2B"

CAMK2B <- log(tumor_Data[,"CAMK2B"]+1)

temp_summary<-cbind(mean(CAMK2B[GROUP_X_PART],na.rm=TRUE),-sd(CAMK2B[GROUP_X_PART],na.rm=TRUE),mean(CAMK2B[GROUP_Y_PART],na.rm=TRUE),-sd(CAMK2B[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(CAMK2B[GROUP_X_PART],na.rm=TRUE),"",median(CAMK2B[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(CAMK2B,na.rm=TRUE),-sd(CAMK2B,na.rm=TRUE) )
temp_total_2<-cbind( median(CAMK2B,na.rm=TRUE),"" )

summary_group_1 <-cbind("CAMK2B","mean (SD)",temp_summary,t.test(CAMK2B[GROUP_X_PART],CAMK2B[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("CAMK2B","median",temp_summary_2,wilcox.test(CAMK2B[GROUP_X_PART],CAMK2B[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "APOL6"

APOL6 <- log(tumor_Data[,"APOL6"]+1)

temp_summary<-cbind(mean(APOL6[GROUP_X_PART],na.rm=TRUE),-sd(APOL6[GROUP_X_PART],na.rm=TRUE),mean(APOL6[GROUP_Y_PART],na.rm=TRUE),-sd(APOL6[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(APOL6[GROUP_X_PART],na.rm=TRUE),"",median(APOL6[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(APOL6,na.rm=TRUE),-sd(APOL6,na.rm=TRUE) )
temp_total_2<-cbind( median(APOL6,na.rm=TRUE),"" )

summary_group_1 <-cbind("APOL6","mean (SD)",temp_summary,t.test(APOL6[GROUP_X_PART],APOL6[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("APOL6","median",temp_summary_2,wilcox.test(APOL6[GROUP_X_PART],APOL6[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "ARHGDIG"

ARHGDIG <- log(tumor_Data[,"ARHGDIG"]+1)

temp_summary<-cbind(mean(ARHGDIG[GROUP_X_PART],na.rm=TRUE),-sd(ARHGDIG[GROUP_X_PART],na.rm=TRUE),mean(ARHGDIG[GROUP_Y_PART],na.rm=TRUE),-sd(ARHGDIG[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(ARHGDIG[GROUP_X_PART],na.rm=TRUE),"",median(ARHGDIG[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(ARHGDIG,na.rm=TRUE),-sd(ARHGDIG,na.rm=TRUE) )
temp_total_2<-cbind( median(ARHGDIG,na.rm=TRUE),"" )

summary_group_1 <-cbind("ARHGDIG","mean (SD)",temp_summary,t.test(ARHGDIG[GROUP_X_PART],ARHGDIG[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("ARHGDIG","median",temp_summary_2,wilcox.test(ARHGDIG[GROUP_X_PART],ARHGDIG[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "WIPI2"

WIPI2 <- log(tumor_Data[,"WIPI2"]+1)

temp_summary<-cbind(mean(WIPI2[GROUP_X_PART],na.rm=TRUE),-sd(WIPI2[GROUP_X_PART],na.rm=TRUE),mean(WIPI2[GROUP_Y_PART],na.rm=TRUE),-sd(WIPI2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(WIPI2[GROUP_X_PART],na.rm=TRUE),"",median(WIPI2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(WIPI2,na.rm=TRUE),-sd(WIPI2,na.rm=TRUE) )
temp_total_2<-cbind( median(WIPI2,na.rm=TRUE),"" )

summary_group_1 <-cbind("WIPI2","mean (SD)",temp_summary,t.test(WIPI2[GROUP_X_PART],WIPI2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("WIPI2","median",temp_summary_2,wilcox.test(WIPI2[GROUP_X_PART],WIPI2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "TRIM39"

TRIM39 <- log(tumor_Data[,"TRIM39"]+1)

temp_summary<-cbind(mean(TRIM39[GROUP_X_PART],na.rm=TRUE),-sd(TRIM39[GROUP_X_PART],na.rm=TRUE),mean(TRIM39[GROUP_Y_PART],na.rm=TRUE),-sd(TRIM39[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(TRIM39[GROUP_X_PART],na.rm=TRUE),"",median(TRIM39[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(TRIM39,na.rm=TRUE),-sd(TRIM39,na.rm=TRUE) )
temp_total_2<-cbind( median(TRIM39,na.rm=TRUE),"" )

summary_group_1 <-cbind("TRIM39","mean (SD)",temp_summary,t.test(TRIM39[GROUP_X_PART],TRIM39[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("TRIM39","median",temp_summary_2,wilcox.test(TRIM39[GROUP_X_PART],TRIM39[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "BEND3"

BEND3 <- log(tumor_Data[,"BEND3"]+1)

temp_summary<-cbind(mean(BEND3[GROUP_X_PART],na.rm=TRUE),-sd(BEND3[GROUP_X_PART],na.rm=TRUE),mean(BEND3[GROUP_Y_PART],na.rm=TRUE),-sd(BEND3[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(BEND3[GROUP_X_PART],na.rm=TRUE),"",median(BEND3[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(BEND3,na.rm=TRUE),-sd(BEND3,na.rm=TRUE) )
temp_total_2<-cbind( median(BEND3,na.rm=TRUE),"" )

summary_group_1 <-cbind("BEND3","mean (SD)",temp_summary,t.test(BEND3[GROUP_X_PART],BEND3[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("BEND3","median",temp_summary_2,wilcox.test(BEND3[GROUP_X_PART],BEND3[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "MYST1"

MYST1 <- log(tumor_Data[,"MYST1"]+1)

temp_summary<-cbind(mean(MYST1[GROUP_X_PART],na.rm=TRUE),-sd(MYST1[GROUP_X_PART],na.rm=TRUE),mean(MYST1[GROUP_Y_PART],na.rm=TRUE),-sd(MYST1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(MYST1[GROUP_X_PART],na.rm=TRUE),"",median(MYST1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(MYST1,na.rm=TRUE),-sd(MYST1,na.rm=TRUE) )
temp_total_2<-cbind( median(MYST1,na.rm=TRUE),"" )

summary_group_1 <-cbind("MYST1","mean (SD)",temp_summary,t.test(MYST1[GROUP_X_PART],MYST1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("MYST1","median",temp_summary_2,wilcox.test(MYST1[GROUP_X_PART],MYST1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "TYMS"

TYMS <- log(tumor_Data[,"TYMS"]+1)

temp_summary<-cbind(mean(TYMS[GROUP_X_PART],na.rm=TRUE),-sd(TYMS[GROUP_X_PART],na.rm=TRUE),mean(TYMS[GROUP_Y_PART],na.rm=TRUE),-sd(TYMS[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(TYMS[GROUP_X_PART],na.rm=TRUE),"",median(TYMS[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(TYMS,na.rm=TRUE),-sd(TYMS,na.rm=TRUE) )
temp_total_2<-cbind( median(TYMS,na.rm=TRUE),"" )

summary_group_1 <-cbind("TYMS","mean (SD)",temp_summary,t.test(TYMS[GROUP_X_PART],TYMS[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("TYMS","median",temp_summary_2,wilcox.test(TYMS[GROUP_X_PART],TYMS[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "CHAC2"

CHAC2 <- log(tumor_Data[,"CHAC2"]+1)

temp_summary<-cbind(mean(CHAC2[GROUP_X_PART],na.rm=TRUE),-sd(CHAC2[GROUP_X_PART],na.rm=TRUE),mean(CHAC2[GROUP_Y_PART],na.rm=TRUE),-sd(CHAC2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(CHAC2[GROUP_X_PART],na.rm=TRUE),"",median(CHAC2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(CHAC2,na.rm=TRUE),-sd(CHAC2,na.rm=TRUE) )
temp_total_2<-cbind( median(CHAC2,na.rm=TRUE),"" )

summary_group_1 <-cbind("CHAC2","mean (SD)",temp_summary,t.test(CHAC2[GROUP_X_PART],CHAC2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("CHAC2","median",temp_summary_2,wilcox.test(CHAC2[GROUP_X_PART],CHAC2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "FADS6"

FADS6 <- log(tumor_Data[,"FADS6"]+1)

temp_summary<-cbind(mean(FADS6[GROUP_X_PART],na.rm=TRUE),-sd(FADS6[GROUP_X_PART],na.rm=TRUE),mean(FADS6[GROUP_Y_PART],na.rm=TRUE),-sd(FADS6[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(FADS6[GROUP_X_PART],na.rm=TRUE),"",median(FADS6[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(FADS6,na.rm=TRUE),-sd(FADS6,na.rm=TRUE) )
temp_total_2<-cbind( median(FADS6,na.rm=TRUE),"" )

summary_group_1 <-cbind("FADS6","mean (SD)",temp_summary,t.test(FADS6[GROUP_X_PART],FADS6[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("FADS6","median",temp_summary_2,wilcox.test(FADS6[GROUP_X_PART],FADS6[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "FUT10"

FUT10 <- log(tumor_Data[,"FUT10"]+1)

temp_summary<-cbind(mean(FUT10[GROUP_X_PART],na.rm=TRUE),-sd(FUT10[GROUP_X_PART],na.rm=TRUE),mean(FUT10[GROUP_Y_PART],na.rm=TRUE),-sd(FUT10[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(FUT10[GROUP_X_PART],na.rm=TRUE),"",median(FUT10[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(FUT10,na.rm=TRUE),-sd(FUT10,na.rm=TRUE) )
temp_total_2<-cbind( median(FUT10,na.rm=TRUE),"" )

summary_group_1 <-cbind("FUT10","mean (SD)",temp_summary,t.test(FUT10[GROUP_X_PART],FUT10[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("FUT10","median",temp_summary_2,wilcox.test(FUT10[GROUP_X_PART],FUT10[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "PLAG1"

PLAG1 <- log(tumor_Data[,"PLAG1"]+1)

temp_summary<-cbind(mean(PLAG1[GROUP_X_PART],na.rm=TRUE),-sd(PLAG1[GROUP_X_PART],na.rm=TRUE),mean(PLAG1[GROUP_Y_PART],na.rm=TRUE),-sd(PLAG1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(PLAG1[GROUP_X_PART],na.rm=TRUE),"",median(PLAG1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(PLAG1,na.rm=TRUE),-sd(PLAG1,na.rm=TRUE) )
temp_total_2<-cbind( median(PLAG1,na.rm=TRUE),"" )

summary_group_1 <-cbind("PLAG1","mean (SD)",temp_summary,t.test(PLAG1[GROUP_X_PART],PLAG1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("PLAG1","median",temp_summary_2,wilcox.test(PLAG1[GROUP_X_PART],PLAG1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "RIMS4"

RIMS4 <- log(tumor_Data[,"RIMS4"]+1)

temp_summary<-cbind(mean(RIMS4[GROUP_X_PART],na.rm=TRUE),-sd(RIMS4[GROUP_X_PART],na.rm=TRUE),mean(RIMS4[GROUP_Y_PART],na.rm=TRUE),-sd(RIMS4[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(RIMS4[GROUP_X_PART],na.rm=TRUE),"",median(RIMS4[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(RIMS4,na.rm=TRUE),-sd(RIMS4,na.rm=TRUE) )
temp_total_2<-cbind( median(RIMS4,na.rm=TRUE),"" )

summary_group_1 <-cbind("RIMS4","mean (SD)",temp_summary,t.test(RIMS4[GROUP_X_PART],RIMS4[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("RIMS4","median",temp_summary_2,wilcox.test(RIMS4[GROUP_X_PART],RIMS4[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "C13orf15"

C13orf15 <- log(tumor_Data[,"C13orf15"]+1)

temp_summary<-cbind(mean(C13orf15[GROUP_X_PART],na.rm=TRUE),-sd(C13orf15[GROUP_X_PART],na.rm=TRUE),mean(C13orf15[GROUP_Y_PART],na.rm=TRUE),-sd(C13orf15[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(C13orf15[GROUP_X_PART],na.rm=TRUE),"",median(C13orf15[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(C13orf15,na.rm=TRUE),-sd(C13orf15,na.rm=TRUE) )
temp_total_2<-cbind( median(C13orf15,na.rm=TRUE),"" )

summary_group_1 <-cbind("C13orf15","mean (SD)",temp_summary,t.test(C13orf15[GROUP_X_PART],C13orf15[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("C13orf15","median",temp_summary_2,wilcox.test(C13orf15[GROUP_X_PART],C13orf15[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "C13orf15"

C13orf15 <- log(tumor_Data[,"C13orf15"]+1)

temp_summary<-cbind(mean(C13orf15[GROUP_X_PART],na.rm=TRUE),-sd(C13orf15[GROUP_X_PART],na.rm=TRUE),mean(C13orf15[GROUP_Y_PART],na.rm=TRUE),-sd(C13orf15[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(C13orf15[GROUP_X_PART],na.rm=TRUE),"",median(C13orf15[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(C13orf15,na.rm=TRUE),-sd(C13orf15,na.rm=TRUE) )
temp_total_2<-cbind( median(C13orf15,na.rm=TRUE),"" )

summary_group_1 <-cbind("C13orf15","mean (SD)",temp_summary,t.test(C13orf15[GROUP_X_PART],C13orf15[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("C13orf15","median",temp_summary_2,wilcox.test(C13orf15[GROUP_X_PART],C13orf15[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "SCN5A"

SCN5A <- log(tumor_Data[,"SCN5A"]+1)

temp_summary<-cbind(mean(SCN5A[GROUP_X_PART],na.rm=TRUE),-sd(SCN5A[GROUP_X_PART],na.rm=TRUE),mean(SCN5A[GROUP_Y_PART],na.rm=TRUE),-sd(SCN5A[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(SCN5A[GROUP_X_PART],na.rm=TRUE),"",median(SCN5A[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(SCN5A,na.rm=TRUE),-sd(SCN5A,na.rm=TRUE) )
temp_total_2<-cbind( median(SCN5A,na.rm=TRUE),"" )

summary_group_1 <-cbind("SCN5A","mean (SD)",temp_summary,t.test(SCN5A[GROUP_X_PART],SCN5A[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("SCN5A","median",temp_summary_2,wilcox.test(SCN5A[GROUP_X_PART],SCN5A[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "GTF2E2"

GTF2E2 <- log(tumor_Data[,"GTF2E2"]+1)

temp_summary<-cbind(mean(GTF2E2[GROUP_X_PART],na.rm=TRUE),-sd(GTF2E2[GROUP_X_PART],na.rm=TRUE),mean(GTF2E2[GROUP_Y_PART],na.rm=TRUE),-sd(GTF2E2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(GTF2E2[GROUP_X_PART],na.rm=TRUE),"",median(GTF2E2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(GTF2E2,na.rm=TRUE),-sd(GTF2E2,na.rm=TRUE) )
temp_total_2<-cbind( median(GTF2E2,na.rm=TRUE),"" )

summary_group_1 <-cbind("GTF2E2","mean (SD)",temp_summary,t.test(GTF2E2[GROUP_X_PART],GTF2E2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("GTF2E2","median",temp_summary_2,wilcox.test(GTF2E2[GROUP_X_PART],GTF2E2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "DDR1"

DDR1 <- log(tumor_Data[,"DDR1"]+1)

temp_summary<-cbind(mean(DDR1[GROUP_X_PART],na.rm=TRUE),-sd(DDR1[GROUP_X_PART],na.rm=TRUE),mean(DDR1[GROUP_Y_PART],na.rm=TRUE),-sd(DDR1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(DDR1[GROUP_X_PART],na.rm=TRUE),"",median(DDR1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(DDR1,na.rm=TRUE),-sd(DDR1,na.rm=TRUE) )
temp_total_2<-cbind( median(DDR1,na.rm=TRUE),"" )

summary_group_1 <-cbind("DDR1","mean (SD)",temp_summary,t.test(DDR1[GROUP_X_PART],DDR1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("DDR1","median",temp_summary_2,wilcox.test(DDR1[GROUP_X_PART],DDR1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "VPS37A"

VPS37A <- log(tumor_Data[,"VPS37A"]+1)

temp_summary<-cbind(mean(VPS37A[GROUP_X_PART],na.rm=TRUE),-sd(VPS37A[GROUP_X_PART],na.rm=TRUE),mean(VPS37A[GROUP_Y_PART],na.rm=TRUE),-sd(VPS37A[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(VPS37A[GROUP_X_PART],na.rm=TRUE),"",median(VPS37A[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(VPS37A,na.rm=TRUE),-sd(VPS37A,na.rm=TRUE) )
temp_total_2<-cbind( median(VPS37A,na.rm=TRUE),"" )

summary_group_1 <-cbind("VPS37A","mean (SD)",temp_summary,t.test(VPS37A[GROUP_X_PART],VPS37A[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("VPS37A","median",temp_summary_2,wilcox.test(VPS37A[GROUP_X_PART],VPS37A[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "KCNT1"

KCNT1 <- log(tumor_Data[,"KCNT1"]+1)

temp_summary<-cbind(mean(KCNT1[GROUP_X_PART],na.rm=TRUE),-sd(KCNT1[GROUP_X_PART],na.rm=TRUE),mean(KCNT1[GROUP_Y_PART],na.rm=TRUE),-sd(KCNT1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(KCNT1[GROUP_X_PART],na.rm=TRUE),"",median(KCNT1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(KCNT1,na.rm=TRUE),-sd(KCNT1,na.rm=TRUE) )
temp_total_2<-cbind( median(KCNT1,na.rm=TRUE),"" )

summary_group_1 <-cbind("KCNT1","mean (SD)",temp_summary,t.test(KCNT1[GROUP_X_PART],KCNT1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("KCNT1","median",temp_summary_2,wilcox.test(KCNT1[GROUP_X_PART],KCNT1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "HAP1"

HAP1 <- log(tumor_Data[,"HAP1"]+1)

temp_summary<-cbind(mean(HAP1[GROUP_X_PART],na.rm=TRUE),-sd(HAP1[GROUP_X_PART],na.rm=TRUE),mean(HAP1[GROUP_Y_PART],na.rm=TRUE),-sd(HAP1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(HAP1[GROUP_X_PART],na.rm=TRUE),"",median(HAP1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(HAP1,na.rm=TRUE),-sd(HAP1,na.rm=TRUE) )
temp_total_2<-cbind( median(HAP1,na.rm=TRUE),"" )

summary_group_1 <-cbind("HAP1","mean (SD)",temp_summary,t.test(HAP1[GROUP_X_PART],HAP1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("HAP1","median",temp_summary_2,wilcox.test(HAP1[GROUP_X_PART],HAP1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "FKBP9"

FKBP9 <- log(tumor_Data[,"FKBP9"]+1)

temp_summary<-cbind(mean(FKBP9[GROUP_X_PART],na.rm=TRUE),-sd(FKBP9[GROUP_X_PART],na.rm=TRUE),mean(FKBP9[GROUP_Y_PART],na.rm=TRUE),-sd(FKBP9[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(FKBP9[GROUP_X_PART],na.rm=TRUE),"",median(FKBP9[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(FKBP9,na.rm=TRUE),-sd(FKBP9,na.rm=TRUE) )
temp_total_2<-cbind( median(FKBP9,na.rm=TRUE),"" )

summary_group_1 <-cbind("FKBP9","mean (SD)",temp_summary,t.test(FKBP9[GROUP_X_PART],FKBP9[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("FKBP9","median",temp_summary_2,wilcox.test(FKBP9[GROUP_X_PART],FKBP9[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "DDAH2"

DDAH2 <- log(tumor_Data[,"DDAH2"]+1)

temp_summary<-cbind(mean(DDAH2[GROUP_X_PART],na.rm=TRUE),-sd(DDAH2[GROUP_X_PART],na.rm=TRUE),mean(DDAH2[GROUP_Y_PART],na.rm=TRUE),-sd(DDAH2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(DDAH2[GROUP_X_PART],na.rm=TRUE),"",median(DDAH2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(DDAH2,na.rm=TRUE),-sd(DDAH2,na.rm=TRUE) )
temp_total_2<-cbind( median(DDAH2,na.rm=TRUE),"" )

summary_group_1 <-cbind("DDAH2","mean (SD)",temp_summary,t.test(DDAH2[GROUP_X_PART],DDAH2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("DDAH2","median",temp_summary_2,wilcox.test(DDAH2[GROUP_X_PART],DDAH2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "PAX7"

PAX7 <- log(tumor_Data[,"PAX7"]+1)

temp_summary<-cbind(mean(PAX7[GROUP_X_PART],na.rm=TRUE),-sd(PAX7[GROUP_X_PART],na.rm=TRUE),mean(PAX7[GROUP_Y_PART],na.rm=TRUE),-sd(PAX7[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(PAX7[GROUP_X_PART],na.rm=TRUE),"",median(PAX7[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(PAX7,na.rm=TRUE),-sd(PAX7,na.rm=TRUE) )
temp_total_2<-cbind( median(PAX7,na.rm=TRUE),"" )

summary_group_1 <-cbind("PAX7","mean (SD)",temp_summary,t.test(PAX7[GROUP_X_PART],PAX7[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("PAX7","median",temp_summary_2,wilcox.test(PAX7[GROUP_X_PART],PAX7[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "MC1R"

MC1R <- log(tumor_Data[,"MC1R"]+1)

temp_summary<-cbind(mean(MC1R[GROUP_X_PART],na.rm=TRUE),-sd(MC1R[GROUP_X_PART],na.rm=TRUE),mean(MC1R[GROUP_Y_PART],na.rm=TRUE),-sd(MC1R[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(MC1R[GROUP_X_PART],na.rm=TRUE),"",median(MC1R[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(MC1R,na.rm=TRUE),-sd(MC1R,na.rm=TRUE) )
temp_total_2<-cbind( median(MC1R,na.rm=TRUE),"" )

summary_group_1 <-cbind("MC1R","mean (SD)",temp_summary,t.test(MC1R[GROUP_X_PART],MC1R[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("MC1R","median",temp_summary_2,wilcox.test(MC1R[GROUP_X_PART],MC1R[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "TNFRSF10A"

TNFRSF10A <- log(tumor_Data[,"TNFRSF10A"]+1)

temp_summary<-cbind(mean(TNFRSF10A[GROUP_X_PART],na.rm=TRUE),-sd(TNFRSF10A[GROUP_X_PART],na.rm=TRUE),mean(TNFRSF10A[GROUP_Y_PART],na.rm=TRUE),-sd(TNFRSF10A[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(TNFRSF10A[GROUP_X_PART],na.rm=TRUE),"",median(TNFRSF10A[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(TNFRSF10A,na.rm=TRUE),-sd(TNFRSF10A,na.rm=TRUE) )
temp_total_2<-cbind( median(TNFRSF10A,na.rm=TRUE),"" )

summary_group_1 <-cbind("TNFRSF10A","mean (SD)",temp_summary,t.test(TNFRSF10A[GROUP_X_PART],TNFRSF10A[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("TNFRSF10A","median",temp_summary_2,wilcox.test(TNFRSF10A[GROUP_X_PART],TNFRSF10A[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "CXCL10"

CXCL10 <- log(tumor_Data[,"CXCL10"]+1)

temp_summary<-cbind(mean(CXCL10[GROUP_X_PART],na.rm=TRUE),-sd(CXCL10[GROUP_X_PART],na.rm=TRUE),mean(CXCL10[GROUP_Y_PART],na.rm=TRUE),-sd(CXCL10[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(CXCL10[GROUP_X_PART],na.rm=TRUE),"",median(CXCL10[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(CXCL10,na.rm=TRUE),-sd(CXCL10,na.rm=TRUE) )
temp_total_2<-cbind( median(CXCL10,na.rm=TRUE),"" )

summary_group_1 <-cbind("CXCL10","mean (SD)",temp_summary,t.test(CXCL10[GROUP_X_PART],CXCL10[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("CXCL10","median",temp_summary_2,wilcox.test(CXCL10[GROUP_X_PART],CXCL10[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "LRCH4"

LRCH4 <- log(tumor_Data[,"LRCH4"]+1)

temp_summary<-cbind(mean(LRCH4[GROUP_X_PART],na.rm=TRUE),-sd(LRCH4[GROUP_X_PART],na.rm=TRUE),mean(LRCH4[GROUP_Y_PART],na.rm=TRUE),-sd(LRCH4[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(LRCH4[GROUP_X_PART],na.rm=TRUE),"",median(LRCH4[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(LRCH4,na.rm=TRUE),-sd(LRCH4,na.rm=TRUE) )
temp_total_2<-cbind( median(LRCH4,na.rm=TRUE),"" )

summary_group_1 <-cbind("LRCH4","mean (SD)",temp_summary,t.test(LRCH4[GROUP_X_PART],LRCH4[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("LRCH4","median",temp_summary_2,wilcox.test(LRCH4[GROUP_X_PART],LRCH4[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "C18orf55"

C18orf55 <- log(tumor_Data[,"C18orf55"]+1)

temp_summary<-cbind(mean(C18orf55[GROUP_X_PART],na.rm=TRUE),-sd(C18orf55[GROUP_X_PART],na.rm=TRUE),mean(C18orf55[GROUP_Y_PART],na.rm=TRUE),-sd(C18orf55[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(C18orf55[GROUP_X_PART],na.rm=TRUE),"",median(C18orf55[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(C18orf55,na.rm=TRUE),-sd(C18orf55,na.rm=TRUE) )
temp_total_2<-cbind( median(C18orf55,na.rm=TRUE),"" )

summary_group_1 <-cbind("C18orf55","mean (SD)",temp_summary,t.test(C18orf55[GROUP_X_PART],C18orf55[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("C18orf55","median",temp_summary_2,wilcox.test(C18orf55[GROUP_X_PART],C18orf55[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "LEPROTL1"

LEPROTL1 <- log(tumor_Data[,"LEPROTL1"]+1)

temp_summary<-cbind(mean(LEPROTL1[GROUP_X_PART],na.rm=TRUE),-sd(LEPROTL1[GROUP_X_PART],na.rm=TRUE),mean(LEPROTL1[GROUP_Y_PART],na.rm=TRUE),-sd(LEPROTL1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(LEPROTL1[GROUP_X_PART],na.rm=TRUE),"",median(LEPROTL1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(LEPROTL1,na.rm=TRUE),-sd(LEPROTL1,na.rm=TRUE) )
temp_total_2<-cbind( median(LEPROTL1,na.rm=TRUE),"" )

summary_group_1 <-cbind("LEPROTL1","mean (SD)",temp_summary,t.test(LEPROTL1[GROUP_X_PART],LEPROTL1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("LEPROTL1","median",temp_summary_2,wilcox.test(LEPROTL1[GROUP_X_PART],LEPROTL1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "TXNL1"

TXNL1 <- log(tumor_Data[,"TXNL1"]+1)

temp_summary<-cbind(mean(TXNL1[GROUP_X_PART],na.rm=TRUE),-sd(TXNL1[GROUP_X_PART],na.rm=TRUE),mean(TXNL1[GROUP_Y_PART],na.rm=TRUE),-sd(TXNL1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(TXNL1[GROUP_X_PART],na.rm=TRUE),"",median(TXNL1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(TXNL1,na.rm=TRUE),-sd(TXNL1,na.rm=TRUE) )
temp_total_2<-cbind( median(TXNL1,na.rm=TRUE),"" )

summary_group_1 <-cbind("TXNL1","mean (SD)",temp_summary,t.test(TXNL1[GROUP_X_PART],TXNL1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("TXNL1","median",temp_summary_2,wilcox.test(TXNL1[GROUP_X_PART],TXNL1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "COX10"

COX10 <- log(tumor_Data[,"COX10"]+1)

temp_summary<-cbind(mean(COX10[GROUP_X_PART],na.rm=TRUE),-sd(COX10[GROUP_X_PART],na.rm=TRUE),mean(COX10[GROUP_Y_PART],na.rm=TRUE),-sd(COX10[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(COX10[GROUP_X_PART],na.rm=TRUE),"",median(COX10[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(COX10,na.rm=TRUE),-sd(COX10,na.rm=TRUE) )
temp_total_2<-cbind( median(COX10,na.rm=TRUE),"" )

summary_group_1 <-cbind("COX10","mean (SD)",temp_summary,t.test(COX10[GROUP_X_PART],COX10[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("COX10","median",temp_summary_2,wilcox.test(COX10[GROUP_X_PART],COX10[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "C2CD4A"

C2CD4A <- log(tumor_Data[,"C2CD4A"]+1)

temp_summary<-cbind(mean(C2CD4A[GROUP_X_PART],na.rm=TRUE),-sd(C2CD4A[GROUP_X_PART],na.rm=TRUE),mean(C2CD4A[GROUP_Y_PART],na.rm=TRUE),-sd(C2CD4A[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(C2CD4A[GROUP_X_PART],na.rm=TRUE),"",median(C2CD4A[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(C2CD4A,na.rm=TRUE),-sd(C2CD4A,na.rm=TRUE) )
temp_total_2<-cbind( median(C2CD4A,na.rm=TRUE),"" )

summary_group_1 <-cbind("C2CD4A","mean (SD)",temp_summary,t.test(C2CD4A[GROUP_X_PART],C2CD4A[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("C2CD4A","median",temp_summary_2,wilcox.test(C2CD4A[GROUP_X_PART],C2CD4A[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "GLYR1"

GLYR1 <- log(tumor_Data[,"GLYR1"]+1)

temp_summary<-cbind(mean(GLYR1[GROUP_X_PART],na.rm=TRUE),-sd(GLYR1[GROUP_X_PART],na.rm=TRUE),mean(GLYR1[GROUP_Y_PART],na.rm=TRUE),-sd(GLYR1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(GLYR1[GROUP_X_PART],na.rm=TRUE),"",median(GLYR1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(GLYR1,na.rm=TRUE),-sd(GLYR1,na.rm=TRUE) )
temp_total_2<-cbind( median(GLYR1,na.rm=TRUE),"" )

summary_group_1 <-cbind("GLYR1","mean (SD)",temp_summary,t.test(GLYR1[GROUP_X_PART],GLYR1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("GLYR1","median",temp_summary_2,wilcox.test(GLYR1[GROUP_X_PART],GLYR1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "UPK2"

UPK2<- log(tumor_Data[,"UPK2"]+1)

temp_summary<-cbind(mean(UPK2[GROUP_X_PART],na.rm=TRUE),-sd(UPK2[GROUP_X_PART],na.rm=TRUE),mean(UPK2[GROUP_Y_PART],na.rm=TRUE),-sd(UPK2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(UPK2[GROUP_X_PART],na.rm=TRUE),"",median(UPK2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(UPK2,na.rm=TRUE),-sd(UPK2,na.rm=TRUE) )
temp_total_2<-cbind( median(UPK2,na.rm=TRUE),"" )

summary_group_1 <-cbind("UPK2","mean (SD)",temp_summary,t.test(UPK2[GROUP_X_PART],UPK2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("UPK2","median",temp_summary_2,wilcox.test(UPK2[GROUP_X_PART],UPK2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "IRF1"

IRF1 <- log(tumor_Data[,"IRF1"]+1)

temp_summary<-cbind(mean(IRF1[GROUP_X_PART],na.rm=TRUE),-sd(IRF1[GROUP_X_PART],na.rm=TRUE),mean(IRF1[GROUP_Y_PART],na.rm=TRUE),-sd(IRF1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(IRF1[GROUP_X_PART],na.rm=TRUE),"",median(IRF1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(IRF1,na.rm=TRUE),-sd(IRF1,na.rm=TRUE) )
temp_total_2<-cbind( median(IRF1,na.rm=TRUE),"" )

summary_group_1 <-cbind("IRF1","mean (SD)",temp_summary,t.test(IRF1[GROUP_X_PART],IRF1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("IRF1","median",temp_summary_2,wilcox.test(IRF1[GROUP_X_PART],IRF1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "PCGF2"

PCGF2 <- log(tumor_Data[,"PCGF2"]+1)

temp_summary<-cbind(mean(PCGF2[GROUP_X_PART],na.rm=TRUE),-sd(PCGF2[GROUP_X_PART],na.rm=TRUE),mean(PCGF2[GROUP_Y_PART],na.rm=TRUE),-sd(PCGF2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(PCGF2[GROUP_X_PART],na.rm=TRUE),"",median(PCGF2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(PCGF2,na.rm=TRUE),-sd(PCGF2,na.rm=TRUE) )
temp_total_2<-cbind( median(PCGF2,na.rm=TRUE),"" )

summary_group_1 <-cbind("PCGF2","mean (SD)",temp_summary,t.test(PCGF2[GROUP_X_PART],PCGF2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("PCGF2","median",temp_summary_2,wilcox.test(PCGF2[GROUP_X_PART],PCGF2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "IFNG"

IFNG <- log(tumor_Data[,"IFNG"]+1)

temp_summary<-cbind(mean(IFNG[GROUP_X_PART],na.rm=TRUE),-sd(IFNG[GROUP_X_PART],na.rm=TRUE),mean(IFNG[GROUP_Y_PART],na.rm=TRUE),-sd(IFNG[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(IFNG[GROUP_X_PART],na.rm=TRUE),"",median(IFNG[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(IFNG,na.rm=TRUE),-sd(IFNG,na.rm=TRUE) )
temp_total_2<-cbind( median(IFNG,na.rm=TRUE),"" )

summary_group_1 <-cbind("IFNG","mean (SD)",temp_summary,t.test(IFNG[GROUP_X_PART],IFNG[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("IFNG","median",temp_summary_2,wilcox.test(IFNG[GROUP_X_PART],IFNG[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "TMEM105"

TMEM105 <- log(tumor_Data[,"TMEM105"]+1)

temp_summary<-cbind(mean(TMEM105[GROUP_X_PART],na.rm=TRUE),-sd(TMEM105[GROUP_X_PART],na.rm=TRUE),mean(TMEM105[GROUP_Y_PART],na.rm=TRUE),-sd(TMEM105[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(TMEM105[GROUP_X_PART],na.rm=TRUE),"",median(TMEM105[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(TMEM105,na.rm=TRUE),-sd(TMEM105,na.rm=TRUE) )
temp_total_2<-cbind( median(TMEM105,na.rm=TRUE),"" )

summary_group_1 <-cbind("TMEM105","mean (SD)",temp_summary,t.test(TMEM105[GROUP_X_PART],TMEM105[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("TMEM105","median",temp_summary_2,wilcox.test(TMEM105[GROUP_X_PART],TMEM105[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "INTS9"

INTS9 <- log(tumor_Data[,"INTS9"]+1)

temp_summary<-cbind(mean(INTS9[GROUP_X_PART],na.rm=TRUE),-sd(INTS9[GROUP_X_PART],na.rm=TRUE),mean(INTS9[GROUP_Y_PART],na.rm=TRUE),-sd(INTS9[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(INTS9[GROUP_X_PART],na.rm=TRUE),"",median(INTS9[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(INTS9,na.rm=TRUE),-sd(INTS9,na.rm=TRUE) )
temp_total_2<-cbind( median(INTS9,na.rm=TRUE),"" )

summary_group_1 <-cbind("INTS9","mean (SD)",temp_summary,t.test(INTS9[GROUP_X_PART],INTS9[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("INTS9","median",temp_summary_2,wilcox.test(INTS9[GROUP_X_PART],INTS9[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "ZNF821"

ZNF821 <- log(tumor_Data[,"ZNF821"]+1)

temp_summary<-cbind(mean(ZNF821[GROUP_X_PART],na.rm=TRUE),-sd(ZNF821[GROUP_X_PART],na.rm=TRUE),mean(ZNF821[GROUP_Y_PART],na.rm=TRUE),-sd(ZNF821[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(ZNF821[GROUP_X_PART],na.rm=TRUE),"",median(ZNF821[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(ZNF821,na.rm=TRUE),-sd(ZNF821,na.rm=TRUE) )
temp_total_2<-cbind( median(ZNF821,na.rm=TRUE),"" )

summary_group_1 <-cbind("ZNF821","mean (SD)",temp_summary,t.test(ZNF821[GROUP_X_PART],ZNF821[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("ZNF821","median",temp_summary_2,wilcox.test(ZNF821[GROUP_X_PART],ZNF821[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "VPS52"

VPS52 <- log(tumor_Data[,"VPS52"]+1)

temp_summary<-cbind(mean(VPS52[GROUP_X_PART],na.rm=TRUE),-sd(VPS52[GROUP_X_PART],na.rm=TRUE),mean(VPS52[GROUP_Y_PART],na.rm=TRUE),-sd(VPS52[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(VPS52[GROUP_X_PART],na.rm=TRUE),"",median(VPS52[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(VPS52,na.rm=TRUE),-sd(VPS52,na.rm=TRUE) )
temp_total_2<-cbind( median(VPS52,na.rm=TRUE),"" )

summary_group_1 <-cbind("VPS52","mean (SD)",temp_summary,t.test(VPS52[GROUP_X_PART],VPS52[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("VPS52","median",temp_summary_2,wilcox.test(VPS52[GROUP_X_PART],VPS52[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "FOXP4"

FOXP4 <- log(tumor_Data[,"FOXP4"]+1)

temp_summary<-cbind(mean(FOXP4[GROUP_X_PART],na.rm=TRUE),-sd(FOXP4[GROUP_X_PART],na.rm=TRUE),mean(FOXP4[GROUP_Y_PART],na.rm=TRUE),-sd(FOXP4[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(FOXP4[GROUP_X_PART],na.rm=TRUE),"",median(FOXP4[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(FOXP4,na.rm=TRUE),-sd(FOXP4,na.rm=TRUE) )
temp_total_2<-cbind( median(FOXP4,na.rm=TRUE),"" )

summary_group_1 <-cbind("FOXP4","mean (SD)",temp_summary,t.test(FOXP4[GROUP_X_PART],FOXP4[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("FOXP4","median",temp_summary_2,wilcox.test(FOXP4[GROUP_X_PART],FOXP4[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "SH2D4A"

SH2D4A <- log(tumor_Data[,"SH2D4A"]+1)

temp_summary<-cbind(mean(SH2D4A[GROUP_X_PART],na.rm=TRUE),-sd(SH2D4A[GROUP_X_PART],na.rm=TRUE),mean(SH2D4A[GROUP_Y_PART],na.rm=TRUE),-sd(SH2D4A[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(SH2D4A[GROUP_X_PART],na.rm=TRUE),"",median(SH2D4A[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(SH2D4A,na.rm=TRUE),-sd(SH2D4A,na.rm=TRUE) )
temp_total_2<-cbind( median(SH2D4A,na.rm=TRUE),"" )

summary_group_1 <-cbind("SH2D4A","mean (SD)",temp_summary,t.test(SH2D4A[GROUP_X_PART],SH2D4A[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("SH2D4A","median",temp_summary_2,wilcox.test(SH2D4A[GROUP_X_PART],SH2D4A[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "LRP2"

LRP2 <- log(tumor_Data[,"LRP2"]+1)

temp_summary<-cbind(mean(LRP2[GROUP_X_PART],na.rm=TRUE),-sd(LRP2[GROUP_X_PART],na.rm=TRUE),mean(LRP2[GROUP_Y_PART],na.rm=TRUE),-sd(LRP2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(LRP2[GROUP_X_PART],na.rm=TRUE),"",median(LRP2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(LRP2,na.rm=TRUE),-sd(LRP2,na.rm=TRUE) )
temp_total_2<-cbind( median(LRP2,na.rm=TRUE),"" )

summary_group_1 <-cbind("LRP2","mean (SD)",temp_summary,t.test(LRP2[GROUP_X_PART],LRP2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("LRP2","median",temp_summary_2,wilcox.test(LRP2[GROUP_X_PART],LRP2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "EEF1A2"

EEF1A2 <- log(tumor_Data[,"EEF1A2"]+1)

temp_summary<-cbind(mean(EEF1A2[GROUP_X_PART],na.rm=TRUE),-sd(EEF1A2[GROUP_X_PART],na.rm=TRUE),mean(EEF1A2[GROUP_Y_PART],na.rm=TRUE),-sd(EEF1A2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(EEF1A2[GROUP_X_PART],na.rm=TRUE),"",median(EEF1A2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(EEF1A2,na.rm=TRUE),-sd(EEF1A2,na.rm=TRUE) )
temp_total_2<-cbind( median(EEF1A2,na.rm=TRUE),"" )

summary_group_1 <-cbind("EEF1A2","mean (SD)",temp_summary,t.test(EEF1A2[GROUP_X_PART],EEF1A2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("EEF1A2","median",temp_summary_2,wilcox.test(EEF1A2[GROUP_X_PART],EEF1A2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "HDAC5"

HDAC5 <- log(tumor_Data[,"HDAC5"]+1)

temp_summary<-cbind(mean(HDAC5[GROUP_X_PART],na.rm=TRUE),-sd(HDAC5[GROUP_X_PART],na.rm=TRUE),mean(HDAC5[GROUP_Y_PART],na.rm=TRUE),-sd(HDAC5[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(HDAC5[GROUP_X_PART],na.rm=TRUE),"",median(HDAC5[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(HDAC5,na.rm=TRUE),-sd(HDAC5,na.rm=TRUE) )
temp_total_2<-cbind( median(HDAC5,na.rm=TRUE),"" )

summary_group_1 <-cbind("HDAC5","mean (SD)",temp_summary,t.test(HDAC5[GROUP_X_PART],HDAC5[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("HDAC5","median",temp_summary_2,wilcox.test(HDAC5[GROUP_X_PART],HDAC5[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "AXIN1"

AXIN1 <- log(tumor_Data[,"AXIN1"]+1)

temp_summary<-cbind(mean(AXIN1[GROUP_X_PART],na.rm=TRUE),-sd(AXIN1[GROUP_X_PART],na.rm=TRUE),mean(AXIN1[GROUP_Y_PART],na.rm=TRUE),-sd(AXIN1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(AXIN1[GROUP_X_PART],na.rm=TRUE),"",median(AXIN1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(AXIN1,na.rm=TRUE),-sd(AXIN1,na.rm=TRUE) )
temp_total_2<-cbind( median(AXIN1,na.rm=TRUE),"" )

summary_group_1 <-cbind("AXIN1","mean (SD)",temp_summary,t.test(AXIN1[GROUP_X_PART],AXIN1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("AXIN1","median",temp_summary_2,wilcox.test(AXIN1[GROUP_X_PART],AXIN1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "ELP3"

ELP3 <- log(tumor_Data[,"ELP3"]+1)

temp_summary<-cbind(mean(ELP3[GROUP_X_PART],na.rm=TRUE),-sd(ELP3[GROUP_X_PART],na.rm=TRUE),mean(ELP3[GROUP_Y_PART],na.rm=TRUE),-sd(ELP3[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(ELP3[GROUP_X_PART],na.rm=TRUE),"",median(ELP3[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(ELP3,na.rm=TRUE),-sd(ELP3,na.rm=TRUE) )
temp_total_2<-cbind( median(ELP3,na.rm=TRUE),"" )

summary_group_1 <-cbind("ELP3","mean (SD)",temp_summary,t.test(ELP3[GROUP_X_PART],ELP3[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("ELP3","median",temp_summary_2,wilcox.test(ELP3[GROUP_X_PART],ELP3[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "STXBP5L"

STXBP5L <- log(tumor_Data[,"STXBP5L"]+1)

temp_summary<-cbind(mean(STXBP5L[GROUP_X_PART],na.rm=TRUE),-sd(STXBP5L[GROUP_X_PART],na.rm=TRUE),mean(STXBP5L[GROUP_Y_PART],na.rm=TRUE),-sd(STXBP5L[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(STXBP5L[GROUP_X_PART],na.rm=TRUE),"",median(STXBP5L[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(STXBP5L,na.rm=TRUE),-sd(STXBP5L,na.rm=TRUE) )
temp_total_2<-cbind( median(STXBP5L,na.rm=TRUE),"" )

summary_group_1 <-cbind("STXBP5L","mean (SD)",temp_summary,t.test(STXBP5L[GROUP_X_PART],STXBP5L[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("STXBP5L","median",temp_summary_2,wilcox.test(STXBP5L[GROUP_X_PART],STXBP5L[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "EXD2"

EXD2 <- log(tumor_Data[,"EXD2"]+1)

temp_summary<-cbind(mean(EXD2[GROUP_X_PART],na.rm=TRUE),-sd(EXD2[GROUP_X_PART],na.rm=TRUE),mean(EXD2[GROUP_Y_PART],na.rm=TRUE),-sd(EXD2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(EXD2[GROUP_X_PART],na.rm=TRUE),"",median(EXD2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(EXD2,na.rm=TRUE),-sd(EXD2,na.rm=TRUE) )
temp_total_2<-cbind( median(EXD2,na.rm=TRUE),"" )

summary_group_1 <-cbind("EXD2","mean (SD)",temp_summary,t.test(EXD2[GROUP_X_PART],EXD2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("EXD2","median",temp_summary_2,wilcox.test(EXD2[GROUP_X_PART],EXD2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "B3GALT4"

B3GALT4 <- log(tumor_Data[,"B3GALT4"]+1)

temp_summary<-cbind(mean(B3GALT4[GROUP_X_PART],na.rm=TRUE),-sd(B3GALT4[GROUP_X_PART],na.rm=TRUE),mean(B3GALT4[GROUP_Y_PART],na.rm=TRUE),-sd(B3GALT4[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(B3GALT4[GROUP_X_PART],na.rm=TRUE),"",median(B3GALT4[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(B3GALT4,na.rm=TRUE),-sd(B3GALT4,na.rm=TRUE) )
temp_total_2<-cbind( median(B3GALT4,na.rm=TRUE),"" )

summary_group_1 <-cbind("B3GALT4","mean (SD)",temp_summary,t.test(B3GALT4[GROUP_X_PART],B3GALT4[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("B3GALT4","median",temp_summary_2,wilcox.test(B3GALT4[GROUP_X_PART],B3GALT4[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "PPP2R2A"

PPP2R2A <- log(tumor_Data[,"PPP2R2A"]+1)

temp_summary<-cbind(mean(PPP2R2A[GROUP_X_PART],na.rm=TRUE),-sd(PPP2R2A[GROUP_X_PART],na.rm=TRUE),mean(PPP2R2A[GROUP_Y_PART],na.rm=TRUE),-sd(PPP2R2A[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(PPP2R2A[GROUP_X_PART],na.rm=TRUE),"",median(PPP2R2A[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(PPP2R2A,na.rm=TRUE),-sd(PPP2R2A,na.rm=TRUE) )
temp_total_2<-cbind( median(PPP2R2A,na.rm=TRUE),"" )

summary_group_1 <-cbind("PPP2R2A","mean (SD)",temp_summary,t.test(PPP2R2A[GROUP_X_PART],PPP2R2A[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("PPP2R2A","median",temp_summary_2,wilcox.test(PPP2R2A[GROUP_X_PART],PPP2R2A[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "NPR3"

NPR3 <- log(tumor_Data[,"NPR3"]+1)

temp_summary<-cbind(mean(NPR3[GROUP_X_PART],na.rm=TRUE),-sd(NPR3[GROUP_X_PART],na.rm=TRUE),mean(NPR3[GROUP_Y_PART],na.rm=TRUE),-sd(NPR3[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(NPR3[GROUP_X_PART],na.rm=TRUE),"",median(NPR3[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(NPR3,na.rm=TRUE),-sd(NPR3,na.rm=TRUE) )
temp_total_2<-cbind( median(NPR3,na.rm=TRUE),"" )

summary_group_1 <-cbind("NPR3","mean (SD)",temp_summary,t.test(NPR3[GROUP_X_PART],NPR3[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("NPR3","median",temp_summary_2,wilcox.test(NPR3[GROUP_X_PART],NPR3[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "TNFRSF11A"

TNFRSF11A <- log(tumor_Data[,"TNFRSF11A"]+1)

temp_summary<-cbind(mean(TNFRSF11A[GROUP_X_PART],na.rm=TRUE),-sd(TNFRSF11A[GROUP_X_PART],na.rm=TRUE),mean(TNFRSF11A[GROUP_Y_PART],na.rm=TRUE),-sd(TNFRSF11A[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(TNFRSF11A[GROUP_X_PART],na.rm=TRUE),"",median(TNFRSF11A[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(TNFRSF11A,na.rm=TRUE),-sd(TNFRSF11A,na.rm=TRUE) )
temp_total_2<-cbind( median(TNFRSF11A,na.rm=TRUE),"" )

summary_group_1 <-cbind("TNFRSF11A","mean (SD)",temp_summary,t.test(TNFRSF11A[GROUP_X_PART],TNFRSF11A[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("TNFRSF11A","median",temp_summary_2,wilcox.test(TNFRSF11A[GROUP_X_PART],TNFRSF11A[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "FBRS"

FBRS <- log(tumor_Data[,"FBRS"]+1)

temp_summary<-cbind(mean(FBRS[GROUP_X_PART],na.rm=TRUE),-sd(FBRS[GROUP_X_PART],na.rm=TRUE),mean(FBRS[GROUP_Y_PART],na.rm=TRUE),-sd(FBRS[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(FBRS[GROUP_X_PART],na.rm=TRUE),"",median(FBRS[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(FBRS,na.rm=TRUE),-sd(FBRS,na.rm=TRUE) )
temp_total_2<-cbind( median(FBRS,na.rm=TRUE),"" )

summary_group_1 <-cbind("FBRS","mean (SD)",temp_summary,t.test(FBRS[GROUP_X_PART],FBRS[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("FBRS","median",temp_summary_2,wilcox.test(FBRS[GROUP_X_PART],FBRS[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "NR1D1"

NR1D1 <- log(tumor_Data[,"NR1D1"]+1)

temp_summary<-cbind(mean(NR1D1[GROUP_X_PART],na.rm=TRUE),-sd(NR1D1[GROUP_X_PART],na.rm=TRUE),mean(NR1D1[GROUP_Y_PART],na.rm=TRUE),-sd(NR1D1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(NR1D1[GROUP_X_PART],na.rm=TRUE),"",median(NR1D1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(NR1D1,na.rm=TRUE),-sd(NR1D1,na.rm=TRUE) )
temp_total_2<-cbind( median(NR1D1,na.rm=TRUE),"" )

summary_group_1 <-cbind("NR1D1","mean (SD)",temp_summary,t.test(NR1D1[GROUP_X_PART],NR1D1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("NR1D1","median",temp_summary_2,wilcox.test(NR1D1[GROUP_X_PART],NR1D1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "PALM2"

PALM2 <- log(tumor_Data[,"PALM2"]+1)

temp_summary<-cbind(mean(PALM2[GROUP_X_PART],na.rm=TRUE),-sd(PALM2[GROUP_X_PART],na.rm=TRUE),mean(PALM2[GROUP_Y_PART],na.rm=TRUE),-sd(PALM2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(PALM2[GROUP_X_PART],na.rm=TRUE),"",median(PALM2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(PALM2,na.rm=TRUE),-sd(PALM2,na.rm=TRUE) )
temp_total_2<-cbind( median(PALM2,na.rm=TRUE),"" )

summary_group_1 <-cbind("PALM2","mean (SD)",temp_summary,t.test(PALM2[GROUP_X_PART],PALM2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("PALM2","median",temp_summary_2,wilcox.test(PALM2[GROUP_X_PART],PALM2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "HEXB"

HEXB <- log(tumor_Data[,"HEXB"]+1)

temp_summary<-cbind(mean(HEXB[GROUP_X_PART],na.rm=TRUE),-sd(HEXB[GROUP_X_PART],na.rm=TRUE),mean(HEXB[GROUP_Y_PART],na.rm=TRUE),-sd(HEXB[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(HEXB[GROUP_X_PART],na.rm=TRUE),"",median(HEXB[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(HEXB,na.rm=TRUE),-sd(HEXB,na.rm=TRUE) )
temp_total_2<-cbind( median(HEXB,na.rm=TRUE),"" )

summary_group_1 <-cbind("HEXB","mean (SD)",temp_summary,t.test(HEXB[GROUP_X_PART],HEXB[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("HEXB","median",temp_summary_2,wilcox.test(HEXB[GROUP_X_PART],HEXB[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "TMEM198"

TMEM198 <- log(tumor_Data[,"TMEM198"]+1)

temp_summary<-cbind(mean(TMEM198[GROUP_X_PART],na.rm=TRUE),-sd(TMEM198[GROUP_X_PART],na.rm=TRUE),mean(TMEM198[GROUP_Y_PART],na.rm=TRUE),-sd(TMEM198[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(TMEM198[GROUP_X_PART],na.rm=TRUE),"",median(TMEM198[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(TMEM198,na.rm=TRUE),-sd(TMEM198,na.rm=TRUE) )
temp_total_2<-cbind( median(TMEM198,na.rm=TRUE),"" )

summary_group_1 <-cbind("TMEM198","mean (SD)",temp_summary,t.test(TMEM198[GROUP_X_PART],TMEM198[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("TMEM198","median",temp_summary_2,wilcox.test(TMEM198[GROUP_X_PART],TMEM198[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

### "SCO2"

SCO2 <- log(tumor_Data[,"SCO2"]+1)

temp_summary<-cbind(mean(SCO2[GROUP_X_PART],na.rm=TRUE),-sd(SCO2[GROUP_X_PART],na.rm=TRUE),mean(SCO2[GROUP_Y_PART],na.rm=TRUE),-sd(SCO2[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(SCO2[GROUP_X_PART],na.rm=TRUE),"",median(SCO2[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(SCO2,na.rm=TRUE),-sd(SCO2,na.rm=TRUE) )
temp_total_2<-cbind( median(SCO2,na.rm=TRUE),"" )

summary_group_1 <-cbind("SCO2","mean (SD)",temp_summary,t.test(SCO2[GROUP_X_PART],SCO2[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("SCO2","median",temp_summary_2,wilcox.test(SCO2[GROUP_X_PART],SCO2[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "TCF7L1"

TCF7L1 <- log(tumor_Data[,"TCF7L1"]+1)

temp_summary<-cbind(mean(TCF7L1[GROUP_X_PART],na.rm=TRUE),-sd(TCF7L1[GROUP_X_PART],na.rm=TRUE),mean(TCF7L1[GROUP_Y_PART],na.rm=TRUE),-sd(TCF7L1[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(TCF7L1[GROUP_X_PART],na.rm=TRUE),"",median(TCF7L1[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(TCF7L1,na.rm=TRUE),-sd(TCF7L1,na.rm=TRUE) )
temp_total_2<-cbind( median(TCF7L1,na.rm=TRUE),"" )

summary_group_1 <-cbind("TCF7L1","mean (SD)",temp_summary,t.test(TCF7L1[GROUP_X_PART],TCF7L1[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("TCF7L1","median",temp_summary_2,wilcox.test(TCF7L1[GROUP_X_PART],TCF7L1[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)


### "MYLK4"

MYLK4 <- log(tumor_Data[,"MYLK4"]+1)

temp_summary<-cbind(mean(MYLK4[GROUP_X_PART],na.rm=TRUE),-sd(MYLK4[GROUP_X_PART],na.rm=TRUE),mean(MYLK4[GROUP_Y_PART],na.rm=TRUE),-sd(MYLK4[GROUP_Y_PART],na.rm=TRUE))
temp_summary_2<-cbind(median(MYLK4[GROUP_X_PART],na.rm=TRUE),"",median(MYLK4[GROUP_Y_PART],na.rm=TRUE),"")

temp_total<-cbind( mean(MYLK4,na.rm=TRUE),-sd(MYLK4,na.rm=TRUE) )
temp_total_2<-cbind( median(MYLK4,na.rm=TRUE),"" )

summary_group_1 <-cbind("MYLK4","mean (SD)",temp_summary,t.test(MYLK4[GROUP_X_PART],MYLK4[GROUP_Y_PART])$p.value,temp_total)
summary_group_2 <-cbind("MYLK4","median",temp_summary_2,wilcox.test(MYLK4[GROUP_X_PART],MYLK4[GROUP_Y_PART])$p.value,temp_total_2)
summary_table<-rbind(summary_table,summary_group_1,summary_group_2)

file_name=paste("C:/Users/korea/Desktop/DEG_100.csv", sep="")
write.table(summary_table,row.names=FALSE,col.names=FALSE, file=file_name,sep=",")

## Check the wilcoxon test
# rm(list=ls())
# 
# #load data for high and low group in COAD data set
# high_colon <- read.csv("C:/Users/korea/Desktop/COAD/clinical/high_grade_colon_mRNA_standardize.csv", row.names = 1)
# low_colon <- read.csv("C:/Users/korea/Desktop/COAD/clinical/low_grade_colon_mRNA_standardize.csv", row.names = 1)
# 
# 
# wilcoxon_test.p <- vector()
# wilcoxon_test.test <- vector()
# 
# ## wilcoxon test
# for (i in 1:ncol(high_colon)){
#   wilcoxon_test.p[i] <- wilcox.test(high_colon[,i], low_colon[,i])$p.value
# }
# 
# 
# ##DEg result
# DEG <- colnames(high_colon)[wilcoxon_test.p<0.005]
# DEG_result <- wilcoxon_test.p[wilcoxon_test.p<0.005]
# total_DEG_result_0.005 <- cbind(DEG, DEG_result)
# write.csv(total_DEG_result_0.005, "C:/Users/korea/Desktop/COAD/step1.DEG/DEG_result_0.05.csv", row.names=F)





