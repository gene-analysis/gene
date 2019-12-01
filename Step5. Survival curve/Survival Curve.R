rm(list=ls())
library("survival")           
library("glmnet") 

file_name_tumor=paste("file path", sep="")
tumor_Data<-read.table(file_name_tumor, header=TRUE, sep=",")
colnames(tumor_Data)

Day.to.follow.up <- as.numeric(as.character(tumor_Data[,"days_to_last_followup"]))
Day.to.Death <- as.numeric(as.character(tumor_Data[,"days_to_death"]))

#race <- tumor_Data[,"race" ]
#gender <- tumor_Data[,"gender"]
#pathologic_stage <- tumor_Data[,"pathologic_stage"]
APBB1 <- log2(tumor_Data[,"APBB1"]+1)
LYSMD3 <- log2(tumor_Data[,"LYSMD3"]+1)
AHSA2 <- log2(tumor_Data[,"AHSA2"]+1)
SLC7A14 <- log2(tumor_Data[,"SLC7A14"]+1)
SLC22A17 <- log2(tumor_Data[,"SLC22A17"]+1)

TMTC3<- log2(tumor_Data[,"TMTC3"]+1)
TMED7 <- log2(tumor_Data[,"TMED7"]+1)
C12orf48 <- log2(tumor_Data[,"C12orf48"]+1)
C14orf129 <- log2(tumor_Data[,"C14orf129"]+1)
C18orf32 <- log2(tumor_Data[,"C18orf32"]+1)
NUMA1 <- log2(tumor_Data[,"NUMA1"]+1)
GBP4 <- log2(tumor_Data[,"GBP4"]+1)
COPS4 <- log2(tumor_Data[,"COPS4"]+1)
GPRASP1 <- log2(tumor_Data[,"GPRASP1"]+1)
JAK2 <- log2(tumor_Data[,"JAK2"]+1)
SETD1A <- log2(tumor_Data[,"SETD1A"]+1)
RAB27B <- log2(tumor_Data[,"RAB27B"]+1)
TNRC6A <- log2(tumor_Data[,"TNRC6A"]+1)
ZNF767 <- log2(tumor_Data[,"ZNF767"]+1)



temp_TF<-!is.na(Day.to.Death)
Observ<-Day.to.Death
Observ[]<-"Alive"                                 
Observ[temp_TF]<-"Death"
Observ<-as.factor(Observ)

Day.to.Death_0 <- Day.to.Death
Day.to.Death_0[!temp_TF] <- Day.to.follow.up[!temp_TF]


### APBB1
med1 <- median(APBB1, na.rm = TRUE)
APBB11<-APBB1

APBB11[APBB1<med1] <- "Low"
APBB11[APBB1>=med1] <- "High"

#APBB11[APBB1<4.7165] <- "Low"
#APBB11[APBB1>=4.7165] <- "High"



summary(APBB11)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ APBB11 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ APBB11 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("APBB1  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(APBB11)
tmp_2<-table(APBB11[Observ=="Death"])

level_list <- levels(APBB11)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==APBB11],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### LYSMD3

med2 <- median(LYSMD3, na.rm = TRUE)

LYSMD31 <- LYSMD3
LYSMD31[LYSMD3<med2] <- "Low"
LYSMD31[LYSMD3>=med2] <- "High"
#LYSMD31[LYSMD3<6] <- "Low"
#LYSMD31[LYSMD3>=6] <- "High"

summary(LYSMD31)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ LYSMD31 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ LYSMD31 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("LYSMD3  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(LYSMD31)
tmp_2<-table(LYSMD31[Observ=="Death"])

level_list <- levels(LYSMD31)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==LYSMD31],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### AHSA2

med3 <- median(AHSA2, na.rm = TRUE)
AHSA21 <-AHSA2
AHSA21[AHSA2<med3] <- "Low"
AHSA21[AHSA2>=med3] <- "High"

#AHSA21[AHSA2<5.905] <- "Low"
#AHSA21[AHSA2>=5.905] <- "High"


summary(AHSA21)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ AHSA21 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ AHSA21 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("AHSA2  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(AHSA21)
tmp_2<-table(AHSA21[Observ=="Death"])

level_list <- levels(AHSA21)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==AHSA21],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### SLC7A14

med4 <- median(SLC7A14, na.rm = TRUE)
SLC7A141 <- SLC7A14  
SLC7A141[SLC7A14<med4] <- "Low"
SLC7A141[SLC7A14>=med4] <- "High"

#SLC7A141[SLC7A14<0.330] <- "Low"
#SLC7A141[SLC7A14>=0.330] <- "High"

summary(SLC7A141)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ SLC7A141 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ SLC7A141 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("SLC7A14  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(SLC7A141)
tmp_2<-table(SLC7A141[Observ=="Death"])

level_list <- levels(SLC7A141)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==SLC7A141],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### SLC22A17

med5 <- median(SLC22A17, na.rm = TRUE)
SLC22A171 <- SLC22A17  

SLC22A171[SLC22A17<med5] <- "Low"
SLC22A171[SLC22A17>=med5] <- "High"

#SLC22A171[SLC22A17<4.314] <- "Low"
#SLC22A171[SLC22A17>=4.314] <- "High"

summary(SLC22A171)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ SLC22A171 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ SLC22A171 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("SLC22A17  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(SLC22A171)
tmp_2<-table(SLC22A171[Observ=="Death"])

level_list <- levels(SLC22A171)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==SLC22A171],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()
##########positive vs DEG
### TMTC3
med1 <- median(TMTC3, na.rm = TRUE)
TMTC31<-TMTC3

TMTC31[TMTC3<med1] <- "Low"
TMTC31[TMTC3>=med1] <- "High"

#TMTC31[TMTC3<4.7165] <- "Low"
#TMTC31[TMTC3>=4.7165] <- "High"



summary(TMTC31)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ TMTC31 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ TMTC31 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("TMTC3  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(TMTC31)
tmp_2<-table(TMTC31[Observ=="Death"])

level_list <- levels(TMTC31)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==TMTC31],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### TMED7

med2 <- median(TMED7, na.rm = TRUE)

TMED71 <- TMED7
TMED71[TMED7<med2] <- "Low"
TMED71[TMED7>=med2] <- "High"
#TMED71[TMED7<6] <- "Low"
#TMED71[TMED7>=6] <- "High"

summary(TMED71)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ TMED71 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ TMED71 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("TMED7  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(TMED71)
tmp_2<-table(TMED71[Observ=="Death"])

level_list <- levels(TMED71)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==TMED71],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()
##############negative vs DEG
### C12orf48

med3 <- median(C12orf48, na.rm = TRUE)
C12orf481 <-C12orf48
C12orf481[C12orf48<med3] <- "Low"
C12orf481[C12orf48>=med3] <- "High"

#C12orf481[C12orf48<5.905] <- "Low"
#C12orf481[C12orf48>=5.905] <- "High"


summary(C12orf481)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ C12orf481 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ C12orf481 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("C12orf48  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(C12orf481)
tmp_2<-table(C12orf481[Observ=="Death"])

level_list <- levels(C12orf481)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==C12orf481],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### C14orf129

med4 <- median(C14orf129, na.rm = TRUE)
C14orf1291 <- C14orf129  
C14orf1291[C14orf129<med4] <- "Low"
C14orf1291[C14orf129>=med4] <- "High"

#C14orf1291[C14orf129<0.330] <- "Low"
#C14orf1291[C14orf129>=0.330] <- "High"

summary(C14orf1291)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ C14orf1291 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ C14orf1291 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("C14orf129  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(C14orf1291)
tmp_2<-table(C14orf1291[Observ=="Death"])

level_list <- levels(C14orf1291)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==C14orf1291],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### C18orf32

med5 <- median(C18orf32, na.rm = TRUE)
C18orf321 <- C18orf32  

C18orf321[C18orf32<med5] <- "Low"
C18orf321[C18orf32>=med5] <- "High"

#C18orf321[C18orf32<4.314] <- "Low"
#C18orf321[C18orf32>=4.314] <- "High"

summary(C18orf321)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ C18orf321 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ C18orf321 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("C18orf32  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(C18orf321)
tmp_2<-table(C18orf321[Observ=="Death"])

level_list <- levels(C18orf321)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==C18orf321],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### NUMA1
med6 <- median(NUMA1, na.rm = TRUE)
NUMA11<-NUMA1

NUMA11[NUMA1<med6] <- "Low"
NUMA11[NUMA1>=med6] <- "High"

#NUMA11[NUMA1<4.7165] <- "Low"
#NUMA11[NUMA1>=4.7165] <- "High"



summary(NUMA11)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ NUMA11 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ NUMA11 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("NUMA1  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(NUMA11)
tmp_2<-table(NUMA11[Observ=="Death"])

level_list <- levels(NUMA11)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==NUMA11],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### GBP4

med7 <- median(GBP4, na.rm = TRUE)

GBP41 <- GBP4
GBP41[GBP4<med7] <- "Low"
GBP41[GBP4>=med7] <- "High"
#GBP41[GBP4<6] <- "Low"
#GBP41[GBP4>=6] <- "High"

summary(GBP41)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ GBP41 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ GBP41 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("GBP4  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(GBP41)
tmp_2<-table(GBP41[Observ=="Death"])

level_list <- levels(GBP41)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==GBP41],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### COPS4

med8 <- median(COPS4, na.rm = TRUE)
COPS41 <-COPS4
COPS41[COPS4<med8] <- "Low"
COPS41[COPS4>=med8] <- "High"

#COPS41[COPS4<5.905] <- "Low"
#COPS41[COPS4>=5.905] <- "High"


summary(COPS41)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ COPS41 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ COPS41 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("COPS4  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(COPS41)
tmp_2<-table(COPS41[Observ=="Death"])

level_list <- levels(COPS41)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==COPS41],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### GPRASP1

med9 <- median(GPRASP1, na.rm = TRUE)
GPRASP11 <- GPRASP1  
GPRASP11[GPRASP1<med9] <- "Low"
GPRASP11[GPRASP1>=med9] <- "High"

#GPRASP11[GPRASP1<0.330] <- "Low"
#GPRASP11[GPRASP1>=0.330] <- "High"

summary(GPRASP11)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ GPRASP11 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ GPRASP11 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("GPRASP1  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(GPRASP11)
tmp_2<-table(GPRASP11[Observ=="Death"])

level_list <- levels(GPRASP11)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==GPRASP11],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### JAK2

med10 <- median(JAK2, na.rm = TRUE)
JAK21 <- JAK2  

JAK21[JAK2<med10] <- "Low"
JAK21[JAK2>=med10] <- "High"

#JAK21[JAK2<4.314] <- "Low"
#JAK21[JAK2>=4.314] <- "High"

summary(JAK21)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ JAK21 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ JAK21 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("JAK2  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(JAK21)
tmp_2<-table(JAK21[Observ=="Death"])

level_list <- levels(JAK21)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==JAK21],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### SETD1A
med11 <- median(SETD1A, na.rm = TRUE)
SETD1A1<-SETD1A

SETD1A1[SETD1A<med11] <- "Low"
SETD1A1[SETD1A>=med11] <- "High"

#SETD1A1[SETD1A<4.7165] <- "Low"
#SETD1A1[SETD1A>=4.7165] <- "High"



summary(SETD1A1)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ SETD1A1 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ SETD1A1 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("SETD1A  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(SETD1A1)
tmp_2<-table(SETD1A1[Observ=="Death"])

level_list <- levels(SETD1A1)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==SETD1A1],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### RAB27B

med12 <- median(RAB27B, na.rm = TRUE)

RAB27B1 <- RAB27B
RAB27B1[RAB27B<med12] <- "Low"
RAB27B1[RAB27B>=med12] <- "High"
#RAB27B1[RAB27B<6] <- "Low"
#RAB27B1[RAB27B>=6] <- "High"

summary(RAB27B1)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ RAB27B1 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ RAB27B1 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("RAB27B  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(RAB27B1)
tmp_2<-table(RAB27B1[Observ=="Death"])

level_list <- levels(RAB27B1)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==RAB27B1],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### TNRC6A

med13 <- median(TNRC6A, na.rm = TRUE)
TNRC6A1 <-TNRC6A
TNRC6A1[TNRC6A<med13] <- "Low"
TNRC6A1[TNRC6A>=med13] <- "High"

#TNRC6A1[TNRC6A<5.905] <- "Low"
#TNRC6A1[TNRC6A>=5.905] <- "High"


summary(TNRC6A1)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ TNRC6A1 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ TNRC6A1 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("TNRC6A  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(TNRC6A1)
tmp_2<-table(TNRC6A1[Observ=="Death"])

level_list <- levels(TNRC6A1)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==TNRC6A1],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### ZNF767

med14 <- median(ZNF767, na.rm = TRUE)
ZNF7671 <- ZNF767  
ZNF7671[ZNF767<med14] <- "Low"
ZNF7671[ZNF767>=med14] <- "High"

#ZNF7671[ZNF767<0.330] <- "Low"
#ZNF7671[ZNF767>=0.330] <- "High"

summary(ZNF7671)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ ZNF7671 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ ZNF7671 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("ZNF767  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(ZNF7671)
tmp_2<-table(ZNF7671[Observ=="Death"])

level_list <- levels(ZNF7671)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==ZNF7671],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()




########common hub of LN(+), LN(-)
rm(list=ls())
library("survival")           
library("glmnet") 

file_name_tumor=paste("file path", sep="")
tumor_Data<-read.table(file_name_tumor, header=TRUE, sep=",")
colnames(tumor_Data)

Day.to.follow.up <- as.numeric(as.character(tumor_Data[,"days_to_last_followup"]))
Day.to.Death <- as.numeric(as.character(tumor_Data[,"days_to_death"]))

HEG1 <- log2(tumor_Data[,"HEG1"]+1)
SECISBP2L <- log2(tumor_Data[,"SECISBP2L"]+1)
TCF4 <- log2(tumor_Data[,"TCF4"]+1)
CLIP3 <- log2(tumor_Data[,"CLIP3"]+1)
MSRB3 <- log2(tumor_Data[,"MSRB3"]+1)
PCNP<- log2(tumor_Data[,"PCNP"]+1)
VIM <- log2(tumor_Data[,"VIM"]+1)
ATP8B2 <- log2(tumor_Data[,"ATP8B2"]+1)

XPO1 <- log2(tumor_Data[,"XPO1"]+1)
TNS1 <- log2(tumor_Data[,"TNS1"]+1)
PIKFYVE <- log2(tumor_Data[,"PIKFYVE"]+1)
VPS26A <- log2(tumor_Data[,"VPS26A"]+1)
CWC22 <- log2(tumor_Data[,"CWC22"]+1)
CCDC80 <- log2(tumor_Data[,"CCDC80"]+1)
MATR3 <- log2(tumor_Data[,"MATR3"]+1)
ZEB1 <- log2(tumor_Data[,"ZEB1"]+1)
C1S <- log2(tumor_Data[,"C1S"]+1)
AFF4 <- log2(tumor_Data[,"AFF4"]+1)
ZEB2 <- log2(tumor_Data[,"ZEB2"]+1)
LY6G6D <- log2(tumor_Data[,"LY6G6D"]+1)
SSB <- log2(tumor_Data[,"SSB"]+1)


temp_TF<-!is.na(Day.to.Death)
Observ<-Day.to.Death
Observ[]<-"Alive"                                 
Observ[temp_TF]<-"Death"
Observ<-as.factor(Observ)

Day.to.Death_0 <- Day.to.Death
Day.to.Death_0[!temp_TF] <- Day.to.follow.up[!temp_TF]


### HEG1
med1 <- median(HEG1, na.rm = TRUE)
HEG11<-HEG1

HEG11[HEG1<med1] <- "Low"
HEG11[HEG1>=med1] <- "High"

#HEG11[HEG1<4.7165] <- "Low"
#HEG11[HEG1>=4.7165] <- "High"



summary(HEG11)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ HEG11 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ HEG11 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("HEG1  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(HEG11)
tmp_2<-table(HEG11[Observ=="Death"])

level_list <- levels(HEG11)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==HEG11],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### SECISBP2L

med2 <- median(SECISBP2L, na.rm = TRUE)

SECISBP2L1 <- SECISBP2L
SECISBP2L1[SECISBP2L<med2] <- "Low"
SECISBP2L1[SECISBP2L>=med2] <- "High"
#SECISBP2L1[SECISBP2L<6] <- "Low"
#SECISBP2L1[SECISBP2L>=6] <- "High"

summary(SECISBP2L1)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ SECISBP2L1 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ SECISBP2L1 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("SECISBP2L  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(SECISBP2L1)
tmp_2<-table(SECISBP2L1[Observ=="Death"])

level_list <- levels(SECISBP2L1)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==SECISBP2L1],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### TCF4

med3 <- median(TCF4, na.rm = TRUE)
TCF41 <-TCF4
TCF41[TCF4<med3] <- "Low"
TCF41[TCF4>=med3] <- "High"

#TCF41[TCF4<5.905] <- "Low"
#TCF41[TCF4>=5.905] <- "High"


summary(TCF41)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ TCF41 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ TCF41 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("TCF4  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(TCF41)
tmp_2<-table(TCF41[Observ=="Death"])

level_list <- levels(TCF41)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==TCF41],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### CLIP3

med4 <- median(CLIP3, na.rm = TRUE)
CLIP31 <- CLIP3  
CLIP31[CLIP3<med4] <- "Low"
CLIP31[CLIP3>=med4] <- "High"

#CLIP31[CLIP3<0.330] <- "Low"
#CLIP31[CLIP3>=0.330] <- "High"

summary(CLIP31)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ CLIP31 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ CLIP31 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("CLIP3  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(CLIP31)
tmp_2<-table(CLIP31[Observ=="Death"])

level_list <- levels(CLIP31)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==CLIP31],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### MSRB3

med5 <- median(MSRB3, na.rm = TRUE)
MSRB31 <- MSRB3  

MSRB31[MSRB3<med5] <- "Low"
MSRB31[MSRB3>=med5] <- "High"

#MSRB31[MSRB3<4.314] <- "Low"
#MSRB31[MSRB3>=4.314] <- "High"

summary(MSRB31)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ MSRB31 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ MSRB31 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("MSRB3  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(MSRB31)
tmp_2<-table(MSRB31[Observ=="Death"])

level_list <- levels(MSRB31)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==MSRB31],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()
##########positive vs DEG
### PCNP
med1 <- median(PCNP, na.rm = TRUE)
PCNP1<-PCNP

PCNP1[PCNP<med1] <- "Low"
PCNP1[PCNP>=med1] <- "High"

#PCNP1[PCNP<4.7165] <- "Low"
#PCNP1[PCNP>=4.7165] <- "High"



summary(PCNP1)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ PCNP1 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ PCNP1 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("PCNP  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(PCNP1)
tmp_2<-table(PCNP1[Observ=="Death"])

level_list <- levels(PCNP1)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==PCNP1],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### VIM

med2 <- median(VIM, na.rm = TRUE)

VIM1 <- VIM
VIM1[VIM<med2] <- "Low"
VIM1[VIM>=med2] <- "High"
#VIM1[VIM<6] <- "Low"
#VIM1[VIM>=6] <- "High"

summary(VIM1)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ VIM1 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ VIM1 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("VIM  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(VIM1)
tmp_2<-table(VIM1[Observ=="Death"])

level_list <- levels(VIM1)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==VIM1],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### ATP8B2

med3 <- median(ATP8B2, na.rm = TRUE)
ATP8B21 <-ATP8B2
ATP8B21[ATP8B2<med3] <- "Low"
ATP8B21[ATP8B2>=med3] <- "High"

#ATP8B21[ATP8B2<5.905] <- "Low"
#ATP8B21[ATP8B2>=5.905] <- "High"


summary(ATP8B21)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ ATP8B21 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ ATP8B21 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("ATP8B2  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(ATP8B21)
tmp_2<-table(ATP8B21[Observ=="Death"])

level_list <- levels(ATP8B21)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==ATP8B21],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### XPO1

med4 <- median(XPO1, na.rm = TRUE)
XPO11 <- XPO1  
XPO11[XPO1<med4] <- "Low"
XPO11[XPO1>=med4] <- "High"

#XPO11[XPO1<0.330] <- "Low"
#XPO11[XPO1>=0.330] <- "High"

summary(XPO11)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ XPO11 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ XPO11 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("XPO1  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(XPO11)
tmp_2<-table(XPO11[Observ=="Death"])

level_list <- levels(XPO11)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==XPO11],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### TNS1

med5 <- median(TNS1, na.rm = TRUE)
TNS11 <- TNS1  

TNS11[TNS1<med5] <- "Low"
TNS11[TNS1>=med5] <- "High"

#TNS11[TNS1<4.314] <- "Low"
#TNS11[TNS1>=4.314] <- "High"

summary(TNS11)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ TNS11 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ TNS11 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("TNS1  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(TNS11)
tmp_2<-table(TNS11[Observ=="Death"])

level_list <- levels(TNS11)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==TNS11],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### PIKFYVE
med6 <- median(PIKFYVE, na.rm = TRUE)
PIKFYVE1<-PIKFYVE

PIKFYVE1[PIKFYVE<med6] <- "Low"
PIKFYVE1[PIKFYVE>=med6] <- "High"

#PIKFYVE1[PIKFYVE<4.7165] <- "Low"
#PIKFYVE1[PIKFYVE>=4.7165] <- "High"



summary(PIKFYVE1)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ PIKFYVE1 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ PIKFYVE1 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("PIKFYVE  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(PIKFYVE1)
tmp_2<-table(PIKFYVE1[Observ=="Death"])

level_list <- levels(PIKFYVE1)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==PIKFYVE1],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### VPS26A

med7 <- median(VPS26A, na.rm = TRUE)

VPS26A1 <- VPS26A
VPS26A1[VPS26A<med7] <- "Low"
VPS26A1[VPS26A>=med7] <- "High"
#VPS26A1[VPS26A<6] <- "Low"
#VPS26A1[VPS26A>=6] <- "High"

summary(VPS26A1)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ VPS26A1 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ VPS26A1 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("VPS26A  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(VPS26A1)
tmp_2<-table(VPS26A1[Observ=="Death"])

level_list <- levels(VPS26A1)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==VPS26A1],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### CWC22

med8 <- median(CWC22, na.rm = TRUE)
CWC221 <-CWC22
CWC221[CWC22<med8] <- "Low"
CWC221[CWC22>=med8] <- "High"

#CWC221[CWC22<5.905] <- "Low"
#CWC221[CWC22>=5.905] <- "High"


summary(CWC221)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ CWC221 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ CWC221 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("CWC22  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(CWC221)
tmp_2<-table(CWC221[Observ=="Death"])

level_list <- levels(CWC221)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==CWC221],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### CCDC80

med9 <- median(CCDC80, na.rm = TRUE)
CCDC801 <- CCDC80  
CCDC801[CCDC80<med9] <- "Low"
CCDC801[CCDC80>=med9] <- "High"

#CCDC801[CCDC80<0.330] <- "Low"
#CCDC801[CCDC80>=0.330] <- "High"

summary(CCDC801)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ CCDC801 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ CCDC801 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("CCDC80  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(CCDC801)
tmp_2<-table(CCDC801[Observ=="Death"])

level_list <- levels(CCDC801)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==CCDC801],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### MATR3

med10 <- median(MATR3, na.rm = TRUE)
MATR31 <- MATR3  

MATR31[MATR3<med10] <- "Low"
MATR31[MATR3>=med10] <- "High"

#MATR31[MATR3<4.314] <- "Low"
#MATR31[MATR3>=4.314] <- "High"

summary(MATR31)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ MATR31 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ MATR31 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("MATR3  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(MATR31)
tmp_2<-table(MATR31[Observ=="Death"])

level_list <- levels(MATR31)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==MATR31],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### ZEB1
med11 <- median(ZEB1, na.rm = TRUE)
ZEB11<-ZEB1

ZEB11[ZEB1<med11] <- "Low"
ZEB11[ZEB1>=med11] <- "High"

#ZEB11[ZEB1<4.7165] <- "Low"
#ZEB11[ZEB1>=4.7165] <- "High"



summary(ZEB11)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ ZEB11 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ ZEB11 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("ZEB1  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(ZEB11)
tmp_2<-table(ZEB11[Observ=="Death"])

level_list <- levels(ZEB11)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==ZEB11],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### C1S

med12 <- median(C1S, na.rm = TRUE)

C1S1 <- C1S
C1S1[C1S<med12] <- "Low"
C1S1[C1S>=med12] <- "High"
#C1S1[C1S<6] <- "Low"
#C1S1[C1S>=6] <- "High"

summary(C1S1)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ C1S1 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ C1S1 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("C1S  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(C1S1)
tmp_2<-table(C1S1[Observ=="Death"])

level_list <- levels(C1S1)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==C1S1],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### AFF4

med13 <- median(AFF4, na.rm = TRUE)
AFF41 <-AFF4
AFF41[AFF4<med13] <- "Low"
AFF41[AFF4>=med13] <- "High"

#AFF41[AFF4<5.905] <- "Low"
#AFF41[AFF4>=5.905] <- "High"


summary(AFF41)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ AFF41 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ AFF41 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("AFF4  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(AFF41)
tmp_2<-table(AFF41[Observ=="Death"])

level_list <- levels(AFF41)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==AFF41],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### ZEB2

med14 <- median(ZEB2, na.rm = TRUE)
ZEB21 <- ZEB2  
ZEB21[ZEB2<med14] <- "Low"
ZEB21[ZEB2>=med14] <- "High"

#ZEB21[ZEB2<0.330] <- "Low"
#ZEB21[ZEB2>=0.330] <- "High"

summary(ZEB21)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ ZEB21 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ ZEB21 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("ZEB2  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(ZEB21)
tmp_2<-table(ZEB21[Observ=="Death"])

level_list <- levels(ZEB21)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==ZEB21],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### LY6G6D
med1 <- median(LY6G6D, na.rm = TRUE)
LY6G6D1<-LY6G6D

LY6G6D1[LY6G6D<med1] <- "Low"
LY6G6D1[LY6G6D>=med1] <- "High"

#LY6G6D1[LY6G6D<4.7165] <- "Low"
#LY6G6D1[LY6G6D>=4.7165] <- "High"



summary(LY6G6D1)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ LY6G6D1 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ LY6G6D1 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("LY6G6D  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(LY6G6D1)
tmp_2<-table(LY6G6D1[Observ=="Death"])

level_list <- levels(LY6G6D1)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==LY6G6D1],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

### SSB

med2 <- median(SSB, na.rm = TRUE)

SSB1 <- SSB
SSB1[SSB<med2] <- "Low"
SSB1[SSB>=med2] <- "High"
#SSB1[SSB<6] <- "Low"
#SSB1[SSB>=6] <- "High"

summary(SSB1)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ SSB1 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ SSB1 )$chisq,2)

summary(fit.byrisk)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("SSB  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(SSB1)
tmp_2<-table(SSB1[Observ=="Death"])

level_list <- levels(SSB1)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==SSB1],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()

###################################################

race<-as.character(race)
race[race=="[Not Available]"] <- NA
race[race=="[Not Evaluated]"] <- NA
race<-as.factor(as.character(race))

gender<-as.character(gender)
gender[gender=="[Not Available]"] <- NA
gender[gender=="[Not Evaluated]"] <- NA
gender<-as.factor(as.character(gender))

three_tumor_stage <- as.character(pathologic_stage)

three_tumor_stage[three_tumor_stage=="stage i"] <- "Stage_1"
three_tumor_stage[three_tumor_stage=="stage ia"] <- "Stage_1"
three_tumor_stage[three_tumor_stage=="stage ib"] <- "Stage_1"

three_tumor_stage[three_tumor_stage=="stage ii"] <- "Stage_2"
three_tumor_stage[three_tumor_stage=="stage iia"] <- "Stage_2"
three_tumor_stage[three_tumor_stage=="stage iib"] <- "Stage_2"
three_tumor_stage[three_tumor_stage=="stage iic"] <- NA

three_tumor_stage[three_tumor_stage=="stage iii"] <- "Stage_3_more"
three_tumor_stage[three_tumor_stage=="stage iiia"] <- "Stage_3_more"
three_tumor_stage[three_tumor_stage=="stage iiib"] <- "Stage_3_more"
three_tumor_stage[three_tumor_stage=="stage iiic"] <- "Stage_3_more"

three_tumor_stage[three_tumor_stage=="stage iv"] <- "Stage_3_more"
three_tumor_stage[three_tumor_stage=="stage iva"] <- NA
three_tumor_stage[three_tumor_stage=="stage ivb"] <- NA

three_tumor_stage[three_tumor_stage=="[Not Available]"] <- NA
three_tumor_stage[three_tumor_stage=="[Not Evaluated]"] <- NA
three_tumor_stage<-as.factor(as.character(three_tumor_stage))



file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ race )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ race )$chisq,2)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste(" p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Death)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(race)
tmp_2<-table(race[Observ=="Death"])

level_list <- levels(race)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==race],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()



file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ gender )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ gender )$chisq,2)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste(" p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Death)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(gender)
tmp_2<-table(gender[Observ=="Death"])

level_list <- levels(gender)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==gender],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()


file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ three_tumor_stage )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ three_tumor_stage )$chisq,2)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste(" p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Death)", cex.axis=1.05 , cex.lab=1.1 )

tmp<-table(three_tumor_stage)
tmp_2<-table(three_tumor_stage[Observ=="Death"])

level_list <- levels(three_tumor_stage)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==three_tumor_stage],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()



# Update
rm(list=ls())
library("survival")
library("glmnet")


file_name_tumor=paste("file path", sep="")
tumor_Data <- read.table(file_name_tumor, header=TRUE, sep=",")


Day.to.follow.up <- as.numeric(as.character(tumor_Data[,"days_to_last_followup"]))
Day.to.Death <- as.numeric(as.character(tumor_Data[,"days_to_death"]))

##DEG&positive&negative
SLC22A17 <- log2(tumor_Data[,"SLC22A17"]+1)
APBB1 <- log2(tumor_Data[,"APBB1"]+1)
SLC7A14 <- log2(tumor_Data[,"SLC7A14"]+1)
JAM3 <- log2(tumor_Data[,"JAM3"]+1)
PRELP <- log2(tumor_Data[,"PRELP"]+1)
LYSMD3 <- log2(tumor_Data[,"LYSMD3"]+1)
RBPMS2 <- log2(tumor_Data[,"RBPMS2"]+1)
LMOD1 <- log2(tumor_Data[,"LMOD1"]+1)
GEFT <- log2(tumor_Data[,"GEFT"]+1)
SALL2 <- log2(tumor_Data[,"SALL2"]+1)
TNS1 <- log2(tumor_Data[,"TNS1"]+1)
EFEMP2 <- log2(tumor_Data[,"EFEMP2"]+1)
SYDE1 <- log2(tumor_Data[,"SYDE1"]+1)
CLIP3 <- log2(tumor_Data[,"CLIP3"]+1)
MRVI1 <- log2(tumor_Data[,"MRVI1"]+1)
PKN2 <- log2(tumor_Data[,"PKN2"]+1)
AHSA2 <- log2(tumor_Data[,"AHSA2"]+1)
AKAP11 <- log2(tumor_Data[,"AKAP11"]+1)
TIMP2 <- log2(tumor_Data[,"TIMP2"]+1)
CDK1 <- log2(tumor_Data[,"CDK1"]+1)
ABCE1 <- log2(tumor_Data[,"ABCE1"]+1)
PKD1 <- log2(tumor_Data[,"PKD1"]+1)
SGMS2 <- log2(tumor_Data[,"SGMS2"]+1)
MGP <- log2(tumor_Data[,"MGP"]+1)
HSPB8 <- log2(tumor_Data[,"HSPB8"]+1)
BOC <- log2(tumor_Data[,"BOC"]+1)

#DEG&positive
TMTC3 <- log2(tumor_Data[,"TMTC3"]+1)
FXYD6 <- log2(tumor_Data[,"FXYD6"]+1)
PDZD4 <- log2(tumor_Data[,"PDZD4"]+1)
SLC35A3 <- log2(tumor_Data[,"SLC35A3"]+1)
TMED7 <- log2(tumor_Data[,"TMED7"]+1)
SCAF1 <- log2(tumor_Data[,"SCAF1"]+1)
TUB <- log2(tumor_Data[,"TUB"]+1)
MYH11 <- log2(tumor_Data[,"MYH11"]+1)
C14orf132 <- log2(tumor_Data[,"C14orf132"]+1)
SPARCL1 <- log2(tumor_Data[,"SPARCL1"]+1)
TRO <- log2(tumor_Data[,"TRO"]+1)

#DEG&negative
C12orf48 <- log2(tumor_Data[,"C12orf48"]+1)
C14orf129 <- log2(tumor_Data[,"C14orf129"]+1)
C18orf32 <- log2(tumor_Data[,"C18orf32"]+1)
PDLIM7 <- log2(tumor_Data[,"PDLIM7"]+1)
COPS4 <- log2(tumor_Data[,"COPS4"]+1)
ADAMTSL3 <- log2(tumor_Data[,"ADAMTSL3"]+1)
FHL1 <- log2(tumor_Data[,"FHL1"]+1)
GPRASP1 <- log2(tumor_Data[,"GPRASP1"]+1)
HMCN1 <- log2(tumor_Data[,"HMCN1"]+1)
GBP4 <- log2(tumor_Data[,"GBP4"]+1)
JAK2 <- log2(tumor_Data[,"JAK2"]+1)
MXRA8 <- log2(tumor_Data[,"MXRA8"]+1)
SETD1A <- log2(tumor_Data[,"SETD1A"]+1)
RAB27B <- log2(tumor_Data[,"RAB27B"]+1)
TNRC6A <- log2(tumor_Data[,"TNRC6A"]+1)
NUMA1 <- log2(tumor_Data[,"NUMA1"]+1)
MRPL50 <- log2(tumor_Data[,"MRPL50"]+1)
ZNF24 <- log2(tumor_Data[,"ZNF24"]+1)
LONRF2 <- log2(tumor_Data[,"LONRF2"]+1)
ZNF767 <- log2(tumor_Data[,"ZNF767"]+1)
ARFIP1 <- log2(tumor_Data[,"ARFIP1"]+1)
USP33 <- log2(tumor_Data[,"USP33"]+1)
C5orf44 <- log2(tumor_Data[,"C5orf44"]+1)
ZNF720 <- log2(tumor_Data[,"ZNF720"]+1)
UBA3 <- log2(tumor_Data[,"UBA3"]+1)
LDB2 <- log2(tumor_Data[,"LDB2"]+1)
CDK10 <- log2(tumor_Data[,"CDK10"]+1)



temp_TF<-!is.na(Day.to.Death)
Observ<-Day.to.Death
Observ[]<-"Alive"                                 
Observ[temp_TF]<-"Death"
Observ<-as.factor(Observ)

Day.to.Death_0 <- Day.to.Death
Day.to.Death_0[!temp_TF] <- Day.to.follow.up[!temp_TF]


##SLC22A17
med1 <- median(SLC22A17, na.rm = TRUE)
SLC22A171<-SLC22A17

SLC22A171[SLC22A17<med1] <- "Low"
SLC22A171[SLC22A17>=med1] <- "High"

summary(SLC22A171)

file_name_plots=paste("file path", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ SLC22A171 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ SLC22A171 )$chisq,2)

summary(fit.byrisk)

write.csv(fit.byrisk, "file path", row.names=F)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("SLC22A17  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)",cex.main=2.5 ,cex.axis=1.2 , cex.lab=1.2  )

tmp<-table(SLC22A171)
tmp_2<-table(SLC22A171[Observ=="Death"])

level_list <- levels(SLC22A171)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==SLC22A171],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()



