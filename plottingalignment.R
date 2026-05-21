library(ggplot2)
library(data.table)
setwd("/Users/a02499139/Desktop/Gompert_Lab_Research/TimemaSVmethods/progressivecactus_pairwise")

# GSH2_GSR2 ---------------------------------------------------------------

## read synteny dat
dat<-fread("cactusStripe_TcrGSH2_TcrGSR2.psl",header=FALSE)
dfdat<-as.data.frame(dat)

head(dat$V10)
head(dat$V14)

dfdat$V10<- as.numeric(gsub(".*Scaffold_([0-9]+).*", "\\1", dfdat$V10))
dfdat$V14<- as.numeric(gsub(".*Scaffold_([0-9]+).*", "\\1", dfdat$V14))


## keep alignments from major scaffolds
xx_query <- table(dfdat[,10])
xx_target <- table(dfdat[,14])

g1 <- names(xx_query)[xx_query > 100]     #lmelCP
g2 <- names(xx_target)[xx_target > 100]   #lmelOLD  

keep <- (dfdat[,10] %in% g1) & (dfdat[,14] %in% g2)
subDfdat <- dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
## normalize with respect to knulli
ntab<-tab
for(i in 1:12){
  ntab[i,]<-ntab[i,]/sum(ntab[i,],na.rm=TRUE)
}

g1_chr <- as.numeric(gsub("Scaffold_([0-9]+)", "\\1", rownames(ntab)))
g2_chr <- as.numeric(gsub("Scaffold_([0-9]+)", "\\1", colnames(ntab)))



pdf("Syn_GSH2_GSR2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="GSR Hap2",ylab="GSH Hap2",cex.lab=1.4)
axis(1,at=seq(0,12,length.out=12)/12,g1_chr,las=2)
axis(2,at=seq(0,13,length.out=13)/13,g2_chr,las=2)
box()
dev.off()

chtab<-matrix(c(1,1,
                1,2,
                2,5,
                3,9,
                4,11,
                5,6,
                6,10,
                7,8,
                8,13,
                9,4,
                10,7,
                11,3,
                12,12),nrow=13,ncol=2,byrow=TRUE)



## colinearity plots for all homologous chromsomes,
## filling in for forward/reverse reads

pdf("AlnPlots_GSH2_GSR2_fwdrvs.pdf", width=10, height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))

for(i in 1:nrow(chtab)) {
  g1_pattern <- paste0(chtab[i,1])
  g2_pattern <- paste0(chtab[i,2])

  
  g1_chr <- which(subDfdat[,10] == g1_pattern)
  g2_chr <- which(subDfdat[,14] == g2_pattern)
  
  cc <- intersect(g1_chr, g2_chr)
  subd <- subDfdat[cc, ]
  
  # Flip y-axis if mostly reverse strand
  major_strand <- names(which.max(table(subd$V9)))
  if (major_strand %in% c("+-", "-+")) {
    subd$V12 <- max(subd$V12, na.rm=TRUE) - subd$V12
    subd$V13 <- max(subd$V13, na.rm=TRUE) - subd$V13
  }
  
  yub <- max(subd[,13], na.rm=TRUE)
  xub <- max(subd[,17], na.rm=TRUE)
  
  plot(as.numeric(subd[1,16:17]), as.numeric(subd[1,12:13]), type='n',
       xlim=c(0, xub), ylim=c(0, yub),
       cex.lab=1.4, ylab="GSR Hap2", xlab="GSH Hap2")
  
  title(main=paste("Pair", i, "GSR Hap2", chtab[i,1], "vs GSH Hap2", chtab[i,2]), cex.main=1.3)
  
  N<-dim(subd)[1]
  for(j in 1:N){
    x <- as.numeric(subd[j,16:17])   
    y <- as.numeric(subd[j,12:13])   
    
    # Fix reversed block endpoint ordering
    if (subd[j,9] %in% c("+-","-+")) {
      y <- rev(y)
    }
    
    if (x[1] > x[2]) {  
      x <- rev(x)                    
      y <- rev(y)                    
    }
    
    if(subd[j,9] == "++") {
      lines(x, y)                     
    }
    else{
      lines(x, y, col="cadetblue")   
    }
  }
  abline(a=0, b=1, col="gray", lty=2)
}


dev.off()



# GSH2_GSR1 ---------------------------------------------------------------

## read synteny dat
dat<-fread("cactusStripe_TcrGSH2_TcrGSR1.psl",header=FALSE)
dfdat<-as.data.frame(dat)

head(dat$V10)
head(dat$V14)

dfdat$V10<- as.numeric(gsub(".*Scaffold_([0-9]+).*", "\\1", dfdat$V10))
dfdat$V14<- as.numeric(gsub(".*Scaffold_([0-9]+).*", "\\1", dfdat$V14))


## keep alignments from major scaffolds
xx_query <- table(dfdat[,10])
xx_target <- table(dfdat[,14])

g1 <- names(xx_query)[xx_query > 100]     #lmelCP
g2 <- names(xx_target)[xx_target > 100]   #lmelOLD  

keep <- (dfdat[,10] %in% g1) & (dfdat[,14] %in% g2)
subDfdat <- dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
## normalize with respect to knulli
ntab<-tab
for(i in 1:13){
  ntab[i,]<-ntab[i,]/sum(ntab[i,],na.rm=TRUE)
}

g1_chr <- as.numeric(gsub("Scaffold_([0-9]+)", "\\1", rownames(ntab)))
g2_chr <- as.numeric(gsub("Scaffold_([0-9]+)", "\\1", colnames(ntab)))



pdf("Syn_GSH2_GSR1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="GSR Hap1",ylab="GSH Hap2",cex.lab=1.4)
axis(1,at=seq(0,13,length.out=13)/13,g1_chr,las=2)
axis(2,at=seq(0,13,length.out=13)/13,g2_chr,las=2)
box()
dev.off()

chtab<-matrix(c(1,1,
                1,2,
                2,5,
                3,9,
                4,11,
                5,6,
                6,10,
                7,8,
                8,13,
                9,4,
                10,7,
                11,3,
                12,12),nrow=13,ncol=2,byrow=TRUE)



## colinearity plots for all homologous chromsomes,
## filling in for forward/reverse reads

pdf("AlnPlots_GSH2_GSR2_fwdrvs.pdf", width=10, height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))

for(i in 1:nrow(chtab)) {
  g1_pattern <- paste0(chtab[i,1])
  g2_pattern <- paste0(chtab[i,2])
  
  
  g1_chr <- which(subDfdat[,10] == g1_pattern)
  g2_chr <- which(subDfdat[,14] == g2_pattern)
  
  cc <- intersect(g1_chr, g2_chr)
  subd <- subDfdat[cc, ]
  
  # Flip y-axis if mostly reverse strand
  major_strand <- names(which.max(table(subd$V9)))
  if (major_strand %in% c("+-", "-+")) {
    subd$V12 <- max(subd$V12, na.rm=TRUE) - subd$V12
    subd$V13 <- max(subd$V13, na.rm=TRUE) - subd$V13
  }
  
  yub <- max(subd[,13], na.rm=TRUE)
  xub <- max(subd[,17], na.rm=TRUE)
  
  plot(as.numeric(subd[1,16:17]), as.numeric(subd[1,12:13]), type='n',
       xlim=c(0, xub), ylim=c(0, yub),
       cex.lab=1.4, ylab="GSR Hap2", xlab="GSH Hap2")
  
  title(main=paste("Pair", i, "GSR Hap2", chtab[i,1], "vs GSH Hap2", chtab[i,2]), cex.main=1.3)
  
  N<-dim(subd)[1]
  for(j in 1:N){
    x <- as.numeric(subd[j,16:17])   
    y <- as.numeric(subd[j,12:13])   
    
    # Fix reversed block endpoint ordering
    if (subd[j,9] %in% c("+-","-+")) {
      y <- rev(y)
    }
    
    if (x[1] > x[2]) {  
      x <- rev(x)                    
      y <- rev(y)                    
    }
    
    if(subd[j,9] == "++") {
      lines(x, y)                     
    }
    else{
      lines(x, y, col="cadetblue")   
    }
  }
  abline(a=0, b=1, col="gray", lty=2)
}


dev.off()

# GSH2_GSH1 ---------------------------------------------------------------

## read synteny dat
dat<-fread("cactusStripe_TcrGSH2_TcrGSH1.psl",header=FALSE)
dfdat<-as.data.frame(dat)

head(dat$V10)
head(dat$V14)

dfdat$V10<- as.numeric(gsub(".*Scaffold_([0-9]+).*", "\\1", dfdat$V10))
dfdat$V14<- as.numeric(gsub(".*Scaffold_([0-9]+).*", "\\1", dfdat$V14))


## keep alignments from major scaffolds
xx_query <- table(dfdat[,10])
xx_target <- table(dfdat[,14])

g1 <- names(xx_query)[xx_query > 100]     #lmelCP
g2 <- names(xx_target)[xx_target > 100]   #lmelOLD  

keep <- (dfdat[,10] %in% g1) & (dfdat[,14] %in% g2)
subDfdat <- dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
## normalize with respect to knulli
ntab<-tab
for(i in 1:13){
  ntab[i,]<-ntab[i,]/sum(ntab[i,],na.rm=TRUE)
}

g1_chr <- as.numeric(gsub("Scaffold_([0-9]+)", "\\1", rownames(ntab)))
g2_chr <- as.numeric(gsub("Scaffold_([0-9]+)", "\\1", colnames(ntab)))



pdf("Syn_GSH2_GSH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="GSH Hap1",ylab="GSH Hap2",cex.lab=1.4)
axis(1,at=seq(0,13,length.out=13)/13,g1_chr,las=2)
axis(2,at=seq(0,13,length.out=13)/13,g2_chr,las=2)
box()
dev.off()

chtab<-matrix(c(1,1,
                1,2,
                2,5,
                3,9,
                4,11,
                5,6,
                6,10,
                7,8,
                8,13,
                9,4,
                10,7,
                11,3,
                12,12),nrow=13,ncol=2,byrow=TRUE)



## colinearity plots for all homologous chromsomes,
## filling in for forward/reverse reads

pdf("AlnPlots_GSH2_GSR2_fwdrvs.pdf", width=10, height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))

for(i in 1:nrow(chtab)) {
  g1_pattern <- paste0(chtab[i,1])
  g2_pattern <- paste0(chtab[i,2])
  
  
  g1_chr <- which(subDfdat[,10] == g1_pattern)
  g2_chr <- which(subDfdat[,14] == g2_pattern)
  
  cc <- intersect(g1_chr, g2_chr)
  subd <- subDfdat[cc, ]
  
  # Flip y-axis if mostly reverse strand
  major_strand <- names(which.max(table(subd$V9)))
  if (major_strand %in% c("+-", "-+")) {
    subd$V12 <- max(subd$V12, na.rm=TRUE) - subd$V12
    subd$V13 <- max(subd$V13, na.rm=TRUE) - subd$V13
  }
  
  yub <- max(subd[,13], na.rm=TRUE)
  xub <- max(subd[,17], na.rm=TRUE)
  
  plot(as.numeric(subd[1,16:17]), as.numeric(subd[1,12:13]), type='n',
       xlim=c(0, xub), ylim=c(0, yub),
       cex.lab=1.4, ylab="GSR Hap2", xlab="GSH Hap2")
  
  title(main=paste("Pair", i, "GSR Hap2", chtab[i,1], "vs GSH Hap2", chtab[i,2]), cex.main=1.3)
  
  N<-dim(subd)[1]
  for(j in 1:N){
    x <- as.numeric(subd[j,16:17])   
    y <- as.numeric(subd[j,12:13])   
    
    # Fix reversed block endpoint ordering
    if (subd[j,9] %in% c("+-","-+")) {
      y <- rev(y)
    }
    
    if (x[1] > x[2]) {  
      x <- rev(x)                    
      y <- rev(y)                    
    }
    
    if(subd[j,9] == "++") {
      lines(x, y)                     
    }
    else{
      lines(x, y, col="cadetblue")   
    }
  }
  abline(a=0, b=1, col="gray", lty=2)
}


dev.off()

# GSH2_GUSR1 ---------------------------------------------------------------

## read synteny dat
dat<-fread("cactusStripe_TcrGSH2_TcrGUSR1.psl",header=FALSE)
dfdat<-as.data.frame(dat)

head(dat$V10)
head(dat$V14)

dfdat$V10<- as.numeric(gsub(".*Scaffold_([0-9]+).*", "\\1", dfdat$V10))
dfdat$V14<- as.numeric(gsub(".*Scaffold_([0-9]+).*", "\\1", dfdat$V14))


## keep alignments from major scaffolds
xx_query <- table(dfdat[,10])
xx_target <- table(dfdat[,14])

g1 <- names(xx_query)[xx_query > 100]     #lmelCP
g2 <- names(xx_target)[xx_target > 100]   #lmelOLD  

keep <- (dfdat[,10] %in% g1) & (dfdat[,14] %in% g2)
subDfdat <- dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
## normalize with respect to knulli
ntab<-tab
for(i in 1:12){
  ntab[i,]<-ntab[i,]/sum(ntab[i,],na.rm=TRUE)
}

g1_chr <- as.numeric(gsub("Scaffold_([0-9]+)", "\\1", rownames(ntab)))
g2_chr <- as.numeric(gsub("Scaffold_([0-9]+)", "\\1", colnames(ntab)))



pdf("Syn_GSH2_GUSR1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="GUSR Hap1",ylab="GSH Hap2",cex.lab=1.4)
axis(1,at=seq(0,12,length.out=12)/12,g1_chr,las=2)
axis(2,at=seq(0,13,length.out=13)/13,g2_chr,las=2)
box()
dev.off()

chtab<-matrix(c(1,1,
                1,2,
                2,3,
                3,4,
                4,7,
                5,10,
                6,5,
                7,6,
                8,9,
                9,11,
                10,8,
                11,13,
                12,12),nrow=13,ncol=2,byrow=TRUE)



## colinearity plots for all homologous chromsomes,
## filling in for forward/reverse reads

pdf("AlnPlots_GSH2_GUSR1_fwdrvs.pdf", width=10, height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))

for(i in 1:nrow(chtab)) {
  g1_pattern <- paste0(chtab[i,1])
  g2_pattern <- paste0(chtab[i,2])
  
  
  g1_chr <- which(subDfdat[,10] == g1_pattern)
  g2_chr <- which(subDfdat[,14] == g2_pattern)
  
  cc <- intersect(g1_chr, g2_chr)
  subd <- subDfdat[cc, ]
  
  # Flip y-axis if mostly reverse strand
  major_strand <- names(which.max(table(subd$V9)))
  if (major_strand %in% c("+-", "-+")) {
    subd$V12 <- max(subd$V12, na.rm=TRUE) - subd$V12
    subd$V13 <- max(subd$V13, na.rm=TRUE) - subd$V13
  }
  
  yub <- max(subd[,13], na.rm=TRUE)
  xub <- max(subd[,17], na.rm=TRUE)
  
  plot(as.numeric(subd[1,16:17]), as.numeric(subd[1,12:13]), type='n',
       xlim=c(0, xub), ylim=c(0, yub),
       cex.lab=1.4, ylab="GSR Hap2", xlab="GSH Hap2")
  
  title(main=paste("Pair", i, "GUSR Hap1", chtab[i,1], "vs GSH Hap2", chtab[i,2]), cex.main=1.3)
  
  N<-dim(subd)[1]
  for(j in 1:N){
    x <- as.numeric(subd[j,16:17])   
    y <- as.numeric(subd[j,12:13])   
    
    # Fix reversed block endpoint ordering
    if (subd[j,9] %in% c("+-","-+")) {
      y <- rev(y)
    }
    
    if (x[1] > x[2]) {  
      x <- rev(x)                    
      y <- rev(y)                    
    }
    
    if(subd[j,9] == "++") {
      lines(x, y)                     
    }
    else{
      lines(x, y, col="cadetblue")   
    }
  }
  abline(a=0, b=1, col="gray", lty=2)
}


dev.off()

# GSH2_GUSR2 ---------------------------------------------------------------

## read synteny dat
dat<-fread("cactusStripe_TcrGSH2_TcrGUSR2.psl",header=FALSE)
dfdat<-as.data.frame(dat)

head(dat$V10)
head(dat$V14)

dfdat$V10<- as.numeric(gsub(".*Scaffold_([0-9]+).*", "\\1", dfdat$V10))
dfdat$V14<- as.numeric(gsub(".*Scaffold_([0-9]+).*", "\\1", dfdat$V14))


## keep alignments from major scaffolds
xx_query <- table(dfdat[,10])
xx_target <- table(dfdat[,14])

g1 <- names(xx_query)[xx_query > 100]     #lmelCP
g2 <- names(xx_target)[xx_target > 100]   #lmelOLD  

keep <- (dfdat[,10] %in% g1) & (dfdat[,14] %in% g2)
subDfdat <- dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
## normalize with respect to knulli
ntab<-tab
for(i in 1:12){
  ntab[i,]<-ntab[i,]/sum(ntab[i,],na.rm=TRUE)
}

g1_chr <- as.numeric(gsub("Scaffold_([0-9]+)", "\\1", rownames(ntab)))
g2_chr <- as.numeric(gsub("Scaffold_([0-9]+)", "\\1", colnames(ntab)))



pdf("Syn_GSH2_GUSR2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="GUSR Hap2",ylab="GSH Hap2",cex.lab=1.4)
axis(1,at=seq(0,12,length.out=12)/12,g1_chr,las=2)
axis(2,at=seq(0,13,length.out=13)/13,g2_chr,las=2)
box()
dev.off()

chtab<-matrix(c(1,1,
                1,2,
                2,3,
                3,4,
                4,5,
                5,6,
                6,7,
                7,10,
                8,11,
                9,13,
                10,9,
                11,8,
                12,12),nrow=13,ncol=2,byrow=TRUE)



## colinearity plots for all homologous chromsomes,
## filling in for forward/reverse reads

pdf("AlnPlots_GSH2_GUSR2_fwdrvs.pdf", width=10, height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))

for(i in 1:nrow(chtab)) {
  g1_pattern <- paste0(chtab[i,1])
  g2_pattern <- paste0(chtab[i,2])
  
  
  g1_chr <- which(subDfdat[,10] == g1_pattern)
  g2_chr <- which(subDfdat[,14] == g2_pattern)
  
  cc <- intersect(g1_chr, g2_chr)
  subd <- subDfdat[cc, ]
  
  # Flip y-axis if mostly reverse strand
  major_strand <- names(which.max(table(subd$V9)))
  if (major_strand %in% c("+-", "-+")) {
    subd$V12 <- max(subd$V12, na.rm=TRUE) - subd$V12
    subd$V13 <- max(subd$V13, na.rm=TRUE) - subd$V13
  }
  
  yub <- max(subd[,13], na.rm=TRUE)
  xub <- max(subd[,17], na.rm=TRUE)
  
  plot(as.numeric(subd[1,16:17]), as.numeric(subd[1,12:13]), type='n',
       xlim=c(0, xub), ylim=c(0, yub),
       cex.lab=1.4, ylab="GSR Hap2", xlab="GSH Hap2")
  
  title(main=paste("Pair", i, "GUSR Hap2", chtab[i,1], "vs GSH Hap2", chtab[i,2]), cex.main=1.3)
  
  N<-dim(subd)[1]
  for(j in 1:N){
    x <- as.numeric(subd[j,16:17])   
    y <- as.numeric(subd[j,12:13])   
    
    # Fix reversed block endpoint ordering
    if (subd[j,9] %in% c("+-","-+")) {
      y <- rev(y)
    }
    
    if (x[1] > x[2]) {  
      x <- rev(x)                    
      y <- rev(y)                    
    }
    
    if(subd[j,9] == "++") {
      lines(x, y)                     
    }
    else{
      lines(x, y, col="cadetblue")   
    }
  }
  abline(a=0, b=1, col="gray", lty=2)
}


dev.off()
