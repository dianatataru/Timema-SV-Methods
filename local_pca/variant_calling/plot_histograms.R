library(ggplot2)

data<- read.table("FHA_all_oneref_site_metrics.txt", header=F)
colnames(data)<-c("CHROM","POS","QUAL","MQBZ","SCBZ","RPBZ","DP")
data$SCBZ<-as.numeric(data$SCBZ)
data$MQBZ<-as.numeric(data$MQBZ)
data$QUAL<-as.numeric(data$QUAL)
data$RPBZ<-as.numeric(data$RPBZ)
data$DP<-as.numeric(data$DP)
summary(data)

pdf("FHA_all_oneref_histograms.pdf", width=8, height=6)


p2 <- ggplot(data, aes(x = SCBZ)) +
  geom_histogram(binwidth = 1) +
  xlim(-20,20)+ylim(0,4000000)+
  labs(title = "Histogram of SCBZ", x = "SCBZ", y = "Frequency")
print(p2)

p3 <- ggplot(data, aes(x = MQBZ)) +
  geom_histogram(binwidth = 1) +
  xlim(-20,20)+ylim(0,9000000)+
  labs(title = "Histogram of MQBZ", x = "MQBZ", y = "Frequency")
print(p3)

p4 <- ggplot(data, aes(x = QUAL)) +
  geom_histogram(binwidth = 10) +
  xlim(0,1000)+ylim(0,1000000)+
  labs(title = "Histogram of QUAL", x = "QUAL", y = "Frequency")
print(p4)

p6 <- ggplot(data, aes(x = RPBZ)) +
  geom_histogram(binwidth = 1) +
  xlim(-25,25)+ylim(0,740000)+
  labs(title = "Histogram of RPBZ", x = "RPBZ", y = "Frequency")
print(p6)

p7 <- ggplot(data, aes(x = DP )) +
  geom_histogram(binwidth =10) +
  xlim(0,2000)+ylim(0,250000)+
  labs(title = "Histogram of DP", x = "DP", y = "Frequency")
print(p7)
#exclude variants with read depth larger than twice the mean
dev.off()
