## compute depth per individual and SNP to identify SNPs and individuals to drop
## the idea is to get rid of low coverage individuals
## and SNPs with either very high coverage (3SD > mean) or very high variance in coverage
## across individuals

library(data.table)
d<-as.matrix(fread("depth_oneref.txt",header=FALSE))

## mean and SD by SNP
mnc<-apply(d,1,mean)
sdc<-apply(d,1,sd)

## CV 
cvc<-sdc/mnc
meancvc<-mean(cvc)
quantcvc<-quantile(cvc,probs=c(.5,.9,.95,.99,.999,1))

meanmnc3sd<-mean(mnc)+3*sd(mnc)
mean(mnc)
quantmnc<-quantile(mnc,probs=c(.5,.9,.95,.99,.999,1))

## for SNPs, keep if CV < 99.9 percentile and mean < meanmnc3sd
keepSNPs<-as.numeric(mnc < 15 & cvc < 1)
keepSNPs_mean<-mean(keepSNPs)
keepSNPs_sum<-sum(keepSNPs)

## for individuals

mni<-apply(d,2,mean)
plot(sort(mni))
summary(mni)
quantile(mni,probs=c(.025,.99))

## kwwp if mni > 2.5% or <99%
keepInds<-as.numeric(mni > 2.7 & mni < 9)
mean(keepInds)
sum(keepInds)

cat(sprintf("meancvc=%.4f\n", meancvc))
cat(sprintf("quantcvc=%s\n", paste(names(quantcvc), round(quantcvc, 4), sep="=", collapse=", ")))
cat(sprintf("meanmnc3sd=%.4f\n", meanmnc3sd))
cat(sprintf("meanmnc=%.4f\n", mean(mnc)))
cat(sprintf("quantmnc=%s\n", paste(names(quantile(mnc, probs=c(.5,.9,.95,.99,.999,1))), 
                                    round(quantile(mnc, probs=c(.5,.9,.95,.99,.999,1)), 4), 
                                    sep="=", collapse=", ")))
cat(sprintf("keepSNPs_mean=%.4f\n", mean(keepSNPs)))
cat(sprintf("keepSNPs_sum=%d\n", sum(keepSNPs)))
cat(sprintf("meanmni=%s\n", paste(names(summary(mni)), round(summary(mni), 4), sep="=", collapse=", ")))
cat(sprintf("quantmni=%s\n", paste(names(quantile(mni, probs=c(.025,.99))), 
                                    round(quantile(mni, probs=c(.025,.99)), 4), 
                                    sep="=", collapse=", ")))

write.table(file="KeepInds.txt",keepInds,row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(file="KeepSNPs.txt",keepSNPs,row.names=FALSE,col.names=FALSE,quote=FALSE)
