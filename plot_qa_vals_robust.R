library("robustbase");
library("pcaPP");
	
ts19 <- read.table("ts19a_p1_new_mapwise_stats.out",header=T);

ts19.pca <- PCAproj(ts19,k=5,center=colMeans(ts19),scale=sd(ts19))
ts19.pre <- predict(ts19.pca);

ts19.cov <- cov(ts19.pre[,1:5]);
ts19.m <- colMeans(ts19.pre[,1:5]);
ts19.mca <- covMcd(ts19.pre[,1:5]);

ts19.ma <- mahalanobis(ts19.pre[,1:5],ts19.mca$center,ts19.mca$cov);

png(file="ts19_maha_qvals_robust.png",width=800,height=800);
par(ps=14);
plot(ts19.ma,type="n",main="Merck ts19a/b (robust estimator)",xlab="# data set",ylab="Mahalanobis distance");
text(1:54,ts19.ma,as.character(1:54));
abline(h=qchisq(.9,5.0),col="blue");
abline(h=qchisq(.95,5.0),col="green");
legend("topleft",c("90th percentile","95th percentile"),col=c("blue","green"),lwd=3);
dev.off();

# -------------------------------- BSA data set charite 1
bc1 <- read.table("bsa_charite1_mapwise_stats.out",header=T);

bc1.pca <- PCAproj(bc1,k=5,center=colMeans,scale=sd);
bc1.pre <- predict(bc1.pca);

bc1.m <- colMeans(bc1.pre[,1:5]);
bc1.cov <- cov(bc1.pre[,1:5]);
bc1.mca <- covMcd(bc1.pre[,1:5]);

bc1.ma <- mahalanobis(bc1.pre[,1:5],bc1.mca$center,bc1.mca$cov);

png(file="bsa_charite1_maha_qvals_robust.png",width=800,height=800);
par(ps=14);
plot(bc1.ma,type="n",main="Charite bsa1 (robust estimator)",xlab="# data set",ylab="Mahalanobis distance");
text(1:22,bc1.ma,as.character(1:22));
abline(h=qchisq(.9,5.0),col="blue");
abline(h=qchisq(.95,5.0),col="green");
legend("topleft",c("90th percentile","95th percentile"),col=c("blue","green"),lwd=3);
dev.off();

# -------------------------- BSA data set charite 2
bc2 <- read.table("bsa_charite2_mapwise_stats.out",header=T);

bc2.pca <- PCAproj(bc2,k=5,center=colMeans,scale=sd);
bc2.pre <- predict(bc2.pca);

bc2.m <- colMeans(bc2.pre[,1:5]);
bc2.cov <- cov(bc2.pre[,1:5]);
bc2.mca <- covMcd(bc2.pre[,1:5]);

bc2.ma <- mahalanobis(bc2.pre[,1:5],bc2.mca$center,bc2.mca$cov);

png(file="bsa_charite2_maha_qvals_robust.png",width=800,height=800);
par(ps=14);
plot(bc2.ma,type="n",main="Charite bsa2 (robust estimator)",xlab="# data set",ylab="Mahalanobis distance");
text(1:22,bc2.ma,as.character(1:22));
abline(h=qchisq(.9,5.0),col="blue");
abline(h=qchisq(.95,5.0),col="green");
legend("topleft",c("90th percentile","95th percentile"),col=c("blue","green"),lwd=3);
dev.off();


# -------------------- BSA data set Merck
me1 <- read.table("bsa_merck_mapwise_stats.out",header=T);

me1.pca <- PCAproj(me1,k=5,center=colMeans,scale=sd);
me1.pre <- predict(me1.pca);

me1.m <- colMeans(me1.pre[,1:5]);
me1.cov <- cov(me1.pre[,1:5]);
me1.mca <- covMcd(me1.pre[,1:5]);

me1.ma <- mahalanobis(me1.pre[,1:5],me1.mca$center,me1.mca$cov);

png(file="bsa_merck_maha_qvals_robust.png",width=800,height=800);
par(ps=14);
plot(me1.ma,type="n",main="Merck bsa (robust estimator)",xlab="# data set",ylab="Mahalanobis distance");
text(1:22,me1.ma,as.character(1:22));
abline(h=qchisq(.9,5.0),col="blue");
abline(h=qchisq(.95,5.0),col="green");
legend("topleft",c("90th percentile","95th percentile"),col=c("blue","green"),lwd=3);
dev.off();

