	
ts19 <- read.table("ts19a_p1_new_mapwise_stats.out",header=T);

ts19.pca <- prcomp(ts19,scale=T);
ts19.pre <- predict(ts19.pca);

ts19.cov <- cov(ts19.pre[,1:5]);
ts19.m <- colMeans(ts19.pre[,1:5]);

ts19.ma <- mahalanobis(ts19.pre[,1:5],ts19.m,ts19.cov);

png(file="ts19_maha_qvals_normal.png",width=800,height=800);
par(ps=14);
plot(ts19.ma,type="n",main="Merck ts19a/b (normal)",xlab="# data set",ylab="Mahalanobis distance");
text(1:54,ts19.ma,as.character(1:54));
abline(h=qchisq(.9,5.0),col="blue");
abline(h=qchisq(.95,5.0),col="green");
legend("topleft",c("90th percentile","95th percentile"),col=c("blue","green"),lwd=3);
dev.off();

# ---------------------- BSA data set charite 1
bc1 <- read.table("bsa_charite1_mapwise_stats.out",header=T);

bc1.pca <- prcomp(bc1,scale=T);
bc1.pre <- predict(bc1.pca);

bc1.m <- colMeans(bc1.pre[,1:5]);
bc1.cov <- cov(bc1.pre[,1:5]);

bc1.ma <- mahalanobis(bc1.pre[,1:5],bc1.m,bc1.cov);

png(file="bsa_charite1_maha_qvals_normal.png",width=800,height=800);
par(ps=14);
plot(bc1.ma,type="n",main="Charite bsa1 (normal)",xlab="# data set",ylab="Mahalanobis distance");
text(1:22,bc1.ma,as.character(1:22));
abline(h=qchisq(.9,5.0),col="blue");
abline(h=qchisq(.95,5.0),col="green");
legend("topleft",c("90th percentile","95th percentile"),col=c("blue","green"),lwd=3);
dev.off();

# ---------------------- BSA data set charite 2
bc2 <- read.table("bsa_charite2_mapwise_stats.out",header=T);

bc2.pca <- prcomp(bc2,scale=T);
bc2.pre <- predict(bc2.pca);

bc2.m <- colMeans(bc2.pre[,1:5]);
bc2.cov <- cov(bc2.pre[,1:5]);

bc2.ma <- mahalanobis(bc2.pre[,1:5],bc2.m,bc2.cov);

png(file="bsa_charite2_maha_qvals_normal.png",width=800,height=800);
par(ps=14);
plot(bc2.ma,type="n",main="Charite bsa2 (normal)",xlab="# data set",ylab="Mahalanobis distance");
text(1:22,bc2.ma,as.character(1:22));
abline(h=qchisq(.9,5.0),col="blue");
abline(h=qchisq(.95,5.0),col="green");
legend("topleft",c("90th percentile","96th percentile"),col=c("blue","green"),lwd=3);
dev.off();

# -------------------- BSA data set merck
me1 <- read.table("bsa_merck_mapwise_stats.out",header=T);

me1.pca <- prcomp(me1);  #,scale=T);
me1.pre <- predict(me1.pca);

me1.m <- colMeans(me1.pre[,1:5]);
me1.cov <- cov(me1.pre[,1:5]);

me1.ma <- mahalanobis(me1.pre[,1:5],me1.m,me1.cov);

png(file="bsa_merck_maha_qvals_normal.png",width=800,height=800);
par(ps=14);
plot(me1.ma,type="n",main="Merck bsa (normal)",xlab="# data set",ylab="Mahalanobis distance");
text(1:22,me1.ma,as.character(1:22));
abline(h=qchisq(.9,5.0),col="blue");
abline(h=qchisq(.95,5.0),col="green");
legend("topleft",c("90th percentile","95th percentile"),col=c("blue","green"),lwd=3);
dev.off();
