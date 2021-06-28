library('pROC')
library(kernlab)
library(caret)

setwd('/Users/pingyi/Desktop/minor_final/data')
secretome <- read.delim("~/Desktop/minor_final/data/features_secretome6.xls", row.names=1)
PBDGene<- read.delim("~/Desktop/minor_final/data/experimental_characterized_genes7.txt")
colnames(secretome)

secretome[,4:28]=log2(secretome[,4:28])

train=PBDGene
test=merge(secretome,train,by.x="gene",by.y="ids")
write.table(test,"training_set.xls",sep="\t",col.names=NA)

data=test
colnames(data)
data2=data[,-c(1,64:(ncol(data)-1))]
colnames(data2)

set.seed(123)
inTrain=createDataPartition(y=data2$type,p=0.8,list=FALSE)
training=data2[inTrain,]
testing=data2[-inTrain,]

################### 
training_predict_rf <- read.delim("~/Desktop/minor_final/data/training_predict_rf.xls")
training_predict_glmnet <- read.delim("~/Desktop/minor_final/data/training_predict_glmnet.xls")

training_two_predict <- data.frame(obs = training_predict_rf$obs, rf=training_predict_rf$PBD, glmnet=training_predict_glmnet$PBD)
training_two_predict$mean <- rowMeans(training_two_predict[2:3])

#####calculate ROC
rc = roc(training_two_predict$obs, training_two_predict$mean)
auc(rc)
#### visualize ROC
plot(rc,print.auc=TRUE,plot=TRUE,print.thres=TRUE)

#####make confusionMatrix
train_PBD <- training_two_predict[which(training_two_predict$mean >= 0.5),]
l <- length(train_PBD$obs)
train_PBD$predicted <- 'PBD'
train_nonPBD <- training_two_predict[which(training_all_predict$mean < 0.5),]
l2 <- length(train_nonPBD$obs)
train_nonPBD$predicted <- 'nonPBD'
train_PBD2 <- rbind(train_nonPBD,train_PBD)
confusionMatrix(table(train_PBD2$obs,train_PBD2$predicted))
write.table(train_PBD2,"train_PBD_rfglmnet.xls",sep="\t",col.names=NA)


################### testing set
test_predict_rf <- read.delim("~/Desktop/minor_final/data/test_predict_rf.xls")
test_predict_glmnet <- read.delim("~/Desktop/minor_final/data/test_predict_glmnet.xls")

test_two_predict <- data.frame(type=test_predict_rf$testing.type,rf=test_predict_rf$PBD,glmnet=test_predict_glmnet$PBD)
test_two_predict$mean <- rowMeans(test_two_predict[2:3])

#####calculate ROC
rc2 = roc(test_two_predict$type, test_two_predict$mean)
auc(rc2)
#### visualize ROC
plot(rc2,print.auc=TRUE,plot=TRUE,print.thres=TRUE)
#####make confusionMatrix
test_PBD_predicted <- test_two_predict[which(test_two_predict$mean >= 0.5),]
l3 <- length(test_PBD_predicted$type)
test_PBD_predicted$predicted <- 'PBD'
test_nonPBD_predicted <- test_two_predict[which(test_two_predict$mean < 0.5),]
l4 <- length(test_nonPBD_predicted$type)
test_nonPBD_predicted$predicted <- 'nonPBD'
test_PBD_predicted2 <- rbind(test_nonPBD_predicted,test_PBD_predicted)
confusionMatrix(table(test_PBD_predicted2$type,test_PBD_predicted2$predicted))
write.table(test_PBD_predicted2,"test_PBD_rfglmnet.xls",sep="\t",col.names=NA)



############## varify the result
rf_PBD_predict1 <- read.delim("~/Desktop/minor_final/data/rf_PBD_predict1.xls")
glmnet_PBD_predict1 <- read.delim("~/Desktop/minor_final/data/glmnet_PBD_predict1.xls")
rfglmnet_verify <- data.frame(gene=rf_PBD_predict1$gene,rf=rf_PBD_predict1$PBD.1,glmnet=glmnet_PBD_predict1$PBD.1)
rfglmnet_verify$mean=rowMeans(rfglmnet_verify[,2:3])

verification_dataset <- read.delim("~/Desktop/minor_final/data/verification_dataset.xls")
colnames(verification_dataset)
verification_dataset <- verification_dataset[,-c(1)]

verification_dataset2 <- merge(verification_dataset,rfglmnet_verify,x.by='gene',b.by='gene')
rc=roc(verification_dataset2$type,verification_dataset2$mean)
auc(rc)
par(mfrow=c(1,1))
plot(rc,print.auc=TRUE,plot=TRUE,print.thres=TRUE)


a <- verification_dataset2[which(verification_dataset2$mean<0.5),]
a$predicted <- 'nonPBD'
b <- verification_dataset2[which(verification_dataset2$mean>=0.5),]
b$predicted <- 'PBD'
new4 <- rbind(a,b)


write.table(new4,"verify_rfglmnet_PBD.xls",sep="\t",col.names=NA)
as.
confusionMatrix(table(new4$predicted,new4$type))













