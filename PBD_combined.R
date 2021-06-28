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
########################


set.seed(345)
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3,classProbs = TRUE,search='grid',summaryFunction = twoClassSummary)
algorithms_to_use <- c('rf', 'glmnet', 'svmRadial','gbm')
library(caretEnsemble)
stacked_models <- caretList(type ~., data=training, trControl=ctrl, methodList=algorithms_to_use)
stacking_results <- resamples(stacked_models)
summary(stacking_results)

xyplot(resamples(stacked_models))
modelCor(resamples(stacked_models))
########## their predicitons are fairly un-correlated, but their overall accuaracy is similar. We do a simple, linear greedy optimization on AUC using caretEnsemble:
#Make a linear regression ensemble
stackControl <- trainControl(method = "repeatedcv", number = 10, repeats = 3,classProbs = TRUE,search='grid',summaryFunction = twoClassSummary,returnData =TRUE,savePredictions =TRUE )
linear <- caretStack(stacked_models, method='glm', trControl=stackControl)
linear
linear$models



############ get the predictions of the training set
orig_fit=linear
orig_fit
#orig_fit$finalModel
orig_fit$error
stark_training_predict_result <-orig_fit$pred

############ calculate ROC 
l2 = length(stark_training_predict_result$pred)
cvrepeated_time <- 3
l_training <- l2/cvrepeated_time
No1_fold <- stark_training_predict_result[c(1:l_training),]
l_training3 <- l_training*2
l_training2 <- l_training+1
No2_fold <- stark_training_predict_result[c(l_training2:l_training3),]
l_training4 <- l_training3+1
l_training5 <- l_training*3
No3_fold <- stark_training_predict_result[c(l_training4:l_training5),]

No1_fold <- merge(No1_fold,No2_fold,by='rowIndex',all=TRUE)
No1_fold <- merge(No1_fold,No3_fold,by='rowIndex',all=TRUE)
No_fold <- data.frame(PBD1=No1_fold$PBD.x,PBD2=No1_fold$PBD.y,PBD3=No1_fold$PBD,obs=No1_fold$obs)
No_fold$PBD= apply(No_fold[1:3],1,mean)
write.table(No_fold,"training_predict_svm.xls",sep="\t",col.names=NA)

library('pROC')
#row.names(test_prediction_result) <- c(1:40)
rc = roc(No_fold$obs, No_fold$PBD)
auc(rc)

#### visualize ROC
ggroc(rc, alpha = 0.5, colour = "red", linetype = 2, size = 2)
plot(rc,print.auc=TRUE,plot=TRUE,print.thres=TRUE)

############ predict on test set

predictions<-predict(orig_fit,newdata=testing)
orig_fit
predictions
confusionMatrix(predictions,testing$type)

predictions_2 <- predict(orig_fit, newdata = testing, type = "prob") 
test_prediction_result <- data.frame(predictions_2[,1:2],predictions,testing$type)



############## calculate ROC 
library('pROC')
#row.names(test_prediction_result) <- c(1:40)
rc = roc(test_prediction_result$testing.type, test_prediction_result$PBD)
auc(rc)

#### visualize ROC
ggroc(rc, alpha = 0.5, colour = "red", linetype = 2, size = 2)

plot(rc,print.auc=TRUE,plot=TRUE,print.thres=TRUE)

##############
model=orig_fit
new=predict(model, newdata = secretome, type = "prob") 
newPredict<-predict(model,newdata=secretome)
new=cbind(secretome,new) 
new2=cbind(new,newPredict) 
write.table(new2,"svm_PBD_predict.xls",sep="\t",col.names=NA)

############## varify the result
vm_PBD_predict <- read.delim("~/Desktop/minor_final/data/svm_PBD_predict.xls", row.names=1)
predicted_Aspergillus_PBD_genes <- read.delim("~/Desktop/minor_final/data/predicted_Aspergillus_PBD_genes_Ronald (1).txt")

f <- vm_PBD_predict[order(-vm_PBD_predict$PBD.1),]
write.table(f,"svm_PBD_predict2.xls",sep="\t",col.names=NA)
f1 <- data.frame(gene=f$gene,PBD=f$PBD,nonPBD=f$nonPBD,before=f$PBD.1,newPredict=f$newPredict)


f2 <- f1[f1$newPredict=='PBD',]
predicted_PBD_list <- f2$gene


predicted_Aspergillus_PBD_genes2 <- predicted_Aspergillus_PBD_genes[predicted_Aspergillus_PBD_genes$Secretion.signal..SignalP.=='Yes',]
b <- predicted_Aspergillus_PBD_genes2$CBS_ids
c <- secretome$gene
d <- intersect(b,c)
n <- intersect(d,predicted_PBD_list)
length(n)
length(d)
length(n)/length(d)


