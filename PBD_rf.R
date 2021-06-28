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

############ tune the rf
########## ntree=500
set.seed(345) 
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3,allowParallel = TRUE,classProbs = TRUE,search='grid',summaryFunction = twoClassSummary)
set.seed(345)   ## set.seed(345) could reach 83 on training set and 76% accuracy on 70pectinolytic gene
ctrl$sampling <- "down"
tunegrid <- expand.grid(mtry = 1:30)
set.seed(345)
rf_default1 <- train(type~., 
                    data=training, 
                    method='rf', 
                    #metric='Accuracy',  #Metric compare model is Accuracy
                    metric = "ROC",
                    tuneGrid=tunegrid,
                    ntree=500,
                    trControl=ctrl)
rf_default1
plot(rf_default1)
rf_default1$bestTune

########## ntree=1000

set.seed(345)
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3, allowParallel = TRUE,classProbs = TRUE,search='grid',summaryFunction = twoClassSummary)
set.seed(345)   ## set.seed(345) could reach 83 on training set and 76% accuracy on 70pectinolytic gene
ctrl$sampling <- "down"
tunegrid <- expand.grid(mtry = 1:30)
set.seed(345)
rf_default2 <- train(type~.,
                    data=training,
                    method='rf',
                    #metric='Accuracy',  #Metric compare model is Accuracy
                    metric = "ROC",
                    tuneGrid=tunegrid,
                    ntree=1000,
                    trControl=ctrl
)
rf_default2
plot(rf_default2)
rf_default2$bestTune
####plot ntree=500 & ntree=1000 together

par(mfrow=c(1,1))
plot(rf_default1$results$mtry,rf_default1$results$ROC,col='blue',type = 'l',xlab='mtry',ylab = 'ROC(Repeated Cross-Validation)',
        main='Random Forest Hyperparameter Selection')
lines(rf_default2$results$mtry,rf_default2$results$ROC,col='red')
legend("bottomright", legend=c('mtry=500','mtry=1000'),
       fill = c('blue','red'))

############# tune 'rf' with parameter  mtry = 27, ntrr=500
set.seed(345) 
#ctrl = trainControl(method = "LOOCV", classProbs = TRUE,search='grid',summaryFunction = twoClassSummary,returnData =TRUE,savePredictions = TRUE)
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3,classProbs = TRUE,search='grid')
#ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3,classProbs = TRUE,search='grid',summaryFunction = twoClassSummary,returnData =TRUE,savePredictions = TRUE)
#ctrl = trainControl(method = "repeatedcv", number = 4, repeats = 10,classProbs = TRUE,search='grid',summaryFunction = twoClassSummary,returnData =TRUE,savePredictions = TRUE)
#ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3,classProbs = TRUE,search='grid')
set.seed(345)   ## set.seed(345) could reach 83 on training set and 76% accuracy on 70pectinolytic gene
ctrl$sampling <- "down"
tunegrid <- expand.grid(mtry=27)
set.seed(345)
rf_tuned <- train(type~., 
                    data=training, 
                    method='rf', 
                    #metric='Accuracy',  #Metric compare model is Accuracy
                    #metric = "ROC",
                    ntree=500,
                    tuneGrid=tunegrid,
                    trControl=ctrl)
#print(rf_default$finalModel)
#summary(rf_tuned)
print(rf_tuned)

set.seed(345) 
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3,classProbs = TRUE,search='grid',summaryFunction = twoClassSummary,returnData =TRUE,savePredictions = TRUE)
set.seed(345)   ## set.seed(345) could reach 83 on training set and 76% accuracy on 70pectinolytic gene
ctrl$sampling <- "down"
tunegrid <- expand.grid(mtry=27)
set.seed(345)
rf_tuned <- train(type~., 
                  data=training, 
                  method='rf', 
                  #metric='Accuracy',  #Metric compare model is Accuracy
                  metric = "ROC",
                  ntree=500,
                  tuneGrid=tunegrid,
                  trControl=ctrl)
#print(rf_default$finalModel)
#summary(rf_tuned)
print(rf_tuned)

############ get the predictions of the training set
orig_fit=rf_tuned
orig_fit
orig_fit$finalModel
#orig_fit$finalModel$confusion
rf_training_predict_result <-  orig_fit$pred

############ plot ROC 
l2 = length(rf_training_predict_result$pred)
cvrepeated_time <- 3
l_training <- l2/cvrepeated_time
No1_fold <- rf_training_predict_result[c(1:l_training),]
l_training3 <- l_training*2
l_training2 <- l_training+1
No2_fold <- rf_training_predict_result[c(l_training2:l_training3),]
l_training4 <- l_training3+1
l_training5 <- l_training*3
No3_fold <- rf_training_predict_result[c(l_training4:l_training5),]

No1_fold <- merge(No1_fold,No2_fold,by.x='rowIndex',by.y='rowIndex',all=TRUE)
No1_fold <- merge(No1_fold,No3_fold,by.x='rowIndex',by.y='rowIndex',all=TRUE)
No_fold <- data.frame(PBD1=No1_fold$PBD.x,PBD2=No1_fold$PBD.y,PBD3=No1_fold$PBD,obs=No1_fold$obs)
No_fold$PBD= apply(No_fold[1:3],1,mean)
write.table(No_fold,"training_predict_rf.xls",sep="\t",col.names=NA)

library('pROC')
#row.names(test_prediction_result) <- c(1:40)
rc = roc(No_fold$obs, No_fold$PBD)

#### visualize ROC
plot(rc,print.auc=TRUE,plot=TRUE,print.thres=TRUE)

############ predict on test set

predictions<-predict(orig_fit,newdata=testing)
predictions

confusionMatrix(predictions,testing$type)

predictions_2 <- predict(orig_fit, newdata = testing, type = "prob") 
test_prediction_result <- data.frame(predictions_2[,1:2],predictions,testing$type)

write.table(test_prediction_result,"test_predict_rf.xls",sep="\t",col.names=NA)

############## calculate ROC 
library('pROC')
#row.names(test_prediction_result) <- c(1:40)
rc = roc(test_prediction_result$testing.type, test_prediction_result$PBD)
auc(rc)

#### visualize ROC
plot(rc,print.auc=TRUE,plot=TRUE,print.thres=TRUE)

############## predict the result for 'secretome' dataset
model=rf_tuned
new=predict(model, newdata = secretome, type = "prob") 
newPredict<-predict(model,newdata=secretome)
new=cbind(secretome,new) 
new2=cbind(new,newPredict) 
write.table(new2,"rf_PBD_predict1.xls",sep="\t",col.names=NA)

############## varify the result
colnames(secretome)
verification_dataset <- read.delim("~/Desktop/minor_final/data/verification_dataset.xls")
colnames(verification_dataset)
verification_dataset <- verification_dataset[,-c(1)]

new3=predict(model, newdata = verification_dataset, type = "prob") 
rc=roc(verification_dataset$type,new3$PBD)
auc(rc)
par(mfrow=c(1,1))
plot(rc,print.auc=TRUE,plot=TRUE,print.thres=TRUE)

newPredict2<-predict(model,newdata=verification_dataset)
new3 <- cbind(verification_dataset,new3)
new4=cbind(new3,newPredict2)
write.table(new4,"verify_rf_PBD.xls",sep="\t",col.names=NA)
confusionMatrix(new4$newPredict2,new4$type)





