library(readr)
library(randomForest)
load("./models/RF.model.Mar2023.dat")

featureFile = "example/test.txt.gz"
output = "ctcf-insite.persistence.txt"

model = models[["full"]]
dat = read.table(featureFile,header =T, as.is =T)
pred_persistent = predict(model, dat,  type='prob')
cat("Results file: ", output, "\n")
head(pred_persistent [, 2])

write.table(pred_persistent[,2,drop=F], file = output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)








