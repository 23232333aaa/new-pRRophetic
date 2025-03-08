
rm(list = ls())
### 1.加载镜像
options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

###2.安装R包
if(!requireNamespace("car",quietly = TRUE)) install.packages("car",update = F,ask = F)
if(!requireNamespace("ridge",quietly = TRUE)) install.packages("ridge",update = F,ask = F)
if(!requireNamespace("preprocessCore",quietly = TRUE)) BiocManager::install("preprocessCore",update = F,ask = F)
if(!requireNamespace("genefilter",quietly = TRUE)) BiocManager::install("genefilter",update = F,ask = F)
if(!requireNamespace("sva",quietly = TRUE)) BiocManager::install("sva",update = F,ask = F)
install.packages("./resource/pRRophetic_Guozi/", repos = NULL,type = "source")

###3.使用测试
library(pRRophetic)
library(ggplot2)
set.seed(12345)
data("bortezomibData") #exprDataBortezomib, bortIndex, studyResponse and studyIndex
pRRopheticQQplot("Bortezomib")

cvOut <- pRRopheticCV("Bortezomib", cvFold=5, testExprData=exprDataBortezomib)
summary(cvOut)
plot(cvOut)

predictedPtype <- pRRopheticPredict(exprDataBortezomib, "Bortezomib", selection=1)
predictedPtype_blood <- pRRopheticPredict(exprDataBortezomib, "Bortezomib", "blood", selection=1)
predictedPtype_solid <- pRRopheticPredict(exprDataBortezomib, "Bortezomib", "allSolidTumors", selection=1)


t.test(predictedPtype[((studyResponse == "PGx_Responder = NR") & bortIndex)],
       predictedPtype[((studyResponse == "PGx_Responder = R") & bortIndex)],
       alternative="greater")

t.test(predictedPtype_blood[((studyResponse == "PGx_Responder = NR") & bortIndex)],
       predictedPtype_blood[((studyResponse == "PGx_Responder = R") & bortIndex)],
       alternative="greater")

t.test(predictedPtype_solid[((studyResponse == "PGx_Responder = NR") & bortIndex)],
       predictedPtype_solid[((studyResponse == "PGx_Responder = R") & bortIndex)],
       alternative="greater")

### 画图
df <- stack(list(NR=predictedPtype_blood[((studyResponse == "PGx_Responder = NR")& bortIndex)], 
                 R=predictedPtype_blood[((studyResponse == "PGx_Responder = R") & bortIndex)]))
ggplot(data=df, aes(y=values, x=ind)) +
  geom_boxplot(alpha=.3, fill=c("#CC0033", "#006633")) + 
  theme_bw() + 
  ylab("Predicted Bortezomib Sensitivity") + 
  xlab("Clinical Response")

