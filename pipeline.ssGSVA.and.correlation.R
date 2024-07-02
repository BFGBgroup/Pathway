################## Step 4. gene expression —> ssGSVA ####
#Pathway enrichment analysis for each individual. The gene expression profile predicted in the previous step and the gene sets of pathways are required. For the ssGSVA analysis, we have integrated gene sets of 1,649,335 human pathways from multiple databases, including GSEA (MSigDB v7.5.1, URL: https://www.gsea-msigdb.org/gsea/index.jsp), DAVID (v.2022q2, URL: https://david.ncifcrf.gov/), Reactome (v81, URL: https://reactome.org/), NetPath (URL: http://netpath.org/), PANTHER (URL: http://www.pantherdb.org/), WikiPathways (v.20220710, URL: https://classic.wikipathways.org/index.php/), and PathBank (URL: https://pathbank.org/). All gene names or gene IDs were transformed to ensemble gene IDs. Meanwhile, the pathway names were organized in the format of type~ID~term, for example, GOBP~GO:0000002~mitochondrial genome maintenance. Some pathways were in type~term format since they do not have IDs. In addition, the same pathway in different databases was subject to the newer version. The integrated data can be downloaded here in Rdata format (geneSets.Rdata).
#Please use R software to perform this step. Please install dependent R packages GSVA, psych, HHG, doParallel, etc. in advance. The codes and annotations please refer to the R script pipeline.ssGSVA.and.correlation.R.

setwd(your.work.dir)
library(GSVA)
library(psych)
library(HHG)

#4.1 Load the gene expression profiles
geneexpression <- read.table("prediction_OUTPUT_DIR/predict.txt",header = TRUE, row.names = 1)

#4.2 Load the gene sets of pathways
load("geneSets.Rdata")

#4.3 ssGSVA analysis 
gsva_es <- gsva(geneexpression, geneSets, method='gsva', kcdf='Gaussian', abs.ranking=TRUE, mx.diff=TRUE, min.sz=2)
#Please refer to the documentation of the gsva function for parameters.


################## Step 5. Association between disease risk (PRS) and pathways (GSVA) ####
#Analysis of pathways significantly associated with disease risk. This step is connected to the previous part, please continue to use R software.
#Since most ssGSVA scores did not follow a normal distribution (Shapiro-Wilk test p < 0.05), and the correlation pattern between ssGSVA score and disease risk scores was unknown, we employed the HHG method to assess the correlation. If the number of pathways is large, it is recommended to use parallel computing. Please refer to the R package doParallel. It is worth noting that the significant correlation pathways identified by HHG (P < 0.001) aligned with those identified by the permutation procedure. Therefore, if you use HHG P value < 0.001 as the threshold for a pathway to be significantly associated with disease risk, you can ignore the permutation calculation in step 5.6. If using permutation procedure to calculate, please ignore step 5.5. And since the recommended number of calculations for each permutation procedure is 1,000 to10,000 times, which takes a long time, please use parallel computing to produce results as quickly as possible. It is necessary to check whether the parallel computing results are complete, that is, to check whether all pathways have produced permutation results. This is because the calculation process undertaken by some cores may fail, but the program can still continue to run and produce uncomplete results.
#The codes and annotations please refer to the R script pipeline.ssGSVA.and.correlation.R.


#5.1 Load disease riks score
PRS_result <- read.table("OUTPUT_DIR/PRScs.final.result.profile",header = T,row.names = NULL)
PRS_result <- data.frame(SampleID = paste(PRS_result$FID, PRS_result$IID, sep = "_"), PRScore=PRS_result$SCORESUM) 
#Please pay attention to the format of "SampleID" here, which needs to be the same as the "SampleID" of gsva_es_result below.

#5.2 Format conversion of ssGSVA score
gsva_es_result <- as.data.frame(t(gsva_es)) #Convert to data.frame with samples in rows and pathways in columns
gsva_es_result$SampleID <- rownames(gsva_es_result) 
#Please pay attention to the format of "SampleID" here, which needs to be the same as the "SampleID" of PRS_result above.

#5.3 Merge data.frame to calculate the correlation between pathways and disease risk
cor_data <- merge(PRS_result, gsva_es_result, by="SampleID")

#5.4 Remove outliers
cor_data_outlier <- data.frame(SampleID = cor_data[, 1], lapply(cor_data[, -1], function(x){winsor(x)}))

#5.5 Calculate correlation by HHG (optional)
#If there are many pathways, parallel computing can be used to save computation time. You can refer to the use of doParallel in step 5.6 below.

#(1) For sample sizes < 100 , Perform the hhg test:
X = cor_data_outlier[,2]
Dx = as.matrix(dist((X), diag = TRUE, upper = TRUE))

result.table = data.frame(matrix(nrow = ncol(cor_data_outlier)-3,ncol = 9))
colnames(result.table) = c("pathway","sum.chisq","P.sum.chisq","sum.lr","P.sum.lr", "max.chisq","P.max.chisq","max.lr","P.max.lr")

for (i in 3:ncol(cor_data_outlier)) {
  Y = cor_data_outlier[,i]
  Dy = as.matrix(dist((Y), diag = TRUE, upper = TRUE))
  #Set set.seed to maintain the same result for multiple repeated calculations.
  set.seed(1234)
  hhg = hhg.test(Dx, Dy, nr.perm = 1000,nr.threads = 0) 
  #nr.perm: the number of permutations calculated, defaults to 10000. Default nr. threads=0, try using all available kernels.
  result.table[i-2,1] = colnames(cor_data)[i]
  result.table[i-2,2] = hhg$sum.chisq
  result.table[i-2,3] = hhg$perm.pval.hhg.sc
  result.table[i-2,4] = hhg$sum.lr
  result.table[i-2,5] = hhg$perm.pval.hhg.sl
  result.table[i-2,6] = hhg$max.chisq
  result.table[i-2,7] = hhg$perm.pval.hhg.mc
  result.table[i-2,8] = hhg$max.lr
  result.table[i-2,9] = hhg$perm.pval.hhg.ml
}
sig_res = subset(result.table,P.sum.chisq<0.001)

#Save the result to HHG_result.filename (set according to your preference)
write.table(result.table,file = HHG_result.filename ,sep="\t",row.names = FALSE,quote = FALSE)

#(2) For sample sizes>100, Fast.independence.test is the reccomended option:
result.table = data.frame(matrix(nrow = ncol(cor_data_outlier)-3,ncol = 3))
colnames(result.table) = c("pathway","MinP","MinP.pvalue")

N_Large = SAMPLE.NUMBER  #the number of samples in your data
NullTable_for_N_Large_MXL_tables = Fast.independence.test.nulltable(N_Large, variant = 'ADP-EQP',nr.atoms = 40, nr.perm=10000)
#nr.atoms: Brill (2016) suggests a minimum of 40 atoms, with an increase of up to 60 for alternatives which are more difficult to detect (on the expense of computational complexity
for (i in 3:ncol(cor_data_outlier)) {
  X_Large = cor_data_outlier[,2]
  Y_Large = cor_data_outlier[,i]
  
  #Set set.seed to maintain the same result for multiple repeated calculations.
  set.seed(1234)
  ADP_EQP_ML_Result = Fast.independence.test(X_Large,Y_Large, NullTable_for_N_Large_MXL_tables) #HHG
  result.table[i-2,1] = colnames(cor_data)[i]
  result.table[i-2,2] = ADP_EQP_ML_Result$MinP
  result.table[i-2,3] = ADP_EQP_ML_Result$MinP.pvalue
  #MinP: The test statistic when the combining type is "MinP".
  #MinP.pvalue: The p-value when the combining type is "MinP"
}
sig_res = subset(result.table,MinP.pvalue<0.001)

#Save the result to HHG_result.filename (set according to your preference)
write.table(result.table,file = HHG_result.filename ,sep="\t",row.names = FALSE,quote = FALSE)


#5.6 Calculate correlation by HHG and permutation process (optional)
#The permutation calculation takes a long time, and it is strongly recommended to use parallel calculation. Parallel computing relies on the core of the computer. The detailed usage please refer to the documentation of the R package doParallel.

#The permutation process is：Permutations were firstly performed by randomizing the disease risk scores and pathways ssGSVA scores. Between 1000 and 10,000 permutations were carried out for the correlation analysis between each pathway and the disease risk scores, and the stopping condition was set within this range to ensure that at least 15 permuted p-values were smaller than the nominal p-value. The p-threshold of each pathway was determined by the following procedure: Firstly, the empirical p-value of each pathway was defined as the number of permuted p-values smaller than the nominal p-value divided by the total number of permuted p-values (formula n/permutation times). Then, for each brain region, the empirical p-values were corrected to empirical q-values and the threshold parameter was defined as the empirical p-value that corresponds to the empirical q-value closest to 0.05. Finally, the p-threshold value for each pathway is determined as follows: it is the nth value obtained after sorting the permutation p-values from smallest to largest, where n represents the integer resulting from multiplying the threshold parameter of each brain region by the number of permutations for each pathway (ceiling(threshold parameter × permutation times)).

#Take sample sizes>100 as an example for HHG analysis:

library(doParallel)

ALLpath_permutation_p = list()
N_Large = SAMPLE.NUMBER #the number of samples in your data
NullTable_for_N_Large_MXL_tables = Fast.independence.test.nulltable(N_Large, variant = 'ADP-EQP',nr.atoms = 40, nr.perm=10000)
#nr.atoms: Brill (2016) suggests a minimum of 40 atoms, with an increase of up to 60 for alternatives which are more difficult to detect (on the expense of computational complexity

#Parallel computing, set the number of cores according to your computer configuration
num_cores <- 50

#Start parallel computing
registerDoParallel(num_cores)
results <- foreach(i = 3:ncol(cor_data_outlier), .combine = rbind) %dopar% {
  permutation_HHG_result = data.frame(matrix(ncol = 5))
  colnames(permutation_HHG_result) = c("pathway","HHG.test.statistic","HHG.pvalue","permutation.cycles.number","permutation.empirical.p")
  permutation_HHG_result$pathway[1] = colnames(cor_data_outlier)[i]
  
  #The nominal HHG p value
  X_Large = cor_data_outlier[,2]
  Y_Large = cor_data_outlier[,i]
  
  #Set set.seed to maintain the same result for multiple repeated calculations.
  set.seed(1234) 
  ADP_EQP_ML_Result = Fast.independence.test(X_Large,Y_Large, NullTable_for_N_Large_MXL_tables) #HHG
  permutation_HHG_result$HHG.test.statistic[1] = ADP_EQP_ML_Result$MinP
  permutation_HHG_result$HHG.pvalue[1] = ADP_EQP_ML_Result$MinP.pvalue
  #MinP: The test statistic when the combining type is "MinP".
  #MinP.pvalue: The p-value when the combining type is "MinP"
  HHG_res_p = ADP_EQP_ML_Result$MinP.pvalue
  
  #permuted p-value
  result.table = data.frame(matrix(ncol = 3))
  colnames(result.table) = c("pathway","permutation.test.statistic","permutation.p")
  
  for (k in 1:10000){
    if(k<=1000|nrow(subset(result.table,permutation.p<HHG_res_p))<15){  
      #随机打乱，以毫秒时间为随机种子。
      set.seed(as.numeric(gsub("[^0-9]", "", format(Sys.time(), "%H:%M:%OS3"))))
      X_Large = sample(cor_data_outlier[,2])
      Y_Large = sample(cor_data_outlier[,i])
      
      set.seed(1234)
      ADP_EQP_ML_Result = Fast.independence.test(X_Large,Y_Large, NullTable_for_N_Large_MXL_tables) #HHG
      
      hhg_permutation_res = data.frame(colnames(cor_data_outlier)[i],ADP_EQP_ML_Result$MinP,ADP_EQP_ML_Result$MinP.pvalue)
      colnames(hhg_permutation_res) = c("pathway","permutation.test.statistic","permutation.p")
      
      result.table = rbind(result.table,hhg_permutation_res)
    }
  }
  result.table = result.table[-1,]
  ALLpath_permutation_p[[colnames(cor_data_outlier)[i]]] = result.table[,3]
  
  #empirical p-value
  permutation_HHG_result$permutation.cycles.number[1] = nrow(result.table)
  smallerp = subset(result.table,permutation.p<HHG_res_p) 
  permutation_HHG_result$permutation.empirical.p[1] = nrow(smallerp)/nrow(result.table)
  
  return(list(ALLpath_permutation_p, permutation_HHG_result))
}
#End parallel computing
stopImplicitCluster()

#Result organization
ALLpath_permutation_p <- results[, 1]
for(i in 1:length(ALLpath_permutation_p)){
  names(ALLpath_permutation_p)[[i]] <- names(ALLpath_permutation_p[[i]])
}

permutation_HHG_result1 <- results[, 2]
permutation_HHG_result1 <- do.call(rbind,permutation_HHG_result1)

#Calculate the empirical q-value of permutation
permutation_HHG_result1 = permutation_HHG_result1[order(permutation_HHG_result1$permutation.empirical.p),]
permutation_HHG_result1$permutation.empirical.q = p.adjust(permutation_HHG_result1$permutation.empirical.p,method = "BH")

#Calculate the p-value threshold corresponding to each pathway
#The empirical p value of each pathway was identified when the empirical q value is closest to 0.05. Note that the permutation p-values are all greater than the original p-value for some pathways after 10,000 times permutation.
permutation_HHG_result2 = subset(permutation_HHG_result1,permutation.empirical.q!=0) 
empirical_p = permutation_HHG_result2[which.min(abs(permutation_HHG_result2$permutation.empirical.q - 0.05)),"permutation.empirical.p"]
empirical_p = as.numeric(empirical_p)

permutation_HHG_result1$threshold.p = ""
for (m in 1:nrow(permutation_HHG_result1)) {
  cycles_number = as.numeric(permutation_HHG_result1$permutation.cycles.number[m])
  order_p = ceiling(empirical_p*cycles_number)
  
  order_permutation_p = ALLpath_permutation_p[[permutation_HHG_result1$pathway[m]]][[1]]
  order_permutation_p1 = order_permutation_p[order(order_permutation_p)]
  permutation_HHG_result1$threshold.p[m] = order_permutation_p1[order_p]
}
sig_res = subset(permutation_HHG_result1,HHG.pvalue<=threshold.p)

#Save the result to HHG.and.permutation_result.filename (set according to your preference)
write.table(permutation_HHG_result1,file = HHG.and.permutation_result.filename ,sep="\t",row.names = FALSE,quote = FALSE)
