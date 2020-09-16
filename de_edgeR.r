dgList <- DGEList(counts = de, genes=rownames(de))
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion >1
keep <- which(rowSums(countCheck) >=2)
dgList <- dgList[keep, ]
dgList <- calcNormFactors(dgList, method = "TMM")
sampleType <- rep("P", ncol(dgList))
sampleType[4] <- "U"
sampleType[5] <- "U"
sampleType[6] <- "U"
sampleReplicate <- paste("S", rep(1:3, each =2), sep = "")
designMat <- model.matrix(~sampleReplicate+sampleType)
sampleType[grep("U", colnames(dgList))] <- "U"
dgList <- estimateGLMCommonDisp(dgList, design=designMat)
plotBCV(dgList)
fit <- glmFit(dgList, designMat)
ltr <- glmLRT(fit, coef = 4)
edgR_result <- topTags(ltr)
save(edgR_result, file = "fc_results_de.RData")
