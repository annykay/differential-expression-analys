
dgList <- DGEList(counts = de, genes=rownames(de))
countsPerMillion <-cpm(dgList)
countCheck <- countsPerMillion>1
keep <-which(rowSums(countCheck)>=2)
dgList <- dgList[keep,]
dgList <- calcNormFactors(dgList, method = "TMM")
sampleType <- c("Uninvolved", "Uninvolved", "Primary", "primary")
designMat <- model.matrix(~sampleType)
dgList <- estimateGLMCommonDisp(dgList, design=designMat)
fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef=3)
edgeR_result <- topTags(lrt)
deGenes <- decideTestsDGE(lrt, p = 0.005)