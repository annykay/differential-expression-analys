directory <- "D:/Lab/htseq/"
sampleFiles <- grep("counts_htseq", list.files(directory), value = TRUE)
sampleCondition <- c("Uninvolved", "Uninvolved", "Uninvolved", "Primary", "Primary", "Primary")
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design =~condition) 
keep <- rowSums(counts(ddsHTSeq)) >=10
ddsHTSeq <- ddsHTSeq[keep,]
rs <- rowSums( counts(ddsHTSeq) == 0)
idx <-  which.max(colSums(counts(ddsHTSeq) == 0))
ddsHTSeq <- ddsHTSeq[, -idx]
ddsHTSeq <- DESeq(ddsHTSeq)
res <- results(ddsHTSeq)
save(res, file = "htseq_results_de.RData")
