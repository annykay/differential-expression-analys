GSM1401683 <- read.table(file = "second_all_feature_counts.txt", row.names = 1, col.names = c("GeneID", "Length", "GSM1401683"), skip = 1)
GSM1401684 <- read.table(file = "third_all_feature_counts.txt", row.names = 1, col.names = c("GeneID", "Length", "GSM1401684"), skip = 1)
GSM1401781 <- read.table(file = "forth_all_feature_counts.txt", row.names = 1, col.names = c("GeneID", "Length", "GSM1401781"), skip = 1)
GSM1401782 <- read.table(file = "fifth_all_feature_counts.txt", row.names = 1, col.names = c("GeneID", "Length", "GSM1401782"), skip = 1)
GSM1401783 <- read.table(file = "sixth_all_feature_counts.txt", row.names = 1, col.names = c("GeneID", "Length", "GSM1401783"), skip = 1)

condition <- factor(c("GSM1401682","GSM1401683", "GSM1401684","GSM1401781","GSM1401782","GSM1401783"))

de <- merge(GSM1401682, GSM1401683, by="row.names", all = TRUE)
rownames(de) <- de[, 1]
de <- de[, -1]
de <- de[, -3]

de <- merge(de, GSM1401684, by="row.names", all = TRUE)
rownames(de) <- de[, 1]
de <- de[, -1]
de <- de[, -4]

de <- merge(de, GSM1401781, by="row.names", all = TRUE)
rownames(de) <- de[, 1]
de <- de[, -1]
de <- de[, -5]

de <- merge(de, GSM1401782, by="row.names", all = TRUE)
rownames(de) <- de[, 1]
de <- de[, -1]
de <- de[, -6]

de <- merge(de, GSM1401783, by="row.names", all = TRUE)
rownames(de) <- de[, 1]
de <- de[, -1]
de <- de[, -7]

de <-de[,-1]

library("DESeq2")
dds <-DESeqDataSetFromMatrix(de, DataFrame(condition), ~condition)
dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized = TRUE)


dds <-DESeqDataSetFromMatrix(norm, DataFrame(condition), design = designMat)
dds <- deseq2Data[rowSums(counts(dds)) > 5, ]
dds <- DESeq(dds)
dds_results <- results(dds))
save(dds_results, "deseq2_results.RData")
