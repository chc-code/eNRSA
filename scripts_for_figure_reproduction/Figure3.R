library('ggplot2')

#theme_update(text = element_text(size = 30))

counttable <- read.table("c:/Users/wangj52/OneDrive - VUMC/eNRSA/manuscript/batch-April-SWI-SNF-PROseq/count_pp_gb.txt",header=T,as.is=T,sep="\\t")

var.list <- names(counttable)

var.ppc <- var.list[grep("ppc_", var.list)]

var.gbc <- var.list[grep("gbc_", var.list)]

counttable$ppc.sum <- apply(counttable[, var.ppc], 1, sum)

counttable$gbc.sum <- apply(counttable[, var.gbc], 1, sum)

rm(var.ppc); rm(var.gbc)

counttable$SameGene.ppc.sum.max <- ave(counttable$ppc.sum, counttable$Gene, FUN = function(x) ifelse(x %in% x[which(x == max(x))], 1, 0))

counttable <- subset(counttable, SameGene.ppc.sum.max == 1)  

counttable$SameGene.gbc.sum.max <- ave(counttable$gbc.sum, counttable$Gene, FUN = function(x) ifelse(x %in% x[which(x == max(x))], 1, 0))

counttable <- subset(counttable, SameGene.gbc.sum.max == 1)  

counttable$SameGene.Order <- ave(counttable$Gene, counttable$Gene, FUN = function(x) 1:length(x))

counttable <- subset(counttable, SameGene.Order == 1, select = var.list)

rm(var.list)

pp <- counttable[,grep("ppc_",names(counttable))]

gb <- counttable[,grep("gbc_",names(counttable))]

colData <- data.frame(condition=factor(c("DMSO","DMSO","dTAG47","dTAG47"),levels = c("DMSO","dTAG47")),sample=c("DMSO_1","DMSO_2","dTAG47_1","dTAG47_2"),batch=factor(c("b1","b2","b1","b2"),levels = c("b1","b2")))

dds <- DESeqDataSetFromMatrix(gb, colData, formula(~ batch + condition))

dds <- estimateSizeFactors(dds)

sf <- sizeFactors(dds)

vsd<-varianceStabilizingTransformation(dds)

#z <- plotPCA(vsd, intgroup=c("condition"))

nudge <- position_nudge(y = 1)

#z + geom_text(aes(label = vsd$sample), position = nudge, size=2) + ggtitle("Gene body")

pcaData <- plotPCA(vsd, intgroup=c("condition", "batch"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +

  geom_point(size=6) +

  xlab(paste0("PC1: ",percentVar[1],"% variance")) +

  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 

  coord_fixed() + geom_text(aes(label = vsd$sample), position = nudge, size=4) + ggtitle("Gene body") + scale_color_hue(h = c(180, 300)) + theme(text = element_text(size = 20))

dds <- DESeq(dds)

result <- results(dds,contrast=c("condition","dTAG47","DMSO"))

write.table(cbind(counttable[,c(1,2)],result),file="c:/Users/wangj52/OneDrive - VUMC/eNRSA/manuscript/batch-April-SWI-SNF-PROseq/dTAG47vsDMSO-gb-withBatch.txt",quote=F,row.names=F,sep="\\t")

rnk <- cbind(counttable[,c(1,2)],result)

rnk <- na.omit(rnk)

write.table(rnk[,c(2,6)],file="c:/Users/wangj52/OneDrive - VUMC/eNRSA/manuscript/batch-April-SWI-SNF-PROseq/dTAG47vsDMSO-gb-withBatch.rnk",quote=F,row.names=F,sep="\\t")

###

gb1 <- read.table("c:/Users/wangj52/OneDrive - VUMC/eNRSA/manuscript/batch-April-SWI-SNF-PROseq/dTAG47vsDMSO-gb.txt",sep="\\t",header = T)

gb2 <- read.table("c:/Users/wangj52/OneDrive - VUMC/eNRSA/manuscript/batch-April-SWI-SNF-PROseq/dTAG47vsDMSO-gb-withBatch.txt",sep="\\t",header = T)

gb1 <- na.omit((gb1))

gb2 <- na.omit((gb2))

dTAG_up=gb1[gb1$padj<0.05 & gb1$log2FoldChange>0.3,1]\
dTAG_up_batch=gb2[gb2$padj<0.05 & gb2$log2FoldChange>0.3,1]\
dTAG_dn=gb1[gb1$padj<0.05 & gb1$log2FoldChange<(-0.3),1]\
dTAG_dn_batch=gb2[gb2$padj<0.05 & gb2$log2FoldChange<0(-0.3),1])

##

data1 <- read.table("c:/Users/wangj52/OneDrive - VUMC/eNRSA/manuscript/batch-April-SWI-SNF-PROseq/dTAG47vsDMSO-gb-withBatch.txt",sep="\\t",header=T, stringsAsFactors = FALSE)

data1 <- na.omit(data1)

data1$col <- "black"

data1[(data1$padj < 0.05) & (data1$log2FoldChange > 0.3),][, "col"] <- "red"

data1[(data1$padj < 0.05) & (data1$log2FoldChange < (-0.3),][, "col"] <- "blue"

data1$y <- log10(data1$padj)*(-1)

#data2$y <- log10(data2$padj)*(-1)

p1 <- ggplot(data1,aes(x=log2FoldChange,y=y,label=Gene))+geom_point(colour=data1$col)+xlim(-2, 2)+theme_bw()+labs(x = "log2(FoldChange)",y="-log10(FDR)",title="dTAG47 vs DMSO")+theme(legend.position = "none",text = element_text(size = 30))

#

p1 + geom_text(aes(label=ifelse(Gene %in% c("RPS24","RPLP2","RPL35","RPL12","RPL22L1","RPL14","RPL34","RPIA","RPF2","RPL22","RPL32","RPL13"),as.character(Gene),'')),hjust=0,vjust=0)

###

library(ggrepel)

p1 <- ggplot(data1,aes(x=log2FoldChange,y=y))+geom_point(colour=data1$col)+xlim(-2, 2)+theme_bw()+labs(x = "log2(FoldChange)",y="-log10(FDR)",title="dTAG47 vs DMSO")+theme(legend.position = "none",text = element_text(size = 30))

p1 + geom_label_repel(aes(label = ifelse(Gene %in% c("RPS24","RPLP2","RPL35","RPL12","RPL22L1","RPL14","RPL34","RPIA","RPF2","RPL22","RPL32","RPL13"),as.character(Gene),'')),box.padding   = 0.35, point.padding = 0.5,segment.color = 'grey50',max.overlaps = 200) 

##

g1 <- data1[(data1$padj < 0.05) & (data1$log2FoldChange > 0.3),2]

write.table(unique(g1),file="c:/Users/wangj52/OneDrive - VUMC/eNRSA/manuscript/batch-April-SWI-SNF-PROseq/batch_gb_up.gene.txt",sep="\\t",col.names = F,row.names = F,quote = F)

g2 <- data1[(data1$padj < 0.05) & (data1$log2FoldChange < (-0.3)),2]

write.table(unique(g2),file="c:/Users/wangj52/OneDrive - VUMC/eNRSA/manuscript/batch-April-SWI-SNF-PROseq/batch_gb_down.gene.txt",sep="\\t",col.names = F,row.names = F,quote = F)

g <- read.table("c:/Users/wangj52/OneDrive - VUMC/eNRSA/manuscript/batch-April-SWI-SNF-PROseq/batch_gb_down.ORA.txt",sep="\\t",header=T)

g$"adjFDR" <- log10(g$FDR)*(-1)

g <- g[order(g$Ratio),]

temp <- g$Description

temp <- unique(temp)

g$Description <- factor(g$Description, levels = temp)

S1 <- ggplot(g, aes(x= adjFDR, y=Description, size=Ratio)) + geom_point(alpha = 0.8,color="blue") + theme_classic()

#S1 

S1 <- S1+scale_size(range = c(2, 6)) + xlab("-log10 FDR") + ylab("")

S1 + theme(text = element_text(size = 20))

