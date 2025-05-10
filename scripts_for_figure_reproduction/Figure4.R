
#R
data1 <- read.table("c:/Users/wangj52/OneDrive - VUMC/eNRSA/manuscript/V2/revision/RamosvsG401-tss-count-merged-forCMH_with_cmh_pvalue.txt",sep="\t",header = T)
data1 <- na.omit(data1)
length(data1[(data1$FDR < 0.05),])
length(data1$RNAseq_FDR < 0.05),])
length(data1[(data1$FDR < 0.05) & data1$RNAseq_FDR < 0.05,])


data2 <- read.table("c:/Users/wangj52/OneDrive - VUMC/eNRSA/manuscript/V2/revision/RamosvsG401-tts-count-merged-forCMH_with_cmh_pvalue.txt",sep="\t",header = T)
data2 <- na.omit(data2)
length(data2[(data1$FDR < 0.05),])
length(data2$RNAseq_FDR < 0.05),])
length(data2[(data1$FDR < 0.05) & data2$RNAseq_FDR < 0.05,])



#deeptools
computeMatrix scale-regions -S EGFP.bigWig OmoMYC.bigWig -R TRT-TTS-down.bed -o TRT-matrix-down.mat.gz

plotProfile -m RT-matrix-down.mat.gz -out TRT-down-OmoMYCvsEGFP-profile.pdf --plotTitle "Genes with decreased readthrough in OmoMYC vs. EGFP"

computeMatrix scale-regions -S EGFP.bigWig OmoMYC.bigWig -R TRT-TTS-up.bed -o TRT-matrix-up.mat.gz

plotProfile -m RT-matrix-down.mat.gz -out TRT-up-OmoMYCvsEGFP-profile.pdf --plotTitle "Genes with increased readthrough in OmoMYC vs. EGFP"