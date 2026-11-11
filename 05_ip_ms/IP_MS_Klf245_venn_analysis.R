

library(readxl)


setwd("E:/Klf5_paper_submission/NC_revision/Analysis part/K2,4,5 IP-MASS")



klf5_tsc=read.csv("iBAQ_TSC_padj0.1_logFC3_results.csv")

klf5_esc=read.csv("iBAQ_ESC_padj0.1_logFC3_results.csv")
klf2_esc=as.data.frame(read_excel("P&P_2021_Klf2_IP-MASS.xlsx"))
klf4_esc=as.data.frame(read_excel("NCB_2019_Klf4_IP-MASS.xlsx"))


klf5_esc_sig=klf5_esc[klf5_esc$P.Value<0.01,]

klf2_esc_sig=klf2_esc[klf2_esc$`Fold change (FLAGKlf2/control)`>8,]
klf2_esc_sig=klf2_esc_sig[klf2_esc_sig$`p-value by t-test`<0.01,]

klf4_esc_sig=klf4_esc[!is.na(klf4_esc$Filter1_noIGG),]
klf4_esc_sig=klf4_esc_sig[klf4_esc_sig$`Signal/Ctr`>8,]



Venn_