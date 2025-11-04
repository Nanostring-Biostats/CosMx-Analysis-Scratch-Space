


HieraType::markerslist_cd4tminor




data("markerslist_cd4tminor", package="HieraType")


all_markerslists <- 
  grep("markerslist", data(package="HieraType")$results[,"Item"], value=TRUE)

mloutl <- list()
for(mlname in all_markerslists){
  ml <- get(mlname, asNamespace("HieraType"))
  predictors <- lapply(ml, "[[", "predictors")
  index_markers <- lapply(ml, "[[", "index_marker")
  index_marker_txt <- lapply(names(index_markers), function(ct){
    paste0(ct, " = c('",paste0(index_markers[[ct]], collapse="','") , "')")
  })
  index_marker_out <- paste0("index_marker = list(\n"
                             ,paste(index_marker_txt, collapse=",\n")
                             ,"\n)"
                             )
  #cat(index_marker_out)
  
  predictor_txt <- lapply(names(predictors), function(ct){
    paste0(ct, " = c('",paste0(predictors[[ct]], collapse="','") , "')")
  })
  predictor_out <- paste0("predictors = list(\n"
                             ,paste(predictor_txt, collapse=",\n")
                             ,"\n)"
                             )
  #cat(predictor_out)
  
  mlout <- 
  paste0(mlname, " <- HieraType::make_markerslist("
         ,"\n", index_marker_out
         ,"\n,", predictor_out
         ,"\n)"
         ,"\n\nusethis::use_data(",mlname, ", overwrite = TRUE)"
         )
  cat(mlout)
  mloutl[[mlname]] <- mlout
#  browser()  
}

mloutall <- paste0(mloutl, collapse="\n\n")
cat(mloutall)
mlname <- "markerslist_cd4tminor"
cat(mloutl[["gbm_markerslist"]])



lapply(names(ml), function(xx){
  
  
  paste0("HieraType::make_markerslist(")
   
  index_marker = list(
    , function(xx)
  )   
  ("\n)\n")
  })





#
setdiff(names(rna_index_marker), names(rna_predictors))
setdiff(names(rna_index_marker), names(protein_predictors))
multi_index <- lapply(names(rna_predictors), function(xx){
  c(rna_index_marker[[xx]], protein_index_marker[[xx]])
})
multi_pred <- 
lapply(names(rna_predictors), function(xx){
  c(rna_predictors[[xx]], protein_predictors[[xx]])
})
names(multi_index) <- names(multi_pred) <- names(rna_predictors)
gbm_markerlist <- 
HieraType::make_markerslist(
index_marker =   multi_index
,predictors =  multi_pred 
)
lapply(names(rna_predictors), function(xx){
  c(rna_predictors[[xx]], protein_predictors[[xx]])
})
)






rna_index_marker = list(
  # Malignant
  Tumor_cells = c("NES", "VIM", "TOP2A", "MKI67", "SOX2", "CD44"),
  Tumor_OPC_like = c("OLIG2", "SOX10", "PDGFRA", "DLL3", "NKX2-2"),
  Tumor_Astrocyte_like = c("GFAP", "S100B", "ALDH1L1", "AQP4", "SLC1A3"),
  Tumor_Mesenchymal_like = c("CHI3L1", "CD44", "SERPINE1", "COL1A1", "VIM"),
  Tumor_Hypoxic = c("HIF1A", "CA9", "VEGFA", "ENO1", "LDHA"),
  Tumor_Cycling = c("MKI67", "TOP2A", "PCNA", "CENPF", "ASPM"),
  
  # Neural/Glial lineage
  Neuron = c("RBFOX3", "MAP2", "TUBB3", "SYT1", "SYN1"),
  Astrocyte = c("GFAP", "AQP4", "ALDH1L1", "SLC1A3", "S100B"),
  Oligodendrocyte_OPC = c("OLIG1", "OLIG2", "SOX10", "MBP", "MOG", "PDGFRA"),
  
  # Immune
  Microglia = c("CX3CR1", "P2RY12", "TMEM119", "TREM2"),
  Macrophage = c("CD68", "CD163", "CSF1R", "MRC1"),
  Monocyte = c("LYZ", "S100A8", "S100A9", "CCR2", "VCAN"),
  Dendritic = c("ITGAX", "CD1C", "CLEC9A", "LAMP3"),
  Tcell_CD4 = c("CD3D", "CD3E", "CD4", "IL7R"),
  Tcell_CD8 = c("CD3D", "CD3E", "CD8A", "GZMB"),
  Treg = c("CD3D", "CD4", "FOXP3", "IL2RA"),
  NK = c("NCAM1", "NKG7", "GNLY", "PRF1"),
  Mast = c("TPSAB1", "CPA3", "KIT", "MS4A2"),
  
  # Vascular/Stromal
  Endothelial = c("PECAM1", "VWF", "CDH5", "KDR"),
  Endothelial_BBB = c("CLDN5", "SLC2A1", "TJP1", "OCLN"),
  Endothelial_Tip = c("DLL4", "KDR", "ESM1", "ANGPT2"),
  Pericyte = c("PDGFRB", "CSPG4", "RGS5", "MCAM"),
  Fibroblast = c("PDGFRA", "FAP", "COL1A1", "DCN"),
  Ependymal = c("FOXJ1", "DNAH11", "S100A10")
)

rna_predictors = list(
  # Tumor - very broad predictive sets (100+ where appropriate)
  Tumor_cells = c(
    "NES","VIM","MKI67","TOP2A","PCNA","SOX2","PROM1","CD44","FABP7","HES1","HIF1A",
    "CENPF","ASPM","PTTG1","TROAP","MXD3","CCNB1","CCNB2","CCNE1","CDK1","CDK2","MYC",
    "MYCN","LGR5","MSI1","POU3F2","OLIG2","SOX9","PROM1","ID1","ID2","BMI1","EZH2",
    "KDM1A","HDAC1","HDAC2","BRD4","DNMT1","KIF11","AURKA","AURKB","BUB1","BUB1B",
    "MCM2","MCM3","MCM4","MCM5","MCM6","RFC4","RFC5","RPA2","RRM2","TOP2A","TOP2B",
    "PCNA","CDKN2A","CDKN1A","RB1","E2F1","E2F2","E2F3","TP53","MDM2","MDM4","AKT1",
    "MTOR","PIK3CA","PTEN","EGFR","MET","FGFR1","FGFR3","PDGFRA","ERBB2","STAT3","STAT1",
    "NFIB","NFKB1","RELA","RELB","BCL2","BCL2A1","MCL1","SNAI1","SNAI2","ZEB1","ZEB2",
    "TWIST1","TWIST2","VIM","FN1","TNC","SPP1","MMP2","MMP9","MMP14","SERPINE1","LOX",
    "LOXL2","COL1A1","COL3A1","COL4A1","FN1","SLC2A1","PGK1","ALDOA","ENO1","LDHA",
    "CA9","BNIP3","ATF4","DDIT3","HSPA1A","HSPA5","XBP1","ERN1","SQSTM1","NFE2L2",
    "HMOX1","G6PD","PPP1R3C","SLC7A5","SLC1A5","SLC3A2","SLC16A3","SLC2A1",
    "CXCL12","CCL2","CXCL8","CXCL1","CXCL2","IL6","IL8","VEGFA","ANGPT2","ADM","PDGFB",
    "PDGFRB","JUN","JUNB","FOS","FOSL1","SOX2","OLIG2","ASCL1","NEUROD1","POU5F1"
  ),
  
  Tumor_OPC_like = c(
    "OLIG2","OLIG1","SOX10","PDGFRA","DLL3","NKX2-2","CSPG4","PDGFB","PTPRZ1","ASCL1",
    "SOX9","NKX6-2","GPR17","ENPP6","QKI","SOX4","PDGFRA","EGFR","FGFR1","FGFR3",
    "PROM1","SOX11","HES6","TCF4","NKX2-2","MBP","MOG","MOBP","MAG","PLP1","CNP",
    "MYRF","TF","SLC1A2","SLC1A3","VIM","CD44","ID4","TNC","FN1","NOTCH1","NOTCH2",
    "HES1","HES5","DLL1","DLL3","JAG1","EGFR","MET","BCL2","BCL2A1","BMI1","EZH2",
    "ASPM","MKI67","TOP2A","PCNA","CCNB1","CCNB2","CDK1","CDK2","MCM2","MCM3",
    "MCM4","MCM5","ATR","CHEK1","CHEK2","AURKA","AURKB","KIF11","KIF2C","PRC1",
    "NCAPG","NCAPH"
  ),
  
  Tumor_Astrocyte_like = c(
    "GFAP","S100B","ALDH1L1","AQP4","SLC1A3","SLC1A2","SPARCL1","CLU","HES5","FGFR3",
    "SOX9","SERPINA3","C3","MT2A","MT1E","GJA1","TNC","PCDH9","GLUL","GLS","GLUD1",
    "SLC7A10","SLC1A2","SLC4A4","SLC38A1","SLC6A11","VIM","CD44","FN1","SPP1","CHI3L1",
    "MMP2","MMP9","MMP14","HIF1A","CA9","LDHA","SLC2A1","HIF1A","JUN","FOS","ATF3",
    "ATF4","DDIT3","HSPA1A","HSPA5","SOX2","PROM1","NES","OLIG2","ASCL1","NOTCH1",
    "HES1","STAT3","STAT1","IL6","LIF","OSMR","JAK1","JAK2","SOCS3","PDGFRA","PDGFB",
    "ADM","ANGPTL4","VEGFA","CXCL1","CXCL2","CXCL8"
  ),
  
  Tumor_Mesenchymal_like = c(
    "CHI3L1","CD44","SERPINE1","COL1A1","VIM","FN1","MMP9","MMP2","MMP14","TNC","LOX",
    "LOXL2","S100A4","ANXA1","ANXA2","S100A10","SPP1","IGFBP3","FAP","PDGFRB","PDGFRA",
    "PDGFB","PDGFC","TGFB1","TGFB2","ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","TWIST2",
    "NFKB1","RELA","RELB","IL1B","IL6","CXCL8","CXCL1","CXCL2","CCL2","CCL3","CCL4",
    "CCL5","CXCL12","CCR2","ITGAM","ITGB2","ITGB3","ITGA5","ITGB5","COL3A1","COL5A1",
    "COL5A2","TIMP1","SERPINA3","SERPINE2","HIF1A","VEGFA","ANGPT2","ADM"
  ),
  
  Tumor_Hypoxic = c(
    "HIF1A","CA9","VEGFA","ENO1","LDHA","SLC2A1","PGK1","ALDOA","P4HA1","PLOD2","BNIP3",
    "BNIP3L","NDRG1","SLC16A3","SLC7A5","SLC1A5","SLC3A2","SLC2A3","HMOX1","NFE2L2",
    "HSPA1A","HSPA5","XBP1","ATF4","DDIT3","ERRFI1","ADM","ANGPTL4","ANGPT2","VEGFB",
    "VEGFC","VEGFD","IL6","IL8","CXCL1","CXCL2","SOD2","NOS2","NOS3","LDHB"
  ),
  
  Tumor_Cycling = c(
    "MKI67","TOP2A","PCNA","CENPF","ASPM","PTTG1","TROAP","PLK1","AURKA","AURKB",
    "BUB1","BUB1B","CCNB1","CCNB2","CCNE1","CDC20","CDC25A","CDC25B","CDC25C","CDK1",
    "CDK2","CDK4","CDK6","E2F1","E2F2","E2F3","RRM2","RRM1","RPA1","RPA2","MCM2",
    "MCM3","MCM4","MCM5","MCM6","MCM7","ORC1","ORC2","ORC3","ORC4","GMNN","PRC1",
    "KIF11","KIF2C","KIF4A","KIF20A","NUF2","NDC80","NUSAP1","SPC24","SPC25"
  ),
  
  # Neural / glia (broad)
  Neuron = c(
    "RBFOX3","MAP2","TUBB3","SYT1","SYN1","SNAP25","NEFL","NEFM","NEFH","VGF",
    "NRGN","SLC17A7","SLC17A6","SLC6A1","GAD1","GAD2","GRIN1","GRIA1","GRIA2","CAMK2A",
    "CAMK2B","SCN1A","SCN2A","SCN3A","PTPRN","PRKCG","DCX","RELN","BSN","PCLO","RAB3A"
  ),
  
  Astrocyte = c(
    "GFAP","AQP4","ALDH1L1","SLC1A3","S100B","SLC1A2","GJA1","GLUL","SLC4A4","SLC38A1",
    "SOX9","FGFR3","KCNJ10","SLC6A11","GLS","GLUD1","SPARCL1","CLU","SERPINA3","C3",
    "APOE","MT2A","MT1E","VIM","HSPB1","HSPA1A","HSP90AA1","IL6","LIF","OSMR","STAT3"
  ),
  
  Oligodendrocyte_OPC = c(
    "OLIG1","OLIG2","SOX10","PDGFRA","MBP","MOG","MOBP","MAG","PLP1","CNP","CLDN11",
    "NKX2-2","GPR17","ENPP6","QKI","MYRF","TF","SOX8","ERMN","OPALIN","UGT8","FA2H"
  ),
  
  # Immune
  Microglia = c(
    "CX3CR1","P2RY12","TMEM119","TREM2","TYROBP","CSF1R","AIF1","SPI1","C1QA","C1QB",
    "C1QC","APOE","CTSB","CTSD","LYZ","FCGR1A","FCGR1B","FCGR3A","ITGAM","ITGAX",
    "CD68","HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","LILRB4","TYROBP","GPR34","PLP2",
    "SALL1","SALL3","SLC2A5","SLC7A7"
  ),
  
  Macrophage = c(
    "CD68","CD163","MRC1","CSF1R","IL1B","LYZ","S100A8","S100A9","VCAN","CCR2","CCL2",
    "CCL3","CCL4","SPP1","MARCO","MSR1","FCGR3A","FCGR1A","ITGAX","ITGAM","HLA-DRA",
    "HLA-DRB1","CD14","CD74","APOE","TYROBP","CTSD","CTSB","CYBB","NFKB1","RELA","RELB"
  ),
  
  Monocyte = c(
    "LYZ","S100A8","S100A9","CCR2","VCAN","FCN1","CD14","IL1B","IL6","PLBD1","SERPINB2",
    "VCAN","S100A12","CLEC4E","G0S2","FCGR3A","CD93","CD163","SPP1","CXCL2"
  ),
  
  Dendritic = c(
    "ITGAX","CD1C","CLEC9A","LAMP3","BATF3","IRF8","FSCN1","BIRC3","CD83","TRAF1",
    "CCR7","CST3","CD80","CD86","HLA-DRA","HLA-DRB1","XCR1","CLEC10A","PLD4"
  ),
  
  Tcell_CD4 = c(
    "CD3D","CD3E","CD4","IL7R","CCR7","SELL","TCF7","LEF1","CD27","CD28","ICOS","BCL6",
    "GATA3","IFNG","IL4","IL5","IL13","IL2","IL21","STAT3","STAT5A","RUNX1","RUNX3",
    "FOXP1","FOXP3","PDCD1","CTLA4","LAG3","TIGIT","TOX","TOX2","TOX3","CXCL13","CXCR5",
    "CCR4","CCR6","TNFRSF4","TNFRSF9","IL2RA","TNFRSF18"
  ),
  
  Tcell_CD8 = c(
    "CD3D","CD3E","CD8A","CD8B","GZMB","GZMA","PRF1","NKG7","IFNG","TNF","EOMES","TBX21",
    "RUNX3","KLRG1","KLRD1","KLRF1","KLRB1","ZBTB16","TOX","PDCD1","LAG3","TIGIT","HAVCR2",
    "CXCR3","CCR5","GZMK","IFITM1","CX3CR1","CD27","CD28","CD244"
  ),
  
  Treg = c(
    "CD3D","CD4","FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18","TIGIT","CCR8","IL10",
    "ITGAE","GITR","HELIO","FOXP1","BATF","IL1R2","TNFRSF4","TNFRSF9"
  ),
  
  NK = c(
    "NCAM1","NKG7","GNLY","PRF1","KLRD1","KLRF1","KLRB1","KLRG1","GZMB","GZMK","GZMA",
    "IFNG","XCL1","XCL2","FCGR3A","TYROBP","ZBTB16","EOMES","TBX21"
  ),
  
  Mast = c(
    "TPSAB1","TPSB2","CPA3","KIT","MS4A2","HPGDS","IL1RL1","CMKLR1","GATA2","GATA3",
    "CPA3","HDC","SRGN","FCER1A"
  ),
  
  # Vascular / stromal
  Endothelial = c(
    "PECAM1","VWF","CDH5","KDR","FLT1","ENG","CLDN5","SLC2A1","TIE1","TEK","NOS3",
    "EDN1","SELE","ICAM1","VCAM1","APLN","APLN","RAMP2","ESAM","PROX1","PLVAP","NRP1",
    "ADAMTS1","MMP14","ANGPT1","ANGPT2","FOXC2","SOX17","KLF2","KLF4"
  ),
  
  Endothelial_BBB = c(
    "CLDN5","SLC2A1","TJP1","OCLN","ABCB1","MFSD2A","SLC7A5","SLC16A1","SLC7A1","SLC38A5",
    "CAV1","GJA1","RAP1A","RAP1B","SYNJ2BP"
  ),
  
  Endothelial_Tip = c(
    "DLL4","KDR","ESM1","ANGPT2","FLT4","PLXND1","NRP1","ADAMTS1","MMP9","MMP14","CXCL12",
    "PDGFB","PDGFC","SOX17","APLN"
  ),
  
  Pericyte = c(
    "PDGFRB","CSPG4","RGS5","MCAM","ACTA2","TAGLN","COL1A1","COL3A1","COL4A1","LUM",
    "POSTN","PDPN","GLI1","NOTCH3","KCNJ8","ABCC9","DES","MYH11","CALD1","MYLK","THY1"
  ),
  
  Fibroblast = c(
    "PDGFRA","FAP","COL1A1","DCN","COL3A1","COL5A1","COL5A2","COL6A3","PDPN","FN1","LOX",
    "LOXL1","LOXL2","POSTN","MMP2","MMP9","MMP14","TIMP1","SERPINE1","S100A4","PDGFB",
    "TGFBI","TGFB1","TGFB2","TNC","VCAN"
  ),
  
  Ependymal = c(
    "FOXJ1","DNAH11","S100A10","ARL13B","TPPP3","TUBA1A","CFAP44","TEKT1","RSPH1","C9orf116"
  )
)




protein_index_marker = list(
  # Malignant
  Tumor_cells = c("Nestin_protein","Vimentin_protein","pS6_protein","MKI67_protein","SOX2_protein","CD44_protein"),
  Tumor_OPC_like = c("OLIG2_protein","SOX10_protein","PDGFRA_protein","DLL3_protein","CSPG4_protein"),
  Tumor_Astrocyte_like = c("GFAP_protein","Channel-GFAP_protein","S100B_protein","Aldh1l1_protein","AQP4_protein","SLC1A3_protein"),
  Tumor_Mesenchymal_like = c("CHI3L1_protein","CD44_protein","SERPINE1_protein","FN1_protein","VIM_protein"),
  Tumor_Hypoxic = c("HIF1A_protein","CA9_protein","VEGFA_protein","LDHA_protein","GLUT1_protein"),
  Tumor_Cycling = c("MKI67_protein","pHistoneH3_protein","TOP2A_protein","PLK1_protein","ASPM_protein"),
  
  # Neural/Glial lineage
  Neuron = c("NeuN_protein","Channel-NeuN_protein","MAP2_protein","Beta III Tubulin_protein","SNAP25_protein","Synaptophysin_protein"),
  Astrocyte = c("GFAP_protein","Channel-GFAP_protein","Aldh1l1_protein","S100B_protein","EAAT1-GLAST_protein","S100A10_protein"),
  Oligodendrocyte_OPC = c("SOX10_protein","MBP_protein","MOG_protein","PLP1_protein","PDGFRA_protein"),
  
  # Immune
  Microglia = c("P2ry12_protein","TMEM119_protein","Iba1_protein","Channel-Iba1_protein","TREM2_protein"),
  Macrophage = c("CD68_protein","CD11b_protein","CD163_protein","CSF1R_protein","MRC1_protein"),
  Monocyte = c("LYZ_protein","S100A8_protein","S100A9_protein","CCR2_protein","CD14_protein"),
  Dendritic = c("CD11c_protein","MHC II_protein","CD1c_protein","LAMP3_protein"),
  Tcell_CD4 = c("CD3_protein","CD4_protein","IL7R_protein","CD45_protein"),
  Tcell_CD8 = c("CD3_protein","CD8_protein","GZMB_protein","PRF1_protein"),
  Treg = c("CD3_protein","CD4_protein","FOXP3_protein","CD25_protein"),
  NK = c("NKG2A_protein","NKG2D_protein","NKG7_protein","GranzymeB_protein"),
  Mast = c("Tryptase_protein","CPA3_protein","KIT_protein","FcER1_protein"),
  
  # Vascular/Stromal
  Endothelial = c("CD31_protein","VWF_protein","VE-cadherin_protein","KDR_protein"),
  Endothelial_BBB = c("CLDN5_protein","SLC2A1_protein","ZO1_protein","OCLN_protein"),
  Endothelial_Tip = c("DLL4_protein","KDR_protein","ESM1_protein","ANGPT2_protein"),
  Pericyte = c("PDGFRB_protein","CSPG4_protein","RGS5_protein","MCAM_protein"),
  Fibroblast = c("PDGFRA_protein","FAP_protein","Collagen1_protein","Decorin_protein"),
  Ependymal = c("FOXJ1_protein","TPPP3_protein","S100A10_protein")
)

protein_predictors = list(
  Tumor_cells = c(
    "Nestin_protein","Vimentin_protein","MKI67_protein","TOP2A_protein","PCNA_protein",
    "SOX2_protein","PROM1_protein","CD44_protein","FABP7_protein","HES1_protein",
    "HIF1A_protein","CENPF_protein","ASPM_protein","PTTG1_protein","TROAP_protein",
    "MXD3_protein","cJun_protein","gamma-H2AX_protein","p53_protein",
    "pAKT_protein","pERK_protein","pS6_protein","Channel-S6_protein"
  ),
  
  Tumor_OPC_like = c(
    "OLIG2_protein","SOX10_protein","PDGFRA_protein","DLL3_protein","CSPG4_protein",
    "QKI_protein","SOX9_protein","ASCL1_protein","MBP_protein","MOG_protein"
  ),
  
  Tumor_Astrocyte_like = c(
    "GFAP_protein","Channel-GFAP_protein","S100B_protein","Aldh1l1_protein",
    "EAAT1-GLAST_protein","GJA1_protein","GLUL_protein","APOE_protein","C3_protein"
  ),
  
  Tumor_Mesenchymal_like = c(
    "CHI3L1_protein","CD44_protein","SERPINE1_protein","Vimentin_protein",
    "FN1_protein","MMP2_protein","MMP9_protein","S100A4_protein"
  ),
  
  Tumor_Hypoxic = c(
    "HIF1A_protein","CA9_protein","VEGFA_protein","LDHA_protein","SLC2A1_protein","BNIP3_protein"
  ),
  
  Tumor_Cycling = c(
    "MKI67_protein","TOP2A_protein","PCNA_protein","ASPM_protein","PTTG1_protein","TROAP_protein"
  ),
  
  # Neural / glia proteins
  Neuron = c(
    "NeuN_protein","Channel-NeuN_protein","MAP2_protein","Beta III Tubulin_protein",
    "Synaptophysin_protein","Doublecortin_protein","Calbindin_protein","NRGN_protein"
  ),
  
  Astrocyte = c(
    "GFAP_protein","Channel-GFAP_protein","Aldh1l1_protein","S100B_protein","EAAT1-GLAST_protein",
    "S100A10_protein","APOE_protein","C3_protein"
  ),
  
  Oligodendrocyte_OPC = c(
    "SOX10_protein","MBP_protein","MOG_protein","PLP1_protein","MAG_protein"
  ),
  
  # Immune proteins
  Microglia = c(
    "P2ry12_protein","TMEM119_protein","Iba1_protein","Channel-Iba1_protein","TYROBP_protein",
    "APOE_protein","CD68_protein","CD11b_protein"
  ),
  
  Macrophage = c(
    "CD68_protein","CD11b_protein","CD163_protein","MRC1_protein","S100A8_protein","S100A9_protein"
  ),
  
  Monocyte = c(
    "LYZ_protein","S100A8_protein","S100A9_protein","CCR2_protein","CD14_protein"
  ),
  
  Dendritic = c(
    "CD11c_protein","MHC II_protein","LAMP3_protein","BATF3_protein"
  ),
  
  Tcell_CD4 = c(
    "CD3_protein","CD4_protein","IL7R_protein","CD45_protein"
  ),
  
  Tcell_CD8 = c(
    "CD3_protein","CD8_protein","GZMB_protein","PRF1_protein"
  ),
  
  Treg = c(
    "CD3_protein","CD4_protein","FOXP3_protein","CD25_protein"
  ),
  
  NK = c(
    "NKG2A_protein","NKG2D_protein","NKG7_protein","GranzymeB_protein"
  ),
  
  Mast = c(
    "Tryptase_protein","CPA3_protein","KIT_protein","FcER1_protein"
  ),
  
  # Vascular / stromal proteins
  Endothelial = c(
    "CD31_protein","VWF_protein","VE-cadherin_protein","KDR_protein","Laminin_protein"
  ),
  
  Endothelial_BBB = c(
    "CLDN5_protein","SLC2A1_protein","ZO1_protein","OCLN_protein"
  ),
  
  Endothelial_Tip = c(
    "DLL4_protein","KDR_protein","ANGPT2_protein"
  ),
  
  Pericyte = c(
    "PDGFRB_protein","CSPG4_protein","RGS5_protein","MCAM_protein"
  ),
  
  Fibroblast = c(
    "PDGFRA_protein","FAP_protein","Collagen1_protein","Decorin_protein"
  ),
  
  Ependymal = c(
    "FOXJ1_protein","TPPP3_protein","S100A10_protein"
  )
)
#
#
#
#
#
#
