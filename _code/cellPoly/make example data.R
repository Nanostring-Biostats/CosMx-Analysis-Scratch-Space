# run within the SMI0016 folder, under "6.992 analysis with controls"
cols <- readRDS("processed_data/cols.RDS")
clusts = readRDS( "processed_data/clusts.RDS")

## load one FOV's transcript locations:
fov = "16"
tdf = read.csv(paste0("../5.1 Final raw data/SP19_1139_R1080_S3/Run1080_Iter1_PCIter0_TSR0.5_DF1_v5/FOV", fov, "/FOV_Analysis_Summary/Run1080_FOV0", fov, "__complete_code_cell_target_call_coord.csv"))
head(tdf)
tdf = tdf[, c("CellId", "x", "y", "target")]
colnames(tdf)[1] = "cell_id"
tdf = tdf[tdf$cell_id != "0", ]
tdf$y = -tdf$y
tdf$cell_id = paste0("c_5_", fov, "_", tdf$cell_id)


sub = (tdf$x < median(tdf$x)) & (tdf$y < median(tdf$y))
cosmx_nephritis = list(transcript_df = tdf[sub, ],
                       celltype = clusts[is.element(names(clusts), tdf$cell_id)],
                       cellcols= cols)
save(cosmx_nephritis, file = "results/cosmx_nephritis.RData")


