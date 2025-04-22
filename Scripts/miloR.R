library(miloR)


SYN_sce <- as.SingleCellExperiment(SYN_merged)
SYN_milo <- Milo(SYN_sce)

traj_milo <- buildGraph(SYN_milo, k = 10, d = 30)

traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 10, d=30, refined = TRUE)
plotNhoodSizeHist(traj_milo)

traj_milo$sample_id <- SYN_merged$orig.ident   # need multiple samples per condition, I think

traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), samples="sample_id")



