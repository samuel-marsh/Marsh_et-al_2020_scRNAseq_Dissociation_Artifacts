# clusterLouvainJaccard function
    # Now deprecated in LIGER

clusterLouvainJaccard = function(object, resolution = 0.1, k.param=30, n.iter = 10, n.start = 10,
                                 print.output = F, ...)
{
  if (!require("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  temp.seurat = CreateSeuratObject(t(Reduce(rbind,object@scale.data)))
  temp.seurat@scale.data = t(Reduce(rbind,object@scale.data))
  rownames(object@H.norm)=colnames(temp.seurat@scale.data)
  temp.seurat@dr$NMF=new(Class="dim.reduction",cell.embeddings=object@H.norm,key="NMF")
  temp.seurat <- FindClusters(object = temp.seurat, reduction.type = "NMF",
                              dims.use = 1:ncol(object@H.norm),force.recalc=T,
                              save.SNN = T,resolution=resolution,k.param=k.param,
                              n.iter = n.iter, n.start = n.start, print.output = print.output, ...)
  object@clusters = temp.seurat@ident
  return(object)
}
