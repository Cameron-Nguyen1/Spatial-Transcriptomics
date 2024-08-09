deps = c("Seurat","SeuratData","ggplot2","patchwork","dplyr","SeuratDisk","spacexr","grid","gridExtra")
lapply(deps,library,character.only=TRUE)

cluster.cols <- pals::polychrome(n = 36) #Supports up to 36 cluster colors
names(cluster.cols) = as.character(0:35)
cluster.cols["0"] = "#9529DF"
cluster.cols["1"] = "#43CD80"
cluster.cols["4"] = "#97FFFF"
cluster.cols["8"] = "#9AFF9A"
cluster.cols["9"] = "#CDC1C5"

Spatial_Prediction_Plots = function(ret_obj,name_vector,ncols,name,width,height){
    mFP = SpatialFeaturePlot(ret_obj, features = name_vector, pt.size.factor = 2.5, ncol = ncols, crop = TRUE)
    for (i in 1:length(mFP)){
        mFP[[i]]$theme$legend.key.width = unit(1.15,"cm")
    }
    ggsave(plot=mFP,filename=name,width=width,height=height)
}
bdim = function(x){
  DefaultAssay(x) = x[["pca"]]@assay.used
  pct=x[['pca']]@stdev/sum(x[['pca']]@stdev)*100 #Make PCT var for PCA dim selection
  csum = cumsum(pct)
  c1 = which(pct < 5 & csum > 90)[1] #PCA dim must contribute <5 variance and > 90 cumulative variance
  c2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1 #Find dim with less than .1 difference in step
  if (is.na(c1) | is.na(c2)){
    z = c(c1,c2)
    bestdims = z[!is.na(z)]
  }else{
  bestdims = min(c1,c2)
  }
  return(bestdims)
}

buildRCTD = function(spatial_obj,reference_obj,annotation_key,slide,threads){
#Build Reference
    rcounts = reference_obj[["RNA"]]$counts
    rnUMI = reference_obj$nCount_RNA
    z = as.factor(reference_obj[[]][annotation_key][[1]])
    names(z) = rownames(reference_obj[[]][annotation_key])
    annots = z  %>% forcats::fct_lump_min(min=25)
    reference = Reference(rcounts,annots,nUMI=rnUMI)
#Build Query
    count_idx = match(slide,names(spatial_obj@images))
    coords = GetTissueCoordinates(spatial_obj,image=slide) %>% rename(c(x="imagerow",y="imagecol"))
    counts = LayerData(spatial_obj, assay = "Spatial", layer = paste0("counts.",count_idx))
    query = SpatialRNA(coords,counts,nUMI=colSums(counts))
#Assemble RCTD
    RCTD <- create.RCTD(query, reference, max_cores = threads)
    return(RCTD)
}
handleit = function(x,pdfname="Mypdf.pdf",width=20,height=20,paired=FALSE){
    if (paired==FALSE){
        divz=length(x)/4
        if (length(x) %% 4 > 0){row_num = floor(divz)+1}else{row_num = divz}
        g=x+plot_layout(ncol=4,nrow=row_num)
        ggsave(plot=g,file=pdfname,width=width,height=height)
        return("Success")
    }else{return(x)}
}

spatialDEG = function(x,resolution){
    #x = PrepSCTFindMarkers(x)
    clusters = sort(unique(x[[resolution]])[,1])
    for (id1 in clusters){
        if (id1 == (length(clusters)-1)){break
        }else{
            for (id2 in (as.integer(id1)+1):(length(clusters)-1)){
                z = FindMarkers(x,ident.1=id1,ident.2=id2,logfc.threshold=0.5)
                write.csv(z,file=paste0("SpatialDEG_LFC.0.5_Cluster",id1,"_","Cluster",id2,".csv"))
            }
            z2 = FindMarkers(x,ident.1=id1,ident.2=NULL,logfc.threshold=0.5)
            write.csv(z2,file=paste0("SpatialDEG_LFC.0.5_Cluster",id1,"_VsOthers.csv"))
        }
    }
}
highlightByCell = function(x){
    for (cluster in unique(x$SCT_snn_res.0.3)){
        a = SpatialDimPlot(x, cells.highlight = CellsByIdentities(x,idents=x$SCT_snn_res.0.3), cols.highlight=c(cluster.cols[as.integer(cluster)+1],"grey60"))
            #scale_fill_manual(name="Cluster",breaks=names(cluster.cols),values=c(cluster.cols))
        handleit(a,width=20,height=16,pdfname=paste0("Clustering_Highlight_0.3_ClusterNumber",cluster,".pdf"))
    }
}

########## BEGIN SCRIPT ################

slide_vec1 = list.files(pattern="_out$")
mylo = list()
mylo = lapply(X=slide_vec1,FUN=function(x,listo=mylo){
    mylo[x] = Load10X_Spatial(x,slice=x)
    mylo[[x]]$orig.ident = x
    mylo[[x]] = mylo[[x]][,unname(which(colSums(GetAssayData(mylo[[x]]))!=0))]
    mylo[[x]][['percent.mt']] = PercentageFeatureSet(mylo[[x]], pattern = '^MT-|^Mt-|^mt-')
    mylo[[x]] = subset(mylo[[x]], percent.mt < 20 & nCount_Spatial > quantile(mylo[[x]]$nCount_Spatial, 0.05) & nCount_Spatial < quantile(mylo[[x]]$nCount_Spatial, 0.98) )
    mylo[[x]] = SCTransform(mylo[[x]],assay="Spatial",verbose=FALSE)
})

for (i in 1:length(slide_vec1)){
    if (i == 1){
        ret = merge(mylo[[1]],mylo[[2]])
    }
    else if (i == length(slide_vec1)){
        break
    }else{
        ret = merge(ret,mylo[[i+1]])
    }
}
VariableFeatures(ret) = SelectIntegrationFeatures(mylo)

ret = RunPCA(ret, assay="SCT", verbose = TRUE)
m_bdim = bdim(ret)
ret = FindNeighbors(ret, reduction="pca", dims = 1:m_bdim)
ret = FindClusters(ret, verbose = TRUE,resolution=0.3)
ret = RunUMAP(ret, reduction="pca", dims = 1:m_bdim)

umap_c = DimPlot(ret, reduction = "umap", group.by = c("ident", "orig.ident"))
spatial_cluster_c = SpatialDimPlot(ret,label=TRUE,repel=TRUE,label.size=2.5)
Fplot = SpatialFeaturePlot(ret, features = c("Sftpc","Cxcr2"))

ggsave(plot=umap_c,filename="Umap_Comparison.pdf",width=20,height=10)
ggsave(plot=spatial_cluster_c,filename="Spatial_Cluster_Comparison.pdf",width=20,height=12)
ggsave(plot=Fplot,filename="FeaturePlot.pdf",width=20,height=12)

#You'll need a reference.
MouseRef = readRDS("Annotated_reference.rds")

#Annotate using Seurats transfer anchors

anchors = FindTransferAnchors(reference = MouseRef, reference.reduction="pca", query = ret, normalization.method = "SCT")
predictions.immgen.assay = TransferData(anchorset = anchors, refdata = MouseRef$immgen.pred_pruned, prediction.assay = TRUE,
    weight.reduction = ret[["pca"]], dims = 1:m_bdim)
predictions.mrsd.assay = TransferData(anchorset = anchors, refdata = MouseRef$mrsd.pred_pruned, prediction.assay = TRUE,
    weight.reduction = ret[["pca"]], dims = 1:m_bdim)
predictions.azimuth.assay = TransferData(anchorset = anchors, refdata = MouseRef$predicted.celltype_level2, prediction.assay = TRUE,
    weight.reduction = ret[["pca"]], dims = 1:m_bdim)

ret[["predictions.mrsd"]] = predictions.mrsd.assay
ret[["predictions.immgen"]] = predictions.immgen.assay
ret[["predictions.azimuth"]] = predictions.azimuth.assay

immgen_annots = sort(unique(MouseRef$immgen.pred_pruned))
mrsd_annots = sort(unique(MouseRef$mrsd.pred_pruned))
azimuth_annots = sort(unique(MouseRef$predicted.celltype_level2))

DefaultAssay(ret) <- "predictions.immgen"
Spatial_Prediction_Plots(ret,immgen_annots,6,"Immgen_Spatial_Predictions.pdf",25,25)
DefaultAssay(ret) <- "predictions.mrsd"
Spatial_Prediction_Plots(ret,mrsd_annots,6,"MRSD_Spatial_Predictions.pdf",25,25)
DefaultAssay(ret) <- "predictions.azimuth"
Spatial_Prediction_Plots(ret,azimuth_annots,10,"Azimuth_Spatial_Predictions.pdf",30,30)

#Annotate using Spatial Deconvolution

RCTDA1 = buildRCTD(ret,MouseRef,"predicted.celltype_level2","A1",16)
RCTDD1 = buildRCTD(ret,MouseRef,"predicted.celltype_level2","D1",16)

#-- The following commands are computationally intensive, you must run these to proceed --#
#RCTDA1 = runRCTD(RCTDA1,doublet_mode="doublet")
#RCTDD1 = runRCTD(RCTDD1,doublet_mode="doublet")

#saveRDS(object=RCTDA1,file="RCTDA1.rds")
#saveRDS(object=RCTDD1,file="RCTDD1.rds")

#aRCTDA1 = readRDS("RCTDA1.rds")
#aRCTDD1 = readRDS("RCTDD1.rds")


ret$A1_atype1 = aRCTDA1@results$results_df$first_type
ret$A1_atype2 = aRCTDA1@results$results_df$second_type
ret$D1_atype1 = aRCTDD1@results$results_df$first_type
ret$D1_atype2 = aRCTDD1@results$results_df$second_type

ret$A1_itype1 = iRCTDA1@results$results_df$first_type
ret$A1_itype2 = iRCTDA1@results$results_df$second_type
ret$D1_itype1 = iRCTDD1@results$results_df$first_type
ret$D1_itype2 = iRCTDD1@results$results_df$second_type

colors = generate_annot_cols(ret)

A1.1 = SpatialDimPlot(ret,group.by="A1_atype1",images="A1", cols = ,pt.size.factor = 2.2) + 
    guides(fill = guide_legend(override.aes = list(size=3), ncol=2)) +
    scale_fill_manual("Azimuth_A1_T1",values=colors[["A1_atype1"]])
A1.2 = SpatialDimPlot(ret,group.by="A1_atype2",images="A1",pt.size.factor = 2.2,crop=TRUE)+
    guides(fill = guide_legend(override.aes = list(size=3), ncol=2)) +
    scale_fill_manual("Azimuth_A1_T2",values=colors[["A1_atype2"]])
D1.1 = SpatialDimPlot(ret,group.by="D1_atype1",images="D1",pt.size.factor = 2.2,crop=TRUE)+
    guides(fill = guide_legend(override.aes = list(size=3), ncol=2)) +
    scale_fill_manual("Azimuth_D1_T1",values=colors[["D1_atype1"]])
D1.2 = SpatialDimPlot(ret,group.by="D1_atype2",images="D1",pt.size.factor = 2.2,crop=TRUE)+
    guides(fill = guide_legend(override.aes = list(size=3), ncol=2)) +
    scale_fill_manual("Azimuth_D1_T2",values=colors[["D1_atype2"]])
z = A1.1 + A1.2 + D1.1 + D1.2 + patchwork::plot_layout(ncol=2,nrow=2)
ggsave("Azimuth_Deconvolution_2types.pdf",plot=z,width=15,height=15)