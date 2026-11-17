#!/usr/bin/env Rscript

# ==============================================================================
# Human_Mouse_Mphi_Clustering_Resolution_Audit_ProvenanceAnchored_v1.1.0.R
#
# Corrected resolution / classification audit anchored to the ORIGINAL
# macrophage clustering provenance.
#
# HUMAN source of truth:
#   Human_MASH_LT1to8_Myeloid_RPCA_v1.7.rds
#   reduction: integrated.rpca
#   graph:     RNA_snn
#   frozen res0.8 matches RNA_snn_res.0.8 exactly on the 6,720 retained cells.
#
# MOUSE source of truth:
#   Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds
#   reduction: mphi.integrated.rpca
#   graph:     mphi_rpca_snn
#   frozen cluster_res2 matches mphi_rpca_res_2.0 exactly.
#
# IMPORTANT:
# - Uses STORED original resolution columns wherever available.
# - Does not overwrite source RDS.
# - M1/M2 are descriptive functional axes, not clustering ground truth.
# - Foreign-lineage audit uses STRICT lineage anchors; ambient-prone hepatocyte
#   transcripts are scored separately and do not by themselves define contamination.
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260)
set.seed(20261117)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(cluster)
})

VERSION <- "v1.1.0"
ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HUMAN_SOURCE <- file.path(
  ROOT, "Human_MASH_RDS",
  "Human_MASH_LT1to8_Myeloid_RPCA_v1.7",
  "objects", "Human_MASH_LT1to8_Myeloid_RPCA_v1.7.rds"
)

HUMAN_FINAL <- file.path(
  ROOT, "Human_MASH_RDS",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9",
  "objects", "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds"
)

MOUSE_FINAL <- file.path(
  ROOT, "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_CleanB_FINAL_v4.14.5",
  "RDS", "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

OUTDIR <- file.path(
  ROOT, "CrossSpecies_MASH_Mphi",
  "Human_Mouse_Mphi_Clustering_Resolution_Audit_ProvenanceAnchored_v1.1.0"
)
TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")
LOGDIR <- file.path(OUTDIR, "logs")
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(LOGDIR, recursive=TRUE, showWarnings=FALSE)

LOG <- file.path(LOGDIR, "run.log")
cat("", file=LOG)

msg <- function(...) {
  x <- paste0(...)
  cat(x, "\n")
  cat(x, "\n", file=LOG, append=TRUE)
}
write_tsv <- function(x, f) {
  write.table(x, f, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
}
z_safe <- function(x) {
  s <- sd(x, na.rm=TRUE)
  m <- mean(x, na.rm=TRUE)
  if (!is.finite(s) || s == 0) rep(0, length(x)) else (x-m)/s
}
norm_entropy <- function(x) {
  x <- as.numeric(x); x <- x[x>0]
  if (length(x)<=1) return(0)
  p <- x/sum(x)
  -sum(p*log(p))/log(length(p))
}
ari <- function(x,y) {
  ok <- !is.na(x)&!is.na(y)
  x <- as.factor(x[ok]); y <- as.factor(y[ok])
  tab <- table(x,y)
  c2 <- function(z) z*(z-1)/2
  nij <- sum(c2(tab)); a <- rowSums(tab); b <- colSums(tab)
  ai <- sum(c2(a)); bj <- sum(c2(b)); n <- sum(tab); tot <- c2(n)
  if (tot==0) return(NA_real_)
  ex <- ai*bj/tot; mx <- 0.5*(ai+bj); den <- mx-ex
  if (!is.finite(den)||den==0) NA_real_ else (nij-ex)/den
}
safe_silhouette <- function(emb, labels, max_cells=3000L) {
  labels <- as.character(labels)
  ok <- complete.cases(emb) & !is.na(labels)
  emb <- emb[ok,,drop=FALSE]; labels <- labels[ok]
  if (nrow(emb)<20 || length(unique(labels))<2) return(NA_real_)
  spl <- split(seq_len(nrow(emb)), labels)
  per <- max(10L, floor(max_cells/length(spl)))
  idx <- unique(unlist(lapply(spl, function(ix) {
    if (length(ix)<=per) ix else sample(ix, per)
  })))
  if (length(idx)>max_cells) idx <- sample(idx,max_cells)
  sil <- cluster::silhouette(as.integer(factor(labels[idx])),
                             dist(emb[idx,,drop=FALSE]))
  mean(sil[,"sil_width"],na.rm=TRUE)
}
safe_mad_high <- function(x,k=3) {
  med <- median(x,na.rm=TRUE)
  md <- mad(x,center=med,constant=1.4826,na.rm=TRUE)
  if (!is.finite(md)||md==0) rep(FALSE,length(x)) else x>med+k*md
}
first_existing <- function(x, candidates) {
  z <- candidates[candidates %in% x]
  if (length(z)) z[1] else NA_character_
}
normalize_mouse_sample <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("STD","STD_rep1")] <- "STD"
  x[x %in% c("CDHFD","CDHFD_rep1","CDAHFD","CDAHFD_rep1")] <- "CDAHFD"
  x
}

for (f in c(HUMAN_SOURCE,HUMAN_FINAL,MOUSE_FINAL)) {
  if (!file.exists(f)) stop("Missing input: ",f)
}

msg("Loading objects...")
hsrc <- readRDS(HUMAN_SOURCE)
hfin <- readRDS(HUMAN_FINAL)
mouse <- readRDS(MOUSE_FINAL)

# ==============================================================================
# Provenance hard checks
# ==============================================================================

common_h <- intersect(colnames(hfin), colnames(hsrc))
if (length(common_h) != ncol(hfin)) {
  stop("Human final cells are not fully contained in Human v1.7 source.")
}
if (!all(c("integrated.rpca","umap.myeloid.rpca") %in% Reductions(hsrc))) {
  stop("Human source missing expected RPCA/UMAP reductions.")
}
if (!"RNA_snn" %in% names(hsrc@graphs)) {
  stop("Human source missing RNA_snn.")
}
if (!all(c("mphi.integrated.rpca","mphi.umap.rpca") %in% Reductions(mouse))) {
  stop("Mouse final missing expected mphi RPCA/UMAP reductions.")
}
if (!"mphi_rpca_snn" %in% names(mouse@graphs)) {
  stop("Mouse final missing mphi_rpca_snn.")
}

HUMAN_CURRENT <- "res0_8_cluster_frozen_v18"
MOUSE_CURRENT <- "cluster_res2"

stopifnot(HUMAN_CURRENT %in% colnames(hfin@meta.data))
stopifnot(MOUSE_CURRENT %in% colnames(mouse@meta.data))

# Exact provenance validation
h08 <- as.character(hsrc@meta.data[common_h,"RNA_snn_res.0.8"])
hfrozen <- as.character(hfin@meta.data[common_h,HUMAN_CURRENT])
m20 <- as.character(mouse@meta.data[["mphi_rpca_res_2.0"]])
mfrozen <- as.character(mouse@meta.data[[MOUSE_CURRENT]])

prov <- data.frame(
  species=c("Human","Mouse"),
  source_object=c(HUMAN_SOURCE,MOUSE_FINAL),
  source_reduction=c("integrated.rpca","mphi.integrated.rpca"),
  source_graph=c("RNA_snn","mphi_rpca_snn"),
  source_cluster_col=c("RNA_snn_res.0.8","mphi_rpca_res_2.0"),
  frozen_cluster_col=c(HUMAN_CURRENT,MOUSE_CURRENT),
  n_cells=c(length(common_h),ncol(mouse)),
  exact_label_match=c(mean(h08==hfrozen),mean(m20==mfrozen)),
  ARI=c(ari(h08,hfrozen),ari(m20,mfrozen)),
  stringsAsFactors=FALSE
)
write_tsv(prov,file.path(TABDIR,"00_provenance_hardcheck.tsv"))
if (any(prov$exact_label_match != 1) || any(prov$ARI != 1)) {
  stop("Frozen/source cluster provenance did not reproduce exactly.")
}

# ==============================================================================
# Stored resolution columns
# ==============================================================================

HUMAN_RES <- c(
  "0.4"="RNA_snn_res.0.4",
  "0.6"="RNA_snn_res.0.6",
  "0.8"="RNA_snn_res.0.8",
  "1.0"="RNA_snn_res.1",
  "1.2"="RNA_snn_res.1.2",
  "1.5"="RNA_snn_res.1.5",
  "2.0"="RNA_snn_res.2"
)

MOUSE_RES <- c(
  "0.8"="mphi_rpca_res_0.8",
  "1.0"="mphi_rpca_res_1.0",
  "1.2"="mphi_rpca_res_1.2",
  "1.5"="mphi_rpca_res_1.5",
  "2.0"="mphi_rpca_res_2.0"
)

for (cc in HUMAN_RES) {
  if (!cc %in% colnames(hsrc@meta.data)) stop("Human stored resolution missing: ",cc)
}
for (cc in MOUSE_RES) {
  if (!cc %in% colnames(mouse@meta.data)) stop("Mouse stored resolution missing: ",cc)
}

# ==============================================================================
# Functional / identity panels
# ==============================================================================

PAN_H <- list(
  Myeloid=c("LST1","TYROBP","FCER1G","CTSS","LYZ","CSF1R","CD68"),
  M1=c("IL1B","TNF","CCL3","CCL4","CXCL2"),
  M2=c("MRC1","CD163","HMOX1","MAF"),
  Resident=c("MARCO","VSIG4","FOLR2","TIMD4","CD5L","C1QA","C1QB","C1QC"),
  Recruited=c("CCR2","S100A8","S100A9","VCAN","FCN1","SELL"),
  LAM=c("TREM2","GPNMB","CD9","APOE","LPL","FABP5","SPP1"),
  IFN=c("STAT1","CXCL10","ISG15","IFIT1","IFIT2","IFIT3"),
  Resolution=c("MERTK","AXL","GAS6","MFGE8")
)

PAN_M <- list(
  Myeloid=c("Lst1","Tyrobp","Fcer1g","Ctss","Lyz2","Csf1r","Cd68"),
  M1=c("Il1b","Tnf","Ccl3","Ccl4","Cxcl2"),
  M2=c("Mrc1","Cd163","Hmox1","Maf"),
  Resident=c("Marco","Vsig4","Folr2","Timd4","Cd5l","C1qa","C1qb","C1qc"),
  Recruited=c("Ccr2","S100a8","S100a9","Vcan","Fcna","Sell"),
  LAM=c("Trem2","Gpnmb","Cd9","Apoe","Lpl","Fabp5","Spp1"),
  IFN=c("Stat1","Cxcl10","Isg15","Ifit1","Ifit2","Ifit3"),
  Resolution=c("Mertk","Axl","Gas6","Mfge8")
)

# Strict non-myeloid lineage anchors.
# Ambient-prone hepatocyte transcripts are intentionally excluded here.
FOREIGN_H <- list(
  Hepatocyte_strict=c("ASGR1","HNF4A","CPS1","CYP3A4","ALDOB","APOB"),
  Endothelial=c("PECAM1","VWF","KDR","EMCN","RAMP2","CLDN5","EGFL7","PLVAP"),
  HSC_Fibroblast=c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA"),
  Cholangiocyte=c("KRT19","KRT7","EPCAM","SOX9"),
  T_cell=c("CD3D","CD3E","TRAC"),
  NK_cell=c("GNLY","KLRD1","PRF1","NKG7"),
  B_cell=c("MS4A1","CD79A","CD79B","CD19"),
  Plasma_cell=c("MZB1","JCHAIN","SDC1","DERL3"),
  Neutrophil=c("FCGR3B","CSF3R","MPO","ELANE","CEACAM8","OLFM4"),
  Erythroid=c("HBB","HBA1","HBA2","ALAS2","AHSP"),
  Platelet=c("PPBP","PF4","GP9","ITGA2B","GP1BA")
)
AMBIENT_H <- c("ALB","APOA1","APOA2","TTR")

FOREIGN_M <- list(
  Hepatocyte_strict=c("Asgr1","Hnf4a","Cps1","Cyp3a11","Aldob","Apob"),
  Endothelial=c("Pecam1","Vwf","Kdr","Emcn","Ramp2","Cldn5","Egfl7","Plvap"),
  HSC_Fibroblast=c("Col1a1","Col1a2","Col3a1","Dcn","Lum","Pdgfra"),
  Cholangiocyte=c("Krt19","Krt7","Epcam","Sox9"),
  T_cell=c("Cd3d","Cd3e","Trac"),
  NK_cell=c("Klrd1","Prf1","Nkg7","Gzma"),
  B_cell=c("Ms4a1","Cd79a","Cd79b","Cd19"),
  Plasma_cell=c("Mzb1","Jchain","Sdc1","Derl3"),
  Neutrophil=c("Csf3r","Mpo","Elane","Camp","Ltf","Ngp","Retnlg"),
  Erythroid=c("Hbb-bs","Hbb-bt","Hba-a1","Hba-a2","Alas2","Ahsp"),
  Platelet=c("Ppbp","Pf4","Gp9","Itga2b","Gp1ba")
)
AMBIENT_M <- c("Alb","Apoa1","Apoa2","Ttr")

# cDC is related myeloid ambiguity, not automatic non-MΦ contamination.
CDC_H <- c("FCER1A","CD1C","CLEC10A","CLEC9A","XCR1")
CDC_M <- c("Fcer1a","Clec10a","Clec9a","Xcr1","Batf3","Cd209a")

score_cells <- function(obj, panels, foreign, ambient, cdc, sample_vec) {
  counts <- LayerData(obj, assay="RNA", layer="counts")
  genes <- unique(c(unlist(panels),unlist(foreign),ambient,cdc))
  gp <- intersect(genes,rownames(counts))
  mat <- as.matrix(counts[gp,,drop=FALSE])
  lib <- Matrix::colSums(counts)
  lib[lib<=0] <- 1
  logcp10k <- log1p(t(t(mat)/lib)*1e4)

  one <- function(gs) {
    g <- intersect(gs,rownames(logcp10k))
    if (!length(g)) return(list(score=rep(NA_real_,ncol(obj)),z=rep(NA_real_,ncol(obj)),det=rep(0,ncol(obj))))
    sc <- colMeans(logcp10k[g,,drop=FALSE])
    list(score=sc,z=z_safe(sc),det=colSums(mat[g,,drop=FALSE]>0))
  }

  out <- data.frame(cell=colnames(obj),sample=as.character(sample_vec),stringsAsFactors=FALSE)
  for (nm in names(panels)) {
    z <- one(panels[[nm]])
    out[[paste0("z_",nm)]] <- z$z
    out[[paste0("det_",nm)]] <- z$det
  }
  for (nm in names(foreign)) {
    z <- one(foreign[[nm]])
    out[[paste0("foreign_z_",nm)]] <- z$z
    out[[paste0("foreign_det_",nm)]] <- z$det
    out[[paste0("foreign_strict_",nm)]] <- z$det>=2 & z$z>=1.5
  }
  a <- one(ambient)
  out$ambient_hepatocyte_z <- a$z
  out$ambient_hepatocyte_detected <- a$det

  d <- one(cdc)
  out$cDC_z <- d$z
  out$cDC_detected <- d$det
  out$cDC_like_ambiguity <- d$det>=2 & d$z>=1.5

  ncount <- if ("nCount_RNA" %in% colnames(obj@meta.data)) obj$nCount_RNA else Matrix::colSums(counts)
  nfeat <- if ("nFeature_RNA" %in% colnames(obj@meta.data)) obj$nFeature_RNA else Matrix::colSums(counts>0)
  out$nCount_RNA <- as.numeric(ncount)
  out$nFeature_RNA <- as.numeric(nfeat)
  out$high_nCount <- FALSE
  out$high_nFeature <- FALSE
  for (s in unique(out$sample)) {
    ix <- which(out$sample==s)
    out$high_nCount[ix] <- safe_mad_high(out$nCount_RNA[ix])
    out$high_nFeature[ix] <- safe_mad_high(out$nFeature_RNA[ix])
  }

  strict_cols <- grep("^foreign_strict_",colnames(out),value=TRUE)
  out$n_strict_foreign <- rowSums(as.matrix(out[,strict_cols,drop=FALSE]),na.rm=TRUE)
  out$possible_nonMphi_strict <- out$n_strict_foreign>=1 & out$z_Myeloid < -0.5
  out$possible_heterotypic_doublet <- out$n_strict_foreign>=1 & out$z_Myeloid>=0 & (out$high_nCount|out$high_nFeature)
  out
}

hsample <- as.character(hfin$sample_v17)
msample_col <- first_existing(colnames(mouse@meta.data),
  c("sample_display_FIXED2","sample_for_R8plot_FIXED2","sample_for_annotation","sample","orig.ident"))
msample <- normalize_mouse_sample(mouse@meta.data[[msample_col]])

# score final retained Human cells using v1.9 counts
hscore <- score_cells(hfin,PAN_H,FOREIGN_H,AMBIENT_H,CDC_H,hsample)
mscore <- score_cells(mouse,PAN_M,FOREIGN_M,AMBIENT_M,CDC_M,msample)

# ==============================================================================
# Resolution metrics using the correct original embeddings and stored labels
# ==============================================================================

audit_species <- function(species, resolution_map, label_source_obj, final_obj,
                          final_cells, sample_vec, score_df, emb, umap, current_labels) {

  sumrows <- list()
  samprows <- list()
  idrows <- list()
  contamrows <- list()
  lablist <- list()

  for (i in seq_along(resolution_map)) {
    resname <- names(resolution_map)[i]
    col <- resolution_map[[i]]
    labels <- as.character(label_source_obj@meta.data[final_cells,col])
    lablist[[resname]] <- labels
    sizes <- table(labels)
    cs <- table(labels,sample_vec)
    ent <- apply(cs,1,norm_entropy)
    maxfrac <- apply(cs,1,function(x) max(x)/sum(x))

    sumrows[[i]] <- data.frame(
      species=species,
      resolution=as.numeric(resname),
      source_cluster_col=col,
      n_clusters=length(unique(labels)),
      median_cluster_size=median(as.numeric(sizes)),
      min_cluster_size=min(as.numeric(sizes)),
      n_clusters_lt20=sum(sizes<20),
      n_clusters_lt50=sum(sizes<50),
      ARI_vs_current=ari(labels,current_labels),
      mean_PCA_silhouette=safe_silhouette(emb,labels,if(species=="Human")2500L else 3000L),
      median_sample_entropy=median(ent,na.rm=TRUE),
      median_max_sample_fraction=median(maxfrac,na.rm=TRUE),
      max_max_sample_fraction=max(maxfrac,na.rm=TRUE),
      stringsAsFactors=FALSE
    )

    cdf <- as.data.frame(cs,stringsAsFactors=FALSE)
    names(cdf) <- c("cluster","sample","n_cells")
    tot <- aggregate(n_cells~cluster,cdf,sum)
    names(tot)[2] <- "cluster_total"
    cdf <- merge(cdf,tot,by="cluster")
    cdf$fraction_cluster <- cdf$n_cells/cdf$cluster_total
    cdf$species <- species; cdf$resolution <- as.numeric(resname)
    samprows[[i]] <- cdf

    sc <- score_df
    sc$cluster <- labels

    for (cl in sort(unique(labels))) {
      ix <- which(labels==cl)
      ids <- c(
        Resident=mean(sc$z_Resident[ix],na.rm=TRUE),
        Recruited=mean(sc$z_Recruited[ix],na.rm=TRUE),
        LAM=mean(sc$z_LAM[ix],na.rm=TRUE)
      )
      oo <- order(ids,decreasing=TRUE)
      idrows[[length(idrows)+1L]] <- data.frame(
        species=species,resolution=as.numeric(resname),cluster=cl,n_cells=length(ix),
        mean_z_Resident=ids["Resident"],mean_z_Recruited=ids["Recruited"],mean_z_LAM=ids["LAM"],
        dominant_identity=names(ids)[oo[1]],
        identity_margin_top_minus_second=ids[oo[1]]-ids[oo[2]],
        mean_z_M1=mean(sc$z_M1[ix],na.rm=TRUE),
        mean_z_M2=mean(sc$z_M2[ix],na.rm=TRUE),
        mean_z_IFN=mean(sc$z_IFN[ix],na.rm=TRUE),
        mean_z_Resolution=mean(sc$z_Resolution[ix],na.rm=TRUE),
        stringsAsFactors=FALSE
      )

      strict_cols <- grep("^foreign_strict_",colnames(sc),value=TRUE)
      ff <- vapply(strict_cols,function(cc) mean(sc[ix,cc],na.rm=TRUE),numeric(1))
      names(ff) <- sub("^foreign_strict_","",strict_cols)
      top <- names(which.max(ff))
      contamrows[[length(contamrows)+1L]] <- data.frame(
        species=species,resolution=as.numeric(resname),cluster=cl,n_cells=length(ix),
        top_strict_foreign_lineage=top,
        top_strict_foreign_fraction=max(ff,na.rm=TRUE),
        fraction_possible_nonMphi_strict=mean(sc$possible_nonMphi_strict[ix],na.rm=TRUE),
        fraction_possible_heterotypic_doublet=mean(sc$possible_heterotypic_doublet[ix],na.rm=TRUE),
        fraction_cDC_like_ambiguity=mean(sc$cDC_like_ambiguity[ix],na.rm=TRUE),
        mean_ambient_hepatocyte_z=mean(sc$ambient_hepatocyte_z[ix],na.rm=TRUE),
        stringsAsFactors=FALSE
      )
    }

    # UMAP
    d <- data.frame(UMAP_1=umap[,1],UMAP_2=umap[,2],cluster=factor(labels))
    p <- ggplot(d,aes(UMAP_1,UMAP_2,color=cluster))+
      geom_point(size=if(species=="Human")0.45 else 0.28,alpha=0.85)+
      theme_void(base_size=10)+
      guides(color=guide_legend(override.aes=list(size=2.5,alpha=1)))+
      ggtitle(paste0(species," | provenance-anchored resolution ",resname))
    ggsave(file.path(FIGDIR,paste0(species,"_UMAP_res",gsub("\\.","_",resname),"_v1.1.0.pdf")),
           p,width=9,height=7)
  }

  # adjacent ARI
  adj <- list()
  nms <- names(lablist)
  if (length(nms)>=2) {
    for (i in seq_len(length(nms)-1)) {
      adj[[i]] <- data.frame(
        species=species,
        lower_resolution=as.numeric(nms[i]),
        higher_resolution=as.numeric(nms[i+1]),
        ARI_adjacent=ari(lablist[[i]],lablist[[i+1]]),
        stringsAsFactors=FALSE
      )
    }
  }

  # current-to-resolution transitions
  trans <- list()
  for (nm in nms) {
    tb <- as.data.frame(table(current_cluster=current_labels,candidate_cluster=lablist[[nm]]),
                        stringsAsFactors=FALSE)
    tb <- tb[tb$Freq>0,,drop=FALSE]
    tot <- aggregate(Freq~current_cluster,tb,sum); names(tot)[2] <- "current_total"
    tb <- merge(tb,tot,by="current_cluster")
    tb$fraction_of_current_cluster <- tb$Freq/tb$current_total
    tb$species <- species
    tb$resolution <- as.numeric(nm)
    trans[[length(trans)+1L]] <- tb
  }

  list(
    summary=do.call(rbind,sumrows),
    adjacent=do.call(rbind,adj),
    sample=do.call(rbind,samprows),
    identity=do.call(rbind,idrows),
    contamination=do.call(rbind,contamrows),
    transitions=do.call(rbind,trans)
  )
}

# Correct Human geometry from v1.7, restricted to retained v1.9 cells.
hemb <- Embeddings(hsrc,"integrated.rpca")[common_h,1:30,drop=FALSE]
humap <- Embeddings(hsrc,"umap.myeloid.rpca")[common_h,1:2,drop=FALSE]

# Correct Mouse geometry.
memb <- Embeddings(mouse,"mphi.integrated.rpca")[,1:30,drop=FALSE]
mumap <- Embeddings(mouse,"mphi.umap.rpca")[,1:2,drop=FALSE]

ha <- audit_species(
  "Human",HUMAN_RES,hsrc,hfin,common_h,hsample,hscore,hemb,humap,hfrozen
)
ma <- audit_species(
  "Mouse",MOUSE_RES,mouse,mouse,colnames(mouse),msample,mscore,memb,mumap,mfrozen
)

write_tsv(ha$summary,file.path(TABDIR,"01_Human_correct_resolution_summary.tsv"))
write_tsv(ma$summary,file.path(TABDIR,"02_Mouse_correct_resolution_summary.tsv"))
write_tsv(ha$adjacent,file.path(TABDIR,"03_Human_correct_adjacent_ARI.tsv"))
write_tsv(ma$adjacent,file.path(TABDIR,"04_Mouse_correct_adjacent_ARI.tsv"))
write_tsv(ha$sample,file.path(TABDIR,"05_Human_correct_sample_composition.tsv"))
write_tsv(ma$sample,file.path(TABDIR,"06_Mouse_correct_sample_composition.tsv"))
write_tsv(ha$identity,file.path(TABDIR,"07_Human_correct_identity_function.tsv"))
write_tsv(ma$identity,file.path(TABDIR,"08_Mouse_correct_identity_function.tsv"))
write_tsv(ha$contamination,file.path(TABDIR,"09_Human_correct_strict_contamination.tsv"))
write_tsv(ma$contamination,file.path(TABDIR,"10_Mouse_correct_strict_contamination.tsv"))
write_tsv(ha$transitions,file.path(TABDIR,"11_Human_correct_cluster_transitions.tsv"))
write_tsv(ma$transitions,file.path(TABDIR,"12_Mouse_correct_cluster_transitions.tsv"))

# Compact current-cluster strict contamination summary
current_contam <- rbind(
  subset(ha$contamination,resolution==0.8),
  subset(ma$contamination,resolution==2.0)
)
write_tsv(current_contam,file.path(TABDIR,"13_CURRENT_strict_contamination_summary.tsv"))

# Current -> neighboring resolutions of greatest interest
human_focus <- subset(ha$transitions,resolution %in% c(0.6,1.0,1.2))
mouse_focus <- subset(ma$transitions,resolution %in% c(1.2,1.5,2.0))
write_tsv(human_focus,file.path(TABDIR,"14_Human_focus_transitions_0.6_1.0_1.2.tsv"))
write_tsv(mouse_focus,file.path(TABDIR,"15_Mouse_focus_transitions_1.2_1.5_2.0.tsv"))

msg("==============================================================================")
msg("STATUS: SUCCESS")
msg("Correct provenance:")
msg("Human: integrated.rpca -> RNA_snn -> RNA_snn_res.0.8")
msg("Mouse: mphi.integrated.rpca -> mphi_rpca_snn -> mphi_rpca_res_2.0")
msg("Review first:")
msg("  01_Human_correct_resolution_summary.tsv")
msg("  02_Mouse_correct_resolution_summary.tsv")
msg("  07_Human_correct_identity_function.tsv")
msg("  08_Mouse_correct_identity_function.tsv")
msg("  09_Human_correct_strict_contamination.tsv")
msg("  10_Mouse_correct_strict_contamination.tsv")
msg("  14_Human_focus_transitions_0.6_1.0_1.2.tsv")
msg("  15_Mouse_focus_transitions_1.2_1.5_2.0.tsv")
msg("==============================================================================")
