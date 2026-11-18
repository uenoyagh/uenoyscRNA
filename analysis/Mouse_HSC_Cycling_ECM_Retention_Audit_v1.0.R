#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, width = 260)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"
AUDITED_RDS <- file.path(ROOT,"Mouse_MASH_RDS","Mouse_MASH_HSC_PrimaryClean_v7.0","objects","Mouse_MASH_HSC_QCAudited_v7.0.rds")
PRIMARY_RDS <- file.path(ROOT,"Mouse_MASH_RDS","Mouse_MASH_HSC_PrimaryClean_v7.0","objects","Mouse_MASH_HSC_PrimaryClean_v7.0.rds")
SENS_RDS <- file.path(ROOT,"Mouse_MASH_RDS","Mouse_MASH_HSC_PrimaryClean_v7.0","objects","Mouse_MASH_HSC_SensitivityClean_v7.0.rds")
OUTDIR <- file.path(ROOT,"CrossSpecies_MASH_HSC","Mouse_HSC_Cycling_ECM_Retention_Audit_v1.0")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

write_tsv <- function(x,f) write.table(x,f,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
safe_mean <- function(x){x<-x[is.finite(x)]; if(!length(x)) return(NA_real_); mean(x)}
safe_median <- function(x){x<-x[is.finite(x)]; if(!length(x)) return(NA_real_); median(x)}
direction_rule <- function(sham,tx){
  sham<-sham[is.finite(sham)]; tx<-tx[is.finite(tx)]
  if(length(sham)<2 || length(tx)<2) return("insufficient")
  if(min(tx)>max(sham)) return("Tx_above_all_Sham")
  if(max(tx)<min(sham)) return("Tx_below_all_Sham")
  "overlap"
}
get_joined_counts <- function(obj){
  lyr <- Layers(obj[["RNA"]])
  count_layers <- grep("^counts($|\\.)",lyr,value=TRUE)
  if(!length(count_layers)) stop("No RNA counts layer found.")
  if(length(count_layers)>1 || !("counts" %in% count_layers)) obj <- JoinLayers(obj, assay="RNA")
  mat <- LayerData(obj[["RNA"]], layer="counts")
  if(!inherits(mat,"Matrix")) mat <- Matrix(mat,sparse=TRUE)
  mat
}

for(f in c(AUDITED_RDS,PRIMARY_RDS,SENS_RDS)) if(!file.exists(f)) stop("Missing RDS: ",f)
message("Loading audited all-HSC object...")
all_obj <- readRDS(AUDITED_RDS)
message("Loading Primary Clean HSC...")
primary_obj <- readRDS(PRIMARY_RDS)
message("Loading Sensitivity Clean HSC...")
sens_obj <- readRDS(SENS_RDS)

stopifnot(ncol(all_obj)==13166, ncol(primary_obj)==12941, ncol(sens_obj)==13026)
md <- all_obj@meta.data
required <- c("sample_hsc_v610","HSC_state3_v612","HSC_state2_v612","qc_original_frozen_cluster","qc_primary_include","qc_sensitivity_include","qc_final_disposition")
miss <- setdiff(required,colnames(md)); if(length(miss)) stop("Missing metadata: ",paste(miss,collapse=", "))

SAMPLES <- c("Sham1","Sham20","Tx17","Tx5")
GROUP <- c(Sham1="Sham",Sham20="Sham",Tx17="Tx",Tx5="Tx")
CYCLING_GENES <- c("Mki67","Top2a","Cenpf","Ube2c","Ccnb1","Ccnb2")
ECM_GENES <- c("Col1a1","Col1a2","Col3a1","Col5a1","Fn1","Lum","Dcn","Bgn","Timp1","Serpine1")

message("Joining original RNA counts...")
counts <- get_joined_counts(all_obj)
counts <- counts[,colnames(all_obj),drop=FALSE]
cyc_present <- intersect(CYCLING_GENES,rownames(counts))
ecm_present <- intersect(ECM_GENES,rownames(counts))
if(length(cyc_present)<4) stop("Too few cycling genes present")
if(length(ecm_present)<6) stop("Too few ECM genes present")

libsize <- Matrix::colSums(counts)
small_genes <- unique(c(cyc_present,ecm_present))
small_counts <- counts[small_genes,,drop=FALSE]
norm_small <- small_counts %*% Diagonal(x=1e4/libsize)
norm_small@x <- log1p(norm_small@x)

cyc_detect_n <- Matrix::colSums(counts[cyc_present,,drop=FALSE] > 0)
cyc_score <- Matrix::colMeans(norm_small[cyc_present,,drop=FALSE])
ecm_score <- Matrix::colMeans(norm_small[ecm_present,,drop=FALSE])

cyc_lenient <- cyc_detect_n >= 1
cyc_strict <- cyc_detect_n >= 2
ECM_Q75 <- as.numeric(quantile(ecm_score,0.75,na.rm=TRUE))
ECM_Q90 <- as.numeric(quantile(ecm_score,0.90,na.rm=TRUE))
ecm_q75 <- ecm_score >= ECM_Q75
ecm_q90 <- ecm_score >= ECM_Q90

cell_audit <- data.frame(
  cell=colnames(all_obj),
  sample=as.character(md$sample_hsc_v610),
  group=GROUP[as.character(md$sample_hsc_v610)],
  cluster=as.character(md$qc_original_frozen_cluster),
  HSC_state3=as.character(md$HSC_state3_v612),
  HSC_state2=as.character(md$HSC_state2_v612),
  qc_primary_include=md$qc_primary_include,
  qc_sensitivity_include=md$qc_sensitivity_include,
  qc_final_disposition=as.character(md$qc_final_disposition),
  cycling_detected_n=as.numeric(cyc_detect_n),
  cycling_score=as.numeric(cyc_score),
  cycling_lenient=as.logical(cyc_lenient),
  cycling_strict=as.logical(cyc_strict),
  ECM_score=as.numeric(ecm_score),
  ECM_high_Q75=as.logical(ecm_q75),
  ECM_high_Q90=as.logical(ecm_q90),
  cycling_strict_ECM_Q75=as.logical(cyc_strict & ecm_q75),
  cycling_strict_ECM_Q90=as.logical(cyc_strict & ecm_q90),
  stringsAsFactors=FALSE
)

write_tsv(data.frame(
  item=c("cycling_genes_requested","cycling_genes_present","ECM_genes_requested","ECM_genes_present","ECM_score_Q75_source","ECM_score_Q90_source","cycling_lenient_definition","cycling_strict_definition"),
  value=c(paste(CYCLING_GENES,collapse=";"),paste(cyc_present,collapse=";"),paste(ECM_GENES,collapse=";"),paste(ecm_present,collapse=";"),ECM_Q75,ECM_Q90,">=1 detected cycling marker",">=2 detected cycling markers"),
  stringsAsFactors=FALSE
), file.path(OUTDIR,"01_audit_definitions_and_thresholds.tsv"))

cluster_rows <- lapply(sort(unique(cell_audit$cluster)), function(cl){
  d <- cell_audit[cell_audit$cluster==cl,,drop=FALSE]
  data.frame(
    cluster=cl,n_cells=nrow(d),
    primary_retained_n=sum(d$qc_primary_include %in% TRUE),
    sensitivity_retained_n=sum(d$qc_sensitivity_include %in% TRUE),
    cycling_lenient_n=sum(d$cycling_lenient),cycling_lenient_fraction=mean(d$cycling_lenient),
    cycling_strict_n=sum(d$cycling_strict),cycling_strict_fraction=mean(d$cycling_strict),
    mean_cycling_score=safe_mean(d$cycling_score),mean_ECM_score=safe_mean(d$ECM_score),
    ECM_high_Q75_n=sum(d$ECM_high_Q75),ECM_high_Q75_fraction=mean(d$ECM_high_Q75),
    cycling_strict_ECM_Q75_n=sum(d$cycling_strict_ECM_Q75),cycling_strict_ECM_Q75_fraction=mean(d$cycling_strict_ECM_Q75),
    cycling_strict_ECM_Q90_n=sum(d$cycling_strict_ECM_Q90),cycling_strict_ECM_Q90_fraction=mean(d$cycling_strict_ECM_Q90),
    stringsAsFactors=FALSE
  )
})
write_tsv(do.call(rbind,cluster_rows),file.path(OUTDIR,"02_cluster_Cycling_ECM_retention_audit.tsv"))

flagged <- cell_audit[cell_audit$cluster %in% c("14","15","16"),,drop=FALSE]
flagged_summary <- as.data.frame(table(cluster=flagged$cluster,sample=flagged$sample,HSC_state3=flagged$HSC_state3,cycling_strict=flagged$cycling_strict,ECM_high_Q75=flagged$ECM_high_Q75),stringsAsFactors=FALSE)
flagged_summary <- flagged_summary[flagged_summary$Freq>0,,drop=FALSE]
write_tsv(flagged_summary,file.path(OUTDIR,"03_flagged_cluster_state_Cycling_ECM_summary.tsv"))
write_tsv(flagged,file.path(OUTDIR,"04_flagged_cells_Cycling_ECM_inventory.tsv"))

dataset_cells <- list(Source_all_HSC=colnames(all_obj),PrimaryClean=colnames(primary_obj),SensitivityClean=colnames(sens_obj))
sample_rows <- list(); k <- 0L
for(ds in names(dataset_cells)){
  d0 <- cell_audit[cell_audit$cell %in% dataset_cells[[ds]],,drop=FALSE]
  for(s in SAMPLES){
    d <- d0[d0$sample==s,,drop=FALSE]; if(!nrow(d)) next
    subsets <- list(ALL_HSC=rep(TRUE,nrow(d)),ECM_activated_HSC=d$HSC_state3=="ECM-activated HSC",qHSC=d$HSC_state3=="qHSC",Contractile_HSC=d$HSC_state3=="Contractile HSC")
    for(subname in names(subsets)){
      sd <- d[subsets[[subname]],,drop=FALSE]; if(!nrow(sd)) next
      cyc <- sd[sd$cycling_strict,,drop=FALSE]
      k <- k+1L
      sample_rows[[k]] <- data.frame(
        dataset=ds,sample=s,group=GROUP[[s]],subset=subname,n_cells=nrow(sd),
        cycling_lenient_n=sum(sd$cycling_lenient),cycling_lenient_fraction=mean(sd$cycling_lenient),
        cycling_strict_n=sum(sd$cycling_strict),cycling_strict_fraction=mean(sd$cycling_strict),
        ECM_high_Q75_n=sum(sd$ECM_high_Q75),ECM_high_Q75_fraction=mean(sd$ECM_high_Q75),
        cycling_strict_ECM_Q75_n=sum(sd$cycling_strict_ECM_Q75),cycling_strict_ECM_Q75_fraction=mean(sd$cycling_strict_ECM_Q75),
        cycling_strict_ECM_Q90_n=sum(sd$cycling_strict_ECM_Q90),cycling_strict_ECM_Q90_fraction=mean(sd$cycling_strict_ECM_Q90),
        mean_ECM_score_all=safe_mean(sd$ECM_score),mean_ECM_score_cycling_strict=safe_mean(cyc$ECM_score),median_ECM_score_cycling_strict=safe_median(cyc$ECM_score),
        stringsAsFactors=FALSE
      )
    }
  }
}
sample_summary <- do.call(rbind,sample_rows)
write_tsv(sample_summary,file.path(OUTDIR,"05_dataset_sample_Cycling_ECM_summary.tsv"))

metrics <- c("cycling_strict_fraction","cycling_strict_ECM_Q75_fraction","cycling_strict_ECM_Q90_fraction","mean_ECM_score_cycling_strict")
effect_rows <- list(); k <- 0L
for(ds in unique(sample_summary$dataset)){
  for(subname in c("ALL_HSC","ECM_activated_HSC")){
    d <- sample_summary[sample_summary$dataset==ds & sample_summary$subset==subname,,drop=FALSE]
    for(metric in metrics){
      sham <- d[[metric]][d$group=="Sham"]; tx <- d[[metric]][d$group=="Tx"]
      k <- k+1L
      effect_rows[[k]] <- data.frame(dataset=ds,subset=subname,metric=metric,Sham_mean=safe_mean(sham),Tx_mean=safe_mean(tx),delta_Tx_minus_Sham=safe_mean(tx)-safe_mean(sham),direction_concordance=direction_rule(sham,tx),stringsAsFactors=FALSE)
    }
  }
}
write_tsv(do.call(rbind,effect_rows),file.path(OUTDIR,"06_group_effects_Cycling_ECM.tsv"))

# ECM pseudobulk within strict cycling HSC
pb_rows <- list(); k <- 0L
for(ds in names(dataset_cells)){
  for(s in SAMPLES){
    cells <- cell_audit$cell[cell_audit$cell %in% dataset_cells[[ds]] & cell_audit$sample==s & cell_audit$cycling_strict]
    if(!length(cells)) next
    ix <- match(cells,colnames(counts))
    sums <- Matrix::rowSums(counts[ecm_present,ix,drop=FALSE])
    lib <- sum(Matrix::colSums(counts[,ix,drop=FALSE]))
    cp <- as.numeric(sums)/lib*1e4; names(cp)<-ecm_present
    for(g in ecm_present){
      k <- k+1L
      pb_rows[[k]] <- data.frame(dataset=ds,sample=s,group=GROUP[[s]],gene=g,cycling_strict_n=length(cells),cycling_cells_pseudobulk_CP10k=cp[g],stringsAsFactors=FALSE)
    }
  }
}
write_tsv(do.call(rbind,pb_rows),file.path(OUTDIR,"07_ECM_gene_pseudobulk_within_strictCycling_HSC.tsv"))

primary_removed <- !cell_audit$qc_primary_include
sens_removed <- !cell_audit$qc_sensitivity_include
ncyc <- sum(cell_audit$cycling_strict)
nq75 <- sum(cell_audit$cycling_strict_ECM_Q75)
nq90 <- sum(cell_audit$cycling_strict_ECM_Q90)
loss <- data.frame(
  metric=c("all_HSC_cells","Primary_removed_all_cells","Sensitivity_removed_all_cells","strict_cycling_cells_source","strict_cycling_removed_Primary","strict_cycling_removed_Sensitivity","strict_cycling_ECM_Q75_source","strict_cycling_ECM_Q75_removed_Primary","strict_cycling_ECM_Q75_removed_Sensitivity","strict_cycling_ECM_Q90_source","strict_cycling_ECM_Q90_removed_Primary","strict_cycling_ECM_Q90_removed_Sensitivity"),
  n=c(nrow(cell_audit),sum(primary_removed),sum(sens_removed),ncyc,sum(cell_audit$cycling_strict & primary_removed),sum(cell_audit$cycling_strict & sens_removed),nq75,sum(cell_audit$cycling_strict_ECM_Q75 & primary_removed),sum(cell_audit$cycling_strict_ECM_Q75 & sens_removed),nq90,sum(cell_audit$cycling_strict_ECM_Q90 & primary_removed),sum(cell_audit$cycling_strict_ECM_Q90 & sens_removed)),
  denominator=c(nrow(cell_audit),nrow(cell_audit),nrow(cell_audit),nrow(cell_audit),ncyc,ncyc,nrow(cell_audit),nq75,nq75,nrow(cell_audit),nq90,nq90),
  stringsAsFactors=FALSE
)
loss$fraction <- loss$n/loss$denominator
write_tsv(loss,file.path(OUTDIR,"08_explicit_Cycling_ECM_loss_audit.tsv"))

c15 <- cell_audit[cell_audit$cluster=="15",,drop=FALSE]
c15_summary <- data.frame(
  metric=c("c15_cells","c15_strict_cycling_n","c15_strict_cycling_fraction","c15_ECM_high_Q75_n","c15_ECM_high_Q75_fraction","c15_strict_cycling_ECM_Q75_n","c15_strict_cycling_ECM_Q75_fraction","c15_strict_cycling_ECM_Q90_n","c15_strict_cycling_ECM_Q90_fraction","c15_mean_ECM_score"),
  value=c(nrow(c15),sum(c15$cycling_strict),mean(c15$cycling_strict),sum(c15$ECM_high_Q75),mean(c15$ECM_high_Q75),sum(c15$cycling_strict_ECM_Q75),mean(c15$cycling_strict_ECM_Q75),sum(c15$cycling_strict_ECM_Q90),mean(c15$cycling_strict_ECM_Q90),safe_mean(c15$ECM_score)),
  stringsAsFactors=FALSE
)
write_tsv(c15_summary,file.path(OUTDIR,"09_c15_sensitivity_restoration_audit.tsv"))

get_fraction <- function(name) loss$fraction[match(name,loss$metric)]
decision <- data.frame(
  question=c("fraction_all_source_HSC_removed_Primary","fraction_strict_cycling_HSC_removed_Primary","fraction_strict_cycling_ECM_Q75_removed_Primary","fraction_strict_cycling_ECM_Q90_removed_Primary","SensitivityClean_restores_all_c15","QC_based_on_cycling_or_high_complexity_alone"),
  value=c(get_fraction("Primary_removed_all_cells"),get_fraction("strict_cycling_removed_Primary"),get_fraction("strict_cycling_ECM_Q75_removed_Primary"),get_fraction("strict_cycling_ECM_Q90_removed_Primary"),as.character(all(c15$qc_sensitivity_include %in% TRUE)),"FALSE"),
  stringsAsFactors=FALSE
)
write_tsv(decision,file.path(OUTDIR,"10_KEY_retention_decision_table.tsv"))

checks <- data.frame(
  item=c("Audited_source_cells","PrimaryClean_cells","SensitivityClean_cells","all_four_samples_present","original_filtered_counts_used","DecontX_counts_used","reclustering_performed","cycling_based_exclusion_performed","biological_replicate"),
  value=c(ncol(all_obj),ncol(primary_obj),ncol(sens_obj),as.character(all(SAMPLES %in% unique(cell_audit$sample))),"TRUE","FALSE","FALSE","FALSE","sample"),
  stringsAsFactors=FALSE
)
write_tsv(checks,file.path(OUTDIR,"11_analysis_integrity_summary.tsv"))
writeLines(capture.output(sessionInfo()),file.path(OUTDIR,"sessionInfo.txt"),useBytes=TRUE)

message("============================================================")
message("STATUS: SUCCESS")
message("Review first:")
message("  02_cluster_Cycling_ECM_retention_audit.tsv")
message("  05_dataset_sample_Cycling_ECM_summary.tsv")
message("  06_group_effects_Cycling_ECM.tsv")
message("  08_explicit_Cycling_ECM_loss_audit.tsv")
message("  09_c15_sensitivity_restoration_audit.tsv")
message("  10_KEY_retention_decision_table.tsv")
message("  11_analysis_integrity_summary.tsv")
message("============================================================")
