#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,width=240)
suppressPackageStartupMessages({library(Seurat);library(SeuratObject);library(Matrix)})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"
BASE <- file.path(ROOT,"Mouse_MASH_RDS","Mouse_MASH_HSC_PrimaryClean_v7.0","objects")
SRC <- file.path(BASE,"Mouse_MASH_HSC_QCAudited_v7.0.rds")
PRI <- file.path(BASE,"Mouse_MASH_HSC_PrimaryClean_v7.0.rds")
SEN <- file.path(BASE,"Mouse_MASH_HSC_SensitivityClean_v7.0.rds")
OLD <- file.path(ROOT,"Mouse_MASH_Interaction","PDGF_mitogenic_FN1_rewiring_v6.3.5",
                  "Tables","05_ECM_HSC_cycling_like_fraction_by_sample_v6.3.5.csv")
OUT <- file.path(ROOT,"CrossSpecies_MASH_HSC","Mouse_HSC_LegacyCycling_Compatibility_Audit_v1.0")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

w <- function(x,n) write.table(x,file.path(OUT,n),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
gm <- function(x) mean(x[is.finite(x)])
dir4 <- function(sh,tx){
  z <- as.vector(outer(tx,sh,"-"))
  up <- sum(z>0); dn <- sum(z<0); ti <- sum(z==0)
  data.frame(pairwise_up_n=up,pairwise_down_n=dn,pairwise_tie_n=ti,
             pairwise_direction=if(dn>up)"Tx_down" else if(up>dn)"Tx_up" else "mixed",
             pairwise_direction_consistency=max(up,dn,ti)/4)
}
get_counts <- function(o){
  ly <- Layers(o[["RNA"]]); cl <- grep("^counts($|\\.)",ly,value=TRUE)
  if(length(cl)>1 || !("counts"%in%cl)) o <- JoinLayers(o,assay="RNA")
  LayerData(o[["RNA"]],layer="counts")
}

src <- readRDS(SRC); pri <- readRDS(PRI); sen <- readRDS(SEN)
stopifnot(ncol(src)==13166,ncol(pri)==12941,ncol(sen)==13026)

md <- src@meta.data
need <- c("sample_hsc_v610","HSC_state3_v612","qc_primary_include",
          "qc_sensitivity_include","qc_original_frozen_cluster")
stopifnot(all(need %in% colnames(md)))

SAMPLES <- c("Sham1","Sham20","Tx17","Tx5")
GROUP <- c(Sham1="Sham",Sham20="Sham",Tx17="Tx",Tx5="Tx")
CORE <- c("Mki67","Top2a","Pcna","Mcm2","Mcm3","Mcm4","Mcm5","Mcm6","Mcm7",
          "Ccna2","Ccnb1","Cdk1")

cnt <- get_counts(src); cnt <- cnt[,colnames(src),drop=FALSE]
gp <- intersect(CORE,rownames(cnt))
x <- cnt[gp,,drop=FALSE]
nhit <- Matrix::colSums(x>0)
mk <- if("Mki67"%in%rownames(x)) as.logical(x["Mki67",]>0) else rep(FALSE,ncol(x))
tp <- if("Top2a"%in%rownames(x)) as.logical(x["Top2a",]>0) else rep(FALSE,ncol(x))
cyc <- mk | tp | nhit>=2

cell <- data.frame(
  cell=colnames(src),
  sample=as.character(md$sample_hsc_v610),
  condition=GROUP[as.character(md$sample_hsc_v610)],
  HSC_state3=as.character(md$HSC_state3_v612),
  cluster=as.character(md$qc_original_frozen_cluster),
  qc_primary_include=md$qc_primary_include,
  qc_sensitivity_include=md$qc_sensitivity_include,
  n_core_detected=as.numeric(nhit),Mki67_detected=mk,Top2a_detected=tp,
  cycling_like=cyc,stringsAsFactors=FALSE
)

w(data.frame(item=c("legacy_definition","core_genes_requested","core_genes_present"),
             value=c("Mki67+ OR Top2a+ OR >=2 core proliferation genes detected",
                     paste(CORE,collapse=";"),paste(gp,collapse=";"))),
  "01_definition.tsv")

sets <- list(Source_all_HSC=colnames(src),PrimaryClean=colnames(pri),
             SensitivityClean=colnames(sen))
L <- list(); k <- 0
for(ds in names(sets)) for(s in SAMPLES){
  allc <- cell$cell %in% sets[[ds]] & cell$sample==s
  d <- cell[allc & cell$HSC_state3=="ECM-activated HSC",,drop=FALSE]
  nall <- sum(allc); nc <- sum(d$cycling_like)
  k <- k+1
  L[[k]] <- data.frame(dataset=ds,sample=s,condition=GROUP[[s]],
    n_ECM_HSC=nrow(d),total_HSC3=nall,ECM_HSC_fraction=nrow(d)/nall,
    n_cycling_like=nc,cycling_like_fraction=nc/nrow(d),
    Mki67_positive_fraction=mean(d$Mki67_detected),
    Top2a_positive_fraction=mean(d$Top2a_detected),
    mean_core_cycling_genes_detected=mean(d$n_core_detected),
    cycling_ECM_HSC_burden=nc/nall,
    cycling_ECM_HSC_per_1000_HSC3=1000*nc/nall)
}
tab <- do.call(rbind,L); w(tab,"02_by_sample.tsv")

metrics <- c("ECM_HSC_fraction","cycling_like_fraction","cycling_ECM_HSC_burden")
R <- list(); k <- 0
for(ds in names(sets)) for(m in metrics){
  d <- tab[tab$dataset==ds,]
  sh <- d[[m]][d$condition=="Sham"]; tx <- d[[m]][d$condition=="Tx"]
  p <- dir4(sh,tx); k <- k+1
  R[[k]] <- cbind(data.frame(dataset=ds,metric=m,
    Sham1=d[[m]][match("Sham1",d$sample)],
    Sham20=d[[m]][match("Sham20",d$sample)],
    Tx17=d[[m]][match("Tx17",d$sample)],
    Tx5=d[[m]][match("Tx5",d$sample)],
    Sham_mean=gm(sh),Tx_mean=gm(tx),Tx_minus_Sham=gm(tx)-gm(sh),
    both_Tx_below_both_Sham=max(tx)<min(sh),
    both_Tx_above_both_Sham=min(tx)>max(sh)),p)
}
rep <- do.call(rbind,R); w(rep,"03_replicate_summary.tsv")

old <- read.csv(OLD,check.names=FALSE)
new <- tab[tab$dataset=="Source_all_HSC",
           c("sample","n_ECM_HSC","n_cycling_like","cycling_like_fraction",
             "Mki67_positive_fraction","Top2a_positive_fraction",
             "mean_core_cycling_genes_detected")]
cmp <- merge(old,new,by="sample",suffixes=c("_old","_reproduced"))
for(v in c("n_ECM_HSC","n_cycling_like","cycling_like_fraction",
           "Mki67_positive_fraction","Top2a_positive_fraction",
           "mean_core_cycling_genes_detected"))
  cmp[[paste0(v,"_diff")]] <- cmp[[paste0(v,"_reproduced")]]-cmp[[paste0(v,"_old")]]
w(cmp,"04_exact_v635_reproduction.tsv")

isecm <- cell$HSC_state3=="ECM-activated HSC"
iscyc <- isecm & cell$cycling_like
loss <- data.frame(
  metric=c("source_legacy_cycling_ECM_HSC",
           "Primary_removed_legacy_cycling_ECM_HSC",
           "Sensitivity_removed_legacy_cycling_ECM_HSC"),
  n=c(sum(iscyc),sum(iscyc & !cell$qc_primary_include),
      sum(iscyc & !cell$qc_sensitivity_include))
)
loss$fraction_of_legacy_cycling_ECM_HSC <- loss$n/sum(iscyc)
w(loss,"05_QC_loss_audit.tsv")

pri_r <- rep[rep$dataset=="PrimaryClean" & rep$metric%in%c("cycling_like_fraction","cycling_ECM_HSC_burden"),]
sen_r <- rep[rep$dataset=="SensitivityClean" & rep$metric%in%c("cycling_like_fraction","cycling_ECM_HSC_burden"),]
diffcols <- grep("_diff$",colnames(cmp),value=TRUE)
exact <- all(abs(as.matrix(cmp[,diffcols,drop=FALSE]))<1e-12,na.rm=TRUE)

dec <- data.frame(
  item=c("Source_reproduces_v6.3.5_exactly",
         "Primary_cycling_fraction_Tx_minus_Sham",
         "Primary_cycling_fraction_pairwise_direction",
         "Primary_cycling_burden_Tx_minus_Sham",
         "Primary_cycling_burden_pairwise_direction",
         "Sensitivity_cycling_fraction_Tx_minus_Sham",
         "Sensitivity_cycling_fraction_pairwise_direction",
         "Sensitivity_cycling_burden_Tx_minus_Sham",
         "Sensitivity_cycling_burden_pairwise_direction",
         "fraction_legacy_cycling_ECM_HSC_removed_Primary",
         "fraction_legacy_cycling_ECM_HSC_removed_Sensitivity"),
  value=c(exact,
    pri_r$Tx_minus_Sham[pri_r$metric=="cycling_like_fraction"],
    pri_r$pairwise_direction[pri_r$metric=="cycling_like_fraction"],
    pri_r$Tx_minus_Sham[pri_r$metric=="cycling_ECM_HSC_burden"],
    pri_r$pairwise_direction[pri_r$metric=="cycling_ECM_HSC_burden"],
    sen_r$Tx_minus_Sham[sen_r$metric=="cycling_like_fraction"],
    sen_r$pairwise_direction[sen_r$metric=="cycling_like_fraction"],
    sen_r$Tx_minus_Sham[sen_r$metric=="cycling_ECM_HSC_burden"],
    sen_r$pairwise_direction[sen_r$metric=="cycling_ECM_HSC_burden"],
    loss$fraction_of_legacy_cycling_ECM_HSC[2],
    loss$fraction_of_legacy_cycling_ECM_HSC[3]))
w(dec,"06_KEY_decision.tsv")

w(data.frame(item=c("source_cells","primary_cells","sensitivity_cells",
                    "original_filtered_counts_used","DecontX_counts_used",
                    "reclustering_performed","biological_replicate"),
             value=c(ncol(src),ncol(pri),ncol(sen),"TRUE","FALSE","FALSE","sample")),
  "07_integrity.tsv")

writeLines(capture.output(sessionInfo()),file.path(OUT,"sessionInfo.txt"))
message("STATUS: SUCCESS")
message("Review 02_by_sample.tsv, 03_replicate_summary.tsv, 04_exact_v635_reproduction.tsv, 05_QC_loss_audit.tsv, 06_KEY_decision.tsv, 07_integrity.tsv")
