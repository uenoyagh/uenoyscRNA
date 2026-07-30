# v3.4 marker augmentation.
# Source v3.3 first, then apply RDS3-supported additions.

source(file.path(project_root_guess, "R", "02_markers_v3.3.R"))

append_unique <- function(x, add) unique(c(x, add))

# RDS3-supported broad/cell-type markers observed among high-ranking presto markers.
UENO_CELLTYPE_POSITIVE$Cholangiocyte <- append_unique(
  UENO_CELLTYPE_POSITIVE$Cholangiocyte,
  c("Fgfr3","Kifc3","Clu","Mmp7","Kif12","Anxa4","Agrn","Sorbs2")
)

UENO_CELLTYPE_POSITIVE$qHSC <- append_unique(
  UENO_CELLTYPE_POSITIVE$qHSC,
  c("Cxcl12","Hand2","Septin4","C4b","Hgf","Colec11","Raph1","Mxra8")
)

UENO_CELLTYPE_POSITIVE$aHSC <- append_unique(
  UENO_CELLTYPE_POSITIVE$aHSC,
  c("Cxcl12","Ctgf","Timp1","Loxl1","Pdgfrb","Cygb","Col1a2")
)

UENO_CELLTYPE_POSITIVE$Portal_fibroblast <- append_unique(
  UENO_CELLTYPE_POSITIVE$Portal_fibroblast,
  c("Col15a1","Pi16","Dpt","Lum","Colec11","C7","C4b")
)

UENO_CELLTYPE_POSITIVE$Kupffer_macrophage <- append_unique(
  UENO_CELLTYPE_POSITIVE$Kupffer_macrophage,
  c("Csf1r","Csf2ra","Apobec1","Tcirg1","Cd72","Lilrb4a","Psap")
)

UENO_CELLTYPE_POSITIVE$Monocyte_derived_macrophage <- append_unique(
  UENO_CELLTYPE_POSITIVE$Monocyte_derived_macrophage,
  c("Csf1r","Axl","Gpnmb","Adam8","Lilrb4a","Tcirg1","Psap","Ctsb","Ctsd")
)

UENO_CELLTYPE_POSITIVE$LSEC <- append_unique(
  UENO_CELLTYPE_POSITIVE$LSEC,
  c("Egfl7","Jam2","Ptprb","Adgrf5","Adgrl4","Gpihbp1","Cd300lg","Aqp1")
)

UENO_CELLTYPE_POSITIVE$Vascular_endothelial <- append_unique(
  UENO_CELLTYPE_POSITIVE$Vascular_endothelial,
  c("Egfl7","Jam2","Ptprb","Adgrf5","Adgrl4","Cd300lg","Aqp1","Vwf")
)

UENO_CELLTYPE_POSITIVE$cDC1 <- append_unique(
  UENO_CELLTYPE_POSITIVE$cDC1,
  c("Wdfy4","Ciita","Plbd1","Naaa","Cst3")
)

UENO_CELLTYPE_POSITIVE$cDC2 <- append_unique(
  UENO_CELLTYPE_POSITIVE$cDC2,
  c("Ciita","Plbd1","Naaa","Cst3","Il4i1","Relb")
)

# Stronger mutual exclusions.
UENO_CELLTYPE_NEGATIVE$Mature_hepatocyte <- append_unique(
  UENO_CELLTYPE_NEGATIVE$Mature_hepatocyte,
  c("Tyrobp","Lyz2","Csf1r","Cdh5","Pecam1","Col1a2","Dcn")
)
UENO_CELLTYPE_NEGATIVE$Hepatic_progenitor <- append_unique(
  UENO_CELLTYPE_NEGATIVE$Hepatic_progenitor,
  c("Tyrobp","Lyz2","Csf1r","Col1a1","Col1a2")
)
UENO_CELLTYPE_NEGATIVE$Portal_fibroblast <- append_unique(
  UENO_CELLTYPE_NEGATIVE$Portal_fibroblast,
  c("Clec4f","Csf1r","Tyrobp","Cdh5","Pecam1","Alb")
)
UENO_CELLTYPE_NEGATIVE$qHSC <- append_unique(
  UENO_CELLTYPE_NEGATIVE$qHSC,
  c("Clec4f","Csf1r","Tyrobp","Cdh5","Pecam1","Alb")
)
UENO_CELLTYPE_NEGATIVE$aHSC <- append_unique(
  UENO_CELLTYPE_NEGATIVE$aHSC,
  c("Clec4f","Csf1r","Tyrobp","Cdh5","Pecam1","Alb")
)
UENO_CELLTYPE_NEGATIVE$Monocyte <- append_unique(
  UENO_CELLTYPE_NEGATIVE$Monocyte,
  c("C1qa","C1qb","C1qc","Marco","Vsig4","Cd5l")
)
UENO_CELLTYPE_NEGATIVE$cDC2 <- append_unique(
  UENO_CELLTYPE_NEGATIVE$cDC2,
  c("Ly6g","Retnlg","Mpo","Elane","Camp","Ngp","Clec4f","Timd4")
)

# Required-marker gates. Values are marker vectors; thresholds are defined separately.
UENO_REQUIRED_MARKERS <- list(
  lineage = list(
    Hepatocyte = c("Alb","Ttr","Hnf4a","Cps1","Ass1"),
    Biliary = c("Krt19","Krt8","Krt18","Epcam","Sox9"),
    Endothelial = c("Pecam1","Cdh5","Kdr","Eng","Emcn"),
    Mesenchymal = c("Col1a1","Col1a2","Dcn","Pdgfrb","Cygb"),
    Myeloid = c("Ptprc","Lyz2","Tyrobp","Ctss","Aif1","Csf1r"),
    Lymphoid = c("Ptprc","Cd3d","Cd3e","Trac","Nkg7","Cd79a"),
    Erythroid = c("Hbb-bs","Hbb-bt","Hba-a2","Alas2"),
    Megakaryocytic = c("Pf4","Ppbp","Gp9","Itga2b","Tubb1"),
    Cycling = c("Mki67","Top2a","Cenpf","Birc5","Ube2c")
  ),
  celltype = list(
    Mature_hepatocyte = c("Alb","Ttr","Hnf4a","Cps1","Ass1","Tat"),
    Hepatic_progenitor = c("Epcam","Sox9","Krt19","Prom1","Tacstd2"),
    Cholangiocyte = c("Krt19","Krt7","Epcam","Sox9","Mmp7","Kif12","Fgfr3"),
    LSEC = c("Stab1","Stab2","Clec4g","Fcgr2b","Kdr","Pecam1","Cdh5"),
    Vascular_endothelial = c("Pecam1","Cdh5","Kdr","Vwf","Ptprb"),
    qHSC = c("Lrat","Rbp1","Cygb","Dcn","Cxcl12","Pdgfrb"),
    aHSC = c("Acta2","Tagln","Col1a1","Col1a2","Postn","Timp1"),
    Portal_fibroblast = c("Col15a1","Pi16","Dpt","Lum","Colec11"),
    Pericyte_VSMC = c("Rgs5","Cspg4","Mcam","Pdgfrb","Myh11"),
    Kupffer_macrophage = c("Clec4f","Timd4","Marco","Vsig4","Cd5l","C1qa","C1qb","C1qc"),
    Monocyte_derived_macrophage = c("Csf1r","Lyz2","Ctss","Fcgr1","Lgals3","Axl","Gpnmb"),
    Monocyte = c("Ccr2","Ly6c2","S100a8","S100a9","Plac8"),
    Neutrophil = c("Ly6g","Retnlg","Mpo","Elane","Camp","Ngp"),
    cDC1 = c("Xcr1","Clec9a","Batf3","Irf8","Wdfy4"),
    cDC2 = c("Clec10a","Cd209a","Sirpa","Irf4","H2-Ab1","Cd74"),
    B_cell = c("Cd79a","Cd79b","Ms4a1","Cd37"),
    Plasma_cell = c("Jchain","Mzb1","Sdc1","Xbp1"),
    CD4_T_cell = c("Cd3d","Cd3e","Trac","Cd4","Il7r"),
    CD8_T_cell = c("Cd3d","Cd3e","Trac","Cd8a","Cd8b1"),
    NK_cell = c("Nkg7","Klrk1","Klrd1","Prf1"),
    Erythroid = c("Hbb-bs","Hbb-bt","Hba-a2","Alas2"),
    Megakaryocyte = c("Pf4","Ppbp","Gp9","Itga2b","Tubb1"),
    Platelet = c("Pf4","Ppbp","Gp9","Tubb1"),
    Cycling = c("Mki67","Top2a","Cenpf","Birc5","Ube2c")
  )
)

UENO_REQUIRED_MIN <- list(
  lineage = 2L,
  celltype = 2L
)

get_marker_reference_v34 <- function() {
  get_marker_reference_v33()
}
