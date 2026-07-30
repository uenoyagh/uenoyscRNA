# General reference markers and Ueno markers.
# Gene symbols are mouse-style. Missing genes are automatically ignored.

GENERAL_POSITIVE <- list(
  Hepatocyte = c("Alb","Ttr","Apoa1","Apob","Cyp2e1","Cyp3a11","Hnf4a","Ass1","Cps1"),
  Cholangiocyte = c("Krt19","Krt8","Krt18","Epcam","Krt7","Sox9","Klf5"),
  Hepatic_progenitor = c("Epcam","Sox9","Krt19","Krt8","Krt18","Prom1","Tacstd2","Spp1"),
  LSEC = c("Kdr","Klf2","Klf4","Eng","Emcn","Pecam1","Cdh5","Stab1","Stab2","Clec4g","Fcgr2b"),
  Vascular_endothelial = c("Pecam1","Cdh5","Kdr","Eng","Emcn","Esam","Ramp2"),
  qHSC = c("Lrat","Rgs5","Reln","Dcn","Cygb","Rbp1","Pparg","Gpx3","Des"),
  aHSC = c("Acta2","Tagln","Pdgfrb","Col1a1","Col1a2","Col3a1","Timp1","Lox","Postn"),
  Kupffer_macrophage = c("Clec4f","Timd4","Marco","Vsig4","Cd5l","C1qa","C1qb","C1qc","Adgre1"),
  Monocyte_macrophage = c("Lyz2","Ccr2","Ly6c2","S100a8","S100a9","Ctss","Fcgr1","Lgals3"),
  Dendritic_cell = c("Flt3","Itgax","H2-Ab1","Cd74","Clec10a","Xcr1","Ccr7"),
  Neutrophil = c("S100a8","S100a9","Ly6g","Retnlg","Mpo","Elane","Camp","Ngp"),
  T_cell = c("Cd3d","Cd3e","Cd3g","Trac","Lck","Il7r"),
  NK_cell = c("Nkg7","Klrk1","Klrd1","Prf1","Gzmb","Ccl5"),
  B_cell = c("Cd79a","Cd79b","Ms4a1","Cd37","Cd74","H2-Aa"),
  Plasma_cell = c("Jchain","Mzb1","Sdc1","Xbp1","Igha","Ighm"),
  Erythroid = c("Hbb-bs","Hbb-bt","Hba-a1","Hba-a2","Alas2","Gypa"),
  Mesothelial = c("Msln","Wt1","Krt19","Krt8","Krt18","Upk3b"),
  Cycling = c("Mki67","Top2a","Tuba1b","Cenpf","Birc5","Ube2c")
)

GENERAL_NEGATIVE <- list(
  Hepatocyte = c("Ptprc","Pecam1","Col1a1","Krt19"),
  Cholangiocyte = c("Ptprc","Clec4f","Col1a1"),
  Hepatic_progenitor = c("Ptprc","Clec4f","Pecam1"),
  LSEC = c("Ptprc","Alb","Col1a1"),
  Vascular_endothelial = c("Ptprc","Alb"),
  qHSC = c("Ptprc","Pecam1","Alb"),
  aHSC = c("Ptprc","Alb","Clec4f"),
  Kupffer_macrophage = c("Alb","Pecam1","Col1a1","S100a8","S100a9"),
  Monocyte_macrophage = c("Alb","Pecam1","Col1a1"),
  Dendritic_cell = c("Alb","Clec4f"),
  Neutrophil = c("Alb","Clec4f","C1qc"),
  T_cell = c("Alb","Pecam1","Col1a1"),
  NK_cell = c("Alb","Pecam1","Col1a1"),
  B_cell = c("Alb","Pecam1","Col1a1"),
  Plasma_cell = c("Alb","Pecam1"),
  Erythroid = c("Ptprc","Pecam1"),
  Mesothelial = c("Ptprc","Alb"),
  Cycling = character()
)

UENO_POSITIVE <- list(
  `Resident Kupffer-like` = c("Clec4f","Timd4","Marco","Vsig4","Cd5l","C1qa","C1qb","C1qc"),
  `Monocyte-like` = c("Ccr2","Ly6c2","S100a8","S100a9","Lyz2","Ctss","Fcgr1"),
  `Inflammatory M1-like` = c("Tnf","Il1b","Il6","Il12b","Il23a","Ccl2","Ccl3","Ccl4",
                             "Cxcl9","Cxcl10","Nos2","Ptgs2","Cd80","Cd86"),
  `Pro-resolution M2-like` = c("Il10","Mrc1","Arg1","Retnla","Chil3","Chil4","Ccl17",
                               "Ccl22","Ccl24","Maf","Cd163"),
  `SPP1/TREM2 MASH-associated` = c("Spp1","Trem2","Gpnmb","Lgals3","Cd9","Lpl","Fabp5","Ctsb","Ctsd"),
  `Efferocytosis/phagocytosis-high` = c("Mertk","Axl","Gas6","Mfge8","Timd4","Marco","Cd36",
                                        "Msr1","Lrp1","C1qa","C1qb","C1qc"),
  `IL10-response-high Mphi` = c("Stat3","Socs3","Bcl3","Sbno2","Dusp1","Dusp2","Il4ra","Maf"),
  `Fibrosis-associated Mphi` = c("Spp1","Trem2","Lgals3","Tgfb1","Pdgfb","Mmp12","Mmp14","Ctsb","Ctsk")
)

UENO_NEGATIVE <- list(
  `Resident Kupffer-like` = c("S100a8","S100a9","Ly6c2","Ccr2"),
  `Monocyte-like` = c("Clec4f","Timd4","Vsig4"),
  `Inflammatory M1-like` = c("Alb","Pecam1","Col1a1"),
  `Pro-resolution M2-like` = c("Alb","Pecam1","S100a8","S100a9"),
  `SPP1/TREM2 MASH-associated` = c("Alb","Pecam1"),
  `Efferocytosis/phagocytosis-high` = c("Alb","Pecam1"),
  `IL10-response-high Mphi` = c("Alb","Pecam1"),
  `Fibrosis-associated Mphi` = c("Alb","Pecam1")
)

make_marker_table <- function(marker_list, source, direction) {

  marker_tables <- lapply(names(marker_list), function(celltype) {

    genes <- unique(stats::na.omit(marker_list[[celltype]]))
    genes <- genes[nzchar(genes)]

    if (!length(genes)) {
      return(NULL)
    }

    data.frame(
      source = rep(source, length(genes)),
      direction = rep(direction, length(genes)),
      label = rep(celltype, length(genes)),
      gene = genes,
      stringsAsFactors = FALSE
    )
  })

  marker_tables <- Filter(Negate(is.null), marker_tables)

  do.call(rbind, marker_tables)
}

get_marker_reference <- function() {
  rbind(
    make_marker_table(GENERAL_POSITIVE, "General", "positive"),
    make_marker_table(GENERAL_NEGATIVE, "General", "negative"),
    make_marker_table(UENO_POSITIVE, "Ueno", "positive"),
    make_marker_table(UENO_NEGATIVE, "Ueno", "negative")
  )
}
