# Hierarchical General/Ueno marker reference for RDS3 annotation validation v3.3.
# Mouse-style gene symbols. Missing genes are ignored automatically.
#
# Levels:
#   General       : conventional broad cell-type comparison
#   Ueno_lineage  : broad lineage
#   Ueno_celltype : lineage-restricted cell type
#   Ueno_subtype  : cell-type-restricted functional or zonation subtype

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

UENO_LINEAGE_POSITIVE <- list(
  Hepatocyte = c("Alb","Ttr","Apoa1","Apob","Hnf4a","Cps1","Ass1"),
  Biliary = c("Krt19","Krt8","Krt18","Krt7","Epcam","Sox9"),
  Endothelial = c("Pecam1","Cdh5","Kdr","Eng","Emcn","Esam","Ramp2"),
  Mesenchymal = c("Col1a1","Col1a2","Dcn","Col3a1","Pdgfra","Pdgfrb","Des"),
  Myeloid = c("Ptprc","Lyz2","Tyrobp","Fcerg1","Ctss","LST1","Aif1"),
  Lymphoid = c("Ptprc","Cd3d","Cd3e","Trac","Nkg7","Cd79a","Ms4a1"),
  Erythroid = c("Hbb-bs","Hbb-bt","Hba-a1","Hba-a2","Alas2","Gypa"),
  Megakaryocytic = c("Pf4","Ppbp","Gp9","Itga2b","Nrgp","Tubb1"),
  Cycling = c("Mki67","Top2a","Cenpf","Birc5","Ube2c","Tuba1b")
)

UENO_LINEAGE_NEGATIVE <- list(
  Hepatocyte = c("Ptprc","Pecam1","Col1a1","Krt19"),
  Biliary = c("Ptprc","Clec4f","Pecam1"),
  Endothelial = c("Ptprc","Alb","Col1a1"),
  Mesenchymal = c("Ptprc","Alb","Pecam1"),
  Myeloid = c("Alb","Krt19","Pecam1"),
  Lymphoid = c("Alb","Col1a1","Pecam1"),
  Erythroid = c("Ptprc","Pecam1"),
  Megakaryocytic = c("Alb","Krt19"),
  Cycling = character()
)

UENO_CELLTYPE_POSITIVE <- list(
  Mature_hepatocyte = c("Alb","Ttr","Apoa1","Cps1","Ass1","Hnf4a","Tat"),
  Hepatic_progenitor = c("Epcam","Sox9","Krt19","Prom1","Tacstd2","Spp1"),
  Cholangiocyte = c("Krt19","Krt7","Krt8","Krt18","Epcam","Klf5"),
  LSEC = c("Stab1","Stab2","Clec4g","Fcgr2b","Kdr","Klf2","Emcn"),
  Vascular_endothelial = c("Pecam1","Cdh5","Kdr","Esam","Ramp2","Rgcc"),
  Lymphatic_endothelial = c("Prox1","Lyve1","Pdpn","Flt4","Ccl21a"),
  qHSC = c("Lrat","Rbp1","Reln","Cygb","Pparg","Gpx3","Dcn"),
  aHSC = c("Acta2","Tagln","Col1a1","Col1a2","Col3a1","Postn","Lox"),
  Portal_fibroblast = c("Col15a1","Pi16","Dpt","Dcn","Col1a1","Lum"),
  Pericyte_VSMC = c("Rgs5","Cspg4","Mcam","Pdgfrb","Acta2","Tagln","Myh11"),
  Kupffer_macrophage = c("Clec4f","Timd4","Marco","Vsig4","Cd5l","C1qa","C1qb","C1qc"),
  Monocyte_derived_macrophage = c("Lyz2","Ccr2","Ly6c2","Ctss","Fcgr1","Lgals3","Spp1"),
  Monocyte = c("Ccr2","Ly6c2","S100a8","S100a9","Lyz2","Plac8"),
  Neutrophil = c("Ly6g","S100a8","S100a9","Retnlg","Mpo","Elane","Camp","Ngp"),
  cDC1 = c("Xcr1","Clec9a","Batf3","Irf8","Cadm1"),
  cDC2 = c("Clec10a","Cd209a","Sirpa","Irf4","H2-Ab1","Cd74"),
  pDC = c("Gzmb","Siglech","Bst2","Tcf4","Irf7"),
  Mast_cell = c("Kit","Fcer1a","Cpa3","Tpsb2","Mcpt4"),
  Basophil = c("Mcpt8","Prss34","Il3ra","Fcer1a","Ccl3"),
  B_cell = c("Cd79a","Cd79b","Ms4a1","Cd37","Cd74","H2-Aa"),
  Plasma_cell = c("Jchain","Mzb1","Sdc1","Xbp1","Igha","Ighm"),
  CD4_T_cell = c("Cd3d","Cd3e","Trac","Cd4","Il7r","Ltb"),
  CD8_T_cell = c("Cd3d","Cd3e","Trac","Cd8a","Cd8b1","Ccl5"),
  Treg = c("Foxp3","Il2ra","Ctla4","Ikzf2","Tnfrsf18"),
  Gamma_delta_T = c("Trdc","Trgc1","Trgc2","Cd3d","Cd3e"),
  NK_cell = c("Nkg7","Klrk1","Klrd1","Prf1","Gzmb","Ccl5"),
  NKT_cell = c("Cd3d","Cd3e","Trac","Nkg7","Klrd1","Klrk1"),
  ILC = c("Il7r","Id2","Rora","Gata3","Tox","Kit"),
  Erythroid = c("Hbb-bs","Hbb-bt","Hba-a1","Hba-a2","Alas2","Gypa"),
  Megakaryocyte = c("Pf4","Ppbp","Gp9","Itga2b","Nrgp","Tubb1"),
  Platelet = c("Pf4","Ppbp","Gp9","Tubb1","Nrgp"),
  Cycling = c("Mki67","Top2a","Cenpf","Birc5","Ube2c","Tuba1b")
)

UENO_CELLTYPE_NEGATIVE <- list(
  Mature_hepatocyte = c("Ptprc","Pecam1","Krt19","Col1a1"),
  Hepatic_progenitor = c("Ptprc","Clec4f","Pecam1"),
  Cholangiocyte = c("Ptprc","Clec4f","Col1a1"),
  LSEC = c("Ptprc","Alb","Col1a1"),
  Vascular_endothelial = c("Ptprc","Alb"),
  Lymphatic_endothelial = c("Ptprc","Alb"),
  qHSC = c("Ptprc","Alb","Pecam1"),
  aHSC = c("Ptprc","Alb","Pecam1"),
  Portal_fibroblast = c("Ptprc","Alb","Pecam1"),
  Pericyte_VSMC = c("Ptprc","Alb"),
  Kupffer_macrophage = c("Alb","Pecam1","S100a8","S100a9","Ly6c2"),
  Monocyte_derived_macrophage = c("Alb","Pecam1","Clec4f"),
  Monocyte = c("Alb","Pecam1","Clec4f","Timd4"),
  Neutrophil = c("Alb","Clec4f","C1qc"),
  cDC1 = c("Alb","Clec4f"),
  cDC2 = c("Alb","Clec4f"),
  pDC = c("Alb","Clec4f"),
  Mast_cell = c("Alb","Pecam1"),
  Basophil = c("Alb","Pecam1"),
  B_cell = c("Alb","Pecam1","Col1a1"),
  Plasma_cell = c("Alb","Pecam1"),
  CD4_T_cell = c("Alb","Pecam1","Col1a1"),
  CD8_T_cell = c("Alb","Pecam1","Col1a1"),
  Treg = c("Alb","Pecam1","Col1a1"),
  Gamma_delta_T = c("Alb","Pecam1","Col1a1"),
  NK_cell = c("Alb","Pecam1","Col1a1"),
  NKT_cell = c("Alb","Pecam1","Col1a1"),
  ILC = c("Alb","Pecam1","Col1a1"),
  Erythroid = c("Ptprc","Pecam1"),
  Megakaryocyte = c("Alb","Krt19"),
  Platelet = c("Alb","Krt19"),
  Cycling = character()
)

# Parent lineage for each Ueno cell type.
UENO_CELLTYPE_PARENT <- c(
  Mature_hepatocyte = "Hepatocyte",
  Hepatic_progenitor = "Hepatocyte",
  Cholangiocyte = "Biliary",
  LSEC = "Endothelial",
  Vascular_endothelial = "Endothelial",
  Lymphatic_endothelial = "Endothelial",
  qHSC = "Mesenchymal",
  aHSC = "Mesenchymal",
  Portal_fibroblast = "Mesenchymal",
  Pericyte_VSMC = "Mesenchymal",
  Kupffer_macrophage = "Myeloid",
  Monocyte_derived_macrophage = "Myeloid",
  Monocyte = "Myeloid",
  Neutrophil = "Myeloid",
  cDC1 = "Myeloid",
  cDC2 = "Myeloid",
  pDC = "Myeloid",
  Mast_cell = "Myeloid",
  Basophil = "Myeloid",
  B_cell = "Lymphoid",
  Plasma_cell = "Lymphoid",
  CD4_T_cell = "Lymphoid",
  CD8_T_cell = "Lymphoid",
  Treg = "Lymphoid",
  Gamma_delta_T = "Lymphoid",
  NK_cell = "Lymphoid",
  NKT_cell = "Lymphoid",
  ILC = "Lymphoid",
  Erythroid = "Erythroid",
  Megakaryocyte = "Megakaryocytic",
  Platelet = "Megakaryocytic",
  Cycling = "Cycling"
)

UENO_SUBTYPE_POSITIVE <- list(
  Periportal_hepatocyte = c("Cps1","Ass1","Asl","Arg1","Pck1","G6pc","Hal"),
  Midzonal_hepatocyte = c("Hamp","Igfbp2","Cyp8b1","Gstm1"),
  Pericentral_hepatocyte = c("Glul","Cyp2e1","Cyp1a2","Axin2","Lect2","Oat"),
  Stress_response_hepatocyte = c("Fos","Jun","Atf3","Ddit3","Hspa1a","Hspa1b"),
  IFN_response_hepatocyte = c("Isg15","Ifit1","Ifit2","Ifit3","Mx1","Oas1a"),
  Cycling_hepatocyte = c("Mki67","Top2a","Cenpf","Birc5","Ube2c"),
  Reactive_cholangiocyte = c("Krt19","Sox9","Spp1","Krt8","Krt18","Lgals3"),
  Periportal_LSEC = c("Efnb2","Dll4","Sox17","Gja5","Kdr"),
  Pericentral_LSEC = c("Wnt2","Wnt9b","Rspo3","Clec4g","Stab2"),
  Capillarized_LSEC = c("Pecam1","Cdh5","Kdr","Vwf","Klf2"),
  Angiogenic_LSEC = c("Kdr","Esm1","Apln","Kcne3","Ramp2"),
  Inflammatory_LSEC = c("Icam1","Vcam1","Sele","Cxcl10","Ccl2"),
  Early_activated_HSC = c("Pdgfrb","Tagln","Timp1","Col1a1","Ctgf"),
  Myofibroblastic_aHSC = c("Acta2","Tagln","Myh11","Cnn1","Lox"),
  Fibrogenic_aHSC = c("Col1a1","Col1a2","Col3a1","Postn","Timp1","Loxl1"),
  Inflammatory_aHSC = c("Il6","Ccl2","Cxcl10","Icam1","Tnfaip3"),
  Resident_Kupffer_like = c("Clec4f","Timd4","Marco","Vsig4","Cd5l","C1qa","C1qb","C1qc"),
  Monocyte_like = c("Ccr2","Ly6c2","S100a8","S100a9","Lyz2","Ctss","Fcgr1"),
  Inflammatory_M1_like = c("Tnf","Il1b","Il6","Il12b","Il23a","Ccl2","Ccl3","Ccl4","Cxcl9","Cxcl10","Nos2","Ptgs2","Cd80","Cd86"),
  Pro_resolution_M2_like = c("Il10","Mrc1","Arg1","Retnla","Chil3","Chil4","Ccl17","Ccl22","Ccl24","Maf","Cd163"),
  SPP1_TREM2_MASH_associated = c("Spp1","Trem2","Gpnmb","Lgals3","Cd9","Lpl","Fabp5","Ctsb","Ctsd"),
  Lipid_associated_macrophage = c("Trem2","Lpl","Cd9","Gpnmb","Fabp5","Apoe","Ctsd"),
  Efferocytosis_phagocytosis_high = c("Mertk","Axl","Gas6","Mfge8","Timd4","Marco","Cd36","Msr1","Lrp1","C1qa","C1qb","C1qc"),
  IL10_response_high_Mphi = c("Stat3","Socs3","Bcl3","Sbno2","Dusp1","Dusp2","Il4ra","Maf"),
  Fibrosis_associated_Mphi = c("Spp1","Trem2","Lgals3","Tgfb1","Pdgfb","Mmp12","Mmp14","Ctsb","Ctsk"),
  IFN_response_macrophage = c("Isg15","Ifit1","Ifit2","Ifit3","Cxcl10","Mx1"),
  Classical_monocyte = c("Ccr2","Ly6c2","S100a8","S100a9","Plac8"),
  Nonclassical_monocyte = c("Nr4a1","Cx3cr1","LST1","Fcgr3","Lair1"),
  Inflammatory_neutrophil = c("S100a8","S100a9","Il1b","Cxcl2","Ccl3"),
  Aged_neutrophil = c("Cxcr4","Sell","Icam1","Bcl2a1b","Socs3"),
  Naive_CD4_T = c("Il7r","Ltb","Malat1","Ccr7","Lef1","Tcf7"),
  Activated_CD4_T = c("Cd69","Icos","Il2ra","Fos","Jun","Tnfrsf4"),
  Naive_CD8_T = c("Ccr7","Lef1","Tcf7","Il7r","Ltb"),
  Cytotoxic_CD8_T = c("Nkg7","Ccl5","Gzmb","Prf1","Gzmk"),
  Exhausted_CD8_T = c("Pdcd1","Lag3","Tigit","Havcr2","Ctla4"),
  Activated_NK = c("Nkg7","Ccl5","Gzmb","Prf1","Ifng"),
  Naive_B = c("Ms4a1","Cd79a","Cd74","H2-Aa","Ighm"),
  Memory_B = c("Cd37","Cd79a","Cd74","Cd44","Cd69"),
  Plasma_cell_subtype = c("Jchain","Mzb1","Xbp1","Sdc1","Igha")
)

UENO_SUBTYPE_NEGATIVE <- lapply(UENO_SUBTYPE_POSITIVE, function(x) character())

UENO_SUBTYPE_PARENT <- c(
  Periportal_hepatocyte = "Mature_hepatocyte",
  Midzonal_hepatocyte = "Mature_hepatocyte",
  Pericentral_hepatocyte = "Mature_hepatocyte",
  Stress_response_hepatocyte = "Mature_hepatocyte",
  IFN_response_hepatocyte = "Mature_hepatocyte",
  Cycling_hepatocyte = "Mature_hepatocyte",
  Reactive_cholangiocyte = "Cholangiocyte",
  Periportal_LSEC = "LSEC",
  Pericentral_LSEC = "LSEC",
  Capillarized_LSEC = "LSEC",
  Angiogenic_LSEC = "LSEC",
  Inflammatory_LSEC = "LSEC",
  Early_activated_HSC = "aHSC",
  Myofibroblastic_aHSC = "aHSC",
  Fibrogenic_aHSC = "aHSC",
  Inflammatory_aHSC = "aHSC",
  Resident_Kupffer_like = "Kupffer_macrophage",
  Monocyte_like = "Monocyte_derived_macrophage",
  Inflammatory_M1_like = "Monocyte_derived_macrophage",
  Pro_resolution_M2_like = "Monocyte_derived_macrophage",
  SPP1_TREM2_MASH_associated = "Monocyte_derived_macrophage",
  Lipid_associated_macrophage = "Monocyte_derived_macrophage",
  Efferocytosis_phagocytosis_high = "Monocyte_derived_macrophage",
  IL10_response_high_Mphi = "Monocyte_derived_macrophage",
  Fibrosis_associated_Mphi = "Monocyte_derived_macrophage",
  IFN_response_macrophage = "Monocyte_derived_macrophage",
  Classical_monocyte = "Monocyte",
  Nonclassical_monocyte = "Monocyte",
  Inflammatory_neutrophil = "Neutrophil",
  Aged_neutrophil = "Neutrophil",
  Naive_CD4_T = "CD4_T_cell",
  Activated_CD4_T = "CD4_T_cell",
  Naive_CD8_T = "CD8_T_cell",
  Cytotoxic_CD8_T = "CD8_T_cell",
  Exhausted_CD8_T = "CD8_T_cell",
  Activated_NK = "NK_cell",
  Naive_B = "B_cell",
  Memory_B = "B_cell",
  Plasma_cell_subtype = "Plasma_cell"
)

make_marker_table_v33 <- function(marker_list, source, level, direction,
                                  parent_lineage = NULL, parent_celltype = NULL) {
  out <- lapply(names(marker_list), function(label) {
    genes <- unique(stats::na.omit(marker_list[[label]]))
    genes <- genes[nzchar(genes)]
    if (!length(genes)) return(NULL)

    pl <- if (is.null(parent_lineage)) NA_character_ else unname(parent_lineage[label])
    pc <- if (is.null(parent_celltype)) NA_character_ else unname(parent_celltype[label])

    data.frame(
      source = source,
      level = level,
      parent_lineage = pl,
      parent_celltype = pc,
      direction = direction,
      label = label,
      gene = genes,
      stringsAsFactors = FALSE
    )
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) {
    return(data.frame(
      source=character(), level=character(), parent_lineage=character(),
      parent_celltype=character(), direction=character(), label=character(),
      gene=character(), stringsAsFactors=FALSE
    ))
  }
  do.call(rbind, out)
}

get_marker_reference_v33 <- function() {
  rbind(
    make_marker_table_v33(GENERAL_POSITIVE, "General", "general", "positive"),
    make_marker_table_v33(GENERAL_NEGATIVE, "General", "general", "negative"),
    make_marker_table_v33(UENO_LINEAGE_POSITIVE, "Ueno", "lineage", "positive"),
    make_marker_table_v33(UENO_LINEAGE_NEGATIVE, "Ueno", "lineage", "negative"),
    make_marker_table_v33(
      UENO_CELLTYPE_POSITIVE, "Ueno", "celltype", "positive",
      parent_lineage = UENO_CELLTYPE_PARENT
    ),
    make_marker_table_v33(
      UENO_CELLTYPE_NEGATIVE, "Ueno", "celltype", "negative",
      parent_lineage = UENO_CELLTYPE_PARENT
    ),
    make_marker_table_v33(
      UENO_SUBTYPE_POSITIVE, "Ueno", "subtype", "positive",
      parent_celltype = UENO_SUBTYPE_PARENT
    ),
    make_marker_table_v33(
      UENO_SUBTYPE_NEGATIVE, "Ueno", "subtype", "negative",
      parent_celltype = UENO_SUBTYPE_PARENT
    )
  )
}
