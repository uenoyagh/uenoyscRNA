# ============================================================
# Mouse macrophage functional signatures
# ============================================================

mouse_macrophage_signatures <- list(
  Resident_Kupffer_like = c(
    "Clec4f", "Timd4", "Vsig4", "Marco", "Cd5l", "C1qa", "C1qb", "C1qc"
  ),
  Monocyte_like = c(
    "Ly6c2", "Ccr2", "S100a8", "S100a9", "Lyz2", "Ctss", "Fcgr3", "Ms4a7"
  ),
  Inflammatory_M1_like = c(
    "Tnf", "Il1b", "Il6", "Nos2", "Ptgs2", "Cd80", "Cd86", "Cxcl10"
  ),
  Pro_resolution_M2_like = c(
    "Mrc1", "Arg1", "Retnla", "Chil3", "Il10", "Ccl24", "Maf", "Klf4"
  ),
  IL10_response = c(
    "Il10ra", "Il10rb", "Stat3", "Socs3", "Bcl3", "Dusp1", "Nfkbia", "Sbno2"
  ),
  Efferocytosis_phagocytosis = c(
    "Mertk", "Axl", "Tyro3", "Gas6", "Mfge8", "Lrp1", "Cd36", "Msr1"
  ),
  SPP1_TREM2_MASH_assoc = c(
    "Spp1", "Trem2", "Gpnmb", "Lgals3", "Fabp5", "Ctsb", "Ctsd", "Lpl"
  ),
  Fibrosis_assoc_macrophage = c(
    "Spp1", "Tgfb1", "Pdgfb", "Mmp12", "Mmp14", "Lgals3", "Fn1", "Thbs1"
  ),
  Complement = c(
    "C1qa", "C1qb", "C1qc", "C3", "C4b", "Cfb", "Serping1", "Cfh"
  ),
  Cycling = c(
    "Mki67", "Top2a", "Pcna", "Cdk1", "Ccnb1", "Ccnb2", "Ube2c", "Tpx2"
  )
)
