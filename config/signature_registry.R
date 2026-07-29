# ============================================================
# Signature registry
# ============================================================

source(file.path(
  "/Users/uenoya/Projects/uenoyscRNA",
  "config",
  "signatures_mouse_liver_allcell.R"
))

source(file.path(
  "/Users/uenoya/Projects/uenoyscRNA",
  "config",
  "signatures_mouse_macrophage.R"
))

signature_registry <- list(
  Mouse_Liver_AllCell = mouse_liver_allcell_signatures,
  Mouse_Macrophage = mouse_macrophage_signatures
)
