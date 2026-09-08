#!/bin/bash
set -euo pipefail

BIOLOGICAL_VERSION="v1.4.2B.3"
TECHNICAL_VERSION="v1.4.2B.3.3.1"

SRC="/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_iPSC_derived_v1.4.2B.3.2/analysis/v1.4.2B.3.2_IL10_probe_integrated_audit"

DST="/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_iPSC_derived_v1.4.2B.3.3.1/analysis/v1.4.2B.3.3.1_IL10_probe_integrated_audit"

TAB="$DST/tables"
FIG="$DST/figures"

echo "============================================================"
echo "Human MASH ${TECHNICAL_VERSION}"
echo "IL10 post-audit figure repair + harmonized sensitivity"
echo "Biological analysis version: ${BIOLOGICAL_VERSION}"
echo "============================================================"
echo "SOURCE: $SRC"
echo "DEST:   $DST"
echo

if [ ! -d "$SRC" ]; then
    echo "ERROR: validated v1.4.2B.3.2 source directory not found:"
    echo "$SRC"
    exit 1
fi

REQ=(
    "$SRC/tables/IL10_cell_level_all_RDS.tsv"
    "$SRC/tables/IL10_cell_level_probe_or_gene_positive.tsv"
    "$SRC/tables/IL10_probe_filter_status_by_sample.tsv"
    "$SRC/tables/IL10_sample_summary.tsv"
    "$SRC/tables/IL10_gene_vs_probe_validation.tsv"
)

for F in "${REQ[@]}"; do
    if [ ! -f "$F" ]; then
        echo "ERROR: required source file missing:"
        echo "$F"
        exit 1
    fi
done

mkdir -p "$TAB" "$FIG"

cp -f "$SRC"/tables/* "$TAB"/
cp -f "$SRC"/figures/* "$FIG"/

if [ -f "$SRC/IL10_probe_integrated_summary.md" ]; then
    cp -f "$SRC/IL10_probe_integrated_summary.md" \
          "$DST/IL10_probe_integrated_summary_v3.2_original.md"
fi

export SRC_DIR="$SRC"
export DST_DIR="$DST"
export TABLE_DIR="$TAB"
export FIG_DIR="$FIG"

Rscript --vanilla - <<'RS'
suppressPackageStartupMessages({
    library(ggplot2)
})

SRC <- Sys.getenv("SRC_DIR")
DST <- Sys.getenv("DST_DIR")
TAB <- Sys.getenv("TABLE_DIR")
FIG <- Sys.getenv("FIG_DIR")

SAMPLES <- c("SfLB1", "SfLB2", "Tx17", "Tx5")
PROBES <- c("6b6892b", "b1bd8e2", "e99f652")

read_tsv <- function(path) {
    read.delim(
        path,
        stringsAsFactors=FALSE,
        check.names=FALSE
    )
}

to_bool <- function(x) {
    tolower(as.character(x)) %in% c(
        "true", "1", "t", "yes"
    )
}

cell <- read_tsv(
    file.path(
        SRC,
        "tables",
        "IL10_cell_level_all_RDS.tsv"
    )
)

probe_meta <- read_tsv(
    file.path(
        SRC,
        "tables",
        "IL10_probe_filter_status_by_sample.tsv"
    )
)

validation <- read_tsv(
    file.path(
        SRC,
        "tables",
        "IL10_gene_vs_probe_validation.tsv"
    )
)

required_cell_cols <- c(
    "sample",
    "rds_cell",
    "barcode_core",
    "gene_IL10_UMI_RDS",
    "lineage",
    "state",
    "audit_class",
    "IL10_6b6892b_UMI",
    "IL10_b1bd8e2_UMI",
    "IL10_e99f652_UMI",
    "accepted_6b6892b",
    "accepted_b1bd8e2",
    "accepted_e99f652",
    "accepted_IL10_UMI",
    "filteredout_IL10_UMI",
    "all3_IL10_probe_UMI",
    "n_positive_probe_accepted",
    "gene_vs_accepted_concordant"
)

missing_cell <- setdiff(
    required_cell_cols,
    colnames(cell)
)

if (length(missing_cell) > 0) {
    stop(
        "Missing required cell columns: ",
        paste(missing_cell, collapse=", ")
    )
}

if (!setequal(unique(cell$sample), SAMPLES)) {
    stop(
        "Unexpected sample set in cell table: ",
        paste(unique(cell$sample), collapse=", ")
    )
}

if (!all(to_bool(validation$exact_total_match))) {
    stop(
        "Source v3.2 validation is not fully exact."
    )
}

if (any(validation$n_cellwise_mismatches != 0)) {
    stop(
        "Source v3.2 has cellwise gene/probe mismatches."
    )
}

if (!all(to_bool(cell$gene_vs_accepted_concordant))) {
    stop(
        "Source cell table contains gene/probe discordance."
    )
}

probe_meta$filtered_probes <- to_bool(
    probe_meta$filtered_probes
)

probe_meta$included <- to_bool(
    probe_meta$included
)

meta_expected_n <- length(SAMPLES) * length(PROBES)

if (nrow(probe_meta) != meta_expected_n) {
    stop(
        "Unexpected probe metadata row count: ",
        nrow(probe_meta),
        " expected ",
        meta_expected_n
    )
}

probe_accept_matrix <- matrix(
    FALSE,
    nrow=length(PROBES),
    ncol=length(SAMPLES),
    dimnames=list(PROBES, SAMPLES)
)

for (P in PROBES) {
    for (S in SAMPLES) {
        z <- probe_meta[
            probe_meta$short_probe == P &
            probe_meta$rds_sample == S,
            ,
            drop=FALSE
        ]

        if (nrow(z) != 1) {
            stop(
                "Probe metadata ambiguity for ",
                S,
                " / ",
                P
            )
        }

        probe_accept_matrix[P, S] <-
            isTRUE(z$filtered_probes[1])
    }
}

shared_accepted_probes <- rownames(
    probe_accept_matrix
)[
    apply(
        probe_accept_matrix,
        1,
        all
    )
]

if (length(shared_accepted_probes) == 0) {
    stop(
        "No IL10 probe is accepted in all four samples."
    )
}

cat("\n============================================================\n")
cat("SHARED ACCEPTED IL10 PROBES\n")
cat("============================================================\n")
print(probe_accept_matrix)
cat(
    "\nCommon accepted probe(s): ",
    paste(shared_accepted_probes, collapse=", "),
    "\n",
    sep=""
)

write.table(
    data.frame(
        probe=rownames(probe_accept_matrix),
        probe_accept_matrix,
        shared_accepted_across_all_samples=
            apply(probe_accept_matrix, 1, all),
        check.names=FALSE
    ),
    file=file.path(
        TAB,
        "IL10_probe_acceptance_matrix.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

umi_col_for_probe <- function(p) {
    paste0(
        "IL10_",
        p,
        "_UMI"
    )
}

shared_umi_cols <- vapply(
    shared_accepted_probes,
    umi_col_for_probe,
    character(1)
)

missing_umi <- setdiff(
    shared_umi_cols,
    colnames(cell)
)

if (length(missing_umi) > 0) {
    stop(
        "Missing UMI columns for shared probe(s): ",
        paste(missing_umi, collapse=", ")
    )
}

cell$strict_harmonized_IL10_UMI <- rowSums(
    cell[
        ,
        shared_umi_cols,
        drop=FALSE
    ]
)

cell$strict_harmonized_positive <-
    cell$strict_harmonized_IL10_UMI > 0

cell$n_positive_shared_accepted_probes <- rowSums(
    cell[
        ,
        shared_umi_cols,
        drop=FALSE
    ] > 0
)

cell$shared_accepted_probe_set <- paste(
    shared_accepted_probes,
    collapse=";"
)

write.table(
    cell,
    file=file.path(
        TAB,
        "IL10_cell_level_all_RDS_v3.3_with_strict_harmonized.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

write.table(
    cell[
        cell$strict_harmonized_positive,
        ,
        drop=FALSE
    ],
    file=file.path(
        TAB,
        "IL10_strict_harmonized_positive_cells.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

sample_summary <- do.call(
    rbind,
    lapply(
        SAMPLES,
        function(S) {
            x <- cell[
                cell$sample == S,
                ,
                drop=FALSE
            ]

            strongest <- if (
                max(x$strict_harmonized_IL10_UMI) > 0
            ) {
                i <- which.max(
                    x$strict_harmonized_IL10_UMI
                )
                list(
                    cell=x$rds_cell[i],
                    lineage=x$lineage[i],
                    state=x$state[i]
                )
            } else {
                list(
                    cell=NA_character_,
                    lineage=NA_character_,
                    state=NA_character_
                )
            }

            data.frame(
                sample=S,
                n_RDS_cells=nrow(x),
                shared_accepted_probe_set=
                    paste(
                        shared_accepted_probes,
                        collapse=";"
                    ),
                strict_harmonized_IL10_total_UMI=
                    sum(
                        x$strict_harmonized_IL10_UMI
                    ),
                strict_harmonized_IL10_positive_cells=
                    sum(
                        x$strict_harmonized_positive
                    ),
                max_strict_harmonized_IL10_UMI=
                    max(
                        x$strict_harmonized_IL10_UMI
                    ),
                strongest_cell=strongest$cell,
                strongest_lineage=strongest$lineage,
                strongest_state=strongest$state,
                stringsAsFactors=FALSE
            )
        }
    )
)

write.table(
    sample_summary,
    file=file.path(
        TAB,
        "IL10_strict_harmonized_sample_summary.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

lineage_key <- interaction(
    cell$sample,
    cell$lineage,
    drop=TRUE,
    lex.order=TRUE
)

lineage_summary <- do.call(
    rbind,
    lapply(
        split(cell, lineage_key),
        function(x) {
            data.frame(
                sample=unique(x$sample),
                lineage=unique(x$lineage),
                n_cells=nrow(x),
                strict_harmonized_IL10_total_UMI=
                    sum(
                        x$strict_harmonized_IL10_UMI
                    ),
                strict_harmonized_IL10_positive_cells=
                    sum(
                        x$strict_harmonized_positive
                    ),
                max_strict_harmonized_IL10_UMI=
                    max(
                        x$strict_harmonized_IL10_UMI
                    ),
                stringsAsFactors=FALSE
            )
        }
    )
)

write.table(
    lineage_summary,
    file=file.path(
        TAB,
        "IL10_strict_harmonized_lineage_summary.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ------------------------------------------------------------
# Corrected Fig03
# IMPORTANT:
# - factor identity = full RDS cell ID (unique)
# - display label = shortened barcode for readability only
# This prevents different cells from collapsing into one row.
# ------------------------------------------------------------

heat <- cell[
    cell$all3_IL10_probe_UMI > 0 |
    cell$gene_IL10_UMI_RDS > 0,
    ,
    drop=FALSE
]

heat$sample <- factor(
    heat$sample,
    levels=SAMPLES
)

heat <- heat[
    order(
        heat$sample,
        -heat$accepted_IL10_UMI,
        -heat$gene_IL10_UMI_RDS,
        heat$lineage,
        heat$barcode_core
    ),
    ,
    drop=FALSE
]

heat$cell_key <- heat$rds_cell

heat$display_label <- paste0(
    heat$sample,
    " | ",
    substr(
        heat$barcode_core,
        1,
        10
    ),
    " | ",
    heat$lineage,
    " | gene=",
    heat$gene_IL10_UMI_RDS
)

# If shortened display labels collide, append a deterministic suffix.
dup_groups <- split(
    seq_len(nrow(heat)),
    heat$display_label
)

for (idx in dup_groups) {
    if (length(idx) > 1) {
        heat$display_label[idx] <- paste0(
            heat$display_label[idx],
            " | #",
            seq_along(idx)
        )
    }
}

if (anyDuplicated(heat$cell_key)) {
    stop(
        "Full RDS cell IDs are not unique in corrected Fig03 input."
    )
}

if (anyDuplicated(heat$display_label)) {
    stop(
        "Display labels still collide after repair."
    )
}

label_map <- setNames(
    heat$display_label,
    heat$cell_key
)

heat_long <- do.call(
    rbind,
    lapply(
        PROBES,
        function(P) {
            data.frame(
                sample=heat$sample,
                cell_key=heat$cell_key,
                probe=P,
                UMI=heat[[
                    paste0(
                        "IL10_",
                        P,
                        "_UMI"
                    )
                ]],
                accepted=to_bool(
                    heat[[
                        paste0(
                            "accepted_",
                            P
                        )
                    ]]
                ),
                stringsAsFactors=FALSE
            )
        }
    )
)

heat_long$sample <- factor(
    heat_long$sample,
    levels=SAMPLES
)

heat_long$probe <- factor(
    heat_long$probe,
    levels=PROBES
)

heat_long$cell_key <- factor(
    heat_long$cell_key,
    levels=rev(
        heat$cell_key
    )
)

p3 <- ggplot(
    heat_long,
    aes(
        x=probe,
        y=cell_key,
        fill=UMI,
        alpha=accepted
    )
) +
    geom_tile() +
    geom_text(
        aes(
            label=ifelse(
                UMI > 0,
                UMI,
                ""
            )
        ),
        size=3
    ) +
    scale_y_discrete(
        labels=label_map
    ) +
    scale_alpha_manual(
        values=c(
            `TRUE`=1,
            `FALSE`=0.35
        ),
        name="Cell Ranger accepted"
    ) +
    facet_grid(
        sample ~ .,
        scales="free_y",
        space="free_y"
    ) +
    labs(
        title="IL10 probe composition of IL10-audited cells",
        subtitle=paste0(
            "Corrected: one heatmap row per unique RDS cell (n=",
            nrow(heat),
            ")"
        ),
        x="IL10 probe",
        y=NULL,
        fill="UMI"
    ) +
    theme_bw(
        base_size=9
    ) +
    theme(
        legend.position="bottom",
        axis.text.y=element_text(
            size=7
        )
    )

heat_h <- max(
    8,
    0.36 * nrow(heat) + 3
)

ggsave(
    file.path(
        FIG,
        "Fig03_IL10_cell_probe_heatmap.pdf"
    ),
    p3,
    width=11.5,
    height=heat_h,
    limitsize=FALSE
)

ggsave(
    file.path(
        FIG,
        "Fig03_IL10_cell_probe_heatmap.png"
    ),
    p3,
    width=11.5,
    height=heat_h,
    dpi=300,
    limitsize=FALSE
)

# ------------------------------------------------------------
# Fig05: strict harmonized sensitivity analysis
# ------------------------------------------------------------

sample_summary$sample <- factor(
    sample_summary$sample,
    levels=SAMPLES
)

p5 <- ggplot(
    sample_summary,
    aes(
        x=sample,
        y=strict_harmonized_IL10_total_UMI
    )
) +
    geom_col() +
    geom_text(
        aes(
            label=paste0(
                "UMI=",
                strict_harmonized_IL10_total_UMI,
                "\n",
                "cells=",
                strict_harmonized_IL10_positive_cells
            )
        ),
        vjust=-0.25,
        size=3.4
    ) +
    labs(
        title="Strict harmonized IL10 sensitivity analysis",
        subtitle=paste0(
            "Shared Cell Ranger-accepted probe(s) across all samples: ",
            paste(
                shared_accepted_probes,
                collapse=", "
            )
        ),
        x=NULL,
        y="Harmonized IL10 UMI"
    ) +
    expand_limits(
        y=max(
            sample_summary$strict_harmonized_IL10_total_UMI
        ) + 1
    ) +
    theme_bw(
        base_size=11
    )

ggsave(
    file.path(
        FIG,
        "Fig05_IL10_strict_harmonized_shared_probe_by_sample.pdf"
    ),
    p5,
    width=8,
    height=5
)

ggsave(
    file.path(
        FIG,
        "Fig05_IL10_strict_harmonized_shared_probe_by_sample.png"
    ),
    p5,
    width=8,
    height=5,
    dpi=300
)

# ------------------------------------------------------------
# Fig06: strict harmonized signal by lineage
# ------------------------------------------------------------

lin_plot <- lineage_summary[
    lineage_summary$strict_harmonized_IL10_total_UMI > 0,
    ,
    drop=FALSE
]

if (nrow(lin_plot) > 0) {
    lin_plot$sample <- factor(
        lin_plot$sample,
        levels=SAMPLES
    )

    p6 <- ggplot(
        lin_plot,
        aes(
            x=lineage,
            y=strict_harmonized_IL10_total_UMI,
            fill=sample
        )
    ) +
        geom_col(
            position="dodge"
        ) +
        labs(
            title="Strict harmonized IL10 signal by lineage",
            subtitle=paste0(
                "Shared accepted probe(s): ",
                paste(
                    shared_accepted_probes,
                    collapse=", "
                )
            ),
            x=NULL,
            y="Harmonized IL10 UMI",
            fill="Sample"
        ) +
        theme_bw(
            base_size=10
        ) +
        theme(
            axis.text.x=element_text(
                angle=35,
                hjust=1
            ),
            legend.position="bottom"
        )

    ggsave(
        file.path(
            FIG,
            "Fig06_IL10_strict_harmonized_by_lineage.pdf"
        ),
        p6,
        width=9,
        height=5
    )

    ggsave(
        file.path(
            FIG,
            "Fig06_IL10_strict_harmonized_by_lineage.png"
        ),
        p6,
        width=9,
        height=5,
        dpi=300
    )
}

# ------------------------------------------------------------
# Audit report
# ------------------------------------------------------------

audit_path <- file.path(
    DST,
    "IL10_probe_integrated_summary_v1.4.2B.3.3.md"
)

con <- file(
    audit_path,
    open="wt"
)

writeLines(
    c(
        "# Human MASH v1.4.2B.3 IL10 probe audit",
        "",
        "Technical post-audit implementation: v1.4.2B.3.3.1",
        "",
        "## Source validation",
        "",
        "- Source: validated v1.4.2B.3.2 integrated IL10 probe audit.",
        "- Gene-level IL10 and sample-specific Cell Ranger accepted-probe IL10 were exactly concordant in SfLB1, SfLB2, Tx17, and Tx5.",
        "- Cellwise mismatches in the validated source were zero for all four samples.",
        "",
        "## Post-audit correction",
        "",
        "- Fig03 in v1.4.2B.3.2 used a shortened barcode as the factor identity.",
        "- Multiple Tx cells shared that shortened suffix and were collapsed visually.",
        "- v1.4.2B.3.3 uses the full RDS cell ID as the factor identity and a shortened barcode only as a display label.",
        paste0(
            "- Corrected Fig03 contains one row per unique IL10-audited RDS cell (n=",
            nrow(heat),
            ")."
        ),
        "",
        "## Strict harmonized sensitivity analysis",
        "",
        paste0(
            "- Shared Cell Ranger-accepted IL10 probe(s) across all four samples: ",
            paste(
                shared_accepted_probes,
                collapse=", "
            ),
            "."
        ),
        "- This sensitivity analysis avoids comparing sample-specific probe combinations across the v1.0.1 SfLB and v2.0.0 Tx probe panels.",
        ""
    ),
    con
)

for (
    i in seq_len(
        nrow(sample_summary)
    )
) {
    x <- sample_summary[
        i,
        ,
        drop=FALSE
    ]

    writeLines(
        paste0(
            "- ",
            as.character(x$sample),
            ": ",
            x$strict_harmonized_IL10_total_UMI,
            " UMI / ",
            x$strict_harmonized_IL10_positive_cells,
            " positive cells; max=",
            x$max_strict_harmonized_IL10_UMI,
            "."
        ),
        con
    )
}

writeLines(
    c(
        "",
        "## Interpretation boundary",
        "",
        "- The primary analysis remains the sample-specific Cell Ranger accepted-probe reconstruction because it exactly reproduces the gene-level IL10 matrix.",
        "- The strict harmonized analysis is a sensitivity analysis based only on probe(s) accepted in all four samples.",
        "- The IL10 probes are not transgene-specific and cannot distinguish endogenous IL10 from transgene-derived IL10 by sequence alone.",
        "- The data therefore support concordance with the mesenchymal IL10-transgene design but do not directly prove transcript origin from the transgene.",
        ""
    ),
    con
)

close(con)

cat("\n============================================================\n")
cat("STRICT HARMONIZED SAMPLE SUMMARY\n")
cat("============================================================\n")
print(
    sample_summary,
    row.names=FALSE
)

cat("\n============================================================\n")
cat("CORRECTED FIG03 AUDIT\n")
cat("============================================================\n")
cat(
    "Unique audited RDS cells: ",
    nrow(heat),
    "\n",
    sep=""
)
cat(
    "Unique factor identities: ",
    length(unique(heat$cell_key)),
    "\n",
    sep=""
)
cat(
    "Unique display labels: ",
    length(unique(heat$display_label)),
    "\n",
    sep=""
)

if (
    nrow(heat) !=
    length(unique(heat$cell_key))
) {
    stop(
        "Fig03 cell identity collision remains."
    )
}

if (
    nrow(heat) !=
    length(unique(heat$display_label))
) {
    stop(
        "Fig03 display-label collision remains."
    )
}

cat("\n============================================================\n")
cat("v1.4.2B.3.3.1 POST-AUDIT REPAIR COMPLETE\n")
cat("============================================================\n")
RS

cat > "$DST/VERSION.txt" <<EOF
Biological analysis version: v1.4.2B.3
Validated integrated implementation source: v1.4.2B.3.2
Post-audit technical implementation: v1.4.2B.3.3.1
Changes:
- corrected Fig03 unique-cell identity
- added strict harmonized shared-accepted-probe sensitivity analysis
- primary Cell Ranger accepted-probe reconstruction unchanged
EOF

cat > "$DST/POST_AUDIT_CHANGELOG.txt" <<EOF
v1.4.2B.3.3.1

1. No biological annotation, RDS, Cell Ranger counts, or probe acceptance statuses were changed.
2. v1.4.2B.3.2 remained the validated source for gene-level vs accepted-probe exact concordance.
3. Fig03 was repaired because shortened barcode suffixes visually collapsed distinct Tx cells.
4. Full RDS cell ID is now the heatmap row identity.
5. Strict harmonized sensitivity analysis was added using only IL10 probe(s) with filtered_probes=TRUE in all four samples.
EOF

echo
echo "============================================================"
echo "OUTPUT"
echo "============================================================"
echo "$DST"
echo
echo "Figures:"
ls -lh "$FIG"
echo
echo "New/updated tables:"
ls -lh \
    "$TAB/IL10_probe_acceptance_matrix.tsv" \
    "$TAB/IL10_cell_level_all_RDS_v3.3_with_strict_harmonized.tsv" \
    "$TAB/IL10_strict_harmonized_positive_cells.tsv" \
    "$TAB/IL10_strict_harmonized_sample_summary.tsv" \
    "$TAB/IL10_strict_harmonized_lineage_summary.tsv"
