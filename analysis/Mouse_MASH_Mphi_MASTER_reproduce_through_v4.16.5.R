#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH scRNA-seq macrophage analysis
# MASTER reproducibility pipeline through v4.16.5
#
# File:
#   Mouse_MASH_Mphi_MASTER_reproduce_through_v4.16.5.R
#
# PURPOSE
#   Reproduce the FINAL Mouse MASH macrophage analysis and manuscript figures
#   from the frozen Clean-B FINAL object through the latest v4.16.5 outputs.
#
# REPRODUCIBILITY CHECKPOINTS
#   v4.14.5  FINAL Clean-B annotation
#   v4.15.5  FINAL publication UMAP / summary figures
#   v4.16.1  Extended final figure collection
#   v4.16.2  Manuscript-gap completion analysis
#   v4.16.3  Program-level remodeling heatmap
#   v4.16.5  Manuscript visual refinement
#
# IMPORTANT
#   - This master script does NOT regenerate the v4.14.5 FINAL RDS from raw
#     matrices. It treats the frozen FINAL RDS as the reproducible parent input.
#   - The parent classification remains Res2.0.
#   - FINAL annotation column:
#         macrophage_class_Res2_FINAL_v4145
#   - Child scripts are executed in a fresh R process using Rscript so that
#     one script's global environment does not contaminate the next.
#   - By default, child scripts are read from:
#         /Users/uenoya/Projects/uenoyscRNA/analysis/
#
# GITHUB CHECKPOINT RECOMMENDATION
#   Tag:
#       mphi-cleanB-master-v4.16.5
# ==============================================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 0. User-adjustable configuration
# ------------------------------------------------------------------------------

PROJECT_DIR <- "/Users/uenoya/Projects/uenoyscRNA"

ANALYSIS_DIR <- file.path(
    PROJECT_DIR,
    "analysis"
)

DATA_ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

FINAL_RDS <- file.path(
    DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "RDS",
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

MASTER_OUTPUT_DIR <- file.path(
    DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "MASTER_Reproducibility_v4.16.5"
)

LOG_DIR <- file.path(
    MASTER_OUTPUT_DIR,
    "Logs"
)

AUDIT_DIR <- file.path(
    MASTER_OUTPUT_DIR,
    "Audit"
)

dir.create(
    LOG_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    AUDIT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

# Set TRUE to run each analysis.
RUN_V4155 <- TRUE
RUN_V4161 <- TRUE
RUN_V4162 <- TRUE
RUN_V4163 <- TRUE
RUN_V4165 <- TRUE

# Set TRUE if an error in one child script should stop the master pipeline.
STOP_ON_CHILD_ERROR <- TRUE

# ------------------------------------------------------------------------------
# 1. Child scripts
# ------------------------------------------------------------------------------

PIPELINE <- data.frame(
    order = c(
        1L,
        2L,
        3L,
        4L,
        5L
    ),

    version = c(
        "v4.15.5",
        "v4.16.1",
        "v4.16.2",
        "v4.16.3",
        "v4.16.5"
    ),

    purpose = c(
        "FINAL publication UMAP and summary figures",
        "Extended FINAL figure collection",
        "Manuscript-gap completion analysis",
        "Program-level functional remodeling heatmap",
        "Manuscript visual refinement: CLR and dense DotPlot"
    ),

    script = c(
        "Mouse_MASH_Mphi_FINAL_figures_v4.15.5_DENSER_DOTS.R",
        "Mouse_MASH_Mphi_FINAL_extended_figures_v4.16.1.R",
        "Mouse_MASH_Mphi_manuscript_gap_completion_v4.16.2.R",
        "Mouse_MASH_Mphi_program_remodeling_heatmap_v4.16.3.R",
        "Mouse_MASH_Mphi_manuscript_visual_refinement_v4.16.5.R"
    ),

    run = c(
        RUN_V4155,
        RUN_V4161,
        RUN_V4162,
        RUN_V4163,
        RUN_V4165
    ),

    stringsAsFactors = FALSE
)

PIPELINE$path <- file.path(
    ANALYSIS_DIR,
    PIPELINE$script
)

# ------------------------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------------------------

msg <- function(...) {

    message(
        "[",
        format(
            Sys.time(),
            "%Y-%m-%d %H:%M:%S"
        ),
        "] ",
        paste0(...)
    )
}


timestamp_id <- function() {

    format(
        Sys.time(),
        "%Y%m%d_%H%M%S"
    )
}


write_lines_safe <- function(
    x,
    path
) {

    writeLines(
        x,
        con = path,
        useBytes = TRUE
    )
}


file_md5_safe <- function(
    path
) {

    if (!file.exists(
        path
    )) {

        return(
            NA_character_
        )
    }

    unname(
        tools::md5sum(
            path
        )
    )
}


run_child_script <- function(
    script_path,
    version,
    purpose
) {

    run_id <- timestamp_id()

    log_path <- file.path(
        LOG_DIR,
        paste0(
            version,
            "_",
            run_id,
            ".log"
        )
    )

    msg(
        "Starting ",
        version,
        ": ",
        purpose
    )

    msg(
        "Script: ",
        script_path
    )

    rscript_bin <- file.path(
        R.home(
            "bin"
        ),
        "Rscript"
    )

    if (!file.exists(
        rscript_bin
    )) {

        rscript_bin <- Sys.which(
            "Rscript"
        )
    }

    if (!nzchar(
        rscript_bin
    )) {

        stop(
            "Rscript executable not found."
        )
    }

    start_time <- Sys.time()

    output <- tryCatch(
        system2(
            command = rscript_bin,
            args = shQuote(
                script_path
            ),
            stdout = TRUE,
            stderr = TRUE
        ),
        error = function(e) {
            structure(
                paste0(
                    "MASTER system2 ERROR: ",
                    conditionMessage(
                        e
                    )
                ),
                status = 999L
            )
        }
    )

    status <- attr(
        output,
        "status"
    )

    if (is.null(
        status
    )) {

        status <- 0L
    }

    end_time <- Sys.time()

    write_lines_safe(
        c(
            paste0(
                "MASTER CHILD RUN: ",
                version
            ),
            paste0(
                "Purpose: ",
                purpose
            ),
            paste0(
                "Script: ",
                script_path
            ),
            paste0(
                "Start: ",
                start_time
            ),
            paste0(
                "End: ",
                end_time
            ),
            paste0(
                "Elapsed_sec: ",
                round(
                    as.numeric(
                        difftime(
                            end_time,
                            start_time,
                            units = "secs"
                        )
                    ),
                    2
                )
            ),
            paste0(
                "Exit_status: ",
                status
            ),
            "",
            "----- CHILD OUTPUT -----",
            output
        ),
        log_path
    )

    msg(
        "Completed ",
        version,
        " with exit status ",
        status
    )

    list(
        version = version,
        script = script_path,
        purpose = purpose,
        start = start_time,
        end = end_time,
        elapsed_sec = as.numeric(
            difftime(
                end_time,
                start_time,
                units = "secs"
            )
        ),
        status = as.integer(
            status
        ),
        log = log_path
    )
}

# ------------------------------------------------------------------------------
# 3. Parent-input validation
# ------------------------------------------------------------------------------

msg(
    "Mouse MASH MΦ MASTER reproducibility pipeline v4.16.5"
)

msg(
    "Project: ",
    PROJECT_DIR
)

msg(
    "Analysis directory: ",
    ANALYSIS_DIR
)

msg(
    "FINAL RDS: ",
    FINAL_RDS
)

if (!dir.exists(
    PROJECT_DIR
)) {

    stop(
        "Project directory not found:\n",
        PROJECT_DIR
    )
}

if (!dir.exists(
    ANALYSIS_DIR
)) {

    stop(
        "Analysis directory not found:\n",
        ANALYSIS_DIR
    )
}

if (!file.exists(
    FINAL_RDS
)) {

    stop(
        "Frozen FINAL RDS not found:\n",
        FINAL_RDS
    )
}

parent_info <- file.info(
    FINAL_RDS
)

parent_audit <- data.frame(
    item = c(
        "FINAL_RDS",
        "FINAL_RDS_size_bytes",
        "FINAL_RDS_modified",
        "FINAL_RDS_md5",
        "FINAL_annotation_column",
        "Parent_resolution"
    ),

    value = c(
        FINAL_RDS,
        as.character(
            parent_info$size
        ),
        as.character(
            parent_info$mtime
        ),
        file_md5_safe(
            FINAL_RDS
        ),
        "macrophage_class_Res2_FINAL_v4145",
        "Res2.0"
    ),

    stringsAsFactors = FALSE
)

write.csv(
    parent_audit,
    file.path(
        AUDIT_DIR,
        "01_parent_FINAL_RDS_audit_v4.16.5.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Child-script existence / checksum audit
# ------------------------------------------------------------------------------

PIPELINE$exists <- file.exists(
    PIPELINE$path
)

PIPELINE$md5 <- vapply(
    PIPELINE$path,
    file_md5_safe,
    character(
        1
    )
)

write.csv(
    PIPELINE,
    file.path(
        AUDIT_DIR,
        "02_child_script_audit_v4.16.5.csv"
    ),
    row.names = FALSE
)

missing_required <- PIPELINE[
    PIPELINE$run &
        !PIPELINE$exists,
    ,
    drop = FALSE
]

if (nrow(
    missing_required
) > 0L) {

    stop(
        paste0(
            "Required child script(s) missing:\n",
            paste(
                missing_required$path,
                collapse = "\n"
            )
        )
    )
}

# ------------------------------------------------------------------------------
# 5. Environment audit
# ------------------------------------------------------------------------------

env_audit <- c(
    paste0(
        "MASTER version: v4.16.5"
    ),
    paste0(
        "Date: ",
        Sys.time()
    ),
    paste0(
        "R version: ",
        R.version.string
    ),
    paste0(
        "Platform: ",
        R.version$platform
    ),
    paste0(
        "Project directory: ",
        PROJECT_DIR
    ),
    paste0(
        "Analysis directory: ",
        ANALYSIS_DIR
    ),
    paste0(
        "FINAL RDS: ",
        FINAL_RDS
    )
)

write_lines_safe(
    env_audit,
    file.path(
        AUDIT_DIR,
        "03_environment_audit_v4.16.5.txt"
    )
)

# ------------------------------------------------------------------------------
# 6. Execute pipeline
# ------------------------------------------------------------------------------

run_results <- list()

selected_pipeline <- PIPELINE[
    PIPELINE$run,
    ,
    drop = FALSE
]

if (nrow(
    selected_pipeline
) == 0L) {

    stop(
        "No child analyses selected."
    )
}

master_start <- Sys.time()

for (i in seq_len(
    nrow(
        selected_pipeline
    )
)) {

    row_now <- selected_pipeline[
        i,
        ,
        drop = FALSE
    ]

    result_now <- run_child_script(
        script_path = row_now$path,
        version = row_now$version,
        purpose = row_now$purpose
    )

    run_results[[
        row_now$version
    ]] <- result_now

    if (
        result_now$status != 0L &&
        STOP_ON_CHILD_ERROR
    ) {

        stop(
            paste0(
                "Child script failed: ",
                row_now$version,
                "\nSee log:\n",
                result_now$log
            )
        )
    }
}

master_end <- Sys.time()

# ------------------------------------------------------------------------------
# 7. Run summary
# ------------------------------------------------------------------------------

run_summary <- do.call(
    rbind,
    lapply(
        run_results,
        function(x) {

            data.frame(
                version = x$version,
                purpose = x$purpose,
                script = x$script,
                start = as.character(
                    x$start
                ),
                end = as.character(
                    x$end
                ),
                elapsed_sec = x$elapsed_sec,
                exit_status = x$status,
                log = x$log,
                stringsAsFactors = FALSE
            )
        }
    )
)

write.csv(
    run_summary,
    file.path(
        MASTER_OUTPUT_DIR,
        "MASTER_RUN_SUMMARY_v4.16.5.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Expected output audit
# ------------------------------------------------------------------------------

expected_outputs <- data.frame(
    version = c(
        "v4.15.5",
        "v4.16.1",
        "v4.16.2",
        "v4.16.3",
        "v4.16.5"
    ),

    key_output = c(
        file.path(
            DATA_ROOT,
            "Mouse_MASH_Mphi_RDS",
            "Mphi_Res2_CleanB_FINAL_v4.14.5",
            "Final_Figures_v4.15.5",
            "Figures",
            "01_FINAL_Mphi_subtype_UMAP_v4.15.5.pdf"
        ),

        file.path(
            DATA_ROOT,
            "Mouse_MASH_Mphi_RDS",
            "Mphi_Res2_CleanB_FINAL_v4.14.5",
            "FINAL_Extended_Figures_v4.16.1",
            "03_M1_M2_programs",
            "03_M1_M2_marker_programs_sample_level_Sham_vs_Tx_v4.16.1.pdf"
        ),

        file.path(
            DATA_ROOT,
            "Mouse_MASH_Mphi_RDS",
            "Mphi_Res2_CleanB_FINAL_v4.14.5",
            "FINAL_Manuscript_Gap_Analysis_v4.16.2",
            "Figures",
            "03_FINAL_subtype_marker_gene_DotPlot_v4.16.2.pdf"
        ),

        file.path(
            DATA_ROOT,
            "Mouse_MASH_Mphi_RDS",
            "Mphi_Res2_CleanB_FINAL_v4.14.5",
            "FINAL_Program_Remodeling_v4.16.3",
            "Figures",
            "01_FINAL_Mphi_program_remodeling_heatmap_v4.16.3.pdf"
        ),

        file.path(
            DATA_ROOT,
            "Mouse_MASH_Mphi_RDS",
            "Mphi_Res2_CleanB_FINAL_v4.14.5",
            "FINAL_Manuscript_Visual_Refinement_v4.16.5",
            "Figures",
            "04_FINAL_subtype_marker_DotPlot_dense_unclipped_v4.16.5.pdf"
        )
    ),

    stringsAsFactors = FALSE
)

expected_outputs$exists <- file.exists(
    expected_outputs$key_output
)

expected_outputs$size_bytes <- ifelse(
    expected_outputs$exists,
    file.info(
        expected_outputs$key_output
    )$size,
    NA_real_
)

expected_outputs$md5 <- vapply(
    expected_outputs$key_output,
    file_md5_safe,
    character(
        1
    )
)

write.csv(
    expected_outputs,
    file.path(
        MASTER_OUTPUT_DIR,
        "MASTER_EXPECTED_OUTPUT_AUDIT_v4.16.5.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Git checkpoint audit
# ------------------------------------------------------------------------------

git_log <- tryCatch(
    system2(
        "git",
        c(
            "-C",
            shQuote(
                PROJECT_DIR
            ),
            "log",
            "--oneline",
            "--decorate",
            "-10"
        ),
        stdout = TRUE,
        stderr = TRUE
    ),
    error = function(e) {
        paste0(
            "git log failed: ",
            conditionMessage(
                e
            )
        )
    }
)

git_status <- tryCatch(
    system2(
        "git",
        c(
            "-C",
            shQuote(
                PROJECT_DIR
            ),
            "status",
            "--short"
        ),
        stdout = TRUE,
        stderr = TRUE
    ),
    error = function(e) {
        paste0(
            "git status failed: ",
            conditionMessage(
                e
            )
        )
    }
)

write_lines_safe(
    c(
        "----- git log -----",
        git_log,
        "",
        "----- git status --short -----",
        git_status
    ),
    file.path(
        AUDIT_DIR,
        "04_git_checkpoint_audit_v4.16.5.txt"
    )
)

# ------------------------------------------------------------------------------
# 10. MASTER README
# ------------------------------------------------------------------------------

readme <- c(
    "Mouse MASH MΦ MASTER reproducibility pipeline through v4.16.5",
    "",
    paste0(
        "Generated: ",
        Sys.time()
    ),
    "",
    "Frozen parent input:",
    paste0(
        "  ",
        FINAL_RDS
    ),
    "",
    "Parent annotation:",
    "  macrophage_class_Res2_FINAL_v4145",
    "",
    "Parent clustering:",
    "  Res2.0",
    "",
    "Pipeline order:",
    "  1. v4.15.5 FINAL UMAP / summary figures",
    "  2. v4.16.1 extended final figures",
    "  3. v4.16.2 manuscript-gap completion",
    "  4. v4.16.3 functional-program remodeling",
    "  5. v4.16.5 CLR / DotPlot visual refinement",
    "",
    "Important interpretation constraints:",
    "  STD vs CDAHFD: n=1 vs n=1; descriptive only.",
    "  Sham vs Tx: n=2 vs n=2 biological samples.",
    "  Cell-level measurements must not be treated as biological replicates.",
    "  CLR / Aitchison analyses are composition-aware.",
    "  v4.16.5 observed range bars are min-max, not SEM or 95% CI.",
    "",
    "Recommended Git checkpoint:",
    "  tag = mphi-cleanB-master-v4.16.5",
    "",
    paste0(
        "Master elapsed seconds: ",
        round(
            as.numeric(
                difftime(
                    master_end,
                    master_start,
                    units = "secs"
                )
            ),
            2
        )
    )
)

write_lines_safe(
    readme,
    file.path(
        MASTER_OUTPUT_DIR,
        "README_MASTER_REPRODUCIBILITY_v4.16.5.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        AUDIT_DIR,
        "05_MASTER_sessionInfo_v4.16.5.txt"
    )
)

# ------------------------------------------------------------------------------
# 11. Final status
# ------------------------------------------------------------------------------

all_child_success <- all(
    run_summary$exit_status ==
        0L
)

all_key_outputs_exist <- all(
    expected_outputs$exists
)

msg(
    "MASTER pipeline completed."
)

msg(
    "All child scripts successful: ",
    all_child_success
)

msg(
    "All key outputs present: ",
    all_key_outputs_exist
)

msg(
    "MASTER output: ",
    MASTER_OUTPUT_DIR
)

print(
    run_summary
)

print(
    expected_outputs
)

if (!all_child_success) {

    warning(
        "One or more child scripts returned non-zero exit status."
    )
}

if (!all_key_outputs_exist) {

    warning(
        "One or more expected key outputs are missing. Check child logs and output paths."
    )
}
