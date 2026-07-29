#!/usr/bin/env bash

# ============================================================
# uenoy scRNAseq Framework
# GitHub initial backup script
#
# Intended release:
#   v3.0.0
#
# Usage:
#   chmod +x scripts/01_publish_github_v3.0.sh
#
#   GITHUB_REPO_URL="git@github.com:YOUR_ACCOUNT/uenoyscRNA.git" \
#   ./scripts/01_publish_github_v3.0.sh
#
# HTTPS example:
#   GITHUB_REPO_URL="https://github.com/YOUR_ACCOUNT/uenoyscRNA.git" \
#   ./scripts/01_publish_github_v3.0.sh
# ============================================================

set -Eeuo pipefail

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

PROJECT_ROOT="/Users/uenoya/Projects/uenoyscRNA"

RELEASE_VERSION="v3.0.0"

COMMIT_MESSAGE="Release ${RELEASE_VERSION}: framework baseline before metadata engine v3.1"

DEFAULT_BRANCH="main"

REMOTE_NAME="origin"

# Environment variable supplied at execution time.
GITHUB_REPO_URL="${GITHUB_REPO_URL:-}"

# ------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------

log() {
  printf '\n[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
  printf '\nERROR: %s\n' "$*" >&2
  exit 1
}

command_exists() {
  command -v "$1" >/dev/null 2>&1
}

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------

command_exists git || die "git command was not found."

[[ -d "${PROJECT_ROOT}" ]] || die \
  "Project directory was not found: ${PROJECT_ROOT}"

[[ -n "${GITHUB_REPO_URL}" ]] || die \
  "GITHUB_REPO_URL is not set.

Example:
GITHUB_REPO_URL=\"git@github.com:YOUR_ACCOUNT/uenoyscRNA.git\" \\
  \"$0\""

cd "${PROJECT_ROOT}"

log "Project root: ${PROJECT_ROOT}"
log "GitHub repository: ${GITHUB_REPO_URL}"
log "Release tag: ${RELEASE_VERSION}"

# ------------------------------------------------------------
# Create or update .gitignore
# ------------------------------------------------------------

log "Creating .gitignore"

cat > .gitignore <<'EOF'
# ============================================================
# macOS
# ============================================================
.DS_Store
.AppleDouble
.LSOverride
._*

# ============================================================
# R / RStudio
# ============================================================
.Rhistory
.RData
.Ruserdata
.Rapp.history
.Rproj.user/
*.Rproj.user
.Renviron
.Rprofile.local

# Local R package libraries and caches
renv/library/
renv/staging/
renv/python/
renv/cache/
.Rcache/
.cache/

# Keep reproducibility files
!renv.lock
!renv/activate.R
!renv/settings.json

# ============================================================
# Local configuration and credentials
# ============================================================
config/local_config.R
config/*_local.R
config/secrets.R
config/credentials.R
.env
.env.*
*.pem
*.key
*.p12
*.pfx
credentials.json
token.json
secrets.json

# Example/template local configuration may be committed
!config/local_config.example.R
!config/example_local_config.R

# ============================================================
# scRNA-seq data and large binary objects
# ============================================================
*.rds
*.RDS
*.rda
*.RData
*.h5
*.h5ad
*.h5seurat
*.loom
*.mtx
*.mtx.gz
*.bam
*.bai
*.cram
*.crai
*.fastq
*.fastq.gz
*.fq
*.fq.gz
*.tar
*.tar.gz
*.tgz
*.zip
*.7z

# Matrix directories and raw data
filtered_feature_bc_matrix/
raw_feature_bc_matrix/
outs/
data/
raw_data/
external_data/
input_data/

# ============================================================
# Analysis outputs
# ============================================================
results/
analysis_results/
output/
outputs/
figures/
plots/
tables/
reports/
logs/
tmp/
temp/
cache/

# Common generated formats
*.pdf
*.png
*.jpg
*.jpeg
*.tiff
*.tif
*.svg
*.html
*.csv
*.tsv
*.xlsx
*.xls
*.parquet
*.feather

# Permit intentionally versioned documentation assets
!README/**/*.png
!README/**/*.jpg
!README/**/*.jpeg
!README/**/*.svg
!docs/**/*.png
!docs/**/*.jpg
!docs/**/*.jpeg
!docs/**/*.svg

# ============================================================
# IDE and editor files
# ============================================================
.vscode/
.idea/
*.swp
*.swo
*~
EOF

# ------------------------------------------------------------
# Create local configuration template when possible
# ------------------------------------------------------------

if [[ -f "config/local_config.R" ]] &&
   [[ ! -f "config/local_config.example.R" ]]; then

  log "Creating config/local_config.example.R"

  cat > config/local_config.example.R <<'EOF'
# ============================================================
# Local configuration template
#
# Copy this file to:
#   config/local_config.R
#
# local_config.R is excluded from Git.
# Do not store passwords, access tokens, or private paths here
# unless the file remains local.
# ============================================================

overwrite_existing <- FALSE

# Example:
# data_root <- "/Volumes/YOUR_EXTERNAL_DISK/scRNA_data"
EOF
fi

# ------------------------------------------------------------
# Inspect files that should not be committed
# ------------------------------------------------------------

log "Checking for large files and sensitive filenames"

LARGE_FILES="$(
  find . \
    -type f \
    -size +50M \
    -not -path './.git/*' \
    -print 2>/dev/null || true
)"

if [[ -n "${LARGE_FILES}" ]]; then
  printf '\nFiles larger than 50 MB were found:\n%s\n' "${LARGE_FILES}"
  printf '\nThese files should normally remain excluded from Git.\n'
fi

SENSITIVE_FILES="$(
  find . \
    -type f \
    \( \
      -iname '*password*' -o \
      -iname '*passwd*' -o \
      -iname '*secret*' -o \
      -iname '*credential*' -o \
      -iname '*token*' -o \
      -iname '*.pem' -o \
      -iname '*.key' \
    \) \
    -not -path './.git/*' \
    -print 2>/dev/null || true
)"

if [[ -n "${SENSITIVE_FILES}" ]]; then
  printf '\nPotentially sensitive files were found:\n%s\n' \
    "${SENSITIVE_FILES}"
  printf '\nReview these files before continuing.\n'
  read -r -p "Continue? [y/N]: " ANSWER
  case "${ANSWER}" in
    y|Y|yes|YES)
      ;;
    *)
      die "Stopped before Git commit."
      ;;
  esac
fi

# ------------------------------------------------------------
# Initialize Git repository
# ------------------------------------------------------------

if [[ ! -d ".git" ]]; then
  log "Initializing Git repository"
  git init
fi

# Configure default branch.
CURRENT_BRANCH="$(git branch --show-current 2>/dev/null || true)"

if [[ -z "${CURRENT_BRANCH}" ]]; then
  git checkout -b "${DEFAULT_BRANCH}"
elif [[ "${CURRENT_BRANCH}" != "${DEFAULT_BRANCH}" ]]; then
  log "Renaming current branch to ${DEFAULT_BRANCH}"
  git branch -M "${DEFAULT_BRANCH}"
fi

# ------------------------------------------------------------
# Configure remote
# ------------------------------------------------------------

if git remote get-url "${REMOTE_NAME}" >/dev/null 2>&1; then
  EXISTING_REMOTE="$(git remote get-url "${REMOTE_NAME}")"

  if [[ "${EXISTING_REMOTE}" != "${GITHUB_REPO_URL}" ]]; then
    log "Updating existing remote URL"
    git remote set-url "${REMOTE_NAME}" "${GITHUB_REPO_URL}"
  fi
else
  log "Adding GitHub remote"
  git remote add "${REMOTE_NAME}" "${GITHUB_REPO_URL}"
fi

# ------------------------------------------------------------
# Show ignored and staged status
# ------------------------------------------------------------

log "Checking Git status"

git status --short

log "Adding version-controlled files"

git add .gitignore
git add --all

# Explicitly remove local-only files from the index, if previously tracked.
git rm --cached --ignore-unmatch \
  config/local_config.R \
  .Renviron \
  .env \
  credentials.json \
  token.json \
  secrets.json \
  >/dev/null 2>&1 || true

# ------------------------------------------------------------
# Safety check: prohibit tracked large files
# ------------------------------------------------------------

TRACKED_LARGE_FILES="$(
  git diff --cached --name-only --diff-filter=ACM |
  while IFS= read -r FILE; do
    [[ -f "${FILE}" ]] || continue
    SIZE="$(stat -f '%z' "${FILE}" 2>/dev/null || stat -c '%s' "${FILE}")"
    if [[ "${SIZE}" -gt 50000000 ]]; then
      printf '%s\n' "${FILE}"
    fi
  done
)"

if [[ -n "${TRACKED_LARGE_FILES}" ]]; then
  printf '\nThe following staged files exceed 50 MB:\n%s\n' \
    "${TRACKED_LARGE_FILES}"
  die "Remove them from Git or configure Git LFS before continuing."
fi

log "Files prepared for commit"

git status --short

printf '\nReview the files listed above.\n'
read -r -p "Create the ${RELEASE_VERSION} commit? [y/N]: " ANSWER

case "${ANSWER}" in
  y|Y|yes|YES)
    ;;
  *)
    die "Stopped before commit."
    ;;
esac

# ------------------------------------------------------------
# Commit
# ------------------------------------------------------------

if git diff --cached --quiet; then
  log "No staged changes were found"
else
  log "Creating commit"
  git commit -m "${COMMIT_MESSAGE}"
fi

# ------------------------------------------------------------
# Create release tag
# ------------------------------------------------------------

if git rev-parse "${RELEASE_VERSION}" >/dev/null 2>&1; then
  log "Tag already exists: ${RELEASE_VERSION}"
else
  log "Creating annotated tag: ${RELEASE_VERSION}"
  git tag -a "${RELEASE_VERSION}" \
    -m "uenoy scRNAseq Framework ${RELEASE_VERSION}"
fi

# ------------------------------------------------------------
# Push
# ------------------------------------------------------------

log "Pushing branch to GitHub"

git push -u "${REMOTE_NAME}" "${DEFAULT_BRANCH}"

log "Pushing release tag"

git push "${REMOTE_NAME}" "${RELEASE_VERSION}"

# ------------------------------------------------------------
# Create development branches
# ------------------------------------------------------------

if git show-ref --verify --quiet refs/heads/develop; then
  log "Local develop branch already exists"
else
  log "Creating develop branch"
  git branch develop "${DEFAULT_BRANCH}"
fi

if git show-ref --verify --quiet \
  refs/heads/feature/metadata-engine; then
  log "Local feature/metadata-engine branch already exists"
else
  log "Creating feature/metadata-engine branch"
  git branch feature/metadata-engine develop
fi

git push -u "${REMOTE_NAME}" develop
git push -u "${REMOTE_NAME}" feature/metadata-engine

log "GitHub backup completed"

printf '\nRepository: %s\n' "${GITHUB_REPO_URL}"
printf 'Stable branch: %s\n' "${DEFAULT_BRANCH}"
printf 'Development branch: develop\n'
printf 'v3.1 work branch: feature/metadata-engine\n'
printf 'Release tag: %s\n' "${RELEASE_VERSION}"

