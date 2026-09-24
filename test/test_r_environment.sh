#!/usr/bin/env bash
set -euo pipefail

readonly PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
readonly TEMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/dhole-r-environment.XXXXXX")"
trap 'rm -rf -- "$TEMP_DIR"' EXIT

mkdir -p "$TEMP_DIR/foreign library"
printf 'stop("A foreign R profile was loaded")\n' > "$TEMP_DIR/Rprofile"
printf 'R_LIBS_USER="%s/foreign library"\n' "$TEMP_DIR" > "$TEMP_DIR/Renviron"

R_LIBS="$TEMP_DIR/foreign library" \
R_LIBS_USER="$TEMP_DIR/foreign library" \
R_LIBS_SITE="$TEMP_DIR/foreign library" \
R_PROFILE="$TEMP_DIR/Rprofile" \
R_PROFILE_USER="$TEMP_DIR/Rprofile" \
R_ENVIRON="$TEMP_DIR/Renviron" \
R_ENVIRON_USER="$TEMP_DIR/Renviron" \
bash "$PROJECT_DIR/tools/run-r" - "argument with spaces" <<'RSCRIPT'
stopifnot(identical(commandArgs(trailingOnly = TRUE), "argument with spaces"))
stopifnot(identical(.libPaths(), normalizePath(.Library)))
stopifnot(startsWith(normalizePath(R.home()), paste0(Sys.getenv("CONDA_PREFIX"), "/")))
invisible(loadNamespace("rlang"))
invisible(loadNamespace("dplyr"))
stopifnot(identical(normalizePath(dirname(find.package("rlang"))), normalizePath(.Library)))

# Child R processes must inherit the same isolation (e.g. package installers).
status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c("--vanilla", "-e", shQuote("stopifnot(identical(.libPaths(), normalizePath(.Library)))"))
)
stopifnot(status == 0L)
cat("TEST PASSED: Conda R ignores external libraries and startup files.\n")
RSCRIPT

status=0
bash "$PROJECT_DIR/tools/run-r" -e 'quit(status = 23L)' >/dev/null 2>&1 || status=$?
[[ "$status" -eq 23 ]] || {
  printf 'TEST FAILED: run-r did not preserve the R exit status\n' >&2
  exit 1
}
