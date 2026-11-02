#!/usr/bin/env bash
set -euo pipefail

readonly SOURCE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
readonly TEMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/dhole-install-log.XXXXXX")"
trap 'rm -rf -- "$TEMP_DIR"' EXIT

# Load only the function under test; never start Conda or download packages.
source <(sed -n '/^install_r_dependencies() {$/,/^verify_primer_qc_dependencies() {$/p' \
  "$SOURCE_DIR/install.sh" | sed '$d')

readonly R_RUNNER="$TEMP_DIR/mock-run-r"
cat > "$R_RUNNER" <<'MOCK'
set -eu
[[ "$PKG_SYSREQS" == false && "$USE_BUNDLED_LIBUV" == 1 ]]
[[ "$PKGCACHE_CONNECTTIMEOUT" == 15 && "$PKGCACHE_TIMEOUT" == "$EXPECTED_TIMEOUT" ]]
[[ "$PKGCACHE_LOW_SPEED_LIMIT" == 1 && "$PKGCACHE_LOW_SPEED_TIME" == 60 ]]
[[ "$PKGCACHE_RETRY" == false && "$PKG_HTTP_RETRY" == false ]]
[[ -s "$1" ]]
printf 'mock R stdout\n'
printf 'mock R stderr\n' >&2
exit "$MOCK_R_STATUS"
MOCK

unset PKGCACHE_CONNECTTIMEOUT PKGCACHE_TIMEOUT PKGCACHE_LOW_SPEED_LIMIT
unset PKGCACHE_LOW_SPEED_TIME PKGCACHE_RETRY PKG_HTTP_RETRY
export EXPECTED_TIMEOUT=180

for expected in 0 23; do
  PROJECT_DIR="$TEMP_DIR/project with spaces $expected"
  export MOCK_R_STATUS=$expected
  status=0
  install_r_dependencies > "$TEMP_DIR/console.log" 2>&1 || status=$?
  [[ "$status" -eq "$expected" ]]
  logs=("$PROJECT_DIR"/local/install-logs/r-dependencies.*/install.log)
  [[ ${#logs[@]} -eq 1 && -s "${logs[0]}" ]]
  [[ -s "${logs[0]%install.log}install.R" ]]
  for output in 'mock R stdout' 'mock R stderr' "exit status: $expected"; do
    grep -Fq "$output" "${logs[0]}"
    grep -Fq "$output" "$TEMP_DIR/console.log"
  done
done

PROJECT_DIR="$TEMP_DIR/custom timeout"
export PKGCACHE_TIMEOUT=600 EXPECTED_TIMEOUT=600 MOCK_R_STATUS=0
install_r_dependencies > "$TEMP_DIR/console.log" 2>&1

printf 'TEST PASSED: R installation logs, exit status, network defaults and timeout override.\n'
