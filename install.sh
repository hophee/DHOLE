#!/usr/bin/env bash
set -uo pipefail

readonly PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly ENV_FILE="$PROJECT_DIR/env.yml"
readonly CONDA_ENV_NAME="oligo_design"

readonly CHOPCHOP_ENV_FILE="$PROJECT_DIR/env_chopchop.yml"
readonly CHOPCHOP_ENV_NAME="oligo_design_chopchop"

readonly VIENNARNA_ENV_FILE="$PROJECT_DIR/env_viennarna.yml"
readonly VIENNARNA_ENV_NAME="oligo_design_viennarna"


readonly TEST_DIR="$PROJECT_DIR/test"
readonly MELTING_WRAPPER="$PROJECT_DIR/tools/melting-batch"
readonly R_RUNNER="$PROJECT_DIR/tools/run-r"

readonly MG1655_REFSEQ_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"
readonly MG1655_GENOME_URL="$MG1655_REFSEQ_URL/GCF_000005845.2_ASM584v2_genomic.fna.gz"
readonly MG1655_GFF_URL="$MG1655_REFSEQ_URL/GCF_000005845.2_ASM584v2_genomic.gff.gz"

declare -a FAILED_COMPONENTS=()


run_step() {
  local component=$1
  shift

  if "$@"; then
    printf 'OK: %s\n' "$component"
  else
    printf 'FAILED: %s\n' "$component" >&2
    FAILED_COMPONENTS+=("$component")
  fi
}


conda_environment_prefix() {
  local environment_name=$1

  conda env list |
    awk -v name="$environment_name" '$1 == name { print $NF; exit }'
}


clone_if_missing() {
  local repository=$1
  local destination=$2

  if [[ -e "$destination" && ! -d "$destination/.git" ]]; then
    printf 'ERROR: %s already exists but is not a Git checkout.\n' \
      "$destination" >&2
    return 1
  fi

  if [[ ! -e "$destination" ]]; then
    GIT_TERMINAL_PROMPT=0 \
      git clone --depth 1 "$repository" "$destination" || return 1
  fi
}


install_conda_environment() {
  local environment_name=$1
  local environment_file=$2

  [[ -f "$environment_file" ]] || {
    printf 'ERROR: Conda environment file not found: %s\n' \
      "$environment_file" >&2
    return 1
  }

  if conda env list |
    awk 'NR > 2 { print $1 }' |
    grep -Fqx "$environment_name"; then

    conda env update \
      --name "$environment_name" \
      --file "$environment_file" \
      --prune
  else
    conda env create \
      --name "$environment_name" \
      --file "$environment_file" \
      --yes
  fi
}


install_melting_wrapper() {
  local environment_prefix

  environment_prefix="$(
    conda_environment_prefix "$CONDA_ENV_NAME"
  )"

  [[ -n "$environment_prefix" && -d "$environment_prefix/bin" ]] || {
    printf 'ERROR: cannot find Conda environment prefix for %s.\n' \
      "$CONDA_ENV_NAME" >&2
    return 1
  }

  [[ -f "$MELTING_WRAPPER" ]] || {
    printf 'ERROR: MELTING wrapper not found: %s\n' \
      "$MELTING_WRAPPER" >&2
    return 1
  }

  [[ -x "$MELTING_WRAPPER" ]] ||
    chmod +x "$MELTING_WRAPPER" ||
    return 1

  ln -sfn \
    "$MELTING_WRAPPER" \
    "$environment_prefix/bin/melting-batch"
}


install_viennarna_commands() {
  local main_prefix
  local viennarna_prefix
  local command_name

  main_prefix="$(
    conda_environment_prefix "$CONDA_ENV_NAME"
  )"

  viennarna_prefix="$(
    conda_environment_prefix "$VIENNARNA_ENV_NAME"
  )"

  [[ -n "$main_prefix" &&
     -n "$viennarna_prefix" &&
     -d "$main_prefix/bin" ]] || {
    printf 'ERROR: cannot locate the main or ViennaRNA environment.\n' >&2
    return 1
  }

  for command_name in RNAfold RNAduplex RNAsubopt; do
    [[ -x "$viennarna_prefix/bin/$command_name" ]] || {
      printf 'ERROR: ViennaRNA command not found: %s\n' \
        "$command_name" >&2
      return 1
    }

    ln -sfn \
      "$viennarna_prefix/bin/$command_name" \
      "$main_prefix/bin/$command_name" ||
      return 1
  done
}


install_chopchop_python() {
  local main_prefix
  local chopchop_prefix

  main_prefix="$(
    conda_environment_prefix "$CONDA_ENV_NAME"
  )"

  chopchop_prefix="$(
    conda_environment_prefix "$CHOPCHOP_ENV_NAME"
  )"

  [[ -n "$main_prefix" &&
     -n "$chopchop_prefix" &&
     -d "$main_prefix/bin" &&
     -x "$chopchop_prefix/bin/python" ]] || {
    printf 'ERROR: cannot locate the isolated CHOPCHOP Python 2 runtime.\n' >&2
    return 1
  }

  ln -sfn \
    "$chopchop_prefix/bin/python" \
    "$main_prefix/bin/chopchop-python"
}


install_r_dependencies() {
  local temporary_script

  temporary_script="$(
    mktemp "${TMPDIR:-/tmp}/install-r-dependencies.XXXXXX.R"
  )" || return 1

  cat > "$temporary_script" <<'RSCRIPT'
options(
  pkg.sysreqs = FALSE,
  pkg.use_bioconductor = TRUE,

  # Verbose terminal progress
  cli.progress_show_after = 0,
  cli.progress_clear = FALSE,
  rlib_interactive = TRUE
)

if (!requireNamespace("pak", quietly = TRUE)) {
  stop("pak is not installed", call. = FALSE)
}

packages <- c(
  "cran::argparser",
  "cran::dplyr",
  "cran::readr",
  "cran::ape",
  "cran::janitor",
  "cran::digest",
  "cran::foreach",
  "cran::doParallel",

  "bioc::Biostrings@2.78.0",
  "bioc::DECIPHER@3.6.0",
  "bioc::openPrimeR@1.32.0",
  "bioc::rmelting@1.26.0"
)

target_library <- normalizePath(.Library, mustWork = TRUE)
.libPaths(target_library)

message("R version: ", R.version.string)
message("Conda R library: ", target_library)
message("Installing R dependencies with pak...")
message(
  "pak may report missing OS packages even with PKG_SYSREQS=false. ",
  "Java is supplied by Conda; the final runtime check verifies the JVM."
)

result <- pak::pkg_install(
  packages,
  lib = target_library,
  upgrade = FALSE,
  ask = FALSE,
  dependencies = NA
)

message("pak installation completed.")

print(result)
RSCRIPT

  PKG_SYSREQS=false \
  USE_BUNDLED_LIBUV=1 \
  bash "$R_RUNNER" "$temporary_script"

  local status=$?
  rm -f "$temporary_script"
  return "$status"
}


verify_primer_qc_dependencies() {
  local project_dir=$PROJECT_DIR

  TWO_PAC_PROJECT_DIR="$project_dir" \
  bash "$R_RUNNER" -e '

    #
    # OligoArrayAux data files
    #
    oligoarray_data <- file.path(
      Sys.getenv("CONDA_PREFIX"),
      "share",
      "oligoarrayaux"
    )

    if (file.exists(
      file.path(oligoarray_data, "stack.DAT")
    )) {
      Sys.setenv(
        UNAFOLDDAT = oligoarray_data
      )
    }


    target_library <- normalizePath(.Library, mustWork = TRUE)
    .libPaths(target_library)
    packages <- c(
      "argparser",
      "dplyr",
      "readr",
      "ape",
      "janitor",
      "digest",
      "foreach",
      "doParallel",
      "Biostrings",
      "DECIPHER",
      "openPrimeR",
      "rJava",
      "rmelting"
    )

    for (package_name in packages) {
      message(
        "Checking R package: ",
        package_name
      )

      package_path <- tryCatch(
        find.package(package_name, lib.loc = target_library),

        error = function(e) {
          stop(
            sprintf(
              "Cannot load R package %s: %s",
              package_name,
              conditionMessage(e)
            ),
            call. = FALSE
          )
        }
      )

      if (normalizePath(dirname(package_path)) != target_library) {
        stop(
          sprintf(
            "R package %s is outside the Conda library: %s",
            package_name,
            package_path
          ),
          call. = FALSE
        )
      }

      loadNamespace(package_name)
    }

    message("Checking Java runtime through rJava...")
    if (rJava::.jinit() < 0L) {
      stop("Cannot initialize the Java VM through rJava", call. = FALSE)
    }
    message("Java VM initialized successfully.")

    expected_versions <- c(
      Biostrings = "2.78.0",
      DECIPHER = "3.6.0",
      openPrimeR = "1.32.0",
      rmelting = "1.26.0"
    )
    installed_versions <- vapply(
      names(expected_versions),
      function(package_name) {
        as.character(packageVersion(package_name, lib.loc = target_library))
      },
      character(1)
    )
    if (!identical(unname(installed_versions), unname(expected_versions))) {
      stop(
        sprintf(
          "Unexpected R package versions: %s",
          paste(names(installed_versions), installed_versions, collapse = ", ")
        ),
        call. = FALSE
      )
    }


    #
    # External commands required by primer QC.
    #
    commands <- c(
      "hybrid-min",
      "RNAfold",
      "RNAduplex",
      "RNAsubopt",
      "melting-batch",
      "java",
      "javac",
      "chopchop-python"
    )

    command_paths <- Sys.which(commands)

    missing_commands <- names(
      command_paths
    )[!nzchar(command_paths)]

    if (length(missing_commands) > 0L) {
      stop(
        sprintf(
          "Required commands are missing from PATH: %s",
          paste(
            missing_commands,
            collapse = ", "
          )
        ),
        call. = FALSE
      )
    }


    project_dir <- Sys.getenv("TWO_PAC_PROJECT_DIR", unset = NA_character_)
    if (is.na(project_dir) || !dir.exists(project_dir)) {
      stop("TWO_PAC_PROJECT_DIR is invalid", call. = FALSE)
    }
    setwd(project_dir)
    source("oligo_designer.R")
    configure_openprimer_environment()
    load_openprimer_settings()
  '
}


validate_call_primer3() {
  [[ -s "$PROJECT_DIR/callPrimer3.R" ]] || {
    printf 'ERROR: tracked callPrimer3.R is missing; restore it from the project source.\n' >&2
    return 1
  }
}


download_gzip_if_missing() {
  local url=$1
  local destination=$2
  local temporary_file

  [[ -s "$destination" ]] && return 0

  temporary_file="$(
    mktemp "$TEST_DIR/.download.XXXXXX.gz"
  )" || return 1


  if ! curl \
    --fail \
    --location \
    --silent \
    --show-error \
    "$url" \
    --output "$temporary_file"; then

    rm -f "$temporary_file"
    return 1
  fi


  if ! gzip -t "$temporary_file" ||
     ! gzip -dc "$temporary_file" > "$destination"; then

    rm -f \
      "$temporary_file" \
      "$destination"

    return 1
  fi


  rm -f "$temporary_file"
}


prepare_test_data() {
  mkdir -p "$TEST_DIR" || return 1

  download_gzip_if_missing \
    "$MG1655_GENOME_URL" \
    "$TEST_DIR/MG1655.fna" ||
    return 1

  download_gzip_if_missing \
    "$MG1655_GFF_URL" \
    "$TEST_DIR/MG1655.gff" ||
    return 1


  if [[ ! -f "$TEST_DIR/test_run.sh" ]]; then
    printf 'ERROR: test runner was not found: %s\n' \
      "$TEST_DIR/test_run.sh" >&2
    return 1
  fi

  chmod +x "$TEST_DIR/test_run.sh"
}


patch_chopchop_max_offtargets_type() {
  local chopchop_file="$PROJECT_DIR/chopchop/chopchop.py"

  local old_argument='parser.add_argument("-m", "--maxOffTargets", metavar="MAX_HITS", help="The maximum number of off targets allowed.")'
  local new_argument='parser.add_argument("-m", "--maxOffTargets", type=int, metavar="MAX_HITS", help="The maximum number of off targets allowed.")'


  [[ -f "$chopchop_file" ]] || {
    printf 'ERROR: CHOPCHOP source file not found: %s\n' \
      "$chopchop_file" >&2
    return 1
  }


  if grep -Fq \
    "$new_argument" \
    "$chopchop_file"; then

    return 0
  fi


  if ! grep -Fq \
    "$old_argument" \
    "$chopchop_file"; then

    printf 'ERROR: CHOPCHOP maxOffTargets argument was not found.\n' >&2
    return 1
  fi


  sed -i \
    "s/${old_argument}/${new_argument}/" \
    "$chopchop_file"
}


cd "$PROJECT_DIR" || exit 1


#
# Conda environments
#
if command -v conda >/dev/null 2>&1; then

  run_step \
    "Conda environment $CONDA_ENV_NAME" \
    install_conda_environment \
    "$CONDA_ENV_NAME" \
    "$ENV_FILE"


  run_step \
    "Conda environment $CHOPCHOP_ENV_NAME" \
    install_conda_environment \
    "$CHOPCHOP_ENV_NAME" \
    "$CHOPCHOP_ENV_FILE"


  run_step \
    "Conda environment $VIENNARNA_ENV_NAME" \
    install_conda_environment \
    "$VIENNARNA_ENV_NAME" \
    "$VIENNARNA_ENV_FILE"


  #
  # Make isolated tools visible from the main environment.
  #
  run_step \
    "MELTING command wrapper" \
    install_melting_wrapper


  run_step \
    "isolated CHOPCHOP Python 2" \
    install_chopchop_python


  run_step \
    "isolated ViennaRNA commands" \
    install_viennarna_commands


  #
  # Install CRAN + Bioconductor packages with pak.
  #
  run_step \
    "R dependencies via pak" \
    install_r_dependencies


else

  printf 'FAILED: Conda environment %s (conda was not found in PATH)\n' \
    "$CONDA_ENV_NAME" >&2

  FAILED_COMPONENTS+=(
    "Conda environment $CONDA_ENV_NAME"
  )
fi


#
# primer3
#
run_step \
  "primer3 source" \
  clone_if_missing \
  https://github.com/primer3-org/primer3.git \
  primer3


run_step \
  "primer3 build" \
  make -C primer3/src


#
# CHOPCHOP
#
run_step \
  "CHOPCHOP source" \
  clone_if_missing \
  https://github.com/JokingHero/chopchop.git \
  chopchop


run_step \
  "CHOPCHOP maxOffTargets integer compatibility" \
  patch_chopchop_max_offtargets_type


#
# Primer3 R wrapper
#
run_step \
  "callPrimer3.R" \
  validate_call_primer3


#
# Test data
#
run_step \
  "MG1655 test data and test runner" \
  prepare_test_data


#
# Runtime verification after all project tools are present.
#
if command -v conda >/dev/null 2>&1; then
  run_step \
    "Primer QC dependencies" \
    verify_primer_qc_dependencies
fi


#
# Final result
#
if (( ${#FAILED_COMPONENTS[@]} == 0 )); then
  printf 'Installation completed successfully.\n'
else
  printf '\nInstallation completed with missing components:\n' >&2

  printf '  - %s\n' \
    "${FAILED_COMPONENTS[@]}" >&2

  exit 1
fi
