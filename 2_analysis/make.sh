#!/bin/bash   

# Trap to handle shell script errors 
trap 'error_handler' ERR
error_handler() {
    error_time=$(date '+%Y-%m-%d %H:%M:%S')
    echo -e "\n\033[0;31mWarning\033[0m: make.sh failed at ${error_time}. Check above for details." # display warning in terminal
    exit 1 # early exit with error code
}

# Set paths
# (Make sure REPO_ROOT is set to point to the root of the repository!)
MAKE_SCRIPT_DIR="$(cd "$(dirname -- "$0")" && pwd -P)"
REPO_ROOT="$(git rev-parse --show-toplevel)"
MODULE=$(basename "$MAKE_SCRIPT_DIR")
LOGFILE="${MAKE_SCRIPT_DIR}/output/make.log"

# Check setup
source "${REPO_ROOT}/lib/shell/check_setup.sh"

# Tell user what we're doing
echo -e "\n\nMaking module \033[35m${MODULE}\033[0m with shell ${SHELL}"

# Load settings & tools
source "${REPO_ROOT}/local_env.sh"
source "${REPO_ROOT}/lib/shell/run_python.sh"
source "${REPO_ROOT}/lib/shell/run_stata.sh"
source "${REPO_ROOT}/lib/shell/run_R.sh"

# Clear output directory
# (Guarantees that all output is produced from a clean run of the code)
rm -rf "${MAKE_SCRIPT_DIR}/output"
mkdir -p "${MAKE_SCRIPT_DIR}/output"

# Add symlink input files to local /input/ directory
# (Make sure get_inputs.sh is updated to pull in all needed input files!)
(   cd "${MAKE_SCRIPT_DIR}"
    source "${MAKE_SCRIPT_DIR}/get_inputs.sh"
)

# Run scripts
# (Do this in a subshell so we return to the original working directory
# after scripts are run)
 echo -e "\nmake.sh started at $(date '+%Y-%m-%d %H:%M:%S')"

(
cd "${MAKE_SCRIPT_DIR}/source"

# Phase 1: Core estimation (produces markups consumed by later scripts)
run_python paper_results.py "${LOGFILE}" || exit 1
run_python premium_timeseries.py "${LOGFILE}" || exit 1
run_python summary_stats.py "${LOGFILE}" || exit 1
run_python specification_sensitivity_table.py "${LOGFILE}" || exit 1

# Phase 2: Analysis scripts
run_python bmy_czech_analysis.py "${LOGFILE}" || exit 1
run_python dls_markup_comparison.py "${LOGFILE}" || exit 1
run_python dls_table2_replication.py "${LOGFILE}" || exit 1
run_python dlw_treatment_eval.py "${LOGFILE}" || exit 1
run_python cwdl_robustness.py "${LOGFILE}" || exit 1
run_python raval_test.py "${LOGFILE}" || exit 1
run_python klms_analysis.py "${LOGFILE}" || exit 1
run_python panel_treatment.py "${LOGFILE}" || exit 1
run_R fect_estimation.R "${LOGFILE}" || exit 1
run_python aggregate_markup_trends.py "${LOGFILE}" || exit 1
run_python favoritism_decomposition.py "${LOGFILE}" || exit 1
run_python orbis_acf_estimation.py "${LOGFILE}" || exit 1
run_python strong_exclusion_diagnostic.py "${LOGFILE}" || exit 1
run_python acf_strong_exclusion.py "${LOGFILE}" || exit 1
run_python borusyak_hull_randomization.py "${LOGFILE}" || exit 1
run_python ags_twostep_identification.py "${LOGFILE}" || exit 1
run_python misspecification_diagnostics.py "${LOGFILE}" || exit 1
run_python adl_instrument_comparison.py "${LOGFILE}" || exit 1
run_R lalonde_estimation.R "${LOGFILE}" || exit 1

# Phase 3: Specification curve (reads all prior outputs)
run_python specification_curve.py "${LOGFILE}" || exit 1

# Phase 4: Stata replication (DGM-style ACF pipeline — parallel to Python)
# Run from stata/ subdirectory so relative paths work
echo "Running Stata replication pipeline..." | tee -a "${LOGFILE}"
(cd source/stata && /Applications/StataNow/StataMP.app/Contents/MacOS/stata-mp -e do launcher.do) \
    >> "${LOGFILE}" 2>&1 || echo "Warning: Stata replication had errors" | tee -a "${LOGFILE}"

# Phase 5: Stata table formatting
run_stata paper_tables.do "${LOGFILE}" || exit 1
) || false

echo -e "\nmake.sh finished at $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${LOGFILE}"
