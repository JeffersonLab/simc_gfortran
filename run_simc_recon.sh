#!/bin/bash
# run_simc_recon.sh
#
# Automate:
# 1) Run SIMC (convert input, build, run SIMC)
# 2) Run recon_hcana on the generated SIMC ROOT/HIST files
#
# Usage:
#   ./run_simc_recon.sh <stem> <reaction> [hadron_type] [Earm_HMS]
#
# where:
#   <stem>        = base name of your SIMC input/output (without .inp)
#   <reaction>    = "heep", "sidis", "rho", "delta", "exclusive", ...
#   [hadron_type] = "mpi" (default) or "mk". For heep, proton mass is hard-coded
#   [Earm_HMS]    = 1 (default: electron in HMS) or 0 (electron in SHMS)
#
# Example (SIDIS):
#   ./run_simc_recon.sh coin_7p87deg_3p632gev_hyd_rsidis sidis mpi 1
#
# Example (heep with electron in SHMS):
#   ./run_simc_recon.sh heep_rsidis_4pass_elastic1 heep mpi 0

set -e  # exit on first error

############################
# CONFIGURE THESE PATHS
############################

# Top-level SIMC directory
SIMC_DIR="$(pwd)"

# Path to SIMC input files
SIMC_INFILES_DIR="${SIMC_DIR}/infiles"

# Path to SIMC ROOT tree util
SIMC_ROOT_TREE_DIR="${SIMC_DIR}/util/root_tree"

# Where SIMC writes .hist files
SIMC_OUTFILES_DIR="${SIMC_DIR}/outfiles"

# Where SIMC writes .root tree files
SIMC_WORK_DIR="${SIMC_DIR}/worksim"

# Directory where recon_hcana.C / recon_hcana.h live
RECON_DIR="${SIMC_DIR}/util/recon_hcana"

############################
# PARSE ARGUMENTS
############################

if [ $# -lt 2 ]; then
  echo "Usage: $0 <stem> <reaction> [hadron_type] [Earm_HMS]"
  echo "  <stem>        = SIMC basename (no .inp)"
  echo "  <reaction>    = heep | sidis | rho | delta | exclusive | ..."
  echo "  [hadron_type] = mpi (default) or mk"
  echo "  [Earm_HMS]    = 1 (default: electron in HMS) or 0 (electron in SHMS)"
  exit 1
fi

STEM="$1"
REACTION="$2"
HADRON_TYPE="${3:-mpi}"
EARM_HMS_INT="${4:-1}"

# Map 1/0 → kTRUE/kFALSE for ROOT
if [ "$EARM_HMS_INT" -eq 0 ]; then
  EARM_FLAG="kFALSE"
else
  EARM_FLAG="kTRUE"
fi

echo "========================================="
echo " Running SIMC + recon_hcana"
echo "-----------------------------------------"
echo "  STEM        = ${STEM}"
echo "  REACTION    = ${REACTION}"
echo "  HADRON_TYPE = ${HADRON_TYPE}"
echo "  Earm_HMS    = ${EARM_FLAG}  (1=HMS, 0=SHMS)"
echo "========================================="

############################
# 0. CHECK INPUT .INP FILE
############################

INP_FILE="${SIMC_INFILES_DIR}/${STEM}.inp"

if [ ! -f "${INP_FILE}" ]; then
  echo "ERROR: Input file not found: ${INP_FILE}"
  echo "Please create/write ${STEM}.inp in ${SIMC_INFILES_DIR} first."
  exit 1
fi

############################
# 1. convert_inputfile.sh
############################

echo
echo ">>> [1] Converting SIMC input file"
cd "${SIMC_INFILES_DIR}"

if [ ! -x "./convert_inputfile.sh" ]; then
  echo "ERROR: convert_inputfiles.sh not found or not executable in ${SIMC_INFILES_DIR}"
  exit 1
fi

./convert_inputfile.sh "${STEM}.inp"

############################
# 2. Build ROOT tree code
############################

echo
echo ">>> [2] Building SIMC ROOT tree (make clean; make)"
cd "${SIMC_ROOT_TREE_DIR}"

make clean
make

############################
# 3. Run SIMC to make .hist/.root
############################

echo
echo ">>> [3] Running SIMC (./run_simc_tree ${STEM})"
cd "${SIMC_DIR}"

if [ ! -x "./run_simc_tree" ]; then
  echo "ERROR: run_simc_tree not found or not executable in ${SIMC_DIR}"
  exit 1
fi

./run_simc_tree "${STEM}"

# Expect .hist and .root now:
HIST_FILE="${SIMC_OUTFILES_DIR}/${STEM}.hist"
ROOT_FILE="${SIMC_WORK_DIR}/${STEM}.root"

echo
echo "Checking SIMC outputs:"
echo "  HIST: ${HIST_FILE}"
echo "  ROOT: ${ROOT_FILE}"

if [ ! -f "${HIST_FILE}" ]; then
  echo "ERROR: SIMC .hist file not found: ${HIST_FILE}"
  exit 1
fi

if [ ! -f "${ROOT_FILE}" ]; then
  echo "ERROR: SIMC .root file not found: ${ROOT_FILE}"
  exit 1
fi

############################
# 4. Run recon_hcana
############################

echo
echo ">>> [4] Running recon_hcana"

cd "${RECON_DIR}"

# Make sure recon_hcana.C is visible
if [ ! -f "recon_hcana.C" ]; then
  echo "ERROR: recon_hcana.C not found in ${RECON_DIR}"
  exit 1
fi

root -l -b -q "recon_hcana.C+(\"${STEM}\",\"${REACTION}\",\"${HADRON_TYPE}\",${EARM_FLAG})"

RETVAL=$?

if [ ${RETVAL} -ne 0 ]; then
  echo "ERROR: recon_hcana.C failed with exit code ${RETVAL}"
  exit ${RETVAL}
fi

echo
echo "================================================================="
echo " All done!"
echo " - SIMC input:    ${INP_FILE}"
echo " - SIMC ROOT:     ${ROOT_FILE}"
echo " - SIMC HIST:     ${HIST_FILE}"
echo " - recon_hcana output: ${SIMC_WORK_DIR}/recon_hcana_${STEM}.root"
echo "================================================================="
