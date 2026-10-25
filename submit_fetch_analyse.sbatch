#!/usr/bin/env bash
# Routine: fetch new Microsynth result mails (01), identify the reads and email the PDF report (02).
# Output: output/sanger-microsynth/ (needs ~/.config/sanger-microsynth/mail.env; 02 SENDS a mail unless --no-send)
#SBATCH --job-name=sanger-microsynth
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/sanger-microsynth/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/sanger-microsynth
OUTPUT=$OUTPUT_ROOT/sanger-microsynth
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_sanger_ms

python "$REPO/scripts/01_fetch_mail.py"
python "$REPO/scripts/02_analyse_report.py" "$@"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT" "$REPO" submit_fetch_analyse.sbatch "$DATA/sanger/microsynth_mail" "$GTDB_REPRESENTATIVES" "$HSP65_AMPLICONS"
