#!/usr/bin/env bash
#SBATCH --job-name=sanger-microsynth
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/sanger-microsynth/logs/%x-%j.log

set -euo pipefail

REPO_ROOT="/shares/sander.imm.uzh/MM/kansasii/repos/sanger-microsynth"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_sanger_ms

python "$REPO_ROOT/main1_fetch_mail/run_fetch_mail.py"
python "$REPO_ROOT/main2_analyse_report/run_analyse_report.py"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-docs/bin/write_manifest.sh \
  "/shares/sander.imm.uzh/MM/kansasii/output/sanger-microsynth" \
  /shares/sander.imm.uzh/MM/kansasii/repos/sanger-microsynth submit_fetch_analyse.sbatch
