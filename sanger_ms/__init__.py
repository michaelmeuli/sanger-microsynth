from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SHARES = REPO_ROOT.parent.parent
# Attachments of the fetched Microsynth mails, one subfolder per mail.
MAIL_DIR = SHARES / "data" / "sanger" / "microsynth_mail"
OUTPUT = SHARES / "output" / "sanger-microsynth"
STATE_FILE = MAIL_DIR / "processed_messages.json"
# Credentials live outside the repo and the shares.
MAIL_ENV = Path.home() / ".config" / "sanger-microsynth" / "mail.env"
# One representative genome per kansasii-complex species, for the read-to-reference alignment.
REFERENCES = SHARES / "data" / "gtdb_genomes" / "Mycobacteriaceae" / "kansasii_complex_gtdb_representatives"
# Order of the species in the alignment/matrix (same as mlsa-kansasii's SPECIES).
SPECIES = ["kansasii", "persicum", "pseudokansasii", "innocens", "attenuatum", "ostraviense", "gastri"]
# The in-silico hsp65 amplicon of every GTDB genome of the complex (mlsa-kansasii scripts/make_hsp65_references.py):
# reference set for hsp65 reads, so within-species diversity is covered and the atypical-hsp65 kansasii is explicit.
REFERENCES_HSP65 = SHARES / "data" / "gtdb_genomes" / "Mycobacteriaceae" / "kansasii_complex_hsp65_amplicons"
