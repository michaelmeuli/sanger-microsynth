# Species test panel (batch `20261016_species_panel`)

Ten hsp65 Sanger traces of isolates with a known species (`data/imm/screening_map_results.csv`, column `species`),
copied from `data/sanger/seq_kansasii/` into `data/sanger/microsynth_mail/20261016_species_panel/` as `<TNR> hsp65_Tb-11.ab1`.
Analysed with `main2_analyse_report/run_analyse_report.py --batch <dir> --no-send --no-blast`.

| Species | Isolate | Trace (TNR) | Closest reference |
|---|---|---|---|
| ostraviense | Mkan329-181 | 2022500313 | ostraviense 99.74% |
| innocens | Mkan329-008 | 2025500758 | innocens 100% |
| attenuatum | Mkan329-073 | 2018500044, `hsp65_11_wdh` | attenuatum 100% |
| attenuatum | Mkan329-073 | 2018500044, `hsp65_11_12w` | attenuatum 100% (only 141 bp) |
| pseudokansasii | Mkan329-119 | 2014500061 | pseudokansasii 99.74% |
| pseudokansasii | Mkan329-084 | 2017500086 | pseudokansasii 99.73% |
| persicum | Mkan329-082 | 2017500195 | persicum 99.74% |
| persicum | Mkan329-002 | 2026500191 | persicum 99.74% |
| kansasii | Mkan329-017 | 2025500342 | kansasii 99.49% |
| kansasii | Mkan329-011 | 2025500620 | kansasii 99.74% |

Notes:
- Only two attenuatum isolates exist (Mkan329-073, Mkan329-030), and -030 has no hsp65 trace, so both attenuatum reads are from -073.
- Other candidate isolates were dropped when their trace gave a short or divergent read (e.g. Mkan329-003, -004, -006, -007).
