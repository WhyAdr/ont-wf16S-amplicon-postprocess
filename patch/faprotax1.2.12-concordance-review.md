# FAPROTAX 1.2.12 engine concordance review

Date: 2026-09-06

This review compares the tracked Ambar Ayunda fixture after the same input
filter used by the opt-in module: 1,836 positive-count Bacteria/Archaea rows,
with the canonical unclassified row, classified non-prokaryotes, and zero-count
rows excluded. The comparison uses the seven fields passed to microeco as
`Kingdom;Phylum;Class;Order;Family;Genus;Species`; the pipeline's NCBI
`kingdom` field is not passed to either engine.

The official reference was FAPROTAX 1.2.12's complete package and its
`collapse_table.py` script, downloaded from the [Louca Lab download page](https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/lib/php/index.php?section=Download).
The command used raw counts and no normalization:

```text
collapse_table.py -i official_input.tsv -g FAPROTAX.txt \
  -o collapsed.tsv --out_groups2records_table groups2records.tsv \
  -c '#' -d taxonomy --omit_columns 0 \
  --column_names_are_in last_comment_line -n none --force
```

The official workflow loaded 92 group definitions, represented 77 groups,
and established 2,524 taxon-function assignments across 827 input taxa. Its
function-count total was 120,352 reads and the union of reads belonging to at
least one function was 50,233. The microeco 2.3.0 run used by this integration
returned 93 binary function columns, 78 functions with positive counts, 3,133
positive taxon-function rows across 924 feature IDs, a function-count total of
158,498 reads, and 52,403 reads belonging to at least one function. Both
engines' mapped-read values are distinct from the sum of overlapping function
counts.

Thirty-four function counts differed. The largest difference was the
microeco-only derived `anaerobic_chemoheterotrophy` column (33,543 reads;
official output 0). Other differences were:

| Function | Official reads | microeco reads |
|---|---:|---:|
| aerobic_chemoheterotrophy | 14,966 | 17,119 |
| aliphatic_non_methane_hydrocarbon_degradation | 3 | 6 |
| animal_parasites_or_symbionts | 1,911 | 1,958 |
| aromatic_compound_degradation | 2,237 | 2,238 |
| aromatic_hydrocarbon_degradation | 372 | 373 |
| chemoheterotrophy | 48,499 | 50,662 |
| chitinolysis | 139 | 154 |
| dark_hydrogen_oxidation | 622 | 633 |
| fermentation | 32,589 | 32,615 |
| fumarate_respiration | 84 | 96 |
| human_associated | 1,846 | 1,892 |
| human_gut | 1,764 | 1,757 |
| human_pathogens_all | 91 | 144 |
| human_pathogens_meningitis | 1 | 11 |
| human_pathogens_pneumonia | 10 | 20 |
| human_pathogens_septicemia | 4 | 2 |
| hydrocarbon_degradation | 378 | 388 |
| iron_respiration | 179 | 181 |
| ligninolysis | 41 | 42 |
| mammal_gut | 1,764 | 1,757 |
| manganese_respiration | 64 | 75 |
| methanotrophy | 4 | 10 |
| methylotrophy | 218 | 222 |
| nitrate_reduction | 1,557 | 1,621 |
| nitrate_respiration | 1,106 | 1,116 |
| nitrogen_fixation | 248 | 259 |
| nitrogen_respiration | 1,274 | 1,289 |
| nitrite_ammonification | 231 | 236 |
| nitrite_respiration | 890 | 895 |
| plant_pathogen | 122 | 5 |
| reductive_acetogenesis | 60 | 66 |
| sulfur_respiration | 48 | 59 |
| ureolysis | 590 | 614 |

The difference is an engine/database-representation difference, not a failure
of the pipeline's read accounting. `collapse_table.py` applies the official
text database's group matcher and set operations. `microeco::trans_func$cal_func`
uses its embedded regex, additive, and subtractive representation and adds the
derived anaerobic-chemoheterotrophy column. This is consistent with the
[microeco issue documenting differences between the two engines](https://github.com/ChiLiubio/microeco/issues/534).
Accordingly, this integration uses the microeco result by design, records the
engine and database versions, and makes no byte-for-byte concordance claim.
