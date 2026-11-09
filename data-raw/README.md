# MagnaGenes data

The app never reads these raw files directly. `prepare_data.R` turns them into the compact
`.rds` files in `app/data/` that are bundled with the Shinylive site.

```sh
# from the repository root
Rscript data-raw/prepare_data.R
```

## Layout

```
data-raw/
  prepare_data.R                  builds app/data/*.rds
  build_fin_ortho.R               builds source/orthology/fin_ortho.tsv from OrthoFinder output
  metadata/
    rnaseq_experiments.csv        RNA-seq experiments: label, sample order, paper link
  source/                         raw data (not in git: too large for GitHub)
    annotation/fin_data.tsv       gene annotation, long format (Gene, variable, value)
    transcriptomics/*.tsv         RNA-seq TPM, long format (Gene, sample, TPM, Experiment)
    phosphoproteomics/fin_phosmag.tsv
    orthology/                    OrthoFinder 3 results (50 species):
      Phylogenetic_Hierarchical_Orthogroups/N0.tsv   hierarchical orthogroups (used)
      Orthogroups/Orthogroups_UnassignedGenes.tsv    genes without orthologs
      Species_Tree/SpeciesTree_rooted_node_labels.txt  species tree, rooted on Rhizophagus
      SpeciesTree_rooted_Rhizophagus.txt             the rooted tree given to OrthoFinder (-s)
    orthology/Lifestyle_extended.tsv  species name, lifestyle and appressorium type; only
                                  species listed here are used (Candida albicans
                                  "Candi_albic" was removed; the JGI genome Canalb1 is kept)
    orthology/fin_ortho.tsv       orthogroup counts per species (built by build_fin_ortho.R)
    gene_ontology/go_universe.tsv gene to GO term mapping
    downloads/                    cached Ensembl FASTA (downloaded automatically)
```

## Adding a new RNA-seq experiment

1. Save the data as a tab-separated file in `source/transcriptomics/`, for example
   `source/transcriptomics/smith_2027_septins.tsv`, with these columns:

   | Gene | sample | TPM | Experiment |
   |---|---|---|---|
   | MGG_00005 | 0hpi | 12.4 | Guy11 - septin mutant |

   `Experiment` is the name that links the rows to the registry. Every `.tsv` file in the
   folder is read and combined, so existing files do not need to change.

   To show error bars, add an optional `Replicate` column and give one row per
   biological replicate (for example `rep1`, `rep2`, `rep3`). The app then plots the mean
   with the standard deviation. Without it, each value is taken as the averaged TPM of the
   replicates and the plots say so.

   | Gene | sample | TPM | Experiment | Replicate |
   |---|---|---|---|---|
   | MGG_00005 | 0hpi | 11.9 | Guy11 - septin mutant | rep1 |

2. Add one row to `metadata/rnaseq_experiments.csv`:

   | column | meaning |
   |---|---|
   | experiment | exactly the value used in the `Experiment` column |
   | label | short name shown in the app |
   | samples | sample order, separated by `;` (samples not listed are added at the end) |
   | reference | citation shown next to the plot, e.g. `Smith et al. (2027) Nature` |
   | url | link to the paper (a DOI link is best) |
   | description | one sentence describing the experiment |

   The row order in this file is the order of the plots in the app.

3. Run `Rscript data-raw/prepare_data.R`, then `Rscript tests/smoke_test.R`.

The Transcriptomics tab picks the new experiment up automatically: it gets its own
check box, paper link and framed plot.

## Other data types

- **Annotation:** add rows to `fin_data.tsv` (`Gene`, `variable`, `value`). New variables
  appear in downloads; to show one in the gene profile, add it in `app/app.R`.
- **Orthology:** replace the OrthoFinder result folders and `Lifestyle_extended.tsv` (every species
  in the OrthoFinder run needs a row), then run `Rscript data-raw/build_fin_ortho.R` and
  `Rscript data-raw/prepare_data.R`. A new lifestyle also needs a colour in `COL_LIFESTYLE`
  in `app/app.R`.
- **Phosphoproteomics, GO:** replace the file in its folder with the same columns
  and rerun the script.

## Keeping the raw data

`source/` is listed in `.gitignore` because the orthology files (over 100 MB in total) is larger than
GitHub's 100 MB file limit. Keep a copy of `source/` somewhere safe (lab storage, Zenodo,
or Git LFS) so the data can be rebuilt.
