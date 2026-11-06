# MagnaGenes v2.0

An integrated multi-omics explorer for *Magnaporthe oryzae* (*Pyricularia oryzae*), the rice blast fungus. It runs entirely in the web browser with [Shinylive](https://posit-dev.github.io/r-shinylive/) (R compiled to WebAssembly), so it needs no server, no API keys and no AI services, and no data leave the user's computer.

MagnaGenes v2.0 is an extension of MagnaGenes v1.0 by Andrew Foster: <https://zenodo.org/records/4647766>.

## Features

| Tab | What it does |
|---|---|
| Gene Search | Search by MGG ID, gene name (e.g. `PMK1`) or keyword in any annotation field. Gene profiles show genome location, InterPro/Pfam domains, mutant phenotypes, literature, an expression chart and a protein architecture diagram, with links to UniProt, AlphaFold DB, Ensembl Fungi, InterPro, Pfam and PubMed. Exports annotation TSV and protein/CDS FASTA. |
| Transcriptomics | TPM line plots per experiment (error bars when replicate values are available, optional box plot of all genes), clustered row z-score heatmaps with optional TPM labels, gene-gene Pearson correlation on a range fitted to the selected genes, and a data table across 31 conditions in 5 experiments. |
| Phosphoproteomics | One plot per phosphopeptide with Guy11 above Δpmk1 (bars are biological replicates), abundance heatmaps, log2 fold-change heatmaps with Welch t-tests (BH-corrected), and a protein site map with labelled sites and Pfam domains. |
| Phylogenomics | Ortholog counts across 42 fungi aligned to the species tree (phylogram or cladogram), annotated by lifestyle and appressorium type; conservation by lifestyle; M. oryzae orthogroup members with FASTA export. |
| GO Enrichment | Fisher's exact test (hypergeometric) on propagated GO annotations, BH or Bonferroni correction, optional custom background, dot or lollipop plot, bar chart, treemap, significance-vs-enrichment plot, gene-term matrix and results table with QuickGO links, in a colour of your choice. Replaces topGO, which does not run in WebAssembly. |

The site has a light/dark theme switch, a Reset button and a Tutorial button with a video tour and clickable chapters. Every chart can be saved as a vector SVG from its toolbar (camera icon), and every tab has PDF, TSV and (where relevant) FASTA downloads.

## Repository layout

```
app/                        everything the website needs
  app.R                     the whole application (UI + server in one file)
  data/*.rds                compact data bundle (about 10 MB), built from data-raw/
  www/                      images used by the app
data-raw/                   how the data are built (see data-raw/README.md)
  prepare_data.R            builds app/data/ from the raw files
  metadata/                 rnaseq_experiments.csv: experiment labels, order, paper links
  source/                   raw TSV files by data type (not in git; too large)
site-assets/                loading-screen image, tutorial video and its chapter list
build_site.R                Shinylive export plus the loading screen
tests/smoke_test.R          renders every output and download
tools/record_tutorial.py    re-records the tutorial video (optional; needs Playwright)
.github/workflows/          publishes the site to GitHub Pages on every push
```

## Run locally

```r
install.packages(c("shiny", "bslib", "plotly", "ggplot2", "reactable"))
shiny::runApp("app")
```

## Build and preview the Shinylive site

```r
install.packages("shinylive")
source("build_site.R")             # exports to site/ with the custom loading screen
httpuv::runStaticServer("site")   # then open the printed URL
```

The first load takes about a minute while webR and the R packages download. Later visits load from the browser cache.

## Tutorial video

The Tutorial button plays `site-assets/tutorial.mp4` with chapters from
`site-assets/tutorial_chapters.json`; `build_site.R` copies both next to the app, so the
video is only downloaded when someone opens it. After changing the app, re-record it with:

```sh
pip install playwright imageio-ffmpeg && playwright install ffmpeg
Rscript build_site.R
python -m http.server 8765 --bind 127.0.0.1 --directory site &
python tools/record_tutorial.py
```

## Deploy to GitHub Pages

1. Push this repository to GitHub.
2. In the repository settings, open **Pages** and set **Source** to **GitHub Actions**.
3. Every push to `main` runs `.github/workflows/deploy-shinylive.yml`, which builds the site with `build_site.R` and publishes it at `https://<user>.github.io/<repository>/`.

## Rebuilding the data and adding experiments

The raw TSV files (about 290 MB) live in `data-raw/source/` and are not stored in git.
`data-raw/README.md` explains the formats and how to add a new RNA-seq experiment. To rebuild:

```sh
Rscript data-raw/prepare_data.R
```

The script needs `data.table`, `ape`, `GO.db` and `AnnotationDbi` (local use only), and downloads the *M. oryzae* MG8 protein and CDS FASTA from Ensembl Fungi release 60 on first run.

What it does:

- Annotation: long format converted to one row per gene; UniProt accessions and chromosome coordinates added from Ensembl.
- Expression: a gene x condition matrix of mean TPM (duplicate rows averaged), plus standard deviations when a `Replicate` column is present.
- Phosphoproteomics: normalised replicate abundances only (the "Averaged" rows equal the replicate mean and are recomputed in the app); absolute phosphosite positions derived from peptide coordinates.
- Orthology: only orthogroups that contain M. oryzae genes are kept (10,015 of 78,665); species tree coordinates precomputed.
- GO: annotations propagated to all ancestor terms using GO.db.

## Data sources

- Curated MagnaGenes annotation (BLAST2GO, InterPro, Pfam, KEGG, mutant phenotypes, effectors, literature); WGCNA modules from Yan et al. (2023), *The Plant Cell*.
- Ensembl Fungi release 60, *M. oryzae* MG8 (70-15): sequences, coordinates and UniProt cross-references.
- Gene Ontology term names and structure from the Bioconductor GO.db package.

## Contact

neha.sahu.tsl@gmail.com
