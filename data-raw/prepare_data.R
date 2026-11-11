# =============================================================================
# MagnaGenes v2.0 - data preparation
#
# Converts the raw MagnaGenes TSV files (about 290 MB) into compact .rds files
# that load quickly in the browser (shinylive / webR).
#
# Run from the repository root:
#   Rscript data-raw/prepare_data.R
#
# Inputs (see data-raw/README.md for formats and how to add new data):
#   data-raw/source/annotation/fin_data*.tsv     (latest fin_data_YYYY-MM-DD.tsv; older archived)
#   data-raw/source/transcriptomics/*.tsv        (all files are combined)
#   data-raw/source/phosphoproteomics/fin_phosmag.tsv
#   data-raw/source/orthology/fin_ortho.tsv (from build_fin_ortho.R),
#     Species_Tree/SpeciesTree_rooted_node_labels.txt
#   data-raw/source/gene_ontology/go_universe.tsv
#   data-raw/metadata/rnaseq_experiments.csv     (labels, order, references)
#
# Needs (locally only, never in the app): data.table, ape, GO.db, AnnotationDbi
# Downloads (once, cached in data-raw/source/downloads/):
#   - Ensembl Fungi release 60, Magnaporthe oryzae MG8 protein and CDS FASTA
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(GO.db)
  library(AnnotationDbi)
})

args    <- commandArgs(trailingOnly = TRUE)
raw_dir  <- if (length(args)) args[1] else "data-raw/source"
meta_dir <- "data-raw/metadata"
out_dir  <- "app/data"
dl_dir   <- file.path(raw_dir, "downloads")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dl_dir,  showWarnings = FALSE, recursive = TRUE)

save_rds <- function(x, name) {
  f <- file.path(out_dir, paste0(name, ".rds"))
  saveRDS(x, f, compress = "xz")
  message(sprintf("  wrote %-22s %6.2f MB", basename(f), file.size(f) / 1e6))
}

# -----------------------------------------------------------------------------
# 1. Sequences (Ensembl Fungi MG8)
# -----------------------------------------------------------------------------
message("Sequences")
ens_base <- "https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-60/fasta/magnaporthe_oryzae"
fa_files <- c(pep = "Magnaporthe_oryzae.MG8.pep.all.fa.gz",
              cds = "Magnaporthe_oryzae.MG8.cds.all.fa.gz")
for (k in names(fa_files)) {
  dest <- file.path(dl_dir, fa_files[[k]])
  if (!file.exists(dest)) download.file(file.path(ens_base, k, fa_files[[k]]), dest, mode = "wb")
}

read_fasta <- function(path) {
  lines <- readLines(gzfile(path))
  hdr   <- grepl("^>", lines)
  id    <- cumsum(hdr)
  seqs  <- vapply(split(lines[!hdr], id[!hdr]), paste, "", collapse = "")
  list(header = sub("^>", "", lines[hdr]), seq = unname(seqs))
}

pep <- read_fasta(file.path(dl_dir, fa_files[["pep"]]))
cds <- read_fasta(file.path(dl_dir, fa_files[["cds"]]))

field <- function(h, key) {
  m <- regmatches(h, regexpr(paste0(key, ":[^ ]+"), h))
  out <- rep(NA_character_, length(h))
  out[grepl(paste0(key, ":"), h)] <- sub(paste0(key, ":"), "", m)
  out
}

tx <- data.table(
  transcript  = sub(" .*", "", pep$header),
  gene        = field(pep$header, "gene"),
  location    = field(pep$header, "chromosome"),
  description = sub(" \\[Source.*", "", sub(".*description:", "", pep$header)),
  uniprot_acc = sub("\\]$", "", field(pep$header, "Acc")),
  protein     = pep$seq
)
tx[!grepl("description:", pep$header), description := NA_character_]
tx[, c("chr", "start", "end", "strand") := tstrsplit(sub("^MG8:", "", location), ":")]
tx[, `:=`(start = as.integer(start), end = as.integer(end),
          strand = ifelse(strand == "1", "+", "-"), location = NULL)]
cds_dt <- data.table(transcript = sub(" .*", "", cds$header), cds = cds$seq)
tx <- merge(tx, cds_dt, by = "transcript", all.x = TRUE)
setorder(tx, gene, transcript)

# Saved separately so the app only loads them when a FASTA is requested
seq_source <- "Ensembl Fungi release 60, Magnaporthe oryzae MG8 (70-15)"
save_rds(list(transcript = tx$transcript, gene = tx$gene, seq = tx$protein,
              source = seq_source), "seq_protein")
save_rds(list(transcript = tx$transcript, gene = tx$gene, seq = tx$cds,
              source = seq_source), "seq_cds")

# One canonical transcript (T0 when present) per gene for the gene table
canon <- tx[order(gene, !grepl("T0$", transcript))][, .SD[1], by = gene]

# -----------------------------------------------------------------------------
# 2. Gene annotation table (one row per gene)
# -----------------------------------------------------------------------------
message("Annotations")
source("data-raw/resolve_latest_fin_data.R")
fin_data_path <- resolve_and_archive_fin_data(file.path(raw_dir, "annotation"))
message("  using ", basename(fin_data_path))
ann <- fread(fin_data_path, sep = "\t", quote = "", encoding = "UTF-8")
ann[, value := trimws(value)]
ann <- ann[value != "" & !is.na(value)]
ann <- ann[, .(value = paste(unique(value), collapse = " | ")), by = .(Gene, variable)]
wide <- dcast(ann, Gene ~ variable, value.var = "value")

rename <- c(
  "Gene" = "gene", "Gene.Name" = "gene_name", "BLAST2GO.Annotation" = "blast2go",
  "Size.(aa)" = "length_aa", "Effector" = "effector", "InterPro" = "interpro",
  "Pfam" = "pfam", "KEGG.Pathway" = "kegg", "Mutant.Phenotype" = "phenotype",
  "Process.Affected" = "process", "PMID.or.link" = "reference",
  "Other.links.to.publication" = "other_reference",
  "WGCNA.Co-Expression.Module.(Yan.et.al,.2023,.The.Plant.Cell)" = "wgcna_module",
  "AlphaFoldDB" = "alphafold_url", "Uniprot" = "uniprot_url", "Ensembl.Link" = "ensembl_url"
)
setnames(wide, names(rename)[names(rename) %in% names(wide)],
         rename[names(rename) %in% names(wide)])
wide[, length_aa := as.integer(length_aa)]
wide[, alphafold_acc := sub(".*/entry/", "", alphafold_url)]
wide[!grepl("/entry/", alphafold_url), alphafold_acc := NA_character_]

genes <- merge(wide, canon[, .(gene, transcript, chr, start, end, strand,
                               ensembl_description = description, uniprot_acc)],
               by = "gene", all.x = TRUE)
# UniProt accession: prefer Ensembl cross-reference, fall back to AlphaFold entry
genes[is.na(uniprot_acc), uniprot_acc := alphafold_acc]
genes[, c("alphafold_url", "uniprot_url", "ensembl_url") := NULL]
genes[, effector := !is.na(effector)]
setorder(genes, gene)
genes <- as.data.frame(genes)
save_rds(genes, "genes")

# -----------------------------------------------------------------------------
# 3. Expression (gene x condition matrix of mean TPM)
#    Every TSV in source/transcriptomics/ is read (columns: Gene, sample, TPM,
#    Experiment). Order, labels and paper links come from
#    metadata/rnaseq_experiments.csv; experiments missing there are appended.
# -----------------------------------------------------------------------------
message("Expression")
ex_files <- list.files(file.path(raw_dir, "transcriptomics"), pattern = "\\.tsv$", full.names = TRUE)
ex <- rbindlist(lapply(ex_files, fread, encoding = "UTF-8"), use.names = TRUE, fill = TRUE)
# Optional "Replicate" column: one row per biological replicate. When present,
# the app shows mean +/- SD error bars; otherwise values are already means.
if (!"Replicate" %in% names(ex)) ex[, Replicate := NA_character_]
ex_rep <- ex[!is.na(Replicate), .(TPM = mean(TPM, na.rm = TRUE)), by = .(Gene, sample, Experiment, Replicate)]
ex_sd  <- ex_rep[, .(SD = if (.N > 1) stats::sd(TPM) else NA_real_, N = .N, TPM = mean(TPM)),
                 by = .(Gene, sample, Experiment)]
ex <- rbind(ex_sd[, .(Gene, sample, Experiment, TPM)],
            ex[is.na(Replicate), .(TPM = mean(TPM, na.rm = TRUE)), by = .(Gene, sample, Experiment)])
ex <- ex[!duplicated(ex[, .(Gene, sample, Experiment)])]

reg <- fread(file.path(meta_dir, "rnaseq_experiments.csv"), encoding = "UTF-8")
extra <- setdiff(unique(ex$Experiment), reg$experiment)
if (length(extra)) {
  warning("Experiments not in rnaseq_experiments.csv (added with default labels): ",
          paste(extra, collapse = ", "))
  reg <- rbind(reg, data.table(experiment = extra, label = extra, samples = "", reference = "",
                               url = "", description = ""))
}
reg <- reg[experiment %in% ex$Experiment]

cond <- rbindlist(lapply(seq_len(nrow(reg)), function(i) {
  smp <- unique(ex[Experiment == reg$experiment[i], sample])
  ord <- strsplit(reg$samples[i], ";")[[1]]
  smp <- c(intersect(ord, smp), sort(setdiff(smp, ord)))
  data.table(experiment = reg$experiment[i], sample = smp)
}))
cond[, id := paste(experiment, sample, sep = " | ")]
ex[, id := paste(Experiment, sample, sep = " | ")]
m <- dcast(ex, Gene ~ id, value.var = "TPM")
mat <- as.matrix(m[, cond$id, with = FALSE])
rownames(mat) <- m$Gene
mat <- signif(mat, 4)
reg[, has_replicates := experiment %in% ex_sd[N > 1, unique(Experiment)]]
expression <- list(
  tpm = mat,
  conditions = as.data.frame(cond),
  experiments = as.data.frame(reg),
  sample_levels = unique(cond$sample)
)
if (nrow(ex_sd)) {
  ex_sd[, id := paste(Experiment, sample, sep = " | ")]
  s <- dcast(ex_sd, Gene ~ id, value.var = "SD")
  sdm <- matrix(NA_real_, nrow(mat), ncol(mat), dimnames = dimnames(mat))
  hit <- intersect(names(s), colnames(mat))
  sdm[match(s$Gene, rownames(mat)), hit] <- as.matrix(s[, hit, with = FALSE])
  expression$sd <- signif(sdm, 4)
}
save_rds(expression, "expression")

# -----------------------------------------------------------------------------
# 4. Phosphoproteomics (normalised abundance, peptide x sample/time/replicate)
#    "Averaged" rows equal the replicate mean, so they are recomputed in the app.
# -----------------------------------------------------------------------------
message("Phosphoproteomics")
ph <- fread(file.path(raw_dir, "phosphoproteomics", "fin_phosmag.tsv"), encoding = "UTF-8")[type == "Normalized"]
ph[, pid := .GRP, by = .(Master_Protein_Accessions, Positions_in_Master_Proteins, label, Modifications)]
pep_meta <- unique(ph[, .(pid, protein = Master_Protein_Accessions,
                          position = Positions_in_Master_Proteins,
                          peptide = label, modifications = Modifications)])

# Absolute phosphosite(s) on the protein, e.g. peptide [126-140] + S3 -> S128
site_of <- function(position, mods) {
  start <- suppressWarnings(as.integer(sub("^[^\\[]*\\[([0-9]+)-.*", "\\1", position)))
  ph    <- regmatches(mods, regexpr("Phospho \\[[^]]*\\]", mods))
  if (!length(ph) || is.na(start)) return(NA_character_)
  parts <- strsplit(sub("Phospho \\[(.*)\\]", "\\1", ph), "; ")[[1]]
  out <- vapply(parts, function(s) {
    if (grepl("^[STY][0-9]+", s)) {
      res <- substr(s, 1, 1)
      pos <- as.integer(sub("^[STY]([0-9]+).*", "\\1", s))
      paste0(res, start + pos - 1L)
    } else paste0(s, " (ambiguous)")
  }, "")
  paste(out, collapse = "; ")
}
pep_meta[, sites := mapply(site_of, position, modifications)]
pep_meta[, n_phospho := as.integer(sub("^.*?([0-9]+)xPhospho.*", "\\1", modifications))]
pep_meta[!grepl("xPhospho", modifications), n_phospho := NA_integer_]

ph[, time := factor(time, c("hpi_0", "hpi_1", "hpi_1.5", "hpi_2", "hpi_4", "hpi_6"))]
ph[, sample := factor(sample, c("Guy11", "Δpmk1"))]
ph[, col := paste(sample, time, replicate, sep = "|")]
cols <- unique(ph[, .(sample, time, replicate, col)])
setorder(cols, sample, time, replicate)
pw <- dcast(ph, pid ~ col, value.var = "Abundances")
pm <- as.matrix(pw[, cols$col, with = FALSE])
rownames(pm) <- pw$pid
pm <- signif(pm, 5)
setorder(pep_meta, pid)
phospho <- list(
  meta = as.data.frame(pep_meta),
  abundance = pm[as.character(pep_meta$pid), , drop = FALSE],
  columns = data.frame(sample = as.character(cols$sample),
                       time = as.character(cols$time),
                       hpi = as.numeric(sub("hpi_", "", cols$time)),
                       replicate = cols$replicate, stringsAsFactors = FALSE)
)
save_rds(phospho, "phospho")

# -----------------------------------------------------------------------------
# 5. Orthology (orthogroups that contain M. oryzae genes only)
# -----------------------------------------------------------------------------
message("Orthology")
og <- fread(file.path(raw_dir, "orthology", "fin_ortho.tsv"), encoding = "UTF-8")
members <- unique(og[!is.na(ProtID) & ProtID != "", .(gene = ProtID, og = OG_ID)])
og <- og[OG_ID %in% members$og]
species <- unique(og[, .(code = Species, name = Name, appressorium = Appresor, lifestyle = Lifestyle)])
counts  <- unique(og[, .(OG_ID, Species, Counts)])
cm <- dcast(counts, OG_ID ~ Species, value.var = "Counts", fill = 0L)
cmat <- as.matrix(cm[, -1])
rownames(cmat) <- cm$OG_ID
storage.mode(cmat) <- "integer"

# Species tree: precomputed rectangular-layout coordinates (no ape needed in app)
tree <- read.tree(file.path(raw_dir, "orthology", "Species_Tree", "SpeciesTree_rooted_node_labels.txt"))
tree$node.label <- NULL
tree <- drop.tip(tree, setdiff(tree$tip.label, species$code))   # only species in fin_ortho.tsv
tree <- ladderize(tree, right = FALSE)
tree <- read.tree(text = write.tree(tree))   # renumber nodes in ladderized order
n_tip <- length(tree$tip.label)
edge  <- tree$edge
x_phylo <- node.depth.edgelength(tree)
x_clado <- max(node.depth(tree, method = 2)) - node.depth(tree, method = 2)
y <- numeric(n_tip + tree$Nnode)
tip_order <- edge[edge[, 2] <= n_tip, 2]
y[tip_order] <- rev(seq_along(tip_order))
# Internal node y = midpoint of its children; iterate until stable (tree is small)
repeat {
  changed <- FALSE
  for (nd in unique(edge[, 1])) {
    v <- mean(range(y[edge[edge[, 1] == nd, 2]]))
    if (!isTRUE(all.equal(v, y[nd]))) { y[nd] <- v; changed <- TRUE }
  }
  if (!changed) break
}
tree_layout <- list(
  edge = edge, n_tip = n_tip, tip_code = tree$tip.label,
  x_phylo = x_phylo, x_clado = x_clado, y = y
)
species <- species[match(tree$tip.label, code)]
species[, tree_y := y[seq_len(n_tip)]]

orthology <- list(
  members = as.data.frame(members),
  counts = cmat[, species$code],
  species = as.data.frame(species),
  tree = tree_layout
)
save_rds(orthology, "orthology")

# -----------------------------------------------------------------------------
# 6. Gene Ontology (annotations propagated to all ancestor terms)
# -----------------------------------------------------------------------------
message("Gene Ontology")
gu <- fread(file.path(raw_dir, "gene_ontology", "go_universe.tsv"), sep = "\t", header = TRUE)
direct <- gu[, .(go = trimws(unlist(strsplit(GO_ID, ",")))), by = .(gene = ProtID)]
direct <- unique(direct[grepl("^GO:", go)])

go_ont <- suppressMessages(AnnotationDbi::select(GO.db, keys = unique(direct$go),
                                                 columns = c("TERM", "ONTOLOGY"),
                                                 keytype = "GOID"))
direct <- merge(direct, data.table(go = go_ont$GOID, ont = go_ont$ONTOLOGY), by = "go")
direct <- direct[ont %in% c("BP", "MF", "CC")]

anc_maps <- list(BP = as.list(GOBPANCESTOR), MF = as.list(GOMFANCESTOR), CC = as.list(GOCCANCESTOR))
prop <- rbindlist(lapply(c("BP", "MF", "CC"), function(o) {
  d <- direct[ont == o]
  anc <- anc_maps[[o]][unique(d$go)]
  expd <- rbindlist(lapply(names(anc), function(g)
    data.table(go = g, term = unique(c(g, setdiff(anc[[g]], "all"))))))
  res <- unique(merge(d, expd, by = "go", allow.cartesian = TRUE)[, .(gene, go = term)])
  res[, ont := o]
  res
}))

terms <- unique(prop[, .(go, ont)])
tn <- suppressMessages(AnnotationDbi::select(GO.db, keys = terms$go, columns = "TERM", keytype = "GOID"))
terms[, term := tn$TERM[match(go, tn$GOID)]]
setorder(terms, ont, go)
gene_ids <- sort(unique(prop$gene))
go_data <- list(
  genes = gene_ids,
  terms = as.data.frame(terms),
  # integer pairs: term index x gene index (compact)
  term_idx = match(prop$go, terms$go),
  gene_idx = match(prop$gene, gene_ids),
  direct_n = direct[, .N, by = ont],
  source = paste("GO.db", as.character(packageVersion("GO.db")))
)
save_rds(go_data, "go")

message("Done.")
