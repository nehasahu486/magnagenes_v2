# =============================================================================
# MagnaGenes v2.0 - build orthology/fin_ortho.tsv from OrthoFinder results
#
# Run from the repository root, before prepare_data.R:
#   Rscript data-raw/build_fin_ortho.R
#
# Inputs (data-raw/source/orthology/, OrthoFinder 3 results folder layout):
#   Phylogenetic_Hierarchical_Orthogroups/N0.tsv
#                               hierarchical orthogroups (HOGs) at the root of the
#                               species tree; one row per HOG, genes per species
#   Orthogroups/Orthogroups_UnassignedGenes.tsv
#                               genes without orthologs in any other species
#   Species_Tree/SpeciesTree_rooted_node_labels.txt
#                               species tree rooted on Rhizophagus irregularis
#                               (read by prepare_data.R)
#   Lifestyle_extended.tsv      Species, Full_name, Lifestyle, Appresor
#                               (defines which species are used)
#
# Candida albicans "Candi_albic" (an incomplete assembly) was removed from all
# these files on 2026-11-09; only the JGI genome "Canalb1" is used.
#
# HOGs are used instead of Orthogroups/Orthogroups.tsv because they are split
# with gene trees: large families (kinases, transcription factors, P450s) are
# divided into orthologous groups instead of being lumped into one cluster.
#
# M. oryzae genes that are in no HOG (unassigned genes, and genes that OrthoFinder
# leaves out of N0) are added as single-gene, M. oryzae-specific groups named
# "unassigned_MGG_xxxxx". The gene list comes from the Ensembl MG8 proteome
# downloaded by prepare_data.R (data-raw/source/downloads/).
#
# OrthoFinder was run on 51 proteomes (-M msa -A mafft_memsave -T fasttree_fastest),
# then the last step was rerun with the species tree rooted on Rhizophagus
# irregularis (SpeciesTree_rooted_Rhizophagus.txt): -ft <results> -s <tree>.
#
# Output: fin_ortho.tsv, one row per HOG x species x M. oryzae gene, for HOGs
# that contain M. oryzae genes. Columns:
#   OG_ID, Species, Counts, Appresor, Lifestyle, Name, ProtID
# M. oryzae transcripts (MGG_xxxxxT0, T1, ...) are collapsed to genes, so the
# M. oryzae count is the number of genes, not transcripts.
#
# Needs: data.table
# =============================================================================

suppressPackageStartupMessages(library(data.table))

dir <- "data-raw/source/orthology"
MO  <- "MGG"   # species code of Magnaporthe oryzae in the OrthoFinder run
strip_tx <- function(x) sub("T[0-9]+$", "", x)
n_genes  <- function(x) ifelse(x == "", 0L, lengths(regmatches(x, gregexpr(",", x, fixed = TRUE))) + 1L)

life <- fread(file.path(dir, "Lifestyle_extended.tsv"), encoding = "UTF-8")
life[, (names(life)) := lapply(.SD, trimws, whitespace = "[\\h\\v]")]   # also non-breaking spaces
life[, Full_name := sub("_JGI$", " (JGI)", Full_name)]

hog <- fread(file.path(dir, "Phylogenetic_Hierarchical_Orthogroups", "N0.tsv"), colClasses = "character")
hog[, c("OG", "Gene Tree Parent Clade") := NULL]
setnames(hog, "HOG", "OG_ID")
hog[, OG_ID := sub("^N0\\.", "", OG_ID)]
species_cols <- setdiff(names(hog), "OG_ID")
# Only species listed in Lifestyle_extended.tsv are used. Any other column (for
# example a proteome later judged unreliable) is dropped here.
extra <- setdiff(species_cols, life$Species)
if (length(extra)) {
  message("Ignoring species not in Lifestyle_extended.tsv: ", paste(extra, collapse = ", "))
  hog[, (extra) := NULL]
  species_cols <- setdiff(species_cols, extra)
}
hog <- hog[get(MO) != ""]

members <- hog[, .(ProtID = unique(strip_tx(trimws(strsplit(get(MO), ",", fixed = TRUE)[[1]])))), by = OG_ID]

counts <- melt(hog, id.vars = "OG_ID", variable.name = "Species", value.name = "genes",
               variable.factor = FALSE)
counts[, Counts := n_genes(genes)][, genes := NULL]
mo_genes <- members[, .(n = .N), by = OG_ID]
counts[Species == MO, Counts := mo_genes$n[match(OG_ID, mo_genes$OG_ID)]]

# M. oryzae genes in no HOG
pep <- "data-raw/source/downloads/Magnaporthe_oryzae.MG8.pep.all.fa.gz"
if (!file.exists(pep)) stop("Run prepare_data.R once to download ", basename(pep))
hdr <- grep("^>", readLines(gzfile(pep)), value = TRUE)
all_mo <- unique(strip_tx(sub("^>(\\S+).*", "\\1", hdr)))
un_file <- fread(file.path(dir, "Orthogroups", "Orthogroups_UnassignedGenes.tsv"), select = MO)[[1]]
un_file <- unique(strip_tx(trimws(un_file[un_file != ""])))
un <- setdiff(all_mo, members$ProtID)
message(sprintf("M. oryzae genes in no HOG: %d (%d unassigned by OrthoFinder, %d in an orthogroup but not in N0)",
                length(un), length(intersect(un, un_file)), length(setdiff(un, un_file))))
members <- rbind(members, data.table(OG_ID = paste0("unassigned_", un), ProtID = un))
counts <- rbind(counts, CJ(OG_ID = paste0("unassigned_", un), Species = species_cols)[
  , Counts := ifelse(Species == MO, 1L, 0L)])

missing <- setdiff(life$Species, species_cols)
if (length(missing)) stop("In Lifestyle_extended.tsv but not in N0.tsv: ", paste(missing, collapse = ", "))
counts <- merge(counts, life[, .(Species, Appresor, Lifestyle, Name = Full_name)], by = "Species")

out <- merge(counts, members, by = "OG_ID", allow.cartesian = TRUE)
setcolorder(out, c("OG_ID", "Species", "Counts", "Appresor", "Lifestyle", "Name", "ProtID"))
setorder(out, OG_ID, Species, ProtID)
fwrite(out, file.path(dir, "fin_ortho.tsv"), sep = "\t")

message(sprintf("fin_ortho.tsv: %s rows, %d HOGs with M. oryzae genes, %d M. oryzae genes, %d species",
                format(nrow(out), big.mark = ","), uniqueN(out$OG_ID), uniqueN(out$ProtID),
                uniqueN(out$Species)))
