# =============================================================================
# MagnaGenes v2.0
# Multi-omics explorer for Magnaporthe oryzae (rice blast fungus)
#
# An extension of MagnaGenes v1.0 by Andrew Foster
# (https://zenodo.org/records/4647766). Runs entirely in the browser with Shinylive (webR): no server, no AI services,
# no API keys. All data are bundled in data/ (built by data-raw/prepare_data.R).
# =============================================================================

library(shiny)
library(bslib)
library(plotly)
library(ggplot2)
library(reactable)

# Not called directly. These lines make the Shinylive exporter bundle packages
# that the WebAssembly builds of ggplot2/scales still require.
if (FALSE) {
  library(munsell)
  library(colorspace)
}

APP_VERSION <- "2.0"
V1_URL      <- "https://zenodo.org/records/4647766"   # MagnaGenes v1.0
TSL_URL     <- "https://www.tsl.ac.uk/"
OPENRICEBLAST_URL  <- "http://openriceblast.org"
OPENWHEATBLAST_URL <- "http://openwheatblast.net/"
CONTACT_EMAIL      <- "neha.sahu.tsl@gmail.com"
SITE_URL           <- "https://nehasahu486.github.io/magnagenes_v2/"

# Version history shown on the About page. Add new releases at the top.
VERSION_HISTORY <- list(
  list(version = "2.0", date = "November 2026", current = TRUE,
       notes = c("Redesigned interface with one search box for MGG IDs, gene names and keywords, a light/dark theme and a Reset button.",
                 "Gene profiles with UniProt, AlphaFold DB, Ensembl Fungi, InterPro, Pfam and PubMed links, a per-experiment expression heatmap, a labelled protein architecture diagram and protein/CDS FASTA export.",
                 "Transcriptomics: TPM line plots sized to each experiment, with replicate error bars where available and an optional box plot; row z-score heatmaps with optional TPM labels; gene-gene correlation per experiment on a colour range fitted to the selected genes.",
                 "Phosphoproteomics: one replicate bar chart per phosphopeptide (Guy11 above \u0394pmk1), heatmaps, fold-change statistics and protein site maps with labelled sites and Pfam domains.",
                 "Phylogenomics: ortholog counts across 50 fungi (OrthoFinder) on the species tree, in deeper shades of each lifestyle colour with more copies, with appressorium annotation.",
                 "GO enrichment (Fisher's exact test) with dot/lollipop plots, bar charts, treemaps, significance-vs-enrichment plots and gene-term matrices, in a colour of your choice.",
                 "Video tutorial with chapters, opened from the Tutorial button.",
                 "Runs fully in the browser with Shinylive; the AI chat module was removed.")),
  list(version = "1.0", date = "March 2021", current = FALSE, url = V1_URL,
       notes = c("First release of MagnaGenes by Andrew J. Foster and colleagues: an open science database of gene function studies in the blast fungus Magnaporthe oryzae."))
)

version_history_ui <- function(history = VERSION_HISTORY) {
  tags$ol(class = "version-list",
    lapply(history, function(v) {
      tags$li(class = if (isTRUE(v$current)) "version-item current" else "version-item",
        div(class = "version-head",
            span(class = "version-tag", paste0("v", v$version)),
            span(class = "version-date", v$date),
            if (isTRUE(v$current)) span(class = "chip ok", "current"),
            if (!is.null(v$url)) tags$a(href = v$url, target = "_blank", rel = "noopener", class = "xref", "Zenodo")),
        tags$ul(lapply(v$notes, tags$li)))
    }))
}
DATA_DIR    <- "data"

# Tables follow the light/dark theme through CSS variables
options(reactable.theme = reactableTheme(
  color = "var(--ink)", backgroundColor = "var(--surface)", borderColor = "var(--line)",
  stripedColor = "var(--surface-2)", highlightColor = "var(--chip-ok)",
  inputStyle = list(backgroundColor = "var(--surface)", color = "var(--ink)", borderColor = "var(--line)"),
  selectStyle = list(backgroundColor = "var(--surface)", color = "var(--ink)"),
  pageButtonHoverStyle = list(backgroundColor = "var(--surface-2)"),
  headerStyle = list(color = "var(--ink-2)", fontWeight = 600, fontSize = "0.88rem")))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# DATA
# =============================================================================

genes  <- readRDS(file.path(DATA_DIR, "genes.rds"))
expr   <- readRDS(file.path(DATA_DIR, "expression.rds"))
phos   <- readRDS(file.path(DATA_DIR, "phospho.rds"))
ortho  <- readRDS(file.path(DATA_DIR, "orthology.rds"))
go_db  <- readRDS(file.path(DATA_DIR, "go.rds"))

rownames(genes) <- genes$gene

# Sequences are large, so they are read only when first needed
.seq_cache <- new.env()
get_seqs <- function(type = c("protein", "cds")) {
  type <- match.arg(type)
  if (is.null(.seq_cache[[type]])) {
    .seq_cache[[type]] <- readRDS(file.path(DATA_DIR, paste0("seq_", type, ".rds")))
  }
  .seq_cache[[type]]
}

# ---- Gene labels and input parsing ------------------------------------------

first_name <- sub(" \\|.*$", "", genes$gene_name)
names(first_name) <- genes$gene

gene_short <- function(ids) {
  nm <- unname(first_name[ids])
  ifelse(is.na(nm), ids, paste(ids, nm))
}

gene_desc <- function(ids) {
  d <- genes[ids, "blast2go"]
  e <- genes[ids, "ensembl_description"]
  ifelse(is.na(d), ifelse(is.na(e), "", e), d)
}

name_lookup <- local({
  nm <- toupper(first_name)
  ok <- !is.na(nm) & !duplicated(nm)
  setNames(genes$gene[ok], nm[ok])
})

known_ids <- unique(c(genes$gene, rownames(expr$tpm), phos$meta$protein, ortho$members$gene))

# Accepts MGG IDs (MGG_01234, MGG01234, MGG_01234T0) and gene names (PMK1)
parse_gene_input <- function(txt) {
  tok <- unlist(strsplit(toupper(txt %||% ""), "[[:space:],;]+"))
  tok <- unique(tok[nzchar(tok)])
  tok <- sub("^MGG_?([0-9]{5})(T[0-9]+)?$", "MGG_\\1", tok)
  id  <- ifelse(tok %in% known_ids, tok, unname(name_lookup[tok]))
  list(ids = unique(id[!is.na(id)]), missing = tok[is.na(id)])
}

# ---- Example gene sets --------------------------------------------------------

EX_SIGNALLING <- c("MGG_09565", "MGG_12958", "MGG_14847", "MGG_00800", "MGG_04943",
                   "MGG_01822", "MGG_06368", "MGG_09898", "MGG_05871", "MGG_03860")
EX_PHYLO      <- c("MGG_09565", "MGG_12958", "MGG_06368", "MGG_09898", "MGG_05871",
                   "MGG_07219", "MGG_05059", "MGG_02252", "MGG_10315", "MGG_13863", "MGG_04795")
EX_PHOSPHO    <- "MGG_00800"
EX_GO_MODULE  <- "Module 7"
ex_go_genes   <- function(module = EX_GO_MODULE) genes$gene[which(genes$wgcna_module == module)]

# ---- Expression helpers -------------------------------------------------------

# Experiments, labels and paper links come from data-raw/metadata/rnaseq_experiments.csv
EXP_INFO  <- expr$experiments
EXP_SHORT <- setNames(EXP_INFO$label, EXP_INFO$experiment)
expr$conditions$short <- unname(EXP_SHORT[expr$conditions$experiment])

# Row z-score of log2(TPM + 1): each gene centred and scaled across the given columns
row_z <- function(m) {
  z <- t(scale(t(log2(m + 1))))
  z[!is.finite(z)] <- NA
  dimnames(z) <- dimnames(m)
  z
}

cluster_order <- function(m) {
  if (nrow(m) < 3) return(seq_len(nrow(m)))
  m[is.na(m)] <- 0
  stats::hclust(stats::dist(m), method = "average")$order
}

# ---- Phospho helpers ----------------------------------------------------------

phos$meta$site_label <- make.unique(paste0(
  ifelse(is.na(phos$meta$sites), "site ?", phos$meta$sites), "  ",
  ifelse(nchar(phos$meta$peptide) > 16, paste0(substr(phos$meta$peptide, 1, 15), "…"),
         phos$meta$peptide)), sep = " #")
phos_proteins <- sort(unique(phos$meta$protein))
phos_n <- table(phos$meta$protein)
HPI <- c(0, 1, 1.5, 2, 4, 6)

# ---- Pfam domains (e.g. "AMP-binding enzyme (PF00501) 79->373") ---------------

parse_pfam <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NULL)
  m <- regmatches(x, gregexpr("[^,]+? \\((PF[0-9]+)\\) [0-9]+->[0-9]+", x))[[1]]
  if (!length(m)) return(NULL)
  data.frame(
    name  = trimws(sub(" \\(PF[0-9]+\\).*$", "", m)),
    pfam  = sub(".*\\((PF[0-9]+)\\).*", "\\1", m),
    start = as.integer(sub(".* ([0-9]+)->.*", "\\1", m)),
    end   = as.integer(sub(".*->([0-9]+)$", "\\1", m)),
    stringsAsFactors = FALSE
  )
}

site_positions <- function(sites) {
  pos <- regmatches(sites, gregexpr("[STY][0-9]+", sites))
  lapply(pos, function(p) data.frame(residue = substr(p, 1, 1),
                                     pos = as.integer(substring(p, 2))))
}

# ---- GO helpers ---------------------------------------------------------------

go_terms <- go_db$terms
go_terms$idx <- seq_len(nrow(go_terms))
go_gene_terms <- split(go_db$term_idx, go_db$gene_idx)   # gene index -> term indices
names(go_gene_terms) <- go_db$genes[as.integer(names(go_gene_terms))]
ONT_NAMES <- c(BP = "Biological process", MF = "Molecular function", CC = "Cellular component")
ont_universe <- lapply(c(BP = "BP", MF = "MF", CC = "CC"), function(o) {
  keep <- go_terms$ont[go_db$term_idx] == o
  unique(go_db$genes[go_db$gene_idx[keep]])
})

go_enrich <- function(query, ont, background = NULL, min_size = 5, max_size = 500,
                      adjust = "BH") {
  universe <- ont_universe[[ont]]
  if (length(background)) universe <- intersect(universe, background)
  q <- intersect(query, universe)
  N <- length(universe); n <- length(q)
  if (n == 0) return(list(table = NULL, n_query = 0, n_universe = N, n_tested = 0))
  in_ont <- go_terms$ont == ont
  term_of <- function(ids) {
    t <- unlist(go_gene_terms[ids], use.names = FALSE)
    t[in_ont[t]]
  }
  K <- tabulate(term_of(universe), nbins = nrow(go_terms))
  k <- tabulate(term_of(q), nbins = nrow(go_terms))
  tested <- which(in_ont & K >= min_size & K <= max_size)
  p <- stats::phyper(k[tested] - 1, K[tested], N - K[tested], n, lower.tail = FALSE)
  res <- data.frame(
    go_id = go_terms$go[tested], term = go_terms$term[tested], ontology = ont,
    annotated = K[tested], significant = k[tested],
    expected = round(K[tested] * n / N, 2),
    fold_enrichment = round((k[tested] / n) / (K[tested] / N), 2),
    p_value = p, p_adj = stats::p.adjust(p, method = adjust),
    stringsAsFactors = FALSE
  )
  res <- res[res$significant > 0, ]
  res <- res[order(res$p_value, -res$fold_enrichment), ]
  # Query genes annotated to each reported term
  qt <- lapply(q, function(g) go_terms$go[intersect(go_gene_terms[[g]], tested)])
  names(qt) <- q
  pairs <- data.frame(gene = rep(q, lengths(qt)), go = unlist(qt), stringsAsFactors = FALSE)
  res$genes <- vapply(res$go_id, function(id) paste(pairs$gene[pairs$go == id], collapse = ", "), "")
  rownames(res) <- NULL
  list(table = res, pairs = pairs, n_query = n, n_universe = N, n_tested = length(tested))
}

# ---- External links -----------------------------------------------------------

url_uniprot   <- function(acc) paste0("https://www.uniprot.org/uniprotkb/", acc, "/entry")
url_alphafold <- function(acc) paste0("https://alphafold.ebi.ac.uk/entry/", acc)
url_interpro  <- function(acc) paste0("https://www.ebi.ac.uk/interpro/protein/UniProt/", acc, "/")
url_ensembl   <- function(g) paste0("https://fungi.ensembl.org/Multi/Search/Results?species=all;idx=;q=", g)
url_pfam      <- function(pf) paste0("https://www.ebi.ac.uk/interpro/entry/pfam/", pf, "/")
url_quickgo   <- function(go) paste0("https://www.ebi.ac.uk/QuickGO/term/", go)
url_pubmed    <- function(q) paste0("https://pubmed.ncbi.nlm.nih.gov/?term=", utils::URLencode(q))

xref_html <- function(href, label, cls = "xref") {
  ifelse(is.na(href), "",
         sprintf('<a href="%s" target="_blank" rel="noopener" class="%s">%s</a>', href, cls, label))
}

xref_links_html <- function(g) {
  r <- genes[g, ]
  paste0(
    xref_html(ifelse(is.na(r$uniprot_acc), NA, url_uniprot(r$uniprot_acc)), "UniProt"),
    xref_html(ifelse(is.na(r$alphafold_acc), NA, url_alphafold(r$alphafold_acc)), "AlphaFold"),
    xref_html(url_ensembl(r$gene), "Ensembl"),
    xref_html(ifelse(is.na(r$uniprot_acc), NA, url_interpro(r$uniprot_acc)), "InterPro")
  )
}

linkify <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NULL)
  parts <- strsplit(x, " \\| ")[[1]]
  tagList(lapply(seq_along(parts), function(i) {
    p <- parts[i]
    sep <- if (i < length(parts)) tags$br()
    if (grepl("^https?://", p)) {
      lab <- if (grepl("pubmed.ncbi.nlm.nih.gov/([0-9]+)", p))
        paste0("PMID ", sub(".*pubmed.ncbi.nlm.nih.gov/([0-9]+).*", "\\1", p)) else p
      tagList(tags$a(href = p, target = "_blank", rel = "noopener", lab), sep)
    } else tagList(p, sep)
  }))
}

# ---- FASTA --------------------------------------------------------------------

wrap_seq <- function(s, width = 60) {
  vapply(s, function(x) {
    if (is.na(x) || !nzchar(x)) return("")
    st <- seq(1, nchar(x), by = width)
    paste(substring(x, st, pmin(st + width - 1, nchar(x))), collapse = "\n")
  }, "", USE.NAMES = FALSE)
}

fasta_text <- function(ids, type = "protein", isoforms = FALSE) {
  s <- get_seqs(type)
  keep <- s$gene %in% ids
  if (!isoforms) keep <- keep & !duplicated(s$gene)
  if (!any(keep)) return("")
  g <- s$gene[keep]
  hdr <- paste0(">", s$transcript[keep], " gene=", g,
                ifelse(is.na(first_name[g]), "", paste0(" name=", first_name[g])),
                " ", gene_desc(g), " [", type, "; MG8]")
  paste0(hdr, "\n", wrap_seq(s$seq[keep]), collapse = "\n")
}

write_tsv <- function(df, file) {
  utils::write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

# =============================================================================
# VISUAL SYSTEM
# Categorical, sequential and diverging palettes follow a validated reference
# palette (colour-vision-deficiency safe in the order given).
# =============================================================================

BRAND <- "#1B9E77"
FONT  <- "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif"

# ---- Light / dark plot colours ------------------------------------------------
# INK, INK2, MUTED, GRID, AXIS, SURFACE, EMPTY and BACKBONE are active bindings:
# they return the colour for the current theme, which each plot sets with
# plot_setup() (server) before drawing. PDFs are always drawn in light mode.
THEME_COLS <- list(
  light = list(INK = "#1d2b2a", INK2 = "#52514e", MUTED = "#8a8984", GRID = "#ececea", AXIS = "#c9c8c2",
               SURFACE = "#ffffff", EMPTY = "#f2f2ef", BACKBONE = "#d9d8d3"),
  dark  = list(INK = "#e6ebe9", INK2 = "#b4bdb9", MUTED = "#8d9692", GRID = "#2f3835", AXIS = "#4d5853",
               SURFACE = "#1b2120", EMPTY = "#28302d", BACKBONE = "#4a5450")
)
theme_col <- function(name) THEME_COLS[[if (isTRUE(getOption("mg.dark"))) "dark" else "light"]][[name]]
.app_env <- environment()
for (.nm in names(THEME_COLS$light)) local({
  n <- .nm
  if (exists(n, envir = .app_env, inherits = FALSE)) rm(list = n, envir = .app_env)
  makeActiveBinding(n, function() theme_col(n), .app_env)
})

# Base font size for plots (set per plot from its width by plot_setup())
fsz <- function(delta = 0) getOption("mg.fs", 14) + delta

# ---- Palettes -----------------------------------------------------------------
SERIES   <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
SEQ_RAMP <- c("#f4f8fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
DIV_RAMP <- c("#184f95", "#3987e5", "#9ec5f4", "#f0efec", "#f3a9a8", "#e34948", "#a62b2b")
BWR_RAMP <- c("#184f95", "#3987e5", "#9ec5f4", "#ffffff", "#f3a9a8", "#e34948", "#a62b2b")  # low - white - high
VIRIDIS  <- c("#440154", "#482878", "#3e4a89", "#31688e", "#26828e", "#1f9e89", "#35b779", "#6dcd59",
              "#b4de2c", "#fde725")
PHOS_MARK <- "#eb6834"   # phosphosite markers on protein diagrams

# Phosphoproteomics: greens for Guy11, greys for delta-pmk1 (as in MagnaGenes v1.0)
SAMPLES <- c("Guy11", "Δpmk1")
COL_SAMPLE <- c("Guy11" = "#2e9a4f", "Δpmk1" = "#6b6b6b")
COL_REP <- list("Guy11" = c("#0f6b2f", "#2e9a4f", "#7cc38f"),
                "Δpmk1" = c("#333333", "#6b6b6b", "#a8a8a8"))
RAMP_SAMPLE <- list("Guy11" = c("#f1f9f3", "#c7e9cf", "#8fd0a1", "#4fae6c", "#1f8a45", "#0b5e2c"),
                    "Δpmk1" = c("#f6f6f6", "#d9d9d9", "#b0b0b0", "#7f7f7f", "#4f4f4f", "#222222"))

# Phylogenomics: lifestyle colours from MagnaGenes v1.0 (ColorBrewer Pastel2,
# alphabetical), plus a light cyan for Nematophagous (new in the extended species
# set; Pastel2's last colour is grey, which would read as "no data"). Cells are
# shaded from the pastel (1 copy) to deeper shades of the same hue (more copies).
LIFESTYLES <- c("Biotroph", "Commensal", "Endophyte", "Hemibiotroph", "Necrotroph", "Nematophagous",
                "Saprophyte", "Symbiont")
COL_LIFESTYLE <- setNames(c("#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#A8DDEB", "#FFF2AE", "#F1E2CC"),
                          LIFESTYLES)
# Same hue, rising chroma and falling lightness (HCL), so every shade reads as
# the lifestyle colour rather than drifting towards brown or grey
lifestyle_ramp <- function(col, n = 6) {
  luv <- grDevices::convertColor(t(grDevices::col2rgb(col) / 255), from = "sRGB", to = "Luv")
  h <- atan2(luv[3], luv[2]) * 180 / pi
  c0 <- sqrt(luv[2]^2 + luv[3]^2)
  grDevices::hcl(h, c = seq(c0, max(c0 + 45, 75), length.out = n), l = seq(luv[1], 50, length.out = n),
                 fixup = TRUE)
}
# Appressorium type (not coloured in v1.0): melanin-inspired tones. Species
# without an appressorium are left blank.
COL_APPRESS <- c(Hyaline = "#d8c8b0", Compound = "#9a7752", Melanized = "#4a3424")

colorscale <- function(cols) {
  n <- length(cols)
  lapply(seq_len(n), function(i) list((i - 1) / (n - 1), cols[i]))
}

# Light tint -> chosen colour, for user-coloured GO plots
tint_ramp <- function(col, n = 6) grDevices::colorRampPalette(c("#f4f6f8", col))(n)

ax <- function(title = NULL, ...) {
  modifyList(list(
    title = list(text = title, font = list(size = fsz(0), color = INK2)),
    gridcolor = GRID, zeroline = FALSE, showline = TRUE, linecolor = AXIS,
    ticks = "outside", tickcolor = AXIS, ticklen = 4,
    tickfont = list(size = fsz(-1), color = INK2), automargin = TRUE
  ), list(...))
}

# legend: TRUE/FALSE to show or hide it, or a list of legend settings
finish_plot <- function(p, file = "magnagenes_plot", legend = TRUE, ...) {
  leg <- list(orientation = "h", x = 0, y = 1.02, yanchor = "bottom",
              font = list(size = fsz(-1), color = INK2))
  if (is.list(legend)) { leg <- modifyList(leg, legend); legend <- TRUE }
  p |>
    layout(
      font = list(family = FONT, size = fsz(0), color = INK),
      paper_bgcolor = SURFACE, plot_bgcolor = SURFACE,
      hoverlabel = list(bgcolor = SURFACE, bordercolor = AXIS,
                        font = list(family = FONT, size = fsz(-1), color = INK)),
      showlegend = legend,
      legend = leg,
      ...
    ) |>
    config(displaylogo = FALSE,
           toImageButtonOptions = list(format = "svg", filename = file),
           modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d",
                                      "hoverCompareCartesian", "hoverClosestCartesian"))
}

empty_plot <- function(msg) {
  plot_ly(type = "scatter", mode = "markers") |>
    layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
           annotations = list(list(text = msg, showarrow = FALSE, xref = "paper", yref = "paper",
                                   x = 0.5, y = 0.5, font = list(size = fsz(1), color = MUTED))),
           paper_bgcolor = SURFACE, plot_bgcolor = SURFACE) |>
    config(displayModeBar = FALSE)
}

theme_mg <- function(base = 10) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          panel.border = element_rect(colour = AXIS, linewidth = 0.4, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", hjust = 0, colour = INK),
          strip.clip = "off",
          axis.ticks = element_line(colour = AXIS, linewidth = 0.3),
          axis.text = element_text(colour = INK2),
          legend.position = "bottom", legend.key = element_blank(),
          plot.title = element_text(face = "bold", colour = INK))
}

# Plain-ASCII text for PDF figures (the PDF device lacks subscripts and ellipses)
pdf_safe <- function(x) {
  x <- gsub("₁₀", "10", x)
  x <- gsub("₂", "2", x)
  x <- gsub("−", "-", x)
  x <- gsub("Δ", "delta-", x)
  gsub("…", "...", x)
}

# Labels as plotmath, so PDF strips show a real Delta and an italic gene name
to_plotmath <- function(lbl) vapply(lbl, function(s) {
  q <- function(t) paste0("'", gsub("'", "", t), "'")
  m <- regexpr("Δ([A-Za-z0-9]+)", s)
  if (m < 0) return(q(s))
  g <- sub("Δ", "", regmatches(s, m))
  pre <- substr(s, 1, m - 1); post <- substr(s, m + attr(m, "match.length"), nchar(s))
  paste0(if (nzchar(pre)) paste0(q(pre), "*") else "", "Delta*italic(", q(g), ")",
         if (nzchar(post)) paste0("*", q(post)) else "")
}, "", USE.NAMES = FALSE)
EXP_PARSE <- setNames(to_plotmath(EXP_SHORT), names(EXP_SHORT))

# PDF download handler for a ggplot (or a function that draws one). Light theme.
pdf_download <- function(make_plot, name, width, height) {
  downloadHandler(
    filename = function() paste0(name(), "_", format(Sys.Date(), "%Y%m%d"), ".pdf"),
    content = function(file) {
      old <- options(mg.dark = FALSE, mg.fs = 11)
      on.exit(options(old), add = TRUE)
      g <- make_plot()
      grDevices::pdf(file, width = width(), height = height(), useDingbats = FALSE)
      on.exit(grDevices::dev.off(), add = TRUE)
      if (is.function(g)) g() else print(g)
    }
  )
}

# =============================================================================
# UI HELPERS
# =============================================================================

feature_card <- function(target, icon_name, title, text) {
  tags$button(
    type = "button", class = "feature-card",
    onclick = sprintf("Shiny.setInputValue('nav_to', '%s', {priority: 'event'})", target),
    div(class = "feature-icon", icon(icon_name)),
    div(class = "feature-title", title),
    div(class = "feature-text", text)
  )
}

stat_tile <- function(value, label, sub = NULL) {
  div(class = "stat-tile",
      div(class = "stat-value", value),
      div(class = "stat-label", label),
      if (!is.null(sub)) div(class = "stat-sub", sub))
}

example_link <- function(id, label = "Load example") {
  actionLink(id, tagList(icon("flask"), label), class = "example-link")
}

download_group <- function(...) div(class = "download-group", ...)

pdf_size_inputs <- function(prefix, w = 10, h = 6) {
  div(class = "pdf-size",
      numericInput(paste0(prefix, "_w"), "PDF width (in)", w, min = 3, max = 30, step = 0.5),
      numericInput(paste0(prefix, "_h"), "Height (in)", h, min = 3, max = 30, step = 0.5))
}

section_note <- function(...) div(class = "section-note", ...)

# Tutorial pop-up: loads tutorial.mp4 and tutorial_chapters.json from the
# published site folder (the page that hosts the Shinylive app), and builds the
# chapter list with clickable timestamps.
tutorial_js <- "
(function () {
  var base;
  try { base = window.parent.location.href; } catch (e) { base = window.location.href; }
  var v = document.getElementById('tutorial_video');
  var list = document.getElementById('tutorial_chapters');
  var fail = function () {
    list.innerHTML = '<li class=\"section-note\">The tutorial video is available on the published MagnaGenes website.</li>';
    v.style.display = 'none';
  };
  v.onerror = fail;
  v.src = new URL('tutorial.mp4', base).href;
  var fmt = function (s) { s = Math.round(s); return Math.floor(s / 60) + ':' + ('0' + (s % 60)).slice(-2); };
  fetch(new URL('tutorial_chapters.json', base).href, { cache: 'no-cache' })
    .then(function (r) { if (!r.ok) throw new Error(r.status); return r.json(); })
    .then(function (ch) {
      list.innerHTML = '';
      ch.forEach(function (c) {
        var li = document.createElement('li');
        var a = document.createElement('a');
        a.href = '#'; a.innerHTML = '<span class=\"tutorial-time\">' + fmt(c.time) + '</span> ' + c.title;
        a.onclick = function (e) { e.preventDefault(); v.currentTime = c.time; v.play(); };
        li.appendChild(a); list.appendChild(li);
      });
      v.ontimeupdate = function () {
        var items = list.querySelectorAll('li');
        ch.forEach(function (c, i) {
          var next = ch[i + 1] ? ch[i + 1].time : Infinity;
          items[i].classList.toggle('active', v.currentTime >= c.time && v.currentTime < next);
        });
      };
    })
    .catch(fail);
})();
"

# =============================================================================
# STYLES
# =============================================================================

app_css <- HTML("
:root, [data-bs-theme='light'] {
  --brand: #1B9E77; --brand-dark: #137a5c; --brand-tint: rgba(225, 247, 236, 0.9);
  --brand-tint-solid: #e6f6ee; --brand-line: #cfeadc; --ink: #1d2b2a; --ink-2: #52514e; --muted: #8a8984;
  --line: #e3e3df; --surface: #ffffff; --surface-2: #fafaf8; --page: #f6f8f7; --frame: #1d2b2a;
  --chip-ok: #f1faf5; --chip-ok-line: #bfe6d3; --chip-warn: #fff8ec; --chip-warn-line: #f1c48a;
  --mono: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
}
[data-bs-theme='dark'] {
  --brand-dark: #5fd3ad; --brand-tint: rgba(27, 158, 119, 0.14); --brand-tint-solid: #17241f; --brand-line: #2a4a3e;
  --ink: #e6ebe9; --ink-2: #b4bdb9; --muted: #8d9692; --line: #303a37; --surface: #1b2120; --surface-2: #222a28;
  --page: #121716; --frame: #c9d1ce; --chip-ok: #18302a; --chip-ok-line: #2a5446; --chip-warn: #33291a; --chip-warn-line: #6b5230;
}
html { font-size: 16.5px; }
body { background: var(--page); color: var(--ink); }
.navbar { background: var(--brand) !important; border-bottom: 0; box-shadow: 0 1px 4px rgba(0,0,0,.12); }
.navbar .navbar-brand { color: #fff !important; font-weight: 700; letter-spacing: .01em; display:flex; align-items:center; gap:.5rem; }
.navbar .navbar-brand .brand-mark { width: 30px; height: 30px; border-radius: 7px; background: #fff; color: #137a5c; display:inline-flex; align-items:center; justify-content:center; font-size: 15px; font-weight: 800; }
.navbar .navbar-brand .brand-version { font-weight: 500; color: rgba(255,255,255,.75); font-size: .85em; }
.navbar .nav-link { color: rgba(255,255,255,.88) !important; font-weight: 500; }
.navbar .nav-link:hover { color: #fff !important; }
.navbar .nav-link.active { color: #fff !important; font-weight: 600; box-shadow: inset 0 -3px 0 #fff; }
.navbar .navbar-toggler { border-color: rgba(255,255,255,.6); }
.navbar bslib-input-dark-mode { --text-1: #fff; --text-2: rgba(255,255,255,.85); --vertical-correction: 2px; margin: 0 .35rem; }
.card { border: 1px solid var(--line); border-radius: 10px; box-shadow: none; background: var(--surface); }
.card-header { background: var(--surface); border-bottom: 1px solid var(--line); font-weight: 600; font-size: 1.02rem; }
.bslib-sidebar-layout > .sidebar { background: var(--brand-tint-solid); border-right: 1px solid var(--line); }
.bslib-sidebar-layout > .sidebar .sidebar-title { color: var(--brand-dark); font-weight: 700; font-size: 1.15rem; }
.sidebar .form-label, .sidebar label { font-weight: 600; font-size: .95rem; color: var(--ink); }
.sidebar .form-check-label, .sidebar .checkbox label, .sidebar .radio label { font-weight: 500; }
.sidebar textarea { font-family: var(--mono); font-size: .9rem; }
textarea::placeholder, input::placeholder, .selectize-input input::placeholder { color: var(--muted) !important; opacity: .55 !important; font-style: italic; }
.btn-primary { background: var(--brand); border-color: var(--brand); }
.btn-primary:hover, .btn-primary:focus { background: #137a5c; border-color: #137a5c; }
.btn-outline-primary { color: var(--brand-dark); border-color: var(--brand); }
.btn-outline-primary:hover { background: var(--brand); border-color: var(--brand); color: #fff; }
.download-group { display: flex; flex-direction: column; gap: .4rem; }
.download-group .btn { text-align: left; font-size: .9rem; }
.example-link { font-size: .9rem; color: var(--brand-dark); }
.section-note { font-size: .9rem; color: var(--ink-2); line-height: 1.5; }
.sidebar-heading { font-size: .78rem; text-transform: uppercase; letter-spacing: .08em; color: var(--muted); font-weight: 700; margin: .4rem 0 .2rem; }
.pdf-size { display: grid; grid-template-columns: 1fr 1fr; gap: .5rem; }
.pdf-size .form-group { margin-bottom: .3rem; }
.exp-choice a { font-size: .82rem; font-weight: 500; color: var(--brand-dark); margin-left: .25rem; }
.color-row { display: flex; align-items: center; gap: .6rem; margin-bottom: .8rem; }
.color-row input[type=color] { width: 44px; height: 32px; padding: 0; border: 1px solid var(--line); border-radius: 6px; background: var(--surface); }

/* ---- Home ---- */
.hero { position: relative; border-radius: 14px; overflow: hidden; min-height: 300px; color: #fff;
  background-image: linear-gradient(100deg, rgba(14,61,48,.92) 0%, rgba(19,122,92,.78) 45%, rgba(27,158,119,.15) 100%), url('appressoria_sem.jpg');
  background-size: cover; background-position: center; padding: 3rem 2.5rem; margin-bottom: 1.5rem; }
.hero h1 { font-weight: 800; font-size: 2.7rem; margin-bottom: .3rem; letter-spacing: -.01em; }
.hero h1 .v { font-weight: 400; opacity: .85; }
.hero .lead { max-width: 680px; font-size: 1.15rem; opacity: .95; }
.hero .hero-caption { position: absolute; right: 14px; bottom: 10px; font-size: .75rem; opacity: .8; }
.stat-row { display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: .9rem; margin-bottom: 1.5rem; }
.stat-tile { background: var(--surface); border: 1px solid var(--line); border-radius: 10px; padding: 1rem 1.1rem; }
.stat-value { font-size: 1.8rem; font-weight: 700; color: var(--ink); font-variant-numeric: tabular-nums; line-height: 1.1; }
.stat-label { font-size: .92rem; color: var(--ink-2); margin-top: .2rem; }
.stat-sub { font-size: .8rem; color: var(--muted); }
.feature-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 1rem; margin-bottom: 1.5rem; }
.feature-card { text-align: left; background: var(--brand-tint); border: 1px solid var(--brand-line); border-radius: 10px; padding: 1.1rem 1.2rem; cursor: pointer; transition: box-shadow .15s, transform .15s, background .15s; min-height: 44px; color: var(--ink); }
.feature-card:hover, .feature-card:focus { filter: brightness(0.97); box-shadow: 0 4px 10px rgba(0,0,0,.08); transform: translateY(-1px); outline: none; }
.feature-icon { width: 38px; height: 38px; border-radius: 9px; background: var(--brand); color: #fff; display: flex; align-items: center; justify-content: center; margin-bottom: .6rem; }
.feature-title { font-weight: 700; color: var(--brand-dark); margin-bottom: .25rem; font-size: 1.08rem; }
.feature-text { font-size: .93rem; color: var(--ink-2); line-height: 1.45; }
.home-split { display: grid; grid-template-columns: 1.1fr .9fr; gap: 1rem; margin-bottom: 1.5rem; }
@media (max-width: 900px) { .home-split { grid-template-columns: 1fr; } .hero { padding: 2rem 1.3rem; } .hero h1 { font-size: 2rem; } }
.infobox { background: var(--surface); border: 1px solid var(--line); border-radius: 10px; padding: 1.2rem 1.3rem; }
.infobox h4 { color: var(--brand-dark); font-size: 1.12rem; font-weight: 700; }
.credit-box { border-left: 4px solid var(--brand); }
.version-list { list-style: none; padding-left: 0; margin: 0; border-left: 2px solid var(--line); }
.version-item { position: relative; padding: 0 0 1rem 1.2rem; }
.version-item::before { content: ''; position: absolute; left: -7px; top: 4px; width: 12px; height: 12px; border-radius: 50%; background: var(--surface); border: 2px solid var(--muted); }
.version-item.current::before { background: var(--brand); border-color: var(--brand); }
.version-head { display: flex; flex-wrap: wrap; align-items: center; gap: .5rem; margin-bottom: .3rem; }
.version-tag { font-weight: 700; color: var(--ink); font-size: 1.05rem; }
.version-date { color: var(--muted); font-size: .9rem; }
.version-item ul { margin: 0; padding-left: 1.1rem; font-size: .93rem; color: var(--ink-2); }

/* ---- Results ---- */
.summary-strip { display: flex; flex-wrap: wrap; gap: .5rem; margin-bottom: .8rem; align-items: center; }
.chip { display: inline-flex; align-items: center; gap: .35rem; padding: .25rem .7rem; border-radius: 999px; background: var(--surface); border: 1px solid var(--line); font-size: .9rem; color: var(--ink-2); }
.chip b { color: var(--ink); font-variant-numeric: tabular-nums; }
.chip.warn { border-color: var(--chip-warn-line); background: var(--chip-warn); }
.chip.ok { border-color: var(--chip-ok-line); background: var(--chip-ok); }
.send-bar { display: flex; flex-wrap: wrap; gap: .4rem; align-items: center; margin-left: auto; }
.send-bar .btn { font-size: .85rem; padding: .2rem .6rem; }
.xref { display: inline-block; padding: .08rem .45rem; margin: 0 .25rem .2rem 0; border-radius: 5px; border: 1px solid var(--brand-line); background: var(--chip-ok); color: var(--brand-dark) !important; font-size: .8rem; font-weight: 600; text-decoration: none; }
.xref:hover { background: var(--brand); color: #fff !important; border-color: var(--brand); }
.gene-id { font-family: var(--mono); font-weight: 600; font-size: .9rem; color: var(--ink); }
.badge-effector { background: #fdeee6; color: #a3431b; border: 1px solid #f3c6ad; border-radius: 5px; font-size: .78rem; padding: .05rem .4rem; font-weight: 600; }
.rt-table { font-size: .92rem; }

/* ---- Experiment plots (one framed plot per dataset) ---- */
.exp-plot { border: 1.5px solid var(--frame); border-radius: 8px; padding: .7rem .8rem .3rem; margin-bottom: 1rem; background: var(--surface); max-width: 100%; }
.exp-grid { display: flex; flex-wrap: wrap; gap: 0 1rem; align-items: flex-start; }
.exp-grid .exp-plot { overflow-x: auto; }
.exp-plot-note { font-size: .82rem; color: var(--muted); margin: -.1rem 0 .3rem; }
.exp-plot-head { display: flex; flex-wrap: wrap; align-items: baseline; gap: .6rem; margin-bottom: .2rem; }
.exp-plot-head h5 { margin: 0; font-weight: 700; font-size: 1.08rem; color: var(--ink); }
.exp-plot-head .ref { font-size: .88rem; color: var(--ink-2); }
.exp-plot-head .ref a { color: var(--brand-dark); }
.plot-note { font-size: .88rem; color: var(--ink-2); line-height: 1.5; margin: .2rem 0 .8rem; }

/* ---- Gene profile ---- */
.profile-head { display: flex; flex-wrap: wrap; gap: .6rem; align-items: baseline; margin-bottom: .2rem; }
.profile-head .pid { font-family: var(--mono); font-size: 1.4rem; font-weight: 700; color: var(--ink); }
.profile-head .pname { font-size: 1.3rem; font-weight: 700; color: var(--brand-dark); }
.profile-desc { color: var(--ink-2); margin-bottom: .6rem; }
.facts { display: grid; grid-template-columns: 170px 1fr; gap: .35rem .9rem; font-size: .93rem; margin: 0; }
.facts dt { color: var(--muted); font-weight: 600; }
.facts dd { margin: 0; color: var(--ink); word-break: break-word; }
.profile-actions { display: flex; flex-wrap: wrap; gap: .5rem; margin: .6rem 0 .2rem; }
.seq-block { font-family: var(--mono); font-size: .82rem; background: var(--surface-2); border: 1px solid var(--line); border-radius: 8px; padding: .8rem; max-height: 260px; overflow: auto; white-space: pre; color: var(--ink); }
details > summary { cursor: pointer; color: var(--brand-dark); font-weight: 600; font-size: .95rem; }
.method-box { font-size: .9rem; color: var(--ink-2); background: var(--surface-2); border: 1px solid var(--line); border-radius: 8px; padding: .7rem .9rem; margin-top: .6rem; }
.footer-note { text-align: center; font-size: .82rem; color: var(--muted); padding: 1.2rem 0 1.6rem; }
.footer-note a { color: var(--brand-dark); }
.shiny-output-error-validation { color: var(--ink-2); font-size: 1rem; padding: 1.5rem; text-align: center; }
.peptide-box { max-height: 260px; overflow-y: auto; background: var(--surface); border: 1px solid var(--line); border-radius: 8px; padding: .4rem .6rem; font-family: var(--mono); font-size: .82rem; }
.peptide-box .checkbox { margin: .1rem 0; }
.peptide-seq { font-family: var(--mono); font-size: .9rem; color: var(--ink); word-break: break-all; margin: .1rem 0 .1rem; }
.recalculating { opacity: .45; transition: opacity .25s; }
.tutorial-wrap { display: grid; grid-template-columns: minmax(0, 1fr) 260px; gap: 1rem; align-items: start; }
@media (max-width: 900px) { .tutorial-wrap { grid-template-columns: 1fr; } }
.tutorial-wrap video { width: 100%; border-radius: 8px; border: 1px solid var(--line); background: #000; }
.tutorial-chapters ol { list-style: none; padding: 0; margin: 0; }
.tutorial-chapters li a { display: block; padding: .35rem .5rem; border-radius: 6px; color: var(--ink); text-decoration: none; font-size: .93rem; }
.tutorial-chapters li a:hover { background: var(--surface-2); }
.tutorial-chapters li.active a { background: var(--chip-ok); color: var(--brand-dark); font-weight: 600; }
.tutorial-time { font-family: var(--mono); color: var(--brand-dark); font-size: .85rem; margin-right: .3rem; }
")


# =============================================================================
# UI
# =============================================================================

n_species <- nrow(ortho$species)
n_go_terms <- nrow(go_terms)

home_panel <- nav_panel(
  title = "Home", value = "home", icon = icon("house"),
  div(class = "container-xxl py-3",
    div(class = "hero",
      h1("MagnaGenes ", span(class = "v", paste0("v", APP_VERSION))),
      p(class = "lead",
        "An integrated multi-omics resource for ", tags$em("Magnaporthe oryzae"),
        " (", tags$em("Pyricularia oryzae"), "), the rice blast fungus: gene annotation, ",
        "transcriptomics, phosphoproteomics, comparative genomics and Gene Ontology enrichment."),
      div(class = "d-flex flex-wrap gap-2 mt-3",
        actionButton("home_start", "Search genes", icon = icon("magnifying-glass"),
                     class = "btn-light fw-semibold"),
        actionButton("home_about", "About & credits", icon = icon("circle-info"),
                     class = "btn-outline-light")),
      div(class = "hero-caption", "Scanning electron micrograph: appressoria on a rice leaf")
    ),
    div(class = "stat-row",
      stat_tile(format(nrow(genes), big.mark = ","), "annotated genes", "MG8 / 70-15 reference"),
      stat_tile(nrow(expr$conditions), "RNA-seq conditions", paste(length(EXP_SHORT), "experiments")),
      stat_tile(format(nrow(phos$meta), big.mark = ","), "phosphopeptides",
                paste(format(length(phos_proteins), big.mark = ","), "proteins, Guy11 vs Δpmk1")),
      stat_tile(n_species, "fungal species", "orthology with species tree"),
      stat_tile(format(n_go_terms, big.mark = ","), "GO terms", "BP, MF and CC")
    ),
    div(class = "feature-grid",
      feature_card("search", "magnifying-glass", "Gene Search",
                   "Look up genes by MGG ID, gene name or keyword. Gene profiles with UniProt, AlphaFold DB, Ensembl and InterPro links, and FASTA export."),
      feature_card("transcriptomics", "chart-line", "Transcriptomics",
                   "Expression profiles, clustered heatmaps and gene-gene correlation across infection, surface and Pmk1 experiments."),
      feature_card("phospho", "atom", "Phosphoproteomics",
                   "Phosphosite abundance during appressorium development in Guy11 and Δpmk1, with fold-change heatmaps and protein site maps."),
      feature_card("phylo", "sitemap", "Phylogenomics",
                   paste("Ortholog counts across", n_species, "fungi on the species tree, annotated by lifestyle and appressorium type.")),
      feature_card("go", "diagram-project", "GO Enrichment",
                   "Over-representation analysis (Fisher's exact test, BH correction) with dot plots and gene-term matrices.")
    ),
    div(class = "home-split",
      div(class = "infobox",
        h4("About the organism"),
        p(tags$em("Magnaporthe oryzae"), " causes rice blast, one of the most destructive diseases of rice, and also blast diseases of wheat and other grasses. ",
          "Infection starts when a three-celled conidium lands on the leaf and forms a melanised appressorium. Enormous turgor in the appressorium drives a penetration peg into the plant cell, ",
          "after which invasive hyphae spread from cell to cell. The Pmk1 MAP kinase pathway is essential for appressorium development."),
        p(class = "section-note",
          "MagnaGenes brings together datasets that cover this infection cycle so a gene can be followed from annotation to expression, phosphorylation and evolutionary conservation.")
      ),
      div(class = "infobox text-center",
        tags$img(src = "infection_cycle.png", alt = "Infection cycle of Magnaporthe oryzae",
                 style = "max-width: 100%; max-height: 300px;"),
        div(class = "section-note mt-2", "The ", tags$em("M. oryzae"), " infection cycle on rice."))
    ),
    div(class = "infobox credit-box mb-3",
      h4("Citation"),
      p(class = "mb-0", "If you use MagnaGenes v.2.0 in your research, please cite: ",
        "Sahu N. (2026) MagnaGenes v.2.0: an integrated multi-omics explorer for the rice blast fungus ",
        tags$em("Magnaporthe oryzae", .noWS = "after"), ". ",
        tags$a(SITE_URL, href = SITE_URL, target = "_blank", rel = "noopener"))),
    div(class = "home-split",
      div(class = "infobox credit-box",
        h4("Information"),
        p("MagnaGenes is developed and maintained by the Nick Talbot Lab at ",
          tags$a("The Sainsbury Laboratory", href = TSL_URL, target = "_blank", rel = "noopener", .noWS = "after"),
          ". It aims to provide a one-stop platform for rice blast research and to promote open science."),
        p("MagnaGenes is freely available to the community as part of the ",
          tags$a("OpenRiceBlast", href = OPENRICEBLAST_URL, target = "_blank", rel = "noopener"), " and ",
          tags$a("OpenWheatBlast", href = OPENWHEATBLAST_URL, target = "_blank", rel = "noopener"),
          " initiatives. The database is updated regularly, and additions or corrections from the blast research community are welcome.")),
      div(class = "infobox",
        h4("Help us grow MagnaGenes"),
        p("Join the community of curators and help extend the MagnaGenes database. ",
          tags$a("Share your data and insights", target = "_blank", rel = "noopener",
                 href = "https://docs.google.com/spreadsheets/d/1Zx-qdqlpkeoYlIeAlOVXnHNzLs95kjwhQWP-uXLioP0/edit?gid=0#gid=0",
                 .noWS = "after"),
          "."),
        p("If something does not work as expected, if you spot an error in the data, or if you have a question or a suggestion, please write to ",
          tags$a(CONTACT_EMAIL, href = paste0("mailto:", CONTACT_EMAIL), .noWS = "after"), "."),
        p(class = "section-note mb-0",
          "Please include the gene IDs, the tab you were using and, if possible, a screenshot, so the problem can be reproduced."))
    )
  )
)

search_panel <- nav_panel(
  title = "Gene Search", value = "search", icon = icon("magnifying-glass"),
  layout_sidebar(
    fillable = FALSE,
    sidebar = sidebar(
      width = 330, title = "Query",
      textAreaInput("s_query", "Gene IDs, gene names or keywords", rows = 6, value = "",
                    placeholder = "e.g.\nmelanin\nmembrane\nMGG_09565, PMK1, MST12"),
      div(class = "section-note mb-2",
          "One entry per line or comma-separated. Gene IDs and names are matched exactly; anything else is searched as a keyword (case-insensitive). ",
          "Keywords are matched as written, not as natural language: use single words such as melanin or membrane."),
      example_link("s_example", "Example: pathogenicity signalling genes"),
      selectInput("s_field", "Search keywords in", choices = c(
        "All fields" = "all", "Gene name" = "gene_name", "Description" = "desc",
        "InterPro" = "interpro", "Pfam" = "pfam", "KEGG pathway" = "kegg",
        "Mutant phenotype" = "phenotype", "Process affected" = "process",
        "WGCNA module (exact)" = "wgcna_module")),
      checkboxInput("s_effector", "Effectors only", FALSE),
      checkboxInput("s_pheno", "Genes with a mutant phenotype only", FALSE),
      actionButton("s_go", "Search", icon = icon("magnifying-glass"), class = "btn-primary w-100"),
      div(class = "sidebar-heading mt-3", "Download results"),
      download_group(
        downloadButton("s_dl_tsv", "Annotations (TSV)", class = "btn-outline-primary btn-sm"),
        downloadButton("s_dl_pep", "Protein FASTA", class = "btn-outline-primary btn-sm"),
        downloadButton("s_dl_cds", "CDS FASTA", class = "btn-outline-primary btn-sm")),
      checkboxInput("s_isoforms", "Include all isoforms in FASTA", FALSE)
    ),
    uiOutput("s_summary"),
    navset_card_underline(
      id = "s_view",
      nav_panel("Results", value = "table", reactableOutput("s_table")),
      nav_panel("Gene profile", value = "profile",
        div(class = "d-flex flex-wrap gap-2 align-items-end mb-2",
            div(style = "min-width: 280px;",
                selectizeInput("s_profile_gene", "Gene", choices = NULL, width = "100%"))),
        uiOutput("s_profile"),
        card(card_header("Expression across conditions"),
             plotlyOutput("s_profile_expr", height = "640px")),
        card(card_header("Protein architecture"),
             plotlyOutput("s_profile_arch", height = "260px"))
      )
    )
  )
)

GENE_PLACEHOLDER <- "e.g.\nMGG_09565\nPMK1, MST12\nMGG_06368"

# Experiment check boxes, each with a link to its source paper
exp_choice_names <- lapply(seq_len(nrow(EXP_INFO)), function(i) {
  e <- EXP_INFO[i, ]
  span(class = "exp-choice", e$label,
       if (nzchar(e$url %||% "")) tags$a(href = e$url, target = "_blank", rel = "noopener",
                                         title = e$reference, icon("arrow-up-right-from-square"), "paper"))
})

tx_panel <- nav_panel(
  title = "Transcriptomics", value = "transcriptomics", icon = icon("chart-line"),
  layout_sidebar(
    fillable = FALSE,
    sidebar = sidebar(
      width = 340, title = "Expression query",
      textAreaInput("tx_genes", "Genes (IDs or names)", rows = 6, value = "", placeholder = GENE_PLACEHOLDER),
      example_link("tx_example", "Example: pathogenicity signalling genes"),
      checkboxGroupInput("tx_exps", "Experiments", choiceNames = exp_choice_names,
                         choiceValues = EXP_INFO$experiment, selected = EXP_INFO$experiment),
      div(class = "section-note mb-2",
          "Expression values come from RNA-seq data reanalysed from the publications linked above."),
      div(class = "sidebar-heading mt-2", "Plot options"),
      checkboxInput("tx_box", "Line plots: show box plot of all genes", FALSE),
      checkboxInput("tx_heat_values", "Heatmaps: show TPM values in cells", FALSE),
      checkboxInput("tx_cluster", "Cluster genes (heatmap, correlation)", TRUE),
      actionButton("tx_go", "Plot expression", icon = icon("chart-line"), class = "btn-primary w-100"),
      div(class = "sidebar-heading mt-3", "Download"),
      download_group(
        downloadButton("tx_dl_long", "Data, long format (TSV)", class = "btn-outline-primary btn-sm"),
        downloadButton("tx_dl_wide", "Gene x condition matrix (TSV)", class = "btn-outline-primary btn-sm"),
        downloadButton("tx_dl_pdf", "Current figure (PDF)", class = "btn-outline-primary btn-sm"),
        downloadButton("tx_dl_fasta", "Protein FASTA", class = "btn-outline-primary btn-sm")),
      pdf_size_inputs("tx_pdf", 11, 7),
      section_note("Use the camera icon on any chart to save a vector SVG.")
    ),
    uiOutput("tx_summary"),
    navset_card_underline(
      id = "tx_view",
      nav_panel("Expression profiles", value = "profiles",
        div(class = "plot-note",
            "One plot per experiment. Each line is a gene (TPM, transcripts per million). Error bars show the standard deviation ",
            "of biological replicates where replicate values are available; otherwise the points are averaged TPM from the biological replicates. ",
            "Tick \"show box plot\" in the sidebar to add the spread of all selected genes at each time point."),
        uiOutput("tx_profiles_ui")),
      nav_panel("Heatmap", value = "heatmap",
        div(class = "plot-note",
            "One heatmap per experiment. Values are row z-scores: each gene's log₂(TPM + 1) values are centred and scaled ",
            "within the experiment, so colour shows when a gene is high (red) or low (blue) relative to its own average. ",
            "Genes are in the same order in every plot so rows can be compared."),
        uiOutput("tx_heat_ui")),
      nav_panel("Gene-gene correlation", value = "corr",
        div(class = "plot-note",
            "How similar the expression patterns of each pair of genes are, calculated separately within each experiment. ",
            "Each cell is the Pearson correlation coefficient (r) of log₂(TPM + 1) values: r = 1 means the two genes rise and fall together, ",
            "r = 0 means no linear relationship and r = −1 means opposite patterns. ",
            "Genes are in the same order in every plot. The colour scale is shared and spans only the range of r in the genes you selected ",
            "(purple lowest, yellow highest), so small differences stay visible. Experiments with fewer than three conditions are not shown."),
        uiOutput("tx_corr_ui")),
      nav_panel("Data table", value = "table", reactableOutput("tx_table"))
    )
  )
)

ph_panel <- nav_panel(
  title = "Phosphoproteomics", value = "phospho", icon = icon("atom"),
  layout_sidebar(
    fillable = FALSE,
    sidebar = sidebar(
      width = 340, title = "Phosphoproteome query",
      selectizeInput("ph_protein", "Protein", choices = NULL,
                     options = list(placeholder = "e.g. MST7 or MGG_00800")),
      example_link("ph_example", "Example: MST7 (MAPKK upstream of Pmk1)"),
      div(class = "section-note mb-2 mt-2", textOutput("ph_protein_note", inline = TRUE)),
      uiOutput("ph_peptides_ui"),
      radioButtons("ph_scale", "Abundance scale", inline = TRUE,
                   choices = c("Linear" = "raw", "log₂" = "log"), selected = "raw"),
      checkboxInput("ph_zero_na", "Treat zero abundance as missing", TRUE),
      checkboxInput("ph_rowz", "Heatmap: row z-score", FALSE),
      div(class = "sidebar-heading mt-3", "Download"),
      download_group(
        downloadButton("ph_dl_data", "Replicate data (TSV)", class = "btn-outline-primary btn-sm"),
        downloadButton("ph_dl_stats", "Statistics (TSV)", class = "btn-outline-primary btn-sm"),
        downloadButton("ph_dl_pdf", "Current figure (PDF)", class = "btn-outline-primary btn-sm"),
        downloadButton("ph_dl_fasta", "Protein FASTA", class = "btn-outline-primary btn-sm")),
      pdf_size_inputs("ph_pdf", 11, 7)
    ),
    uiOutput("ph_header"),
    navset_card_underline(
      id = "ph_view",
      nav_panel("Phosphosite abundance", value = "profiles",
        div(class = "plot-note",
            "Normalised abundance of each phosphopeptide in Guy11 (greens) and Δpmk1 (greys) during appressorium development ",
            "on a hydrophobic surface. One plot per phosphopeptide, with Guy11 above Δpmk1 on the same scale; each bar is one biological replicate."),
        uiOutput("ph_profiles_ui")),
      nav_panel("Abundance heatmap", value = "heatmap",
        div(class = "plot-note", "Mean abundance of every phosphopeptide of this protein. Guy11 and Δpmk1 share one value range."),
        uiOutput("ph_heat_ui")),
      nav_panel("Δpmk1 vs Guy11", value = "fc", uiOutput("ph_fc_ui")),
      nav_panel("Protein site map", value = "map",
        div(class = "d-flex flex-wrap gap-3 align-items-end",
            selectInput("ph_map_time", "Fold change at", width = "200px",
                        choices = c("Largest change" = "max", setNames(HPI, paste(HPI, "hpi"))))),
        div(class = "plot-note",
            "Each lollipop is a localised phosphosite, labelled with its residue (for example pS358). Height is the log₂ fold change ",
            "in Δpmk1 relative to Guy11 at the chosen time point; colour shows whether the change is significant. ",
            "The track below shows the protein with its Pfam domains; triangles mark the sites."),
        plotlyOutput("ph_map", height = "520px"),
        uiOutput("ph_map_domains")),
      nav_panel("Statistics", value = "table", reactableOutput("ph_table"),
        div(class = "method-box",
            "Welch two-sample t-test per phosphopeptide and time point on log₂ abundance (n = 3 biological replicates per group), ",
            "Benjamini-Hochberg correction across all tests for this protein. With three replicates the power is limited, so treat p-values as a guide."))
    )
  )
)

og_panel <- nav_panel(
  title = "Phylogenomics", value = "phylo", icon = icon("sitemap"),
  layout_sidebar(
    fillable = FALSE,
    sidebar = sidebar(
      width = 340, title = "Orthology query",
      textAreaInput("og_genes", "Genes (IDs or names, up to 40)", rows = 6, value = "",
                    placeholder = "e.g.\nPMK1\nALB1, RSY1, BUF1\nMGG_13863"),
      example_link("og_example", "Example: signalling, melanin and effector genes"),
      radioButtons("og_tree_mode", "Tree", inline = TRUE,
                   choices = c("Phylogram" = "phylo", "Cladogram" = "clado")),
      checkboxInput("og_show_counts", "Show counts in cells", TRUE),
      actionButton("og_go", "Plot orthologs", icon = icon("sitemap"), class = "btn-primary w-100"),
      div(class = "sidebar-heading mt-3", "Download"),
      download_group(
        downloadButton("og_dl_counts", "Ortholog counts (TSV)", class = "btn-outline-primary btn-sm"),
        downloadButton("og_dl_members", "M. oryzae orthogroup members (TSV)", class = "btn-outline-primary btn-sm"),
        downloadButton("og_dl_fasta", "Orthogroup members, protein FASTA", class = "btn-outline-primary btn-sm"),
        downloadButton("og_dl_pdf", "Current figure (PDF)", class = "btn-outline-primary btn-sm")),
      pdf_size_inputs("og_pdf", 12, 9)
    ),
    uiOutput("og_summary"),
    navset_card_underline(
      id = "og_view",
      nav_panel("Ortholog profile", value = "profile",
        div(class = "plot-note",
            "How many orthologs of each query gene are found in ", n_species, " fungal genomes. Species are arranged on the species tree (left). ",
            "Each column is one ", tags$em("M. oryzae"), " gene and its orthogroup; each cell is the number of genes from that species in the orthogroup. ",
            "Cells take the colour of the species' lifestyle, in deeper shades of that colour with more copies; empty cells mean no ortholog. ",
            "The strips on the right show lifestyle and appressorium type (blank: no appressorium). Based on OrthoFinder hierarchical orthogroups across the ", n_species, " proteomes, at the root of the species tree (rooted on ",
            tags$em("Rhizophagus irregularis", .noWS = "after"), "), both generated for MagnaGenes v2.0. Large gene families, such as PTH11-like receptors or ",
            "polyketide synthases, form one group, so their counts include paralogs."),
        uiOutput("og_plot_ui")),
      nav_panel("Conservation by lifestyle", value = "lifestyle",
        div(class = "plot-note", "Share of species in each lifestyle group with at least one ortholog (",
            tags$em("M. oryzae"), " itself excluded)."),
        uiOutput("og_life_ui")),
      nav_panel("M. oryzae orthogroup members", value = "members", reactableOutput("og_members")),
      nav_panel("Data table", value = "table", reactableOutput("og_table"))
    )
  )
)

go_panel <- nav_panel(
  title = "GO Enrichment", value = "go", icon = icon("diagram-project"),
  layout_sidebar(
    fillable = FALSE,
    sidebar = sidebar(
      width = 350, title = "Enrichment settings",
      textAreaInput("go_genes", "Gene list (IDs or names)", rows = 6, value = "",
                    placeholder = "e.g. paste a list of MGG IDs\nMGG_09565\nMGG_12958\n..."),
      div(class = "d-flex flex-wrap gap-2 align-items-center mb-2",
          example_link("go_example", "Example: WGCNA module"),
          selectInput("go_example_module", NULL, width = "130px",
                      choices = paste("Module", 1:10), selected = EX_GO_MODULE)),
      radioButtons("go_ont", "Ontology", choices = setNames(names(ONT_NAMES), ONT_NAMES)),
      radioButtons("go_bg_mode", "Background", choices = c(
        "All annotated genes" = "all", "Custom background list" = "custom")),
      conditionalPanel("input.go_bg_mode == 'custom'",
        textAreaInput("go_bg", NULL, rows = 4, placeholder = "Paste background genes")),
      layout_columns(col_widths = c(6, 6),
        numericInput("go_min", "Min term size", 5, min = 1, max = 500),
        numericInput("go_max", "Max term size", 500, min = 10, max = 20000)),
      layout_columns(col_widths = c(6, 6),
        selectInput("go_adjust", "Correction", choices = c(
          "Benjamini-Hochberg" = "BH", "Bonferroni" = "bonferroni", "None" = "none")),
        numericInput("go_cutoff", "Adj. p cutoff", 0.05, min = 0, max = 1, step = 0.01)),
      actionButton("go_run", "Run enrichment", icon = icon("play"), class = "btn-primary w-100"),
      div(class = "sidebar-heading mt-3", "Plot options"),
      sliderInput("go_top", "Terms to plot", min = 5, max = 50, value = 20, step = 1),
      tags$label(class = "form-label", "Colour"),
      div(class = "color-row",
          tags$input(type = "color", id = "go_color_picker", value = SERIES[1],
                     oninput = "Shiny.setInputValue('go_color', this.value)",
                     `aria-label` = "Plot colour"),
          checkboxInput("go_viridis", "Use viridis instead", FALSE)),
      checkboxInput("go_lollipop", "Dot plot as lollipop (line + bubble)", FALSE),
      div(class = "sidebar-heading mt-3", "Download"),
      download_group(
        downloadButton("go_dl_table", "Enrichment table (TSV)", class = "btn-outline-primary btn-sm"),
        downloadButton("go_dl_pairs", "Gene-to-term pairs (TSV)", class = "btn-outline-primary btn-sm"),
        downloadButton("go_dl_pdf", "Current figure (PDF)", class = "btn-outline-primary btn-sm")),
      pdf_size_inputs("go_pdf", 10, 7)
    ),
    uiOutput("go_summary"),
    navset_card_underline(
      id = "go_view",
      nav_panel("Dot plot", value = "dot", uiOutput("go_dot_ui")),
      nav_panel("Bar chart", value = "bar", uiOutput("go_bar_ui")),
      nav_panel("Treemap", value = "treemap",
        div(class = "plot-note", "Each tile is a significant term. Tile area is the number of genes from your list in the term; colour is −log₁₀ adjusted p."),
        uiOutput("go_tree_ui")),
      nav_panel("Significance vs enrichment", value = "volcano",
        div(class = "plot-note", "All tested terms. Fold enrichment (x) against −log₁₀ adjusted p (y); bubble area is the number of genes. Terms above the dashed line pass the cutoff."),
        uiOutput("go_volcano_ui")),
      nav_panel("Gene-term matrix", value = "matrix", uiOutput("go_matrix_ui")),
      nav_panel("Results table", value = "table", reactableOutput("go_table"))
    ),
    div(class = "method-box",
        tags$strong("Method. "),
        "One-sided Fisher's exact test (hypergeometric) on GO annotations propagated to all ancestor terms ",
        "(the 'classic' algorithm), followed by multiple-testing correction. Only terms whose size falls within the chosen range are tested. ",
        "Fold enrichment = (genes in term / genes in list) / (term size / background size).")
  )
)

about_panel <- nav_panel(
  title = "About", value = "about", icon = icon("circle-info"),
  div(class = "container-xxl py-3", style = "max-width: 980px;",
    div(class = "infobox credit-box mb-3",
      h4("Credits"),
      p(tags$strong("MagnaGenes v2.0"), " is an extension of MagnaGenes v1.0 by Andrew Foster, available at ",
        tags$a("zenodo.org/records/4647766", href = V1_URL, target = "_blank", rel = "noopener"),
        " (Foster AJ, Were VM, Yan X, Win J, Harant A, Langner T, Kamoun S, Talbot NJ; doi:10.5281/zenodo.4647766)."),
      p("Version 2.0 adds an integrated gene search, gene profiles with ",
        "cross-references to UniProt, AlphaFold DB, Ensembl Fungi and InterPro, protein and CDS FASTA export, clustered heatmaps, ",
        "correlation analysis, phosphosite statistics and site maps, a species-tree ortholog profile, and a lightweight ",
        "GO enrichment engine. It runs entirely in the browser with Shinylive (webR), so no data leave your computer.")
    ),
    div(class = "infobox mb-3",
      h4("Version history"),
      version_history_ui()
    ),
    div(class = "infobox mb-3",
      h4("Data sources"),
      tags$ul(
        tags$li(tags$strong("Gene annotation: "), "MagnaGenes curated annotation (BLAST2GO, InterPro, Pfam, KEGG, mutant phenotypes, effector predictions, literature). ",
                "WGCNA co-expression modules from Yan et al. (2023), ", tags$em("The Plant Cell"), "."),
        tags$li(tags$strong("Genome and sequences: "), "Ensembl Fungi release 60, ", tags$em("M. oryzae"),
                " MG8 (strain 70-15) protein and CDS FASTA, chromosome coordinates and UniProt cross-references."),
        tags$li(tags$strong("Transcriptomics: "), "Mean TPM for Guy11 in planta (conidia to 144 hpi), Guy11 on hydrophobic and hydrophilic surfaces, ",
                "Δpmk1 on a hydrophobic surface, and the Pmk1 analogue-sensitive mutant (active vs inhibited)."),
        tags$li(tags$strong("Phosphoproteomics: "), "Normalised phosphopeptide abundances for Guy11 and Δpmk1 at 0, 1, 1.5, 2, 4 and 6 hpi on a hydrophobic surface (three biological replicates)."),
        tags$li(tags$strong("Orthology: "), "Hierarchical orthogroups (OrthoFinder 3, N0 level) across ", n_species, " fungal species, with the species tree rooted on Rhizophagus irregularis and lifestyle / appressorium annotation, generated for MagnaGenes v2.0 (not previously published)."),
        tags$li(tags$strong("Gene Ontology: "), "MagnaGenes GO annotation, with term names and the ontology graph from ", go_db$source, ".")
      )
    ),
    div(class = "infobox mb-3",
      h4("Methods"),
      tags$ul(
        tags$li("Heatmap clustering: average-linkage hierarchical clustering on Euclidean distances of the displayed values."),
        tags$li("Gene-gene correlation: Pearson correlation of log₂(TPM + 1) across the conditions of each experiment, calculated separately per experiment (experiments with at least three conditions)."),
        tags$li("Phosphosite positions: peptide start in the master protein plus the modified residue's position in the peptide. Ambiguous localisations are listed as ambiguous."),
        tags$li("Phospho statistics: Welch t-test on log₂ abundance per time point, Benjamini-Hochberg correction within the protein."),
        tags$li("GO enrichment: one-sided Fisher's exact test on propagated annotations, BH or Bonferroni correction.")
      )
    ),
    div(class = "infobox",
      h4("Contact"),
      p("For questions, problems or suggestions, please write to ",
        tags$a(CONTACT_EMAIL, href = paste0("mailto:", CONTACT_EMAIL), .noWS = "after"), "."),
      p(class = "section-note mb-0", "Version ", APP_VERSION, ". Built with R, Shiny, bslib, plotly, ggplot2 and reactable.")
    )
  )
)

ui <- page_navbar(
  id = "main_nav",
  title = tags$span("MagnaGenes",
                    tags$span(class = "brand-version", paste0("v", APP_VERSION))),
  window_title = "MagnaGenes v2.0",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = BRAND, success = BRAND,
                   base_font = font_collection("-apple-system", "BlinkMacSystemFont", "Segoe UI",
                                               "Roboto", "Helvetica Neue", "Arial", "sans-serif"),
                   "font-size-base" = "1rem"),
  fillable = FALSE,
  # Load plotly.js with the page so charts drawn inside tabs never race the script
  header = tags$head(htmltools::attachDependencies(
                       tags$span(), c(plot_ly()$dependencies, htmlwidgets::getDependency("plotly", "plotly"))),
                     tags$style(app_css),
                     tags$meta(name = "description", content = "MagnaGenes v2.0: multi-omics explorer for Magnaporthe oryzae")),
  footer = div(class = "footer-note",
               paste0("MagnaGenes v", APP_VERSION)),
  home_panel, search_panel, tx_panel, ph_panel, og_panel, go_panel,
  nav_spacer(),
  nav_item(actionLink("reset_app", tagList(icon("rotate-left"), "Reset"), class = "nav-link",
                      title = "Clear all searches and start again")),
  nav_item(actionLink("tutorial", tagList(icon("circle-play"), "Tutorial"), class = "nav-link",
                      title = "Video tour of MagnaGenes")),
  nav_item(input_dark_mode(id = "dark_mode", mode = "light")),
  about_panel
)

# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {

  # ---------------------------------------------------------------------------
  # Navigation
  # ---------------------------------------------------------------------------
  observeEvent(input$nav_to, nav_select("main_nav", input$nav_to))
  observeEvent(input$home_start, nav_select("main_nav", "search"))
  observeEvent(input$home_about, nav_select("main_nav", "about"))

  # Reset: start a fresh session (packages stay loaded, so this is quick)
  observeEvent(input$reset_app, session$reload())

  # Tutorial: video tour with clickable chapters. The video and its chapter list
  # sit next to the published site (see build_site.R), not in the app bundle, so
  # they are only downloaded when the tutorial is opened.
  observeEvent(input$tutorial, {
    showModal(modalDialog(
      title = tagList(icon("circle-play"), " MagnaGenes video tutorial"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close"),
      div(class = "tutorial-wrap",
          tags$video(id = "tutorial_video", controls = NA, preload = "metadata", playsinline = NA),
          div(class = "tutorial-chapters",
              div(class = "sidebar-heading", "Chapters"),
              tags$ol(id = "tutorial_chapters", tags$li(class = "section-note", "Loading...")))),
      tags$script(HTML(tutorial_js))
    ))
  })

  send_genes <- function(ids, target_input, target_tab, run_input = NULL) {
    updateTextAreaInput(session, target_input, value = paste(ids, collapse = "\n"))
    nav_select("main_nav", target_tab)
  }

  # ===========================================================================
  # GENE SEARCH
  # ===========================================================================
  search_index <- local({
    s <- function(x) ifelse(is.na(x), "", tolower(x))
    list(
      gene_name = s(genes$gene_name),
      desc = paste(s(genes$blast2go), s(genes$ensembl_description)),
      interpro = s(genes$interpro), pfam = s(genes$pfam), kegg = s(genes$kegg),
      phenotype = s(genes$phenotype), process = s(genes$process),
      wgcna_module = s(genes$wgcna_module),
      all = paste(tolower(genes$gene), s(genes$gene_name), s(genes$blast2go),
                  s(genes$ensembl_description), s(genes$interpro), s(genes$pfam),
                  s(genes$kegg), s(genes$phenotype), s(genes$process),
                  s(genes$wgcna_module), s(genes$uniprot_acc), sep = " ; ")
    )
  })

  observeEvent(input$s_example, {
    updateTextAreaInput(session, "s_query", value = paste(EX_SIGNALLING, collapse = "\n"))
  })

  # One box for everything: an entry whose words all resolve to gene IDs/names is
  # an ID lookup; any other entry is a keyword (phrase) search.
  s_result <- eventReactive(input$s_go, {
    entries <- trimws(unlist(strsplit(input$s_query %||% "", "[\n,;]+")))
    entries <- entries[nzchar(entries)]
    validate(need(length(entries) > 0,
                  "Type gene IDs, gene names or keywords in the box on the left, then press Search."))
    ids <- character(0); keywords <- character(0); missing <- character(0)
    for (e in entries) {
      p <- parse_gene_input(e)
      if (length(p$ids) && !length(p$missing)) {
        ids <- c(ids, p$ids)
        missing <- c(missing, setdiff(p$ids, genes$gene))
      } else keywords <- c(keywords, tolower(e))
    }
    hits <- character(0)
    if (length(keywords)) {
      hay <- search_index[[input$s_field %||% "all"]]
      found <- character(0)
      for (k in keywords) {
        # Multi-word keywords match genes that contain every word
        hit <- if (identical(input$s_field, "wgcna_module")) hay == k else
          Reduce(`&`, lapply(strsplit(k, "[[:space:]]+")[[1]], function(w) grepl(w, hay, fixed = TRUE)))
        h <- genes$gene[hit]
        if (length(h)) found <- c(found, k) else missing <- c(missing, k)
        hits <- c(hits, h)
      }
      keywords <- found
    }
    ids <- head(unique(c(intersect(ids, genes$gene), hits)), 5000)
    query <- if (length(keywords)) paste0("“", paste(keywords, collapse = "” or “"), "”")
    df <- genes[ids, , drop = FALSE]
    if (isTRUE(input$s_effector)) df <- df[df$effector, , drop = FALSE]
    if (isTRUE(input$s_pheno)) df <- df[!is.na(df$phenotype), , drop = FALSE]
    list(df = df, missing = missing, query = query)
  }, ignoreNULL = FALSE)

  output$s_summary <- renderUI({
    r <- s_result(); df <- r$df
    div(class = "summary-strip",
      span(class = "chip ok", icon("check"), tags$b(format(nrow(df), big.mark = ",")),
           if (nrow(df) == 1) "gene" else "genes", if (!is.null(r$query)) paste("found, including keyword matches for", r$query)),
      if (length(r$missing))
        span(class = "chip warn", icon("triangle-exclamation"), tags$b(length(r$missing)), "not found:",
             paste(head(r$missing, 6), collapse = ", "), if (length(r$missing) > 6) "…"),
      if (nrow(df)) span(class = "chip", tags$b(sum(df$effector)), "effectors"),
      if (nrow(df)) span(class = "chip", tags$b(sum(!is.na(df$phenotype))), "with mutant phenotype"),
      if (nrow(df)) div(class = "send-bar",
        span(class = "section-note", "Send to:"),
        actionButton("s_send_tx", "Transcriptomics", icon = icon("chart-line"), class = "btn-outline-primary"),
        actionButton("s_send_og", "Phylogenomics", icon = icon("sitemap"), class = "btn-outline-primary"),
        actionButton("s_send_go", "GO enrichment", icon = icon("diagram-project"), class = "btn-outline-primary"))
    )
  })

  observeEvent(input$s_send_tx, send_genes(head(s_result()$df$gene, 200), "tx_genes", "transcriptomics"))
  observeEvent(input$s_send_og, send_genes(head(s_result()$df$gene, 40), "og_genes", "phylo"))
  observeEvent(input$s_send_go, send_genes(s_result()$df$gene, "go_genes", "go"))

  output$s_table <- renderReactable({
    df <- s_result()$df
    validate(need(nrow(df) > 0, "No genes match this query."))
    tab <- data.frame(
      gene = df$gene,
      name = ifelse(is.na(first_name[df$gene]), "", first_name[df$gene]),
      description = gene_desc(df$gene),
      length = df$length_aa,
      effector = ifelse(df$effector, '<span class="badge-effector">effector</span>', ""),
      process = ifelse(is.na(df$process), "", df$process),
      module = ifelse(is.na(df$wgcna_module), "", df$wgcna_module),
      links = vapply(df$gene, xref_links_html, "", USE.NAMES = FALSE),
      stringsAsFactors = FALSE
    )
    rownames(tab) <- NULL
    reactable(
      tab, searchable = TRUE, highlight = TRUE, compact = TRUE, striped = FALSE,
      defaultPageSize = 15, showPageSizeOptions = TRUE, pageSizeOptions = c(15, 30, 100),
      selection = "single", onClick = "select",
      columns = list(
        gene = colDef(name = "Gene", minWidth = 110,
                      cell = function(v) tags$span(class = "gene-id", v)),
        name = colDef(name = "Name", minWidth = 80),
        description = colDef(name = "Description", minWidth = 240),
        length = colDef(name = "Length (aa)", minWidth = 80, align = "right"),
        effector = colDef(name = "", html = TRUE, minWidth = 70, sortable = FALSE),
        process = colDef(name = "Process affected", minWidth = 150),
        module = colDef(name = "WGCNA", minWidth = 80),
        links = colDef(name = "Databases", html = TRUE, minWidth = 230, sortable = FALSE)
      )
    )
  })

  observe({
    df <- s_result()$df
    ch <- if (nrow(df)) setNames(df$gene, gene_short(df$gene)) else character(0)
    updateSelectizeInput(session, "s_profile_gene", choices = ch, selected = ch[1] %||% character(0),
                         server = TRUE)
  })

  observeEvent(getReactableState("s_table", "selected"), {
    sel <- getReactableState("s_table", "selected")
    df <- s_result()$df
    if (length(sel) && sel <= nrow(df)) {
      g <- df$gene[sel]
      updateSelectizeInput(session, "s_profile_gene",
                           choices = setNames(df$gene, gene_short(df$gene)), selected = g, server = TRUE)
      nav_select("s_view", "profile")
    }
  })

  profile_gene <- reactive({
    g <- input$s_profile_gene
    req(g, g %in% genes$gene)
    g
  })

  output$s_profile <- renderUI({
    g <- profile_gene(); r <- genes[g, ]
    nm <- first_name[[g]]
    n_ph <- if (g %in% names(phos_n)) phos_n[[g]] else 0
    og <- ortho$members$og[match(g, ortho$members$gene)]
    n_sp <- if (!is.na(og)) sum(ortho$counts[og, ortho$species$code != "MGG"] > 0) else NA
    n_go <- length(go_gene_terms[[g]] %||% integer(0))
    loc <- if (!is.na(r$chr)) sprintf("Chromosome %s: %s-%s (%s strand)", r$chr,
                                      format(r$start, big.mark = ","), format(r$end, big.mark = ","), r$strand) else "Not available"
    pf <- parse_pfam(r$pfam)
    seqs <- get_seqs("protein")
    aa <- seqs$seq[match(g, seqs$gene)]
    fact <- function(label, value) if (!is.null(value) && length(value) && !all(is.na(value)) && !identical(value, ""))
      tagList(tags$dt(label), tags$dd(value))
    card(
      card_body(
        div(class = "profile-head",
            span(class = "pid", g),
            if (!is.na(nm)) span(class = "pname", nm),
            if (isTRUE(r$effector)) span(class = "badge-effector", "predicted effector")),
        div(class = "profile-desc", gene_desc(g)),
        div(HTML(xref_links_html(g)),
            HTML(xref_html(url_pubmed(paste(g, if (!is.na(nm)) paste("OR", nm), "AND (Magnaporthe OR Pyricularia)")), "PubMed search"))),
        layout_columns(
          col_widths = c(6, 6), class = "mt-3",
          tags$dl(class = "facts",
            fact("Location (MG8)", loc),
            fact("Transcript", r$transcript),
            fact("Protein length", if (!is.na(r$length_aa)) paste(r$length_aa, "aa")),
            fact("UniProt", if (!is.na(r$uniprot_acc)) tags$a(href = url_uniprot(r$uniprot_acc), target = "_blank", r$uniprot_acc)),
            fact("AlphaFold DB", if (!is.na(r$alphafold_acc)) tags$a(href = url_alphafold(r$alphafold_acc), target = "_blank", paste("AF", r$alphafold_acc, "model")) else "No model listed"),
            fact("WGCNA module", r$wgcna_module),
            fact("KEGG pathway", r$kegg)),
          tags$dl(class = "facts",
            fact("InterPro", r$interpro),
            fact("Pfam domains", if (!is.null(pf)) tagList(lapply(seq_len(nrow(pf)), function(i)
              tagList(tags$a(href = url_pfam(pf$pfam[i]), target = "_blank", pf$pfam[i]),
                      sprintf(" %s (%d-%d)", pf$name[i], pf$start[i], pf$end[i]), tags$br())))),
            fact("Mutant phenotype", r$phenotype),
            fact("Process affected", r$process),
            fact("References", linkify(r$reference)),
            fact("Other publications", linkify(r$other_reference)))
        ),
        div(class = "profile-actions",
          actionButton("s_open_tx", "Expression", icon = icon("chart-line"), class = "btn-outline-primary btn-sm"),
          if (n_ph > 0) actionButton("s_open_ph", sprintf("%d phosphopeptides", n_ph), icon = icon("atom"),
                                     class = "btn-outline-primary btn-sm"),
          if (!is.na(n_sp)) actionButton("s_open_og", sprintf("Orthologs in %d of %d species", n_sp, n_species - 1),
                                         icon = icon("sitemap"), class = "btn-outline-primary btn-sm"),
          span(class = "chip", tags$b(n_go), "GO terms (incl. ancestors)"),
          downloadButton("s_dl_gene_pep", "Protein FASTA", class = "btn-outline-primary btn-sm"),
          downloadButton("s_dl_gene_cds", "CDS FASTA", class = "btn-outline-primary btn-sm")),
        if (!is.na(aa)) tags$details(class = "mt-2",
          tags$summary(sprintf("Protein sequence (%d aa)", nchar(aa))),
          div(class = "seq-block mt-2", paste0(">", g, "\n", wrap_seq(aa))))
      )
    )
  })

  observeEvent(input$s_open_tx, send_genes(profile_gene(), "tx_genes", "transcriptomics"))
  observeEvent(input$s_open_og, send_genes(profile_gene(), "og_genes", "phylo"))
  observeEvent(input$s_open_ph, {
    updateSelectizeInput(session, "ph_protein", selected = profile_gene(),
                         choices = setNames(phos_proteins, gene_short(phos_proteins)), server = TRUE)
    nav_select("main_nav", "phospho")
  })

  output$s_dl_gene_pep <- downloadHandler(
    filename = function() paste0(profile_gene(), "_protein.fasta"),
    content = function(file) writeLines(fasta_text(profile_gene(), "protein", TRUE), file))
  output$s_dl_gene_cds <- downloadHandler(
    filename = function() paste0(profile_gene(), "_cds.fasta"),
    content = function(file) writeLines(fasta_text(profile_gene(), "cds", TRUE), file))

  # Font size that follows the width of a plot output (px), within limits
  responsive_font <- function(id, base = 15, min = 12, max = 19) {
    w <- session$clientData[[paste0("output_", id, "_width")]] %||% 900
    max(min, min(max, base * w / 900))
  }

  # Call first in every plot: sets the light/dark colours and the font size
  # fixed: use this font size instead of scaling with width (for plots sized to their content)
  plot_setup <- function(id, base = 15, fixed = NULL) {
    options(mg.dark = identical(input$dark_mode, "dark"),
            mg.fs = if (is.null(fixed)) responsive_font(id, base) else fixed)
    invisible(fsz(0))
  }


  fmt_tpm <- function(x) ifelse(is.na(x), "",
                                ifelse(x >= 1000, format(round(x), big.mark = ",", trim = TRUE),
                                       ifelse(x >= 100, sprintf("%.0f", x),
                                              ifelse(x >= 10, sprintf("%.1f", x), sprintf("%.2f", x)))))

  # One heatmap row per experiment, each coloured on its own scale
  # (log2 TPM rescaled 0-1 within the experiment: blue low, white mid, red high)
  output$s_profile_expr <- renderPlotly({
    plot_setup("s_profile_expr")
    g <- profile_gene()
    if (!g %in% rownames(expr$tpm)) return(empty_plot("No expression data for this gene"))
    fs <- fsz(0)
    cd <- expr$conditions
    exps <- unique(cd$experiment)
    n_max <- max(table(cd$experiment))
    panels <- lapply(seq_along(exps), function(i) {
      k <- cd$experiment == exps[i]
      tpm <- unname(expr$tpm[g, cd$id[k]])
      l <- log2(tpm + 1)
      rng <- range(l, na.rm = TRUE)
      sc <- if (diff(rng) > 0) (l - rng[1]) / diff(rng) else rep(0.5, length(l))
      plot_ly(x = cd$sample[k], y = EXP_SHORT[[exps[i]]], z = matrix(sc, nrow = 1),
              type = "heatmap", zmin = 0, zmax = 1, colorscale = colorscale(BWR_RAMP), showscale = FALSE,
              text = matrix(fmt_tpm(tpm), nrow = 1), texttemplate = "%{text}",
              textfont = list(size = fs, family = FONT),
              xgap = 3, ygap = 0,
              hovertemplate = paste0("<b>", EXP_SHORT[[exps[i]]], "</b><br>%{x}<br>TPM %{text}<extra></extra>")) |>
        layout(xaxis = list(type = "category", range = c(-0.5, n_max - 0.5), showgrid = FALSE, zeroline = FALSE,
                            showline = FALSE, ticks = "", tickfont = list(size = fs, color = INK2), fixedrange = TRUE),
               yaxis = list(visible = FALSE, fixedrange = TRUE))
    })
    sp <- subplot(panels, nrows = length(exps), shareX = FALSE, shareY = FALSE, margin = 0)
    # Equal rows: each band holds a title, the cells and the sample labels
    h <- 1 / length(exps)
    dom <- setNames(lapply(seq_along(exps), function(i)
      list(domain = c(1 - i * h + 0.4 * h, 1 - (i - 1) * h - 0.2 * h), visible = FALSE, fixedrange = TRUE)),
      paste0("yaxis", c("", seq_along(exps)[-1])))
    sp <- do.call(layout, c(list(sp), dom))
    ann <- lapply(seq_along(exps), function(i) {
      sfx <- if (i == 1) "" else i
      list(text = paste0("<b>", EXP_SHORT[[exps[i]]], "</b>"), xref = paste0("x", sfx, " domain"),
           yref = paste0("y", sfx, " domain"), x = 0, y = 1.04, xanchor = "left", yanchor = "bottom",
           showarrow = FALSE, font = list(size = fs + 1, color = INK))
    })
    key <- sprintf(paste0("Cell values: TPM. Colour is scaled within each experiment: ",
                          "<span style='color:%s'>■</span> low  ",
                          "<span style='color:#d0d0cc'>■</span> mid  ",
                          "<span style='color:%s'>■</span> high"), BWR_RAMP[1], BWR_RAMP[7])
    ann[[length(ann) + 1]] <- list(text = key, xref = "paper", yref = "paper", x = 0, y = -0.06,
                                   xanchor = "left", yanchor = "top", showarrow = FALSE,
                                   font = list(size = fs - 1, color = INK2))
    sp |> finish_plot(paste0(g, "_expression_heatmap"), legend = FALSE, annotations = ann,
                      margin = list(t = 30, b = 80, l = 10, r = 10))
  })

  output$s_profile_arch <- renderPlotly({
    plot_setup("s_profile_arch")
    g <- profile_gene()
    protein_architecture_plot(g, width = session$clientData$output_s_profile_arch_width %||% 900, fs = fsz(-1))
  })

  s_ids_for_download <- reactive(s_result()$df$gene)
  output$s_dl_tsv <- downloadHandler(
    filename = function() paste0("MagnaGenes_search_", format(Sys.Date(), "%Y%m%d"), ".tsv"),
    content = function(file) {
      df <- s_result()$df
      df$uniprot_url <- ifelse(is.na(df$uniprot_acc), NA, url_uniprot(df$uniprot_acc))
      df$alphafold_url <- ifelse(is.na(df$alphafold_acc), NA, url_alphafold(df$alphafold_acc))
      df$ensembl_url <- url_ensembl(df$gene)
      write_tsv(df, file)
    })
  output$s_dl_pep <- downloadHandler(
    filename = function() paste0("MagnaGenes_search_protein_", format(Sys.Date(), "%Y%m%d"), ".fasta"),
    content = function(file) writeLines(fasta_text(s_ids_for_download(), "protein", isTRUE(input$s_isoforms)), file))
  output$s_dl_cds <- downloadHandler(
    filename = function() paste0("MagnaGenes_search_cds_", format(Sys.Date(), "%Y%m%d"), ".fasta"),
    content = function(file) writeLines(fasta_text(s_ids_for_download(), "cds", isTRUE(input$s_isoforms)), file))

  # Protein architecture: backbone, Pfam domains labelled with their names,
  # and detected phosphosites as lollipops above the backbone
  protein_architecture_plot <- function(g, width = 900, fs = 13) {
    len <- genes[g, "length_aa"]
    if (is.na(len)) return(empty_plot("Protein length not available"))
    pf <- parse_pfam(genes[g, "pfam"])
    sites <- phos$meta[phos$meta$protein == g, ]
    plot_w <- max(200, width - 60)
    x_pad <- len * 0.04
    px_per_aa <- plot_w / (len + 2 * x_pad)
    shapes <- list(list(type = "rect", x0 = 1, x1 = len, y0 = -0.07, y1 = 0.07, fillcolor = BACKBONE,
                        line = list(width = 0), layer = "below"))
    ann <- list(
      list(text = "<b>N</b>", x = 1, y = 0, xanchor = "right", xshift = -6, showarrow = FALSE,
           font = list(size = fs, color = INK2)),
      list(text = "<b>C</b>", x = len, y = 0, xanchor = "left", xshift = 6, showarrow = FALSE,
           font = list(size = fs, color = INK2)))
    p <- plot_ly()
    if (!is.null(pf)) {
      dom_names <- unique(pf$name)
      dcol <- setNames(if (length(dom_names) <= 8) SERIES[seq_along(dom_names)] else rep(SERIES[1], length(dom_names)),
                       dom_names)
      light <- function(col) { v <- grDevices::col2rgb(col) / 255; sum(c(0.299, 0.587, 0.114) * v) > 0.6 }
      char_px <- fs * 0.6
      lane_end <- c(-Inf, -Inf, -Inf)            # right edge (px) of the last label in each lane
      lane_y <- c(-0.62, -0.88, -1.14)
      to_px <- function(x) (x + x_pad) * px_per_aa
      for (i in seq_len(nrow(pf))) {
        col <- dcol[[pf$name[i]]]
        shapes[[length(shapes) + 1]] <- list(type = "rect", x0 = pf$start[i], x1 = pf$end[i], y0 = -0.3, y1 = 0.3,
                                             fillcolor = col, line = list(color = SURFACE, width = 1.5), layer = "above")
        fit <- floor((pf$end[i] - pf$start[i]) * px_per_aa / char_px) - 1
        nm <- pf$name[i]
        mid <- (pf$start[i] + pf$end[i]) / 2
        if (fit >= 6) {
          lab <- if (nchar(nm) > fit) paste0(substr(nm, 1, fit - 1), "…") else nm
          ann[[length(ann) + 1]] <- list(text = paste0("<b>", lab, "</b>"), x = mid, y = 0, showarrow = FALSE,
                                         font = list(size = fs - 1, color = if (light(col)) INK else "#ffffff"))
        } else {
          # Too narrow: label below with a leader line, in the first lane where it does not overlap
          w <- nchar(nm) * (fs - 2) * 0.6
          x0 <- to_px(mid) - w / 2
          lane <- which(x0 > lane_end + 8)[1]
          if (is.na(lane)) lane <- which.min(lane_end)
          lane_end[lane] <- x0 + w
          yl <- lane_y[lane]
          ann[[length(ann) + 1]] <- list(text = nm, x = mid, y = -0.32, ax = mid, ay = yl, axref = "x", ayref = "y",
                                         showarrow = TRUE, arrowhead = 0, arrowwidth = 1, arrowcolor = MUTED,
                                         yanchor = "top", font = list(size = fs - 2, color = INK2))
        }
      }
      p <- p |> add_trace(x = (pf$start + pf$end) / 2, y = rep(0, nrow(pf)), type = "scatter", mode = "markers",
                          marker = list(size = 22, opacity = 0), showlegend = FALSE, hoverinfo = "text",
                          text = sprintf("<b>%s</b><br>%s<br>%d-%d aa", pf$name, pf$pfam, pf$start, pf$end))
    }
    n_sites <- 0
    if (nrow(sites)) {
      sp <- site_positions(sites$sites)
      sd <- do.call(rbind, lapply(seq_along(sp), function(i)
        if (nrow(sp[[i]])) cbind(sp[[i]], peptide = sites$peptide[i]) else NULL))
      if (!is.null(sd) && nrow(sd)) {
        sd <- sd[!duplicated(sd$pos), ]
        sd <- sd[order(sd$pos), ]
        n_sites <- nrow(sd)
        # Raise every other marker when sites are closer than one marker width
        min_gap <- (fs + 6) / px_per_aa
        sd$y <- 0.85
        if (nrow(sd) > 1) for (j in 2:nrow(sd))
          if (sd$pos[j] - sd$pos[j - 1] < min_gap && sd$y[j - 1] == 0.85) sd$y[j] <- 1.12
        p <- p |>
          add_segments(x = sd$pos, xend = sd$pos, y = 0.3, yend = sd$y - 0.07, showlegend = FALSE,
                       line = list(color = MUTED, width = 1), hoverinfo = "none") |>
          add_trace(x = sd$pos, y = sd$y, type = "scatter", mode = "markers+text",
                    text = "P", textfont = list(size = max(8, fs - 5), color = "#ffffff", family = FONT),
                    textposition = "middle center", showlegend = FALSE,
                    marker = list(size = fs + 4, color = PHOS_MARK, line = list(color = SURFACE, width = 1.5)),
                    hovertext = sprintf("<b>p%s%d</b><br>%s", sd$residue, sd$pos, sd$peptide), hoverinfo = "text")
      }
    }
    if (is.null(pf)) ann[[length(ann) + 1]] <- list(text = "No Pfam domains annotated", x = len / 2, y = -0.5,
                                                    showarrow = FALSE, font = list(size = fs - 1, color = MUTED))
    key <- paste0(g, " · ", format(len, big.mark = ","), " aa · ",
                  if (is.null(pf)) "0" else nrow(pf), " Pfam domain(s) · ", n_sites, " phosphosite(s)")
    ann[[length(ann) + 1]] <- list(text = key, x = 0, y = 1.02, xref = "paper", yref = "paper", xanchor = "left",
                                   yanchor = "bottom", showarrow = FALSE, font = list(size = fs - 1, color = INK2))
    p |> finish_plot(paste0(g, "_architecture"), legend = FALSE, shapes = shapes, annotations = ann,
                     xaxis = ax("Residue", range = c(-x_pad, len + x_pad), showgrid = FALSE, zeroline = FALSE,
                                tickfont = list(size = fs - 1, color = INK2)),
                     yaxis = list(visible = FALSE, range = c(-1.35, 1.3), fixedrange = TRUE),
                     margin = list(t = 30, b = 50, l = 20, r = 20))
  }

  # ===========================================================================
  # TRANSCRIPTOMICS
  # ===========================================================================
  observeEvent(input$tx_example, updateTextAreaInput(session, "tx_genes", value = paste(EX_SIGNALLING, collapse = "\n")))

  tx_query <- eventReactive(input$tx_go, {
    p <- parse_gene_input(input$tx_genes)
    validate(need(length(p$ids) + length(p$missing) > 0,
                  "Type gene IDs or names in the box on the left (or load the example), then press Plot expression."))
    ids <- head(intersect(p$ids, rownames(expr$tpm)), 200)
    validate(need(length(ids) > 0, "None of the genes were found in the expression data."))
    validate(need(length(input$tx_exps) > 0, "Select at least one experiment."))
    exps <- EXP_INFO$experiment[EXP_INFO$experiment %in% input$tx_exps]
    cd <- expr$conditions[expr$conditions$experiment %in% exps, ]
    list(ids = ids, exps = exps, missing = c(p$missing, setdiff(p$ids, rownames(expr$tpm))), cond = cd,
         m = expr$tpm[ids, cd$id, drop = FALSE],
         sd = if (!is.null(expr$sd)) expr$sd[ids, cd$id, drop = FALSE] else NULL)
  }, ignoreNULL = FALSE)

  # Row z-scores of log2(TPM + 1), computed within each experiment
  tx_z <- reactive({
    q <- tx_query()
    v <- q$m
    for (e in q$exps) {
      k <- q$cond$experiment == e
      v[, k] <- row_z(q$m[, k, drop = FALSE])
    }
    v
  })

  tx_long <- reactive({
    q <- tx_query(); z <- tx_z()
    n <- length(q$m)
    data.frame(gene = rep(q$ids, times = ncol(z)),
               name = rep(ifelse(is.na(first_name[q$ids]), "", first_name[q$ids]), times = ncol(z)),
               experiment = rep(q$cond$experiment, each = nrow(z)),
               sample = rep(q$cond$sample, each = nrow(z)),
               TPM = as.vector(q$m),
               SD = if (is.null(q$sd)) rep(NA_real_, n) else as.vector(q$sd),
               log2_TPM_plus_1 = round(log2(as.vector(q$m) + 1), 4),
               row_z_score = round(as.vector(z), 4),
               stringsAsFactors = FALSE)
  })

  output$tx_summary <- renderUI({
    q <- tx_query()
    div(class = "summary-strip",
      span(class = "chip ok", tags$b(length(q$ids)), "genes"),
      span(class = "chip", tags$b(nrow(q$cond)), "conditions in", tags$b(length(q$exps)), "experiments"),
      if (length(q$missing)) span(class = "chip warn", icon("triangle-exclamation"), tags$b(length(q$missing)),
                                  "not found:", paste(head(q$missing, 6), collapse = ", ")),
      if (length(q$ids) > 8) span(class = "chip", icon("circle-info"),
                                  "Line plots show the first 8 genes; box plots and heatmaps use all."))
  })

  exp_has_reps <- function(e) isTRUE(EXP_INFO$has_replicates[match(e, EXP_INFO$experiment)])
  exp_n_samples <- function(e) sum(expr$conditions$experiment == e)

  # Plot widths follow the number of samples, so few time points are not stretched
  tx_prof_width <- function(e, ids) {
    leg <- 40 + 8.5 * max(nchar(gene_short(head(ids, 8))), 12)
    as.integer(max(460, 95 + 78 * exp_n_samples(e) + leg))
  }
  tx_heat_width <- function(e, ids) {
    lab <- 20 + 8.5 * max(nchar(gene_short(ids)))
    as.integer(max(380, lab + 64 * exp_n_samples(e) + 110))
  }

  # A framed block per experiment: title, paper link and the plot
  exp_frame <- function(e, plot_ui, width = NULL, note = NULL) {
    info <- EXP_INFO[EXP_INFO$experiment == e, ]
    div(class = "exp-plot", style = if (!is.null(width)) sprintf("width: %dpx;", width + 30),
        div(class = "exp-plot-head",
            tags$h5(info$label),
            span(class = "ref",
                 if (nzchar(info$description %||% "")) paste0(info$description, ". "),
                 "Reanalysed from ",
                 if (nzchar(info$url %||% "")) tags$a(href = info$url, target = "_blank", rel = "noopener", info$reference,
                                                     .noWS = "after")
                 else info$reference, ".")),
        if (!is.null(note)) div(class = "exp-plot-note", note),
        plot_ui)
  }
  exp_index <- function(e) match(e, EXP_INFO$experiment)

  rep_note <- function(e) {
    if (exp_has_reps(e)) "Points: mean TPM; error bars: standard deviation of biological replicates."
    else "Points: averaged TPM from biological replicates (replicate values not available, so no error bars)."
  }

  output$tx_profiles_ui <- renderUI({
    q <- tx_query()
    div(class = "exp-grid", lapply(q$exps, function(e) {
      w <- tx_prof_width(e, q$ids)
      exp_frame(e, plotlyOutput(paste0("tx_prof_", exp_index(e)), width = paste0(w, "px"), height = "400px"),
                width = w, note = rep_note(e))
    }))
  })

  output$tx_heat_ui <- renderUI({
    q <- tx_query()
    h <- paste0(max(240, 130 + 26 * length(q$ids)), "px")
    div(class = "exp-grid", lapply(q$exps, function(e) {
      w <- tx_heat_width(e, q$ids)
      exp_frame(e, plotlyOutput(paste0("tx_heat_", exp_index(e)), width = paste0(w, "px"), height = h), width = w)
    }))
  })

  tx_heat_order <- reactive({
    v <- tx_z()
    if (isTRUE(input$tx_cluster)) cluster_order(v) else seq_len(nrow(v))
  })

  # One output per registered experiment; only those shown in the UI render
  lapply(seq_len(nrow(EXP_INFO)), function(k) {
    e <- EXP_INFO$experiment[k]
    pid <- paste0("tx_prof_", k)
    hid <- paste0("tx_heat_", k)

    output[[pid]] <- renderPlotly({
      plot_setup(pid, fixed = 15)
      q <- tx_query(); req(e %in% q$exps)
      d <- tx_long(); d <- d[d$experiment == e, ]
      lev <- expr$conditions$sample[expr$conditions$experiment == e]
      ids <- head(q$ids, 8)
      reps <- exp_has_reps(e)
      p <- plot_ly()
      if (isTRUE(input$tx_box)) {
        p <- p |> add_trace(
          x = d$sample, y = d$TPM, type = "box", name = paste0("All genes (n = ", length(q$ids), ")"),
          fillcolor = "rgba(138, 137, 132, 0.16)", line = list(color = MUTED, width = 1.3),
          marker = list(color = MUTED, size = 4, opacity = 0.7), boxpoints = "outliers",
          hoverinfo = "y")
      }
      for (j in seq_along(ids)) {
        dj <- d[d$gene == ids[j], ]
        dj <- dj[match(lev, dj$sample), ]
        p <- p |> add_trace(
          x = dj$sample, y = dj$TPM, type = "scatter", mode = "lines+markers", name = gene_short(ids[j]),
          line = list(color = SERIES[j], width = 2.5),
          marker = list(color = SERIES[j], size = 9, line = list(color = SURFACE, width = 1.5)),
          error_y = if (reps) list(type = "data", array = dj$SD, color = SERIES[j], thickness = 1.3, width = 4),
          customdata = if (reps) round(dj$SD, 2) else rep("", nrow(dj)),
          hovertemplate = paste0("<b>", gene_short(ids[j]), "</b><br>%{x}<br>TPM: %{y:.2f}",
                                 if (reps) "<br>SD: %{customdata}", "<extra></extra>"))
      }
      p |> finish_plot(paste0("expression_profiles_", make.names(EXP_SHORT[[e]])),
                       xaxis = ax(NULL, type = "category", categoryorder = "array", categoryarray = lev,
                                  showgrid = FALSE),
                       yaxis = ax(if (reps) "TPM (mean ± SD)" else "TPM", rangemode = "tozero"),
                       legend = list(orientation = "v", x = 1.02, y = 1, yanchor = "top", xanchor = "left",
                                     font = list(size = fsz(-1), color = INK2)),
                       margin = list(t = 10, r = 10, b = 50, l = 60))
    })

    output[[hid]] <- renderPlotly({
      plot_setup(hid, fixed = 15)
      q <- tx_query(); req(e %in% q$exps)
      k_col <- q$cond$experiment == e
      ord <- tx_heat_order()
      v <- tx_z()[ord, k_col, drop = FALSE]
      tpm <- q$m[ord, k_col, drop = FALSE]
      smp <- q$cond$sample[k_col]
      ylab <- gene_short(rownames(v))
      hover <- matrix(sprintf("<b>%s</b><br>%s<br>TPM: %s<br>Row z-score: %s", rep(ylab, ncol(v)),
                              rep(smp, each = nrow(v)), signif(as.vector(tpm), 4), round(as.vector(v), 2)),
                      nrow = nrow(v))
      lim <- suppressWarnings(max(abs(v), na.rm = TRUE))
      if (!is.finite(lim) || lim == 0) lim <- 1
      show_vals <- isTRUE(input$tx_heat_values)
      plot_ly(x = smp, y = ylab, z = v, type = "heatmap", customdata = hover,
              hovertemplate = "%{customdata}<extra></extra>",
              text = if (show_vals) matrix(fmt_tpm(tpm), nrow = nrow(v)),
              texttemplate = if (show_vals) "%{text}",
              textfont = list(size = fsz(-2), family = FONT),
              colorscale = colorscale(BWR_RAMP), zmin = -lim, zmax = lim, xgap = 2, ygap = 2,
              colorbar = list(title = list(text = "Row<br>z-score", side = "top",
                                           font = list(size = fsz(-1), color = INK2)),
                              tickfont = list(size = fsz(-2), color = INK2), thickness = 14, len = 0.9)) |>
        finish_plot(paste0("expression_heatmap_", make.names(EXP_SHORT[[e]])), legend = FALSE,
                    xaxis = ax(NULL, type = "category", categoryorder = "array", categoryarray = smp,
                               showgrid = FALSE, showline = FALSE),
                    yaxis = ax(NULL, autorange = "reversed", showgrid = FALSE, showline = FALSE,
                               tickfont = list(size = fsz(-1), color = INK)),
                    margin = list(t = 10, r = 10, b = 50, l = 10))
    })
  })

  # Gene-gene correlation, computed separately within each experiment (correlating
  # across different experiments would mix study designs). Experiments with fewer
  # than three conditions are skipped: with two points r is always +1 or -1.
  CORR_MIN_SAMPLES <- 3
  tx_corr_data <- reactive({
    q <- tx_query()
    validate(need(length(q$ids) >= 2, "Correlation needs at least two genes."))
    l <- log2(q$m + 1)
    usable <- Filter(function(e) sum(q$cond$experiment == e) >= CORR_MIN_SAMPLES, q$exps)
    cms <- lapply(setNames(usable, usable), function(e) {
      suppressWarnings(stats::cor(t(l[, q$cond$experiment == e, drop = FALSE]),
                                  use = "pairwise.complete.obs"))
    })
    # One gene order for every experiment, from the mean correlation, so rows line up
    n <- length(q$ids)
    ord <- seq_len(n)
    if (isTRUE(input$tx_cluster) && n > 2 && length(cms)) {
      avg <- Reduce(`+`, lapply(cms, function(m) { m[is.na(m)] <- 0; m })) / length(cms)
      ord <- stats::hclust(stats::as.dist(1 - avg), "average")$order
    }
    list(cms = lapply(cms, function(m) m[ord, ord, drop = FALSE]),
         skipped = setdiff(q$exps, usable))
  })

  # Colour range shared by all experiments: from the lowest r among the selected
  # genes (rounded down to 0.05) up to 1
  corr_range <- function(cms) {
    off <- unlist(lapply(cms, function(cm) cm[upper.tri(cm)]))
    lo <- suppressWarnings(min(off, na.rm = TRUE))
    if (!is.finite(lo)) return(c(-1, 1))
    lo <- max(-1, floor(lo * 20) / 20)
    if (lo > 0.95) lo <- 0.95
    c(lo, 1)
  }
  tx_corr_size <- function(ids) {
    cells <- min(820, max(300, 34 * length(ids)))
    lab <- 20 + 8.5 * max(nchar(gene_short(ids)))
    list(w = as.integer(cells + lab + 120), h = as.integer(cells + lab + 30))
  }

  output$tx_corr_ui <- renderUI({
    q <- tx_query()
    cd <- tx_corr_data()
    sz <- tx_corr_size(q$ids)
    tagList(
      div(class = "exp-grid", lapply(names(cd$cms), function(e) {
        exp_frame(e, plotlyOutput(paste0("tx_corr_", exp_index(e)), width = paste0(sz$w, "px"),
                                  height = paste0(sz$h, "px")),
                  width = sz$w, note = sprintf("Pearson r across the %d conditions of this experiment.",
                                               sum(q$cond$experiment == e)))
      })),
      if (length(cd$skipped))
        div(class = "section-note mt-2",
            "Not shown (fewer than ", CORR_MIN_SAMPLES, " conditions, so a correlation is not meaningful): ",
            paste(EXP_SHORT[cd$skipped], collapse = ", "), "."))
  })
  lapply(seq_len(nrow(EXP_INFO)), function(k) {
    e <- EXP_INFO$experiment[k]
    cid <- paste0("tx_corr_", k)
    output[[cid]] <- renderPlotly({
      plot_setup(cid, fixed = 15)
      cd <- tx_corr_data(); req(e %in% names(cd$cms))
      cm <- cd$cms[[e]]; lab <- gene_short(rownames(cm))
      rng <- corr_range(cd$cms)
      hover <- matrix(sprintf("%s<br>%s<br>Pearson r = %s", rep(lab, ncol(cm)), rep(lab, each = nrow(cm)),
                              round(as.vector(cm), 3)), nrow = nrow(cm))
      plot_ly(x = lab, y = lab, z = cm, type = "heatmap", zmin = rng[1], zmax = rng[2], text = hover,
              hoverinfo = "text", colorscale = colorscale(VIRIDIS), xgap = 2, ygap = 2,
              colorbar = list(title = list(text = "Pearson r", side = "right",
                                           font = list(size = fsz(-1), color = INK2)),
                              tickfont = list(size = fsz(-2), color = INK2), thickness = 14, len = 0.6)) |>
        finish_plot(paste0("expression_correlation_", make.names(EXP_SHORT[[e]])), legend = FALSE,
                    xaxis = ax(NULL, tickangle = -50, showgrid = FALSE, showline = FALSE,
                               tickfont = list(size = fsz(-1), color = INK)),
                    yaxis = ax(NULL, autorange = "reversed", showgrid = FALSE, showline = FALSE,
                               scaleanchor = "x", tickfont = list(size = fsz(-1), color = INK)),
                    margin = list(t = 10, r = 10, b = 10, l = 10))
    })
  })

  output$tx_table <- renderReactable({
    d <- tx_long()
    has_sd <- any(!is.na(d$SD))
    reactable(d, searchable = TRUE, compact = TRUE, highlight = TRUE, defaultPageSize = 20,
              filterable = TRUE,
              columns = list(gene = colDef(name = "Gene", cell = function(v) tags$span(class = "gene-id", v)),
                             name = colDef(name = "Name"),
                             experiment = colDef(name = "Experiment", minWidth = 180),
                             sample = colDef(name = "Sample"),
                             TPM = colDef(name = if (has_sd) "TPM (mean)" else "TPM (averaged)", format = colFormat(digits = 2)),
                             SD = colDef(show = has_sd, format = colFormat(digits = 2)),
                             log2_TPM_plus_1 = colDef(name = "log₂(TPM + 1)", format = colFormat(digits = 2)),
                             row_z_score = colDef(name = "Row z-score", format = colFormat(digits = 2))))
  })

  output$tx_dl_long <- downloadHandler(
    filename = function() paste0("MagnaGenes_expression_long_", format(Sys.Date(), "%Y%m%d"), ".tsv"),
    content = function(file) write_tsv(tx_long(), file))
  output$tx_dl_wide <- downloadHandler(
    filename = function() paste0("MagnaGenes_expression_matrix_TPM_", format(Sys.Date(), "%Y%m%d"), ".tsv"),
    content = function(file) {
      q <- tx_query(); v <- q$m
      colnames(v) <- q$cond$id
      write_tsv(data.frame(gene = rownames(v), name = first_name[rownames(v)], v,
                           check.names = FALSE), file)
    })
  output$tx_dl_fasta <- downloadHandler(
    filename = function() paste0("MagnaGenes_expression_genes_protein_", format(Sys.Date(), "%Y%m%d"), ".fasta"),
    content = function(file) writeLines(fasta_text(tx_query()$ids, "protein"), file))

  tx_gg <- function() {
    view <- input$tx_view %||% "profiles"
    q <- tx_query()
    if (view %in% c("profiles", "table")) {
      d_all <- tx_long()
      d_all$experiment <- factor(EXP_PARSE[d_all$experiment], EXP_PARSE[q$exps])
      d_all$sample <- factor(d_all$sample, expr$sample_levels)
      ids <- head(unique(d_all$gene), 8)
      d <- d_all[d_all$gene %in% ids, ]
      d$label <- factor(gene_short(d$gene), gene_short(ids))
      g <- ggplot(d, aes(sample, TPM, colour = label, group = label))
      if (isTRUE(input$tx_box))
        g <- g + geom_boxplot(data = d_all, aes(sample, TPM, group = sample), inherit.aes = FALSE,
                              fill = "grey90", colour = "grey55", width = 0.6, outlier.size = 0.6)
      if (any(!is.na(d$SD)))
        g <- g + geom_errorbar(aes(ymin = pmax(0, TPM - SD), ymax = TPM + SD), width = 0.2, linewidth = 0.4,
                               na.rm = TRUE)
      g + geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
        facet_wrap(~experiment, scales = "free", labeller = label_parsed) +
        scale_colour_manual(values = SERIES) +
        scale_y_continuous(limits = c(0, NA)) +
        labs(x = NULL, y = if (any(!is.na(d$SD))) "TPM (mean +/- SD)" else "TPM (averaged from biological replicates)",
             colour = NULL) +
        theme_mg() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                           panel.border = element_rect(colour = "black", linewidth = 0.6, fill = NA))
    } else if (view == "heatmap") {
      # One heatmap per experiment, each on its own z-score range, stacked
      ord <- tx_heat_order()
      v_all <- tx_z()[ord, , drop = FALSE]
      tpm_all <- q$m[ord, , drop = FALSE]
      show_vals <- isTRUE(input$tx_heat_values)
      plots <- lapply(q$exps, function(e) {
        k <- q$cond$experiment == e
        v <- v_all[, k, drop = FALSE]
        d <- data.frame(gene = factor(rep(gene_short(rownames(v)), ncol(v)), rev(gene_short(rownames(v)))),
                        sample = factor(rep(q$cond$sample[k], each = nrow(v)), q$cond$sample[k]),
                        value = as.vector(v), tpm = fmt_tpm(as.vector(tpm_all[, k])))
        lim <- suppressWarnings(max(abs(d$value), na.rm = TRUE))
        if (!is.finite(lim) || lim == 0) lim <- 1
        ggplot(d, aes(sample, gene, fill = value)) + geom_tile(colour = "white", linewidth = 0.4) +
          (if (show_vals) geom_text(aes(label = tpm), size = 2.2, colour = "#1d2b2a")) +
          scale_fill_gradientn(colours = BWR_RAMP, limits = c(-lim, lim), name = "Row z-score",
                               na.value = "#f2f2ef") +
          labs(x = NULL, y = NULL, title = parse(text = EXP_PARSE[[e]])[[1]]) + theme_mg(9) +
          theme(panel.grid = element_blank(), legend.position = "right",
                panel.border = element_rect(colour = "black", linewidth = 0.6, fill = NA))
      })
      function() {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(length(plots), 1)))
        for (i in seq_along(plots)) print(plots[[i]], vp = grid::viewport(layout.pos.row = i, layout.pos.col = 1))
      }
    } else {
      cd <- tx_corr_data()
      validate(need(length(cd$cms) > 0, "No selected experiment has enough conditions for a correlation."))
      d <- do.call(rbind, lapply(names(cd$cms), function(e) {
        cm <- cd$cms[[e]]; lab <- gene_short(rownames(cm))
        data.frame(experiment = e, a = factor(rep(lab, ncol(cm)), lab),
                   b = factor(rep(lab, each = nrow(cm)), rev(lab)), r = as.vector(cm))
      }))
      d$experiment <- factor(EXP_PARSE[d$experiment], EXP_PARSE[names(cd$cms)])
      ggplot(d, aes(a, b, fill = r)) + geom_tile(colour = "white", linewidth = 0.3) +
        facet_wrap(~experiment, labeller = label_parsed) +
        scale_fill_gradientn(colours = VIRIDIS, limits = corr_range(cd$cms), oob = scales::squish,
                             name = "Pearson r", na.value = "#f2f2ef") +
        labs(x = NULL, y = NULL,
             caption = "Pearson correlation of log2(TPM + 1) across the conditions of each experiment") +
        coord_equal() + theme_mg() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
              legend.position = "right")
    }
  }
  output$tx_dl_pdf <- pdf_download(tx_gg, function() paste0("MagnaGenes_expression_", input$tx_view %||% "profiles"),
                                   reactive(input$tx_pdf_w %||% 11), reactive(input$tx_pdf_h %||% 7))

  # ===========================================================================
  # PHOSPHOPROTEOMICS
  # ===========================================================================
  updateSelectizeInput(session, "ph_protein", choices = setNames(phos_proteins, gene_short(phos_proteins)),
                       selected = character(0), server = TRUE)
  observeEvent(input$ph_example, {
    updateSelectizeInput(session, "ph_protein", choices = setNames(phos_proteins, gene_short(phos_proteins)),
                         selected = EX_PHOSPHO, server = TRUE)
  })

  ph_protein <- reactive({
    req(input$ph_protein, input$ph_protein %in% phos_proteins)
    input$ph_protein
  })

  ph_meta <- reactive(phos$meta[phos$meta$protein == ph_protein(), ])

  output$ph_protein_note <- renderText({
    if (!isTruthy(input$ph_protein)) return(paste(format(length(phos_proteins), big.mark = ","), "proteins with phosphopeptides."))
    paste0(format(length(phos_proteins), big.mark = ","), " proteins with phosphopeptides. ",
           phos_n[[ph_protein()]], " peptide forms for this protein.")
  })

  output$ph_peptides_ui <- renderUI({
    m <- ph_meta()
    sel <- head(m$pid, 6)
    tagList(
      tags$label(class = "form-label", "Phosphopeptides (site, peptide)"),
      div(class = "peptide-box",
          checkboxGroupInput("ph_pids", NULL, choices = setNames(m$pid, m$site_label), selected = sel)),
      div(class = "d-flex gap-3 mt-1 mb-2",
          actionLink("ph_all", "Select all", class = "example-link"),
          actionLink("ph_none", "Clear", class = "example-link"))
    )
  })
  observeEvent(input$ph_all, updateCheckboxGroupInput(session, "ph_pids", selected = ph_meta()$pid))
  observeEvent(input$ph_none, updateCheckboxGroupInput(session, "ph_pids", selected = character(0)))

  ph_value <- function(a) if (input$ph_scale == "log") log2(a + 1) else a

  ph_long_all <- reactive({
    m <- ph_meta()
    A <- phos$abundance[as.character(m$pid), , drop = FALSE]
    if (isTRUE(input$ph_zero_na)) A[A == 0] <- NA
    cols <- phos$columns
    data.frame(pid = rep(m$pid, times = ncol(A)),
               site = rep(m$site_label, times = ncol(A)),
               peptide = rep(m$peptide, times = ncol(A)),
               modifications = rep(m$modifications, times = ncol(A)),
               sample = rep(cols$sample, each = nrow(A)), hpi = rep(cols$hpi, each = nrow(A)),
               replicate = rep(cols$replicate, each = nrow(A)), abundance = as.vector(A),
               stringsAsFactors = FALSE)
  })

  ph_summary <- function(d) {
    d$value <- ph_value(d$abundance)
    key <- interaction(d$pid, d$sample, d$hpi, drop = TRUE, lex.order = TRUE)
    s <- do.call(rbind, lapply(split(d, key), function(x) {
      v <- x$value[!is.na(x$value)]
      data.frame(pid = x$pid[1], site = x$site[1], sample = x$sample[1], hpi = x$hpi[1],
                 mean = if (length(v)) mean(v) else NA, sd = if (length(v) > 1) stats::sd(v) else NA,
                 n = length(v))
    }))
    rownames(s) <- NULL
    s
  }

  ph_stats <- reactive({
    d <- ph_long_all()
    d$l2 <- log2(d$abundance + 1)
    res <- do.call(rbind, lapply(split(d, interaction(d$pid, d$hpi, drop = TRUE, lex.order = TRUE)), function(x) {
      g <- x$l2[x$sample == SAMPLES[1] & !is.na(x$l2)]
      k <- x$l2[x$sample == SAMPLES[2] & !is.na(x$l2)]
      p <- if (length(g) >= 2 && length(k) >= 2 && (stats::var(g) > 0 || stats::var(k) > 0))
        tryCatch(stats::t.test(k, g)$p.value, error = function(e) NA) else NA
      data.frame(pid = x$pid[1], site = x$site[1], peptide = x$peptide[1], hpi = x$hpi[1],
                 n_guy11 = length(g), n_pmk1 = length(k),
                 mean_log2_guy11 = if (length(g)) mean(g) else NA,
                 mean_log2_pmk1 = if (length(k)) mean(k) else NA,
                 log2FC = if (length(g) && length(k)) mean(k) - mean(g) else NA,
                 p_value = p, stringsAsFactors = FALSE)
    }))
    res$p_adj <- stats::p.adjust(res$p_value, "BH")
    rownames(res) <- NULL
    res[order(res$pid, res$hpi), ]
  })

  ph_selected <- reactive({
    pids <- as.integer(input$ph_pids %||% integer(0))
    pids <- intersect(pids, ph_meta()$pid)
    validate(need(length(pids) > 0, "Select at least one phosphopeptide."))
    pids
  })

  output$ph_header <- renderUI({
    validate(need(isTruthy(input$ph_protein), "Choose a protein in the box on the left (or load the example)."))
    g <- ph_protein(); m <- ph_meta()
    nm <- first_name[[g]] %||% NA
    st <- ph_stats()
    n_sig <- sum(st$p_adj < 0.05, na.rm = TRUE)
    div(class = "summary-strip",
      span(class = "gene-id", style = "font-size: 1.05rem;", g),
      if (!is.na(nm)) span(style = "font-weight: 700; color: var(--brand-dark);", nm),
      span(class = "section-note", gene_desc(g)),
      span(class = "chip", tags$b(nrow(m)), "peptide forms"),
      span(class = "chip", tags$b(length(unique(unlist(regmatches(m$sites, gregexpr("[STY][0-9]+", m$sites)))))), "localised sites"),
      span(class = if (n_sig) "chip ok" else "chip", tags$b(n_sig), "tests with adj. p < 0.05"),
      if (g %in% genes$gene) HTML(xref_links_html(g)))
  })

  ph_x_axis <- function(title = "Hours post inoculation") ax(title, tickvals = HPI, ticktext = HPI, showgrid = FALSE)

  # One framed plot per phosphopeptide: site and full peptide in the header,
  # Guy11 (greens) above delta-pmk1 (greys), one bar per biological replicate
  ph_bar_plot <- function(pid, file) {
    d <- ph_long_all(); d <- d[d$pid == pid, ]
    d$value <- ph_value(d$abundance)
    reps <- sort(unique(d$replicate))
    ylab <- if (input$ph_scale == "log") "log₂ abundance" else "Abundance"
    hlab <- paste(HPI, "hpi")
    ymax <- suppressWarnings(max(d$value, na.rm = TRUE))
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    panels <- lapply(SAMPLES, function(smp) {
      p <- plot_ly()
      for (r in seq_along(reps)) {
        dd <- d[d$sample == smp & d$replicate == reps[r], ]
        dd <- dd[match(HPI, dd$hpi), ]
        rl <- sub("rep_", "replicate ", reps[r])
        p <- p |> add_trace(x = hlab, y = dd$value, type = "bar", name = paste(smp, rl),
                            marker = list(color = COL_REP[[smp]][r], line = list(color = SURFACE, width = 1)),
                            hovertemplate = paste0("<b>", smp, "</b>, %{x}<br>", rl, "<br>", ylab,
                                                   ": %{y:.3s}<extra></extra>"))
      }
      p |> layout(yaxis = ax(ylab, range = c(0, ymax * 1.1)))
    })
    ann <- lapply(1:2, function(j) list(
      text = paste0("<b>", SAMPLES[j], "</b>"), xref = "paper", yref = paste0("y", if (j == 1) "" else 2, " domain"),
      x = 0, y = 1, xanchor = "left", yanchor = "bottom", showarrow = FALSE,
      font = list(size = fsz(0), color = COL_SAMPLE[[j]])))
    subplot(panels, nrows = 2, shareX = TRUE, titleY = TRUE, margin = 0.09) |>
      finish_plot(file, annotations = ann, barmode = "group", bargap = 0.25,
                  xaxis = ax("Hours post inoculation", type = "category", categoryorder = "array",
                             categoryarray = hlab, showgrid = FALSE),
                  legend = list(orientation = "h", x = 0, y = -0.16, yanchor = "top",
                                font = list(size = fsz(-2), color = INK2)),
                  margin = list(t = 30, r = 10, b = 90, l = 60))
  }

  output$ph_profiles_ui <- renderUI({
    pids <- ph_selected(); m <- ph_meta()
    for (pid in pids) local({
      id <- pid; oid <- paste0("ph_prof_", id)
      output[[oid]] <- renderPlotly({
        plot_setup(oid, fixed = 15)
        ph_bar_plot(id, paste0(ph_protein(), "_phosphosite_", id))
      })
    })
    div(class = "exp-grid", lapply(pids, function(id) {
      r <- m[m$pid == id, ]
      div(class = "exp-plot", style = "width: 640px;",
          div(class = "exp-plot-head",
              tags$h5(if (is.na(r$sites)) "Site not localised" else paste("Site", gsub(";", ", ", r$sites))),
              span(class = "ref", r$position)),
          div(class = "peptide-seq", r$peptide),
          div(class = "exp-plot-note", r$modifications),
          plotlyOutput(paste0("ph_prof_", id), width = "610px", height = "470px"))
    }))
  })

  ph_heat_matrix <- reactive({
    d <- ph_long_all()
    s <- ph_summary(d)
    m <- ph_meta()
    cols <- expand.grid(hpi = HPI, sample = SAMPLES, stringsAsFactors = FALSE)
    M <- matrix(NA_real_, nrow(m), nrow(cols), dimnames = list(m$site_label, paste(cols$sample, cols$hpi)))
    idx <- cbind(match(s$pid, m$pid), match(paste(s$sample, s$hpi), colnames(M)))
    M[idx] <- s$mean
    if (isTRUE(input$ph_rowz)) { M <- t(scale(t(M))); M[!is.finite(M)] <- NA }
    list(M = M, cols = cols)
  })

  output$ph_heat_ui <- renderUI({
    plotlyOutput("ph_heat", height = paste0(max(420, 170 + 36 * nrow(ph_meta())), "px"))
  })
  # Two heatmaps on one value range: greens for Guy11, greys for delta-pmk1
  output$ph_heat <- renderPlotly({
    plot_setup("ph_heat")
    h <- ph_heat_matrix(); M <- h$M; cols <- h$cols
    z_mode <- isTRUE(input$ph_rowz)
    lab <- if (z_mode) "Row z-score" else if (input$ph_scale == "log") "log₂ abundance" else "Abundance"
    rng <- range(M, na.rm = TRUE)
    if (!all(is.finite(rng))) rng <- c(0, 1)
    panels <- lapply(seq_along(SAMPLES), function(j) {
      k <- cols$sample == SAMPLES[j]
      Mj <- M[, k, drop = FALSE]
      hover <- matrix(sprintf("<b>%s</b><br>%s, %s hpi<br>%s: %s", rep(rownames(Mj), ncol(Mj)), SAMPLES[j],
                              rep(cols$hpi[k], each = nrow(Mj)), lab, signif(as.vector(Mj), 4)), nrow = nrow(Mj))
      plot_ly(x = paste(cols$hpi[k], "hpi"), y = rownames(Mj), z = Mj, type = "heatmap", text = hover,
              hoverinfo = "text", colorscale = colorscale(RAMP_SAMPLE[[SAMPLES[j]]]), zmin = rng[1], zmax = rng[2],
              xgap = 2, ygap = 2,
              colorbar = list(title = list(text = paste0(SAMPLES[j], "<br>", lab), side = "top",
                                           font = list(size = fsz(-2), color = INK2)),
                              tickfont = list(size = fsz(-2), color = INK2),
                              thickness = 12, len = 0.42, x = 1.02, xanchor = "left",
                              y = if (j == 1) 1 else 0.45, yanchor = "top")) |>
        layout(xaxis = ax("Hours post inoculation", type = "category", showgrid = FALSE, showline = FALSE),
               yaxis = ax(NULL, autorange = "reversed", showgrid = FALSE, showline = FALSE,
                          showticklabels = (j == 1),
                          tickfont = list(size = fsz(-2), color = INK, family = "ui-monospace, Menlo, monospace")))
    })
    subplot(panels, nrows = 1, shareY = TRUE, titleX = TRUE, margin = 0.015) |>
      finish_plot(paste0(ph_protein(), "_phospho_heatmap"), legend = FALSE,
                  annotations = lapply(1:2, function(j) list(
                    text = paste0("<b>", SAMPLES[j], "</b>"), xref = paste0("x", if (j == 1) "" else 2, " domain"),
                    yref = "paper", x = 0.5, y = 1.01, yanchor = "bottom", showarrow = FALSE,
                    font = list(size = fsz(1), color = COL_SAMPLE[[j]]))),
                  margin = list(t = 40, r = 10, b = 50, l = 10))
  })

  ph_fc_matrix <- reactive({
    st <- ph_stats(); m <- ph_meta()
    FC <- matrix(NA_real_, nrow(m), length(HPI), dimnames = list(m$site_label, HPI))
    P <- FC
    idx <- cbind(match(st$pid, m$pid), match(st$hpi, HPI))
    FC[idx] <- st$log2FC; P[idx] <- st$p_adj
    list(FC = FC, P = P, raw_p = { R <- FC; R[idx] <- st$p_value; R })
  })

  output$ph_fc_ui <- renderUI({
    plotlyOutput("ph_fc", height = paste0(max(320, 150 + 20 * nrow(ph_meta())), "px"))
  })
  output$ph_fc <- renderPlotly({
    plot_setup("ph_fc")
    f <- ph_fc_matrix(); FC <- f$FC
    lim <- max(1, max(abs(FC), na.rm = TRUE))
    hover <- matrix(sprintf("<b>%s</b><br>%s hpi<br>log₂FC (Δpmk1 / Guy11): %s<br>p = %s, adj. p = %s",
                            rep(rownames(FC), ncol(FC)), rep(HPI, each = nrow(FC)), round(as.vector(FC), 2),
                            signif(as.vector(f$raw_p), 2), signif(as.vector(f$P), 2)), nrow = nrow(FC))
    sig <- which(!is.na(f$P) & f$P < 0.05, arr.ind = TRUE)
    p <- plot_ly() |>
      add_trace(x = seq_along(HPI), y = rownames(FC), z = FC, type = "heatmap", text = hover, hoverinfo = "text",
                colorscale = colorscale(DIV_RAMP), zmin = -lim, zmax = lim, xgap = 2, ygap = 2,
                colorbar = list(title = list(text = "log₂FC", side = "right"), thickness = 12, len = 0.6))
    if (nrow(sig)) {
      p <- p |> add_trace(x = sig[, 2], y = rownames(FC)[sig[, 1]], type = "scatter", mode = "markers", inherit = FALSE,
                          marker = list(symbol = "circle", size = 7, color = INK, line = list(color = SURFACE, width = 1)),
                          name = "adj. p < 0.05", hoverinfo = "skip", showlegend = TRUE)
    }
    p |> finish_plot(paste0(ph_protein(), "_phospho_log2FC"), legend = nrow(sig) > 0,
                     xaxis = ax("Hours post inoculation", tickvals = seq_along(HPI), ticktext = HPI,
                                showgrid = FALSE, showline = FALSE),
                     yaxis = ax(NULL, autorange = "reversed", showgrid = FALSE, showline = FALSE, type = "category",
                                tickfont = list(size = fsz(-3), color = INK, family = "ui-monospace, Menlo, monospace")),
                     annotations = list(list(text = "Blue: lower in Δpmk1   Red: higher in Δpmk1   Dot: adj. p < 0.05",
                                             xref = "paper", yref = "paper", x = 0, y = 1.01, yanchor = "bottom",
                                             xanchor = "left", showarrow = FALSE, font = list(size = fsz(-2), color = INK2))),
                     margin = list(t = 50, r = 10, b = 50, l = 10))
  })

  ph_map_data <- reactive({
    st <- ph_stats(); m <- ph_meta()
    tsel <- input$ph_map_time %||% "max"
    pick <- do.call(rbind, lapply(split(st, st$pid), function(x) {
      x <- x[!is.na(x$log2FC), ]
      if (!nrow(x)) return(NULL)
      r <- if (tsel == "max") x[which.max(abs(x$log2FC)), ] else x[x$hpi == as.numeric(tsel), ]
      if (!nrow(r)) NULL else r
    }))
    if (is.null(pick)) return(NULL)
    pos <- site_positions(m$sites[match(pick$pid, m$pid)])
    out <- do.call(rbind, lapply(seq_len(nrow(pick)), function(i)
      if (nrow(pos[[i]])) cbind(pick[rep(i, nrow(pos[[i]])), ], pos[[i]]) else NULL))
    out
  })

  # Site map: lollipops of log2FC (delta-pmk1 / Guy11) per localised site, labelled
  # with the residue, above a protein track with the Pfam domains
  output$ph_map <- renderPlotly({
    plot_setup("ph_map")
    g <- ph_protein()
    len <- genes[g, "length_aa"]
    d <- ph_map_data()
    if (is.null(d) || !nrow(d)) return(empty_plot("No localised phosphosites with fold-change estimates for this protein"))
    if (is.na(len)) len <- max(d$pos) + 20
    fs <- fsz(0)
    d$class <- ifelse(!is.na(d$p_adj) & d$p_adj < 0.05, ifelse(d$log2FC > 0, "up", "down"), "ns")
    d$lab <- paste0("p", d$residue, d$pos)
    ccol <- c(up = DIV_RAMP[6], down = DIV_RAMP[2], ns = "#a9a8a2")
    clab <- c(up = "Higher in Δpmk1 (adj. p < 0.05)", down = "Lower in Δpmk1 (adj. p < 0.05)",
              ns = "Not significant")
    ylim <- max(1.5, max(abs(d$log2FC), na.rm = TRUE) * 1.3)

    # Top panel: lollipops with residue labels
    p_top <- plot_ly() |>
      add_segments(x = d$pos, xend = d$pos, y = 0, yend = d$log2FC, showlegend = FALSE, hoverinfo = "skip",
                   line = list(color = "#b5b4ad", width = 1.4))
    for (cl in c("up", "down", "ns")) {
      dd <- d[d$class == cl, ]
      if (!nrow(dd)) next
      p_top <- p_top |> add_trace(
        x = dd$pos, y = dd$log2FC, type = "scatter", mode = "markers+text", name = clab[[cl]],
        text = dd$lab, textposition = ifelse(dd$log2FC >= 0, "top center", "bottom center"),
        textfont = list(size = fs - 1, color = INK),
        marker = list(color = ccol[[cl]], size = 12, line = list(color = SURFACE, width = 1.5)),
        customdata = sprintf("<b>%s</b><br>%s<br>%s hpi: log₂FC %s<br>adj. p = %s",
                             dd$lab, dd$peptide, dd$hpi, round(dd$log2FC, 2), signif(dd$p_adj, 2)),
        hovertemplate = "%{customdata}<extra></extra>")
    }
    p_top <- p_top |> layout(yaxis = ax("log₂FC<br>(Δpmk1 / Guy11)", range = c(-ylim, ylim),
                                        zeroline = TRUE, zerolinecolor = AXIS, zerolinewidth = 1.5))

    # Bottom panel: protein track with domains and site ticks
    pf <- parse_pfam(genes[g, "pfam"])
    p_bot <- plot_ly() |>
      add_trace(x = c(1, len), y = c(0, 0), type = "scatter", mode = "lines", showlegend = FALSE,
                hoverinfo = "skip", line = list(color = BACKBONE, width = 10))
    shapes <- list()
    ann <- list()
    if (!is.null(pf)) {
      dom_names <- unique(pf$name)
      dcol <- setNames(SERIES[(seq_along(dom_names) - 1) %% length(SERIES) + 1], dom_names)
      plot_px <- max(300, (session$clientData$output_ph_map_width %||% 900) - 140)
      px_per_aa <- plot_px / len
      for (i in seq_len(nrow(pf))) {
        shapes[[length(shapes) + 1]] <- list(type = "rect", xref = "x", yref = "y",
                                             x0 = pf$start[i], x1 = pf$end[i], y0 = -0.45, y1 = 0.45,
                                             fillcolor = dcol[[pf$name[i]]], line = list(color = SURFACE, width = 1))
        fit <- floor((pf$end[i] - pf$start[i]) * px_per_aa / (fs * 0.6)) - 1
        nm <- pf$name[i]
        inside <- fit >= 6
        lab <- if (!inside) nm else if (nchar(nm) > fit) paste0(substr(nm, 1, fit - 1), "…") else nm
        ann[[length(ann) + 1]] <- list(text = paste0("<b>", lab, "</b>"), xref = "x", yref = "y2",
                                       x = (pf$start[i] + pf$end[i]) / 2, y = if (inside) 0 else -0.6,
                                       yanchor = if (inside) "middle" else "top", showarrow = FALSE,
                                       font = list(size = fs - 2, color = if (inside) "#ffffff" else INK2))
      }
      p_bot <- p_bot |> add_trace(x = (pf$start + pf$end) / 2, y = rep(0, nrow(pf)), type = "scatter",
                                  mode = "markers", marker = list(size = 20, opacity = 0), showlegend = FALSE,
                                  customdata = sprintf("<b>%s</b> (%s)<br>%d-%d aa", pf$name, pf$pfam, pf$start, pf$end),
                                  hovertemplate = "%{customdata}<extra></extra>")
    } else {
      ann[[1]] <- list(text = "No Pfam domains annotated", xref = "x", yref = "y2", x = len / 2, y = -0.6,
                       yanchor = "top", showarrow = FALSE, font = list(size = fs - 2, color = MUTED))
    }
    u <- d[!duplicated(d$pos), ]
    p_bot <- p_bot |>
      add_trace(x = u$pos, y = rep(0.62, nrow(u)), type = "scatter", mode = "markers", showlegend = FALSE,
                marker = list(symbol = "triangle-down", size = 9, color = PHOS_MARK),
                customdata = u$lab, hovertemplate = "%{customdata}<extra></extra>") |>
      layout(yaxis = list(visible = FALSE, range = c(-1.2, 0.9), fixedrange = TRUE), shapes = shapes)

    # Shapes belong to the bottom panel (subplot maps them onto its axes)
    subplot(p_top, p_bot, nrows = 2, heights = c(0.74, 0.26), shareX = TRUE, titleY = TRUE, margin = 0.02) |>
      finish_plot(paste0(g, "_phosphosite_map"), annotations = ann,
                  xaxis = ax(paste0("Residue (", g, ", ", len, " aa)"), range = c(-len * 0.02, len * 1.02),
                             showgrid = FALSE),
                  legend = list(orientation = "h", x = 0, y = 1.02, yanchor = "bottom",
                                font = list(size = fsz(-1), color = INK2)),
                  margin = list(t = 40, r = 20, b = 50, l = 70))
  })

  # Domain annotation listed under the site map
  output$ph_map_domains <- renderUI({
    g <- ph_protein()
    pf <- parse_pfam(genes[g, "pfam"])
    ipr <- genes[g, "interpro"]
    div(class = "method-box",
        tags$strong("Domains. "),
        if (is.null(pf)) "No Pfam domains with positions are annotated for this protein."
        else tagList(lapply(seq_len(nrow(pf)), function(i) tagList(
          if (i > 1) "; ",
          tags$a(href = url_pfam(pf$pfam[i]), target = "_blank", rel = "noopener", pf$name[i]),
          sprintf(" (%s, %d-%d aa)", pf$pfam[i], pf$start[i], pf$end[i])))),
        if (!is.na(ipr) && nzchar(ipr)) div(class = "mt-1", tags$strong("InterPro: "), ipr))
  })

  output$ph_table <- renderReactable({
    st <- ph_stats()
    st$pid <- NULL
    reactable(st, compact = TRUE, highlight = TRUE, searchable = TRUE, defaultPageSize = 18, filterable = TRUE,
              defaultSorted = list(p_value = "asc"),
              columns = list(
                site = colDef(name = "Site / peptide", minWidth = 200,
                              style = list(fontFamily = "ui-monospace, Menlo, monospace", fontSize = "0.78rem")),
                peptide = colDef(name = "Peptide", minWidth = 200,
                                 style = list(fontFamily = "ui-monospace, Menlo, monospace", fontSize = "0.78rem")),
                hpi = colDef(name = "hpi", maxWidth = 60),
                n_guy11 = colDef(name = "n Guy11", maxWidth = 75), n_pmk1 = colDef(name = "n Δpmk1", maxWidth = 80),
                mean_log2_guy11 = colDef(name = "log₂ Guy11", format = colFormat(digits = 2)),
                mean_log2_pmk1 = colDef(name = "log₂ Δpmk1", format = colFormat(digits = 2)),
                log2FC = colDef(name = "log₂FC", format = colFormat(digits = 2),
                                style = function(v) if (!is.na(v) && abs(v) >= 1) list(fontWeight = 700) else NULL),
                p_value = colDef(name = "p", cell = function(v) if (is.na(v)) "" else formatC(v, format = "g", digits = 2)),
                p_adj = colDef(name = "adj. p", cell = function(v) if (is.na(v)) "" else formatC(v, format = "g", digits = 2))
              ))
  })

  output$ph_dl_data <- downloadHandler(
    filename = function() paste0("MagnaGenes_phospho_", ph_protein(), "_", format(Sys.Date(), "%Y%m%d"), ".tsv"),
    content = function(file) { d <- ph_long_all(); d$pid <- NULL; d$protein <- ph_protein(); write_tsv(d, file) })
  output$ph_dl_stats <- downloadHandler(
    filename = function() paste0("MagnaGenes_phospho_stats_", ph_protein(), "_", format(Sys.Date(), "%Y%m%d"), ".tsv"),
    content = function(file) { d <- ph_stats(); d$pid <- NULL; d$protein <- ph_protein(); write_tsv(d, file) })
  output$ph_dl_fasta <- downloadHandler(
    filename = function() paste0(ph_protein(), "_protein.fasta"),
    content = function(file) writeLines(fasta_text(ph_protein(), "protein", TRUE), file))

  ph_gg <- function() {
    view <- input$ph_view %||% "profiles"
    if (view %in% c("profiles", "table")) {
      pids <- ph_selected()
      d <- ph_long_all(); d <- d[d$pid %in% pids, ]
      d$value <- ph_value(d$abundance)
      lv <- phos$meta$site_label[match(pids, phos$meta$pid)]
      d$site <- factor(pdf_safe(d$site), pdf_safe(lv))
      d$sample <- factor(d$sample, SAMPLES, labels = c("Guy11", "delta-pmk1"))
      d$group <- paste(d$sample, d$replicate)
      reps <- sort(unique(d$replicate))
      fill_vals <- c(setNames(COL_REP[[1]], paste("Guy11", reps)), setNames(COL_REP[[2]], paste("delta-pmk1", reps)))
      ggplot(d, aes(factor(hpi, HPI), value, fill = group)) +
        geom_col(position = position_dodge(width = 0.85), width = 0.8, colour = "white", linewidth = 0.2, na.rm = TRUE) +
        facet_grid(sample ~ site) +
        scale_fill_manual(values = fill_vals, name = NULL, labels = function(x) sub("rep_", "replicate ", x)) +
        labs(x = "Hours post inoculation", y = if (input$ph_scale == "log") "log2 abundance" else "Abundance",
             title = paste(ph_protein(), first_name[[ph_protein()]] %||% "")) +
        theme_mg() + theme(strip.text.x = element_text(family = "mono", size = 6),
                           axis.text.x = element_text(angle = 45, hjust = 1),
                           panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA))
    } else if (view == "heatmap") {
      h <- ph_heat_matrix(); M <- h$M
      rng <- range(M, na.rm = TRUE)
      lab <- if (isTRUE(input$ph_rowz)) "z-score" else if (input$ph_scale == "log") "log2 abundance" else "Abundance"
      plots <- lapply(seq_along(SAMPLES), function(j) {
        k <- h$cols$sample == SAMPLES[j]
        d <- data.frame(site = factor(pdf_safe(rep(rownames(M), sum(k))), rev(pdf_safe(rownames(M)))),
                        hpi = factor(rep(h$cols$hpi[k], each = nrow(M)), HPI), value = as.vector(M[, k]))
        ggplot(d, aes(hpi, site, fill = value)) + geom_tile(colour = "white", linewidth = 0.3) +
          scale_fill_gradientn(colours = RAMP_SAMPLE[[j]], limits = rng, name = lab, na.value = "#f2f2ef") +
          labs(x = "Hours post inoculation", y = NULL, title = c("Guy11", "delta-pmk1")[j]) + theme_mg(9) +
          theme(panel.grid = element_blank(), legend.position = "bottom",
                axis.text.y = if (j == 1) element_text(family = "mono", size = 6) else element_blank())
      })
      function() {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = grid::unit(c(0.58, 0.42), "npc"))))
        for (j in 1:2) print(plots[[j]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = j))
      }
    } else if (view == "fc") {
      f <- ph_fc_matrix()
      d <- data.frame(site = factor(pdf_safe(rep(rownames(f$FC), ncol(f$FC))), rev(pdf_safe(rownames(f$FC)))),
                      hpi = factor(rep(HPI, each = nrow(f$FC)), HPI),
                      fc = as.vector(f$FC), sig = as.vector(!is.na(f$P) & f$P < 0.05))
      lim <- max(1, max(abs(d$fc), na.rm = TRUE))
      ggplot(d, aes(hpi, site, fill = fc)) + geom_tile(colour = SURFACE, linewidth = 0.5) +
        geom_point(data = d[d$sig, ], size = 1.4, colour = INK) +
        scale_fill_gradientn(colours = DIV_RAMP, limits = c(-lim, lim), name = "log2FC\n(delta-pmk1 / Guy11)",
                             na.value = "#f4f4f2") +
        labs(x = "Hours post inoculation", y = NULL, title = ph_protein(),
             caption = "Dots: BH-adjusted p < 0.05 (Welch t-test)") +
        theme_mg() + theme(panel.grid = element_blank(), legend.position = "right",
                           axis.text.y = element_text(family = "mono", size = 7))
    } else {
      d <- ph_map_data(); g <- ph_protein(); len <- genes[g, "length_aa"]
      validate(need(!is.null(d), "No sites to plot"))
      if (is.na(len)) len <- max(d$pos) + 20
      d$class <- ifelse(!is.na(d$p_adj) & d$p_adj < 0.05, ifelse(d$log2FC > 0, "Higher in delta-pmk1", "Lower in delta-pmk1"), "Not significant")
      pf <- parse_pfam(genes[g, "pfam"])
      d$lab <- paste0("p", d$residue, d$pos)
      ylim <- max(1.5, max(abs(d$log2FC), na.rm = TRUE) * 1.3)
      trk <- -ylim - 0.6   # protein track below the lollipops
      ggplot(d) +
        geom_hline(yintercept = 0, colour = "#c9c8c2", linewidth = 0.4) +
        annotate("rect", xmin = 1, xmax = len, ymin = trk - 0.12, ymax = trk + 0.12, fill = "#d9d8d3") +
        (if (!is.null(pf)) annotate("rect", xmin = pf$start, xmax = pf$end, ymin = trk - 0.35, ymax = trk + 0.35,
                                    fill = SERIES[(match(pf$name, unique(pf$name)) - 1) %% length(SERIES) + 1])) +
        (if (!is.null(pf)) annotate("text", x = (pf$start + pf$end) / 2, y = trk - 0.5, label = pdf_safe(pf$name),
                                    vjust = 1, size = 2.6, colour = INK2)) +
        geom_segment(aes(x = pos, xend = pos, y = 0, yend = log2FC), colour = "#b5b4ad") +
        geom_point(aes(pos, log2FC, colour = class), size = 2.5) +
        geom_text(aes(pos, log2FC, label = lab, vjust = ifelse(log2FC >= 0, -1, 2)), size = 2.6, colour = INK) +
        coord_cartesian(ylim = c(trk - 1, ylim), clip = "off") +
        scale_colour_manual(values = c("Higher in delta-pmk1" = DIV_RAMP[6], "Lower in delta-pmk1" = DIV_RAMP[2],
                                       "Not significant" = "#a9a8a2"), name = NULL) +
        labs(x = paste0("Residue (", g, ", ", len, " aa)"), y = "log2FC (delta-pmk1 / Guy11)",
             title = paste(g, "phosphosite map")) + theme_mg()
    }
  }
  output$ph_dl_pdf <- pdf_download(ph_gg, function() paste0("MagnaGenes_phospho_", ph_protein(), "_", input$ph_view %||% "profiles"),
                                   reactive(input$ph_pdf_w %||% 11), reactive(input$ph_pdf_h %||% 7))

  # ===========================================================================
  # PHYLOGENOMICS
  # ===========================================================================
  observeEvent(input$og_example, updateTextAreaInput(session, "og_genes", value = paste(EX_PHYLO, collapse = "\n")))

  og_query <- eventReactive(input$og_go, {
    p <- parse_gene_input(input$og_genes)
    ids <- p$ids[p$ids %in% ortho$members$gene]
    validate(need(length(p$ids) + length(p$missing) > 0,
                  "Type gene IDs or names in the box on the left (or load the example), then press Plot orthologs."))
    validate(need(length(ids) > 0, "None of the genes were found in the orthology data."))
    ids <- head(ids, 40)
    ogs <- ortho$members$og[match(ids, ortho$members$gene)]
    sp <- ortho$species[order(ortho$species$tree_y), ]
    M <- t(ortho$counts[ogs, sp$code, drop = FALSE])   # species x genes
    colnames(M) <- ids
    list(ids = ids, ogs = ogs, M = M, species = sp,
         missing = c(p$missing, setdiff(p$ids, ortho$members$gene)))
  }, ignoreNULL = FALSE)

  output$og_summary <- renderUI({
    q <- og_query()
    other <- q$species$code != "MGG"
    n_cons <- colSums(q$M[other, , drop = FALSE] > 0)
    div(class = "summary-strip",
      span(class = "chip ok", tags$b(length(q$ids)), "genes in", tags$b(length(unique(q$ogs))), "orthogroups"),
      span(class = "chip", "Median", tags$b(stats::median(n_cons)), "of", sum(other), "species with orthologs"),
      if (any(n_cons == 0)) span(class = "chip", tags$b(sum(n_cons == 0)), tags$em("M. oryzae"), "-specific"),
      if (length(q$missing)) span(class = "chip warn", icon("triangle-exclamation"), tags$b(length(q$missing)),
                                  "not found:", paste(head(q$missing, 6), collapse = ", ")))
  })

  tree_segments <- function(mode) {
    tr <- ortho$tree
    x <- if (mode == "phylo") tr$x_phylo else tr$x_clado
    x <- x / max(x)
    y <- tr$y; e <- tr$edge
    hx <- as.vector(rbind(x[e[, 1]], x[e[, 2]], NA)); hy <- as.vector(rbind(y[e[, 2]], y[e[, 2]], NA))
    parents <- unique(e[, 1])
    vx <- as.vector(rbind(x[parents], x[parents], NA))
    vy <- as.vector(rbind(vapply(parents, function(p) min(y[e[e[, 1] == p, 2]]), 0),
                          vapply(parents, function(p) max(y[e[e[, 1] == p, 2]]), 0), NA))
    tips <- seq_len(tr$n_tip)
    list(x = c(hx, vx), y = c(hy, vy), tip_x = x[tips], tip_y = y[tips], tip_code = tr$tip_code)
  }

  output$og_plot_ui <- renderUI({
    plotlyOutput("og_plot", height = paste0(250 + 19 * n_species, "px"))
  })
  output$og_plot <- renderPlotly({
    plot_setup("og_plot")
    q <- og_query(); sp <- q$species; M <- q$M
    tr <- tree_segments(input$og_tree_mode)
    tip_name <- ortho$species$name[match(tr$tip_code, ortho$species$code)]
    is_mo <- tr$tip_code == "MGG"
    p_tree <- plot_ly() |>
      add_trace(x = tr$x, y = tr$y, type = "scatter", mode = "lines", line = list(color = MUTED, width = 1.4),
                hoverinfo = "skip", showlegend = FALSE) |>
      add_trace(x = tr$tip_x + 0.02, y = tr$tip_y, type = "scatter", mode = "text",
                text = ifelse(is_mo, paste0("<b><i>", tip_name, "</i></b>"), paste0("<i>", tip_name, "</i>")),
                textposition = "middle right", textfont = list(size = fsz(-2), color = INK),
                hoverinfo = "skip", showlegend = FALSE) |>
      layout(xaxis = list(visible = FALSE, range = c(-0.02, 1.1 + 0.052 * max(nchar(tip_name))), fixedrange = TRUE))

    collab <- gene_short(colnames(M))
    show_counts <- isTRUE(input$og_show_counts)
    zmax <- max(3, stats::quantile(M[M > 0], 0.95, na.rm = TRUE))
    hover <- matrix(sprintf("<i>%s</i> (%s)<br>%s (%s)<br>%d ortholog(s)", rep(sp$name, ncol(M)),
                            rep(sp$lifestyle, ncol(M)), rep(collab, each = nrow(M)), rep(q$ogs, each = nrow(M)),
                            as.vector(M)), nrow = nrow(M))
    # One heatmap trace per lifestyle: cells take the lifestyle colour, darker with more copies
    p_heat <- plot_ly()
    for (ls in unique(sp$lifestyle)) {
      Z <- pmin(M, zmax); Z[sp$lifestyle != ls, ] <- NA
      ramp <- lifestyle_ramp(COL_LIFESTYLE[[ls]] %||% "#cccccc")
      step <- 1 / zmax
      cs <- c(list(list(0, EMPTY), list(step * 0.999, EMPTY)),
              # seq() ends at exactly 1; plotly rejects a scale that stops at 0.9999...
              Map(function(pos, col) list(pos, col), seq(step, 1, length.out = length(ramp)), ramp))
      p_heat <- p_heat |> add_trace(x = collab, y = sp$tree_y, z = Z, type = "heatmap", text = hover, hoverinfo = "text",
                                    colorscale = cs, zmin = 0, zmax = zmax, xgap = 2, ygap = 2, showscale = FALSE)
    }
    if (show_counts) {
      nz <- which(M > 0, arr.ind = TRUE)
      p_heat <- p_heat |> add_trace(x = collab[nz[, 2]], y = sp$tree_y[nz[, 1]], type = "scatter", mode = "text",
                                    inherit = FALSE, text = M[nz],
                                    textfont = list(size = fsz(-3), color = ifelse(pmin(M[nz], zmax) / zmax > 0.55, "#ffffff", "#1d2b2a")),
                                    hoverinfo = "skip", showlegend = FALSE)
    }
    p_heat <- p_heat |> layout(xaxis = list(side = "top", tickangle = -50, tickfont = list(size = fsz(-2), color = INK),
                                            showgrid = FALSE, type = "category",
                                            categoryorder = "array", categoryarray = collab))
    p_strip <- plot_ly()
    for (ls in LIFESTYLES) {
      k <- sp$lifestyle == ls
      if (!any(k)) next
      p_strip <- p_strip |> add_trace(x = rep("Lifestyle", sum(k)), y = sp$tree_y[k], type = "scatter", mode = "markers",
                                      name = ls, legendgroup = "life",
                                      marker = list(symbol = "square", size = 15, color = COL_LIFESTYLE[[ls]],
                                                    line = list(color = SURFACE, width = 1)),
                                      text = paste0("<i>", sp$name[k], "</i><br>", ls), hoverinfo = "text")
    }
    for (ap in names(COL_APPRESS)) {
      k <- sp$appressorium == ap
      if (!any(k)) next
      p_strip <- p_strip |> add_trace(x = rep("Appressorium", sum(k)), y = sp$tree_y[k], type = "scatter",
                                      mode = "markers", name = paste("Appressorium:", tolower(ap)), legendgroup = "app",
                                      marker = list(symbol = "square", size = 15, color = COL_APPRESS[[ap]],
                                                    line = list(color = SURFACE, width = 1)),
                                      text = paste0("<i>", sp$name[k], "</i><br>Appressorium: ", ap), hoverinfo = "text")
    }
    p_strip <- p_strip |> layout(xaxis = list(side = "top", tickangle = -50, type = "category",
                                              categoryorder = "array", categoryarray = c("Lifestyle", "Appressorium"),
                                              range = c(-0.6, 1.6), tickfont = list(size = fsz(-2), color = INK),
                                              showgrid = FALSE))
    w_heat <- min(0.58, 0.05 + 0.04 * ncol(M))
    subplot(p_tree, p_heat, p_strip, widths = c(0.92 - w_heat, w_heat, 0.08), shareY = TRUE,
            margin = 0.004) |>
      finish_plot("ortholog_profile",
                  yaxis = list(visible = FALSE, range = c(0.4, nrow(sp) + 0.6), fixedrange = FALSE),
                  legend = list(font = list(size = fsz(-2), color = INK2)),
                  margin = list(t = 110, r = 10, b = 70, l = 10))
  })

  og_life <- reactive({
    q <- og_query(); sp <- q$species; M <- q$M
    keep <- sp$code != "MGG"
    groups <- c(LIFESTYLES[LIFESTYLES %in% sp$lifestyle[keep]], "All species")
    F <- sapply(groups, function(g) {
      k <- keep & (if (g == "All species") TRUE else sp$lifestyle == g)
      colMeans(M[k, , drop = FALSE] > 0)
    })
    if (is.null(dim(F))) F <- matrix(F, nrow = 1, dimnames = list(colnames(M), groups))
    n <- sapply(groups, function(g) sum(keep & (if (g == "All species") TRUE else sp$lifestyle == g)))
    list(F = F, n = n)
  })

  output$og_life_ui <- renderUI(plotlyOutput("og_life", height = paste0(max(320, 140 + 22 * length(og_query()$ids)), "px")))
  output$og_life <- renderPlotly({
    plot_setup("og_life")
    l <- og_life(); F <- l$F
    rl <- gene_short(rownames(F)); cl <- paste0(colnames(F), " (n=", l$n, ")")
    hover <- matrix(sprintf("%s<br>%s<br>%s%% of species with an ortholog", rep(rl, ncol(F)), rep(cl, each = nrow(F)),
                            round(100 * as.vector(F))), nrow = nrow(F))
    # Each lifestyle column in its own colour (grey for all species), darker = more conserved
    p <- plot_ly()
    for (j in seq_len(ncol(F))) {
      Z <- F; Z[, -j] <- NA
      base <- COL_LIFESTYLE[colnames(F)[j]]
      ramp <- if (is.na(base)) c("#f2f2f2", "#4f4f4f") else lifestyle_ramp(base)
      p <- p |> add_trace(x = cl, y = rl, z = Z, type = "heatmap", zmin = 0, zmax = 1, text = hover, hoverinfo = "text",
                          colorscale = c(list(list(0, EMPTY)),
                                         Map(function(pos, col) list(pos, col),
                                             seq(0.001, 1, length.out = length(ramp)), ramp)),
                          xgap = 3, ygap = 3, showscale = FALSE)
    }
    p |> add_trace(x = rep(cl, each = nrow(F)), y = rep(rl, ncol(F)), type = "scatter", mode = "text", inherit = FALSE,
                   text = paste0(round(100 * as.vector(F)), "%"), hoverinfo = "skip", showlegend = FALSE,
                   textfont = list(size = fsz(-2), color = ifelse(as.vector(F) > 0.6, "#ffffff", "#1d2b2a"))) |>
      finish_plot("ortholog_lifestyle", legend = FALSE,
                  xaxis = ax(NULL, side = "top", tickangle = -35, showgrid = FALSE, showline = FALSE, type = "category",
                             categoryorder = "array", categoryarray = cl),
                  yaxis = ax(NULL, autorange = "reversed", showgrid = FALSE, showline = FALSE, type = "category",
                             categoryorder = "array", categoryarray = rl),
                  margin = list(t = 100, r = 10, b = 20, l = 10))
  })

  og_members_df <- reactive({
    q <- og_query()
    mem <- ortho$members[ortho$members$og %in% q$ogs, ]
    mem <- mem[order(match(mem$og, q$ogs), mem$gene), ]
    data.frame(orthogroup = mem$og, gene = mem$gene,
               name = ifelse(is.na(first_name[mem$gene]), "", first_name[mem$gene]),
               description = gene_desc(mem$gene),
               in_query = mem$gene %in% q$ids, stringsAsFactors = FALSE)
  })

  output$og_members <- renderReactable({
    d <- og_members_df()
    d$links <- vapply(d$gene, function(g) if (g %in% genes$gene) xref_links_html(g) else "", "", USE.NAMES = FALSE)
    rownames(d) <- NULL
    reactable(d, groupBy = "orthogroup", compact = TRUE, highlight = TRUE, searchable = TRUE, defaultExpanded = TRUE,
              defaultPageSize = 30,
              columns = list(
                orthogroup = colDef(name = "Orthogroup", minWidth = 120),
                gene = colDef(name = "Gene", cell = function(v) tags$span(class = "gene-id", v)),
                in_query = colDef(name = "In query", cell = function(v) if (isTRUE(v)) "✓" else "", maxWidth = 80, align = "center"),
                description = colDef(minWidth = 240),
                links = colDef(name = "Databases", html = TRUE, minWidth = 220, sortable = FALSE)))
  })

  og_counts_df <- reactive({
    q <- og_query()
    data.frame(species = q$species$name, code = q$species$code, lifestyle = q$species$lifestyle,
               appressorium = q$species$appressorium, q$M[, , drop = FALSE], check.names = FALSE)[nrow(q$M):1, ]
  })

  output$og_table <- renderReactable({
    d <- og_counts_df(); rownames(d) <- NULL
    q <- og_query()
    gcols <- setNames(lapply(q$ids, function(i) colDef(name = gene_short(i), align = "center", minWidth = 90,
      style = function(v) if (v == 0) list(color = "#c9c8c2") else list(fontWeight = 600))), q$ids)
    reactable(d, compact = TRUE, highlight = TRUE, pagination = FALSE, searchable = TRUE,
              columns = c(list(species = colDef(name = "Species", minWidth = 220, style = list(fontStyle = "italic"),
                                                sticky = "left"),
                               code = colDef(show = FALSE)), gcols))
  })

  output$og_dl_counts <- downloadHandler(
    filename = function() paste0("MagnaGenes_ortholog_counts_", format(Sys.Date(), "%Y%m%d"), ".tsv"),
    content = function(file) {
      d <- og_counts_df(); q <- og_query()
      hdr <- paste0("# orthogroups: ", paste(paste0(q$ids, "=", q$ogs), collapse = ", "))
      writeLines(hdr, file)
      suppressWarnings(utils::write.table(d, file, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE))
    })
  output$og_dl_members <- downloadHandler(
    filename = function() paste0("MagnaGenes_orthogroup_members_", format(Sys.Date(), "%Y%m%d"), ".tsv"),
    content = function(file) write_tsv(og_members_df(), file))
  output$og_dl_fasta <- downloadHandler(
    filename = function() paste0("MagnaGenes_orthogroup_members_protein_", format(Sys.Date(), "%Y%m%d"), ".fasta"),
    content = function(file) writeLines(fasta_text(og_members_df()$gene, "protein"), file))

  og_gg <- function() {
    view <- input$og_view %||% "profile"
    q <- og_query(); M <- q$M; sp <- q$species
    if (view == "lifestyle") {
      l <- og_life(); F <- l$F
      d <- data.frame(gene = factor(rep(gene_short(rownames(F)), ncol(F)), rev(gene_short(rownames(F)))),
                      group = factor(rep(paste0(colnames(F), "\n(n=", l$n, ")"), each = nrow(F)),
                                     paste0(colnames(F), "\n(n=", l$n, ")")),
                      f = as.vector(F))
      d$life <- rep(colnames(F), each = nrow(F))
      d$fill <- vapply(seq_len(nrow(d)), function(i) {
        base <- COL_LIFESTYLE[d$life[i]]
        ramp <- if (is.na(base)) c("#f2f2f2", "#4f4f4f") else lifestyle_ramp(base)
        if (d$f[i] == 0) "#f2f2ef" else ramp[max(1, ceiling(d$f[i] * length(ramp)))]
      }, "")
      return(ggplot(d, aes(group, gene, fill = fill)) + geom_tile(colour = "white", linewidth = 0.6) +
               geom_text(aes(label = paste0(round(100 * f), "%")), size = 2.6) + scale_fill_identity() +
               labs(x = NULL, y = NULL) + theme_mg() +
               theme(panel.grid = element_blank(), legend.position = "right"))
    }
    tr <- tree_segments(input$og_tree_mode)
    T_W <- 4; L_W <- 6.5; n <- ncol(M)
    seg <- data.frame(x = tr$x * T_W - (T_W + L_W), y = tr$y)
    tips <- data.frame(x = tr$tip_x * T_W - (T_W + L_W) + 0.1, y = tr$tip_y,
                       name = ortho$species$name[match(tr$tip_code, ortho$species$code)])
    cells <- data.frame(x = rep(seq_len(n), each = nrow(M)), y = rep(sp$tree_y, n), count = as.vector(M),
                        lifestyle = rep(sp$lifestyle, n))
    strips <- rbind(data.frame(x = n + 1.2, y = sp$tree_y, fill = sp$lifestyle, kind = "Lifestyle"),
                    data.frame(x = n + 2.2, y = sp$tree_y, fill = paste("Appressorium:", sp$appressorium), kind = "Appressorium"))
    strips <- strips[strips$fill != "Appressorium: None", ]   # no appressorium: left blank
    zmax <- max(3, stats::quantile(M[M > 0], 0.95, na.rm = TRUE))
    cells$fill <- vapply(seq_len(nrow(cells)), function(i) {
      if (cells$count[i] == 0) return("#f2f2ef")
      ramp <- lifestyle_ramp(COL_LIFESTYLE[[cells$lifestyle[i]]] %||% "#cccccc")
      ramp[max(1, ceiling(min(cells$count[i], zmax) / zmax * length(ramp)))]
    }, "")
    fill_vals <- c(COL_LIFESTYLE, setNames(COL_APPRESS, paste("Appressorium:", names(COL_APPRESS))))
    draw <- function() {
      g <- ggplot() +
        geom_path(data = seg, aes(x, y), colour = "#8a8984", linewidth = 0.35, na.rm = TRUE) +
        geom_text(data = tips, aes(x, y, label = name), hjust = 0, size = 2.6, fontface = "italic", colour = INK) +
        geom_tile(data = cells, aes(x, y, fill = fill), colour = "white", linewidth = 0.4,
                  width = 0.95, height = 0.95) +
        scale_fill_identity() +
        (if (isTRUE(input$og_show_counts)) geom_text(data = cells[cells$count > 0, ],
                                                    aes(x, y, label = count), size = 2.2,
                                                    colour = ifelse(cells$count[cells$count > 0] / zmax > 0.55, "white", INK))) +
        annotate("text", x = seq_len(n), y = max(sp$tree_y) + 0.9, label = gene_short(colnames(M)),
                 angle = 50, hjust = 0, size = 2.6, colour = INK) +
        annotate("text", x = c(n + 1.2, n + 2.2), y = max(sp$tree_y) + 0.9, label = c("Lifestyle", "Appressorium"),
                 angle = 50, hjust = 0, size = 2.6, colour = INK) +
        geom_point(data = strips, aes(x, y, colour = fill), shape = 15, size = 3.2) +
        scale_colour_manual(values = fill_vals, name = NULL) +
        coord_cartesian(clip = "off", ylim = c(0.5, max(sp$tree_y) + 0.5)) +
        theme_void(base_size = 9) +
        theme(legend.position = "bottom", plot.margin = margin(60, 10, 10, 10),
              legend.text = element_text(size = 7)) +
        guides(colour = guide_legend(nrow = 3, override.aes = list(size = 3)))
      print(g)
    }
    draw
  }
  output$og_dl_pdf <- pdf_download(og_gg, function() paste0("MagnaGenes_orthologs_", input$og_view %||% "profile"),
                                   reactive(input$og_pdf_w %||% 12), reactive(input$og_pdf_h %||% 9))

  # ===========================================================================
  # GO ENRICHMENT
  # ===========================================================================
  observeEvent(input$go_example, {
    updateTextAreaInput(session, "go_genes", value = paste(ex_go_genes(input$go_example_module), collapse = "\n"))
  })

  go_result <- eventReactive(input$go_run, {
    p <- parse_gene_input(input$go_genes)
    validate(need(length(p$ids) > 0,
                  "Paste a gene list in the box on the left (or load a WGCNA module example), then press Run enrichment."))
    bg <- if (identical(input$go_bg_mode, "custom")) {
      b <- parse_gene_input(input$go_bg)$ids
      validate(need(length(b) > 0, "The custom background is empty."))
      b
    } else NULL
    r <- go_enrich(p$ids, input$go_ont, bg, input$go_min %||% 5, input$go_max %||% 500, input$go_adjust)
    r$ids <- p$ids; r$missing <- p$missing; r$ont <- input$go_ont
    r$cutoff <- input$go_cutoff %||% 0.05
    r
  }, ignoreNULL = FALSE)

  go_sig <- reactive({
    r <- go_result()
    validate(need(!is.null(r$table) && nrow(r$table) > 0,
                  "No genes from the list are annotated in this ontology/background."))
    s <- r$table[r$table$p_adj <= r$cutoff, ]
    validate(need(nrow(s) > 0, sprintf("No terms reach adjusted p ≤ %s. Try a different ontology, relax the cutoff, or check the results table.", r$cutoff)))
    head(s, input$go_top)
  })

  output$go_summary <- renderUI({
    r <- go_result()
    n_sig <- if (!is.null(r$table)) sum(r$table$p_adj <= r$cutoff) else 0
    div(class = "summary-strip",
      span(class = "chip", tags$b(length(r$ids)), "genes in list"),
      span(class = "chip ok", tags$b(r$n_query), "annotated in", ONT_NAMES[[r$ont]]),
      span(class = "chip", "background", tags$b(format(r$n_universe, big.mark = ","))),
      span(class = "chip", tags$b(format(r$n_tested, big.mark = ",")), "terms tested"),
      span(class = if (n_sig) "chip ok" else "chip warn", tags$b(n_sig), "significant terms"),
      if (length(r$missing)) span(class = "chip warn", tags$b(length(r$missing)), "IDs not recognised"))
  })

  go_term_label <- function(t, n = 60) ifelse(nchar(t) > n, paste0(substr(t, 1, n - 1), "…"), t)

  # Colours chosen by the user: one colour (picker) or viridis
  go_main_col <- reactive(if (isTRUE(input$go_viridis)) VIRIDIS[5] else (input$go_color %||% SERIES[1]))
  go_ramp <- reactive(if (isTRUE(input$go_viridis)) VIRIDIS else tint_ramp(go_main_col())[-1])

  go_height <- reactive(paste0(max(380, 140 + 28 * min(input$go_top, nrow(go_result()$table %||% data.frame()))), "px"))
  output$go_dot_ui <- renderUI(plotlyOutput("go_dot", height = go_height()))
  output$go_bar_ui <- renderUI(plotlyOutput("go_bar", height = go_height()))
  output$go_matrix_ui <- renderUI(plotlyOutput("go_matrix", height = go_height()))
  output$go_tree_ui <- renderUI(plotlyOutput("go_tree", height = "620px"))
  output$go_volcano_ui <- renderUI(plotlyOutput("go_volcano", height = "560px"))

  go_hover <- function(s) sprintf("<b>%s</b><br>%s<br>%d of %d genes (expected %s)<br>Fold enrichment %s<br>p = %s<br>adj. p = %s",
                                  s$term, s$go_id, s$significant, s$annotated, s$expected, s$fold_enrichment,
                                  signif(s$p_value, 3), signif(s$p_adj, 3))
  nlp_bar <- function() list(title = list(text = "−log₁₀ adj. p", side = "right",
                                          font = list(size = fsz(-1), color = INK2)),
                             tickfont = list(size = fsz(-2), color = INK2), thickness = 14, len = 0.5)

  output$go_dot <- renderPlotly({
    plot_setup("go_dot")
    s <- go_sig(); s$label <- go_term_label(s$term)
    s$nlp <- -log10(pmax(s$p_adj, 1e-300))
    lev <- rev(s$label)
    p <- plot_ly()
    if (isTRUE(input$go_lollipop)) {
      p <- p |> add_trace(x = as.vector(rbind(0, s$fold_enrichment, NA)), y = as.vector(rbind(s$label, s$label, NA)),
                          type = "scatter", mode = "lines", line = list(color = AXIS, width = 2.5),
                          hoverinfo = "skip", showlegend = FALSE)
    }
    p |> add_trace(data = s, x = ~fold_enrichment, y = ~label, type = "scatter", mode = "markers",
                   text = go_hover(s), hoverinfo = "text", showlegend = FALSE,
                   marker = list(size = ~significant, sizemode = "area", sizeref = 2 * max(s$significant) / 28^2,
                                 sizemin = 6, color = ~nlp, colorscale = colorscale(go_ramp()), showscale = TRUE,
                                 line = list(color = SURFACE, width = 1.5), colorbar = nlp_bar())) |>
      finish_plot(paste0("GO_", go_result()$ont, "_dotplot"), legend = FALSE,
                  xaxis = ax("Fold enrichment", rangemode = "tozero"),
                  yaxis = ax(NULL, type = "category", categoryorder = "array", categoryarray = lev, showgrid = TRUE,
                             tickfont = list(size = fsz(-1), color = INK)),
                  annotations = list(list(text = "Marker area = genes in term", xref = "paper", yref = "paper",
                                          x = 1, y = 1.01, xanchor = "right", yanchor = "bottom", showarrow = FALSE,
                                          font = list(size = fsz(-2), color = INK2))),
                  margin = list(t = 30, r = 10, b = 50, l = 10))
  })

  output$go_bar <- renderPlotly({
    plot_setup("go_bar")
    s <- go_sig(); s$label <- go_term_label(s$term)
    s$nlp <- -log10(pmax(s$p_adj, 1e-300))
    plot_ly(s, x = ~nlp, y = ~label, type = "bar", orientation = "h", text = go_hover(s), hoverinfo = "text",
            textposition = "none", marker = list(color = go_main_col())) |>
      finish_plot(paste0("GO_", go_result()$ont, "_barchart"), legend = FALSE, bargap = 0.35,
                  shapes = list(list(type = "line", x0 = -log10(go_result()$cutoff), x1 = -log10(go_result()$cutoff),
                                     y0 = 0, y1 = 1, yref = "paper", line = list(color = INK2, width = 1))),
                  xaxis = ax("−log₁₀ adjusted p-value", rangemode = "tozero"),
                  yaxis = ax(NULL, categoryorder = "array", categoryarray = rev(s$label), showgrid = FALSE,
                             tickfont = list(size = fsz(-1), color = INK)),
                  margin = list(t = 20, r = 10, b = 50, l = 10))
  })

  output$go_tree <- renderPlotly({
    plot_setup("go_tree")
    s <- go_sig(); s$nlp <- -log10(pmax(s$p_adj, 1e-300))
    plot_ly(type = "treemap", labels = go_term_label(s$term, 45), parents = rep("", nrow(s)), values = s$significant,
            ids = s$go_id, branchvalues = "remainder", text = go_hover(s), hoverinfo = "text",
            textinfo = "label+value", textfont = list(size = fsz(0), family = FONT),
            marker = list(colors = s$nlp, colorscale = colorscale(go_ramp()), showscale = TRUE,
                          line = list(color = SURFACE, width = 2), colorbar = nlp_bar()),
            tiling = list(pad = 2), root = list(color = SURFACE), pathbar = list(visible = FALSE)) |>
      finish_plot(paste0("GO_", go_result()$ont, "_treemap"), legend = FALSE,
                  margin = list(t = 10, r = 10, b = 10, l = 10))
  })

  output$go_volcano <- renderPlotly({
    plot_setup("go_volcano")
    r <- go_result()
    validate(need(!is.null(r$table) && nrow(r$table) > 0, "No results."))
    d <- r$table; d$nlp <- -log10(pmax(d$p_adj, 1e-300)); d$sig <- d$p_adj <= r$cutoff
    top <- head(d[d$sig, ], 8)
    p <- plot_ly()
    if (any(!d$sig)) p <- p |> add_trace(data = d[!d$sig, ], x = ~fold_enrichment, y = ~nlp, type = "scatter", mode = "markers",
                                          name = "Not significant", text = go_hover(d[!d$sig, ]), hoverinfo = "text",
                                          marker = list(size = ~significant, sizemode = "area", sizeref = 2 * max(d$significant) / 26^2,
                                                        sizemin = 4, color = AXIS, opacity = 0.6,
                                                        line = list(color = SURFACE, width = 1)))
    if (any(d$sig)) p <- p |> add_trace(data = d[d$sig, ], x = ~fold_enrichment, y = ~nlp, type = "scatter", mode = "markers",
                                         name = paste("adj. p ≤", r$cutoff), text = go_hover(d[d$sig, ]), hoverinfo = "text",
                                         marker = list(size = ~significant, sizemode = "area", sizeref = 2 * max(d$significant) / 26^2,
                                                       sizemin = 5, color = go_main_col(), opacity = 0.85,
                                                       line = list(color = SURFACE, width = 1.5)))
    p |> finish_plot(paste0("GO_", r$ont, "_significance_vs_enrichment"),
                     shapes = list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = -log10(r$cutoff), y1 = -log10(r$cutoff),
                                        line = list(color = INK2, width = 1, dash = "dash"))),
                     annotations = lapply(seq_len(nrow(top)), function(i) list(
                       x = top$fold_enrichment[i], y = top$nlp[i], text = go_term_label(top$term[i], 34),
                       showarrow = TRUE, arrowhead = 0, arrowcolor = AXIS, ax = 30, ay = -18,
                       font = list(size = fsz(-2), color = INK))),
                     xaxis = ax("Fold enrichment", rangemode = "tozero"),
                     yaxis = ax("−log₁₀ adjusted p-value", rangemode = "tozero"),
                     margin = list(t = 40, r = 20, b = 50, l = 60))
  })

  go_matrix_data <- reactive({
    s <- go_sig(); pr <- go_result()$pairs
    pr <- pr[pr$go %in% s$go_id, ]
    gl <- names(sort(table(pr$gene), decreasing = TRUE))
    gl <- head(gl, 60)
    M <- matrix(0, nrow(s), length(gl), dimnames = list(s$go_id, gl))
    pr <- pr[pr$gene %in% gl, ]
    M[cbind(match(pr$go, s$go_id), match(pr$gene, gl))] <- 1
    list(M = M, s = s)
  })

  output$go_matrix <- renderPlotly({
    plot_setup("go_matrix")
    md <- go_matrix_data(); M <- md$M; s <- md$s
    tl <- go_term_label(s$term, 50); gl <- gene_short(colnames(M))
    hover <- matrix(sprintf("%s<br>%s<br>%s", rep(tl, ncol(M)), rep(gl, each = nrow(M)),
                            ifelse(as.vector(M) == 1, "annotated", "not annotated")), nrow = nrow(M))
    col <- go_main_col()
    plot_ly(x = gl, y = tl, z = M, type = "heatmap", zmin = 0, zmax = 1, text = hover, hoverinfo = "text",
            colorscale = list(list(0, EMPTY), list(0.5, EMPTY), list(0.5, col), list(1, col)),
            showscale = FALSE, xgap = 2, ygap = 2) |>
      finish_plot(paste0("GO_", go_result()$ont, "_gene_term_matrix"), legend = FALSE,
                  xaxis = ax(NULL, side = "top", tickangle = -60, showgrid = FALSE, showline = FALSE,
                             tickfont = list(size = fsz(-3), color = INK)),
                  yaxis = ax(NULL, autorange = "reversed", showgrid = FALSE, showline = FALSE,
                             tickfont = list(size = fsz(-2), color = INK)),
                  margin = list(t = 110, r = 10, b = 20, l = 10))
  })

  output$go_table <- renderReactable({
    r <- go_result()
    validate(need(!is.null(r$table) && nrow(r$table) > 0, "No results."))
    d <- r$table
    d$go_link <- sprintf('<a href="%s" target="_blank" rel="noopener">%s</a>', url_quickgo(d$go_id), d$go_id)
    d <- d[, c("go_link", "term", "annotated", "significant", "expected", "fold_enrichment", "p_value", "p_adj", "genes")]
    reactable(d, compact = TRUE, highlight = TRUE, searchable = TRUE, defaultPageSize = 20, filterable = FALSE,
              rowStyle = function(i) if (d$p_adj[i] <= r$cutoff) list(background = "var(--chip-ok)") else NULL,
              columns = list(
                go_link = colDef(name = "GO ID", html = TRUE, minWidth = 110),
                term = colDef(name = "Term", minWidth = 260),
                annotated = colDef(name = "Size", maxWidth = 70),
                significant = colDef(name = "In list", maxWidth = 70),
                expected = colDef(name = "Expected", maxWidth = 85),
                fold_enrichment = colDef(name = "Fold", maxWidth = 70),
                p_value = colDef(name = "p", cell = function(v) formatC(v, format = "g", digits = 2), maxWidth = 85),
                p_adj = colDef(name = "adj. p", cell = function(v) formatC(v, format = "g", digits = 2), maxWidth = 85),
                genes = colDef(name = "Genes", minWidth = 260, style = list(fontSize = "0.8rem", fontFamily = "ui-monospace, Menlo, monospace"))))
  })

  output$go_dl_table <- downloadHandler(
    filename = function() paste0("MagnaGenes_GO_", go_result()$ont, "_enrichment_", format(Sys.Date(), "%Y%m%d"), ".tsv"),
    content = function(file) write_tsv(go_result()$table, file))
  output$go_dl_pairs <- downloadHandler(
    filename = function() paste0("MagnaGenes_GO_", go_result()$ont, "_gene_term_pairs_", format(Sys.Date(), "%Y%m%d"), ".tsv"),
    content = function(file) {
      pr <- go_result()$pairs
      pr$term <- go_terms$term[match(pr$go, go_terms$go)]
      write_tsv(pr, file)
    })

  # Squarified treemap layout (for the PDF version of the treemap)
  squarify <- function(values, x = 0, y = 0, w = 1, h = 1) {
    a <- values / sum(values) * w * h
    out <- matrix(NA_real_, length(a), 4)
    worst <- function(r, side) { s <- sum(r); max(side^2 * max(r) / s^2, s^2 / (side^2 * min(r))) }
    place <- function(idx) {
      s <- sum(a[idx])
      if (w >= h) {
        cw <- s / h; yy <- y
        for (i in idx) { hh <- a[i] / cw; out[i, ] <<- c(x, x + cw, yy, yy + hh); yy <- yy + hh }
        x <<- x + cw; w <<- w - cw
      } else {
        rh <- s / w; xx <- x
        for (i in idx) { ww <- a[i] / rh; out[i, ] <<- c(xx, xx + ww, y, y + rh); xx <- xx + ww }
        y <<- y + rh; h <<- h - rh
      }
    }
    row <- integer(0); i <- 1
    while (i <= length(a)) {
      side <- min(w, h)
      if (!length(row) || worst(a[c(row, i)], side) <= worst(a[row], side)) { row <- c(row, i); i <- i + 1 }
      else { place(row); row <- integer(0) }
    }
    if (length(row)) place(row)
    data.frame(xmin = out[, 1], xmax = out[, 2], ymin = out[, 3], ymax = out[, 4])
  }

  go_gg <- function() {
    view <- input$go_view %||% "dot"
    ont <- ONT_NAMES[[go_result()$ont]]
    col <- go_main_col(); ramp <- go_ramp()
    if (view == "matrix") {
      md <- go_matrix_data(); M <- md$M
      tl <- pdf_safe(go_term_label(md$s$term, 50))
      d <- data.frame(term = factor(rep(tl, ncol(M)), rev(unique(tl))),
                      gene = factor(rep(gene_short(colnames(M)), each = nrow(M)), gene_short(colnames(M))),
                      v = factor(as.vector(M), c(0, 1)))
      return(ggplot(d, aes(gene, term, fill = v)) + geom_tile(colour = "white", linewidth = 0.5) +
               scale_fill_manual(values = c("0" = "#f2f2ef", "1" = col), guide = "none") +
               scale_x_discrete(position = "top") + labs(x = NULL, y = NULL, title = ont) + theme_mg() +
               theme(axis.text.x = element_text(angle = 60, hjust = 0, size = 6), panel.grid = element_blank()))
    }
    if (view == "volcano") {
      r <- go_result(); d <- r$table; d$nlp <- -log10(pmax(d$p_adj, 1e-300)); d$sig <- d$p_adj <= r$cutoff
      return(ggplot(d, aes(fold_enrichment, nlp, size = significant)) +
               geom_hline(yintercept = -log10(r$cutoff), colour = INK2, linewidth = 0.3, linetype = "dashed") +
               geom_point(aes(colour = sig), alpha = 0.75) +
               scale_colour_manual(values = c("FALSE" = "#c9c8c2", "TRUE" = col), labels = c("Not significant", "Significant"),
                                   name = NULL) +
               scale_size_area(max_size = 7, name = "Genes") +
               labs(x = "Fold enrichment", y = "-log10 adjusted p-value", title = ont) + theme_mg() +
               theme(legend.position = "right"))
    }
    s <- go_sig(); s$label <- factor(pdf_safe(go_term_label(s$term)), rev(unique(pdf_safe(go_term_label(s$term)))))
    s$nlp <- -log10(pmax(s$p_adj, 1e-300))
    if (view == "treemap") {
      tm <- cbind(squarify(s$significant), s)
      tm$lab <- paste0(pdf_safe(go_term_label(s$term, 40)), "\n(", s$significant, ")")
      return(ggplot(tm) +
               geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = nlp), colour = "white", linewidth = 0.8) +
               geom_text(aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = lab), size = 2.4, colour = "white",
                         lineheight = 0.9) +
               scale_fill_gradientn(colours = ramp, name = "-log10 adj. p") +
               scale_y_reverse() + labs(title = ont) + theme_void(base_size = 10) +
               theme(plot.title = element_text(face = "bold"), legend.position = "right"))
    }
    if (view == "bar") {
      ggplot(s, aes(nlp, label)) + geom_col(fill = col, width = 0.65) +
        geom_vline(xintercept = -log10(go_result()$cutoff), colour = INK2, linewidth = 0.3) +
        labs(x = "-log10 adjusted p-value", y = NULL, title = ont) + theme_mg()
    } else {
      ggplot(s, aes(fold_enrichment, label)) +
        (if (isTRUE(input$go_lollipop)) geom_segment(aes(x = 0, xend = fold_enrichment, yend = label),
                                                    colour = "#b5b4ad", linewidth = 0.7)) +
        geom_point(aes(size = significant, colour = nlp)) +
        scale_colour_gradientn(colours = ramp, name = "-log10 adj. p") +
        scale_size_area(max_size = 7, name = "Genes") +
        labs(x = "Fold enrichment", y = NULL, title = ont) + theme_mg() +
        theme(legend.position = "right")
    }
  }
  output$go_dl_pdf <- pdf_download(go_gg, function() paste0("MagnaGenes_GO_", go_result()$ont, "_", input$go_view %||% "dot"),
                                   reactive(input$go_pdf_w %||% 10), reactive(input$go_pdf_h %||% 7))
}

shinyApp(ui, server)
