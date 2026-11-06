# Smoke test: renders every output and download of the app with default inputs.
# Run from the repository root:  Rscript tests/smoke_test.R
suppressPackageStartupMessages(library(shiny))
setwd("app")
app_env <- new.env()
suppressPackageStartupMessages(sys.source("app.R", envir = app_env, toplevel.env = app_env))
server <- app_env$server

check <- function(label, expr) {
  res <- tryCatch({ force(expr); "ok" }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-28s %s\n", label, res))
  invisible(res == "ok")
}

testServer(server, {
  session$setInputs(
    s_mode = "id", s_ids = paste(app_env$EX_SIGNALLING, collapse = "\n"), s_field = "all", s_kw = "",
    s_effector = FALSE, s_pheno = FALSE, s_isoforms = FALSE, s_go = 0, s_profile_gene = "MGG_09565",
    tx_genes = paste(app_env$EX_SIGNALLING, collapse = "\n"), tx_exps = names(app_env$EXP_SHORT),
    tx_scale = "log", tx_cluster = TRUE, tx_go = 0, tx_view = "profiles", tx_pdf_w = 11, tx_pdf_h = 7,
    ph_protein = "MGG_00800", ph_scale = "log", ph_zero_na = TRUE, ph_rowz = FALSE, ph_map_time = "max",
    ph_view = "profiles", ph_pdf_w = 11, ph_pdf_h = 7,
    og_genes = paste(app_env$EX_PHYLO, collapse = "\n"), og_tree_mode = "phylo", og_show_counts = TRUE,
    og_go = 0, og_view = "profile", og_pdf_w = 12, og_pdf_h = 9,
    go_genes = paste(app_env$ex_go_genes(), collapse = "\n"), go_ont = "BP", go_bg_mode = "all",
    go_min = 5, go_max = 500, go_adjust = "BH", go_cutoff = 0.05, go_top = 20, go_run = 0,
    go_view = "dot", go_pdf_w = 10, go_pdf_h = 7)
  session$setInputs(ph_pids = as.character(head(app_env$phos$meta$pid[app_env$phos$meta$protein == "MGG_00800"], 6)))

  outs <- c("s_summary", "s_table", "s_profile", "s_profile_expr", "s_profile_arch",
            "tx_summary", "tx_profiles_ui", "tx_profiles", "tx_heat", "tx_corr", "tx_table",
            "ph_header", "ph_peptides_ui", "ph_profiles", "ph_heat", "ph_fc", "ph_map", "ph_table",
            "og_summary", "og_plot", "og_life", "og_members", "og_table",
            "go_summary", "go_dot", "go_bar", "go_matrix", "go_table")
  for (o in outs) check(o, output[[o]])

  dls <- c("s_dl_tsv", "s_dl_pep", "s_dl_cds", "s_dl_gene_pep", "tx_dl_long", "tx_dl_wide", "tx_dl_fasta",
           "ph_dl_data", "ph_dl_stats", "ph_dl_fasta", "og_dl_counts", "og_dl_members", "og_dl_fasta",
           "go_dl_table", "go_dl_pairs")
  for (d in dls) check(d, { f <- output[[d]]; stopifnot(file.exists(f), file.size(f) > 0) })

  pdfs <- list(tx_dl_pdf = c("tx_view", "profiles", "heatmap", "corr"),
               ph_dl_pdf = c("ph_view", "profiles", "heatmap", "fc", "map"),
               og_dl_pdf = c("og_view", "profile", "lifestyle"),
               go_dl_pdf = c("go_view", "dot", "bar", "matrix"))
  for (d in names(pdfs)) {
    v <- pdfs[[d]]
    for (view in v[-1]) {
      do.call(session$setInputs, setNames(list(view), v[1]))
      check(paste(d, view), { f <- output[[d]]; stopifnot(file.size(f) > 1000) })
    }
  }

  # Keyword search and alternative settings
  session$setInputs(s_mode = "kw", s_kw = "melanin, septin", s_go = 1)
  check("keyword search", stopifnot(nrow(s_result()$df) > 0))
  session$setInputs(tx_scale = "z", tx_go = 1)
  check("tx z-score heatmap", output$tx_heat)
  session$setInputs(og_tree_mode = "clado", og_go = 1)
  check("og cladogram", output$og_plot)
  session$setInputs(go_ont = "MF", go_run = 1)
  check("go MF", output$go_table)
})
print(unique(names(warnings())))
