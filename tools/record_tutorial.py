"""
Record the MagnaGenes video tutorial.

Drives the built site (site/) through every tab with captions, then writes
site-assets/tutorial.mp4 and site-assets/tutorial_chapters.json, which
build_site.R copies next to the published app (Tutorial button in the navbar).

Requirements (not needed by the app itself):
    pip install playwright imageio-ffmpeg     # uses the installed Google Chrome

Usage, from the repository root after `Rscript build_site.R`:
    python -m http.server 8765 --bind 127.0.0.1 --directory site &
    python tools/record_tutorial.py
"""
import json
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import imageio_ffmpeg
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
URL = "http://127.0.0.1:8765/"
OUT_VIDEO = ROOT / "site-assets" / "tutorial.mp4"
OUT_CHAPTERS = ROOT / "site-assets" / "tutorial_chapters.json"
VIEW = {"width": 1440, "height": 900}
LAUNCH = dict(channel="chrome", headless=True, args=["--no-proxy-server"], viewport=VIEW)

CAPTION_JS = """
(text) => {
  let c = document.getElementById('tut-caption');
  if (!c) {
    c = document.createElement('div'); c.id = 'tut-caption';
    c.style.cssText = 'position:fixed;left:50%;bottom:28px;transform:translateX(-50%);z-index:99999;' +
      'max-width:900px;padding:14px 26px;border-radius:12px;background:rgba(20,40,34,.92);color:#fff;' +
      'font:600 21px/1.4 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Arial,sans-serif;' +
      'box-shadow:0 6px 24px rgba(0,0,0,.25);text-align:center;transition:opacity .3s';
    document.body.appendChild(c);
  }
  c.style.opacity = text ? '1' : '0';
  if (text) c.innerHTML = text;
}
"""


def main():
    workdir = Path(tempfile.mkdtemp(prefix="mg_tutorial_"))
    profile = workdir / "profile"
    with sync_playwright() as p:
        # 1. Warm-up visit (not recorded) so packages come from the browser cache
        ctx = p.chromium.launch_persistent_context(str(profile), **LAUNCH)
        pg = ctx.new_page()
        pg.goto(URL)
        pg.frame_locator("#root iframe").locator(".navbar").wait_for(timeout=240000)
        ctx.close()

        # 2. Recorded visit
        ctx = p.chromium.launch_persistent_context(str(profile), record_video_dir=str(workdir),
                                                   record_video_size=VIEW, **LAUNCH)
        page = ctx.new_page()
        t_video = time.time()
        page.goto(URL)
        app = page.frame_locator("#root iframe")
        app.locator(".navbar").wait_for(timeout=240000)
        page.wait_for_selector("#mg-loader", state="detached", timeout=60000)
        time.sleep(0.5)
        t0 = time.time() - t_video          # video time where the tutorial starts
        frame = [f for f in page.frames if f != page.main_frame][0]
        chapters = []

        def now():
            return time.time() - t_video - t0

        def chapter(title):
            chapters.append({"time": round(max(0.0, now()), 1), "title": title})

        def caption(text, hold=0.0):
            page.evaluate(CAPTION_JS, text)
            if hold:
                time.sleep(hold)

        def wait_idle(timeout=25):
            # Shiny marks <html> with "shiny-busy" while R is computing; wait until it
            # has been idle for a moment, then give plotly time to draw
            end = time.time() + timeout
            quiet = 0.0
            time.sleep(0.4)
            while time.time() < end and quiet < 1.0:
                busy = frame.evaluate("document.documentElement.classList.contains('shiny-busy')")
                quiet = 0.0 if busy else quiet + 0.2
                time.sleep(0.2)
            time.sleep(0.8)

        def click(loc, wait=1.2):
            loc.scroll_into_view_if_needed()
            loc.evaluate("el => { el.style.outline = '3px solid #f5a623'; el.style.outlineOffset = '3px'; }")
            time.sleep(0.7)
            loc.click()
            loc.evaluate("el => { el.style.outline = ''; }")
            wait_idle()
            time.sleep(wait)

        def scroll(y, wait=1.6):
            frame.evaluate(f"window.scrollTo({{top: {y}, behavior: 'smooth'}})")
            wait_idle(20)
            time.sleep(wait)

        def scroll_to(loc, wait=1.6):
            loc.evaluate("el => el.scrollIntoView({behavior: 'smooth', block: 'start'})")
            time.sleep(wait)

        def tab(value):
            caption("")
            click(app.locator(f".navbar a[data-value='{value}']"), wait=1.5)
            scroll(0, 0.5)

        def subtab(navset, value, wait=3.0):
            caption("")
            click(app.locator(f"#{navset} a[data-value='{value}']"), wait=wait)

        # ---- Home -----------------------------------------------------------
        chapter("Welcome and home page")
        caption("MagnaGenes v2.0: a multi-omics explorer for <i>Magnaporthe oryzae</i>", 3.5)
        caption("Runs entirely in your browser: no installation and no server", 3)
        scroll(520)
        caption("Each card opens a tool: gene search, transcriptomics, phosphoproteomics, phylogenomics and GO", 4)
        scroll(1400)
        caption("Information, contact details and how to contribute data", 3.5)
        scroll(0, 1)

        # ---- Gene search ----------------------------------------------------
        tab("search")
        chapter("Gene Search")
        caption("Gene Search: type MGG IDs, gene names or keywords", 2)
        box = app.locator("#s_query")
        box.click()
        box.press_sequentially("PMK1, MST12\nmelanin biosynthesis", delay=70)
        time.sleep(0.6)
        caption("IDs and names are matched exactly; other words are searched as keywords", 1)
        click(app.locator("#s_go"), wait=3)
        caption("Results: annotation, phenotypes and links to UniProt, AlphaFold, Ensembl and InterPro", 4)
        scroll(400, 2)
        scroll(0, 0.8)
        subtab("s_view", "profile", wait=3)
        caption("Gene profile: annotation, expression in every experiment and protein architecture", 3)
        scroll(700, 2.5)
        scroll(1300, 2.5)
        caption("Protein and CDS FASTA downloads are in the sidebar", 2.5)

        # ---- Transcriptomics --------------------------------------------------
        tab("transcriptomics")
        chapter("Transcriptomics")
        caption("Transcriptomics: paste genes or load the example list", 1.5)
        click(app.locator("#tx_example"), wait=0.8)
        caption("Choose experiments; each links to its source paper", 2.5)
        click(app.locator("#tx_go"), wait=4)
        caption("One plot per experiment: TPM over time for each gene", 3.5)
        click(app.locator("#tx_box"), wait=3)
        caption("Tick the box plot option to add the spread of all selected genes", 3)
        scroll(650, 3)
        scroll(0, 0.8)
        subtab("tx_view", "heatmap", wait=3)
        caption("Heatmaps show row z-scores: when each gene is high (red) or low (blue)", 3.5)
        click(app.locator("#tx_heat_values"), wait=3)
        caption("Show the TPM values in the cells if needed", 3)
        subtab("tx_view", "corr", wait=3.5)
        caption("Gene-gene correlation; the colour scale fits the genes you selected", 4)
        subtab("tx_view", "table", wait=2)
        caption("Data table and TSV, PDF and FASTA downloads", 3)

        # ---- Phosphoproteomics -------------------------------------------------
        tab("phospho")
        chapter("Phosphoproteomics")
        caption("Phosphoproteomics: choose a protein or load the example (MST7)", 1.5)
        click(app.locator("#ph_example"), wait=5)
        caption("One plot per phosphopeptide: Guy11 (green) above Δpmk1 (grey), bars are replicates", 4)
        scroll(600, 3)
        scroll(0, 0.8)
        subtab("ph_view", "heatmap", wait=3)
        caption("Abundance heatmap of all peptides of this protein", 3)
        subtab("ph_view", "fc", wait=3)
        caption("Fold change in Δpmk1 against Guy11; dots mark significant changes", 3.5)
        subtab("ph_view", "map", wait=3.5)
        caption("Protein site map: labelled sites above the protein and its Pfam domains", 3.5)
        scroll(420, 3)
        scroll(0, 0.8)
        subtab("ph_view", "table", wait=2)
        caption("Welch t-test statistics for every site and time point", 3)

        # ---- Phylogenomics --------------------------------------------------------
        tab("phylo")
        chapter("Phylogenomics")
        caption("Phylogenomics: ortholog counts across 42 fungi", 1.5)
        click(app.locator("#og_example"), wait=0.6)
        click(app.locator("#og_go"), wait=5)
        caption("Cells take the lifestyle colour, deeper with more copies; strips show lifestyle and appressorium", 4.5)
        scroll(600, 3)
        scroll(0, 0.8)
        subtab("og_view", "lifestyle", wait=3)
        caption("Conservation by lifestyle: share of species with an ortholog", 3.5)
        subtab("og_view", "members", wait=2)
        caption("<i>M. oryzae</i> orthogroup members with database links", 3)

        # ---- GO enrichment ----------------------------------------------------------
        tab("go")
        chapter("GO Enrichment")
        caption("GO enrichment: paste a gene list or load an example module", 1.5)
        click(app.locator("#go_example"), wait=0.8)
        click(app.locator("#go_run"), wait=5)
        caption("Dot plot of enriched terms (Fisher's exact test)", 3)
        subtab("go_view", "bar", wait=2.5)
        caption("Bar chart, in a colour of your choice", 2.5)
        subtab("go_view", "treemap", wait=3)
        caption("Treemap of significant terms", 2.5)
        subtab("go_view", "volcano", wait=2.5)
        caption("Significance against enrichment for all tested terms", 2.5)
        subtab("go_view", "matrix", wait=2.5)
        caption("Gene-term matrix: which genes drive each term", 3)

        # ---- Theme, reset, about ----------------------------------------------------------
        chapter("Dark mode, Reset and About")
        caption("Switch between light and dark themes", 1)
        click(app.locator("bslib-input-dark-mode"), wait=2.5)
        click(app.locator("bslib-input-dark-mode"), wait=1.5)
        caption("Reset clears every search and starts again", 3)
        tab("about")
        caption("About: credits, version history, data sources and methods", 3)
        scroll(700, 3)
        caption("Thank you for using MagnaGenes", 3)
        caption("", 0.5)
        duration = now()
        video = page.video.path()
        ctx.close()

    ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
    OUT_VIDEO.parent.mkdir(exist_ok=True)
    subprocess.run([ffmpeg, "-y", "-loglevel", "error", "-ss", f"{t0:.2f}", "-i", str(video),
                    "-t", f"{duration:.2f}", "-c:v", "libx264", "-preset", "slow", "-crf", "28",
                    "-pix_fmt", "yuv420p", "-movflags", "+faststart", "-an", str(OUT_VIDEO)], check=True)
    OUT_CHAPTERS.write_text(json.dumps(chapters, indent=2))
    shutil.rmtree(workdir, ignore_errors=True)
    print(f"Wrote {OUT_VIDEO} ({OUT_VIDEO.stat().st_size / 1e6:.1f} MB, {duration:.0f} s)")
    for c in chapters:
        print(f"  {int(c['time'] // 60)}:{int(c['time'] % 60):02d}  {c['title']}")


if __name__ == "__main__":
    main()
