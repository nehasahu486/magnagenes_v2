"""
Retake the README screenshots (docs/screenshots/01..11_*.png).

Same requirements as record_tutorial.py. From the repository root, after
`Rscript build_site.R` and with site/ served on 127.0.0.1:8765:
    python tools/take_screenshots.py
"""
import tempfile
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

OUT = Path(__file__).resolve().parents[1] / "docs" / "screenshots"
URL = "http://127.0.0.1:8765/"
VIEW = {"width": 1440, "height": 900}


def main():
    with sync_playwright() as p:
        ctx = p.chromium.launch_persistent_context(
            tempfile.mkdtemp(prefix="mg_shots_"), channel="chrome", headless=True,
            args=["--no-proxy-server"], viewport=VIEW)
        page = ctx.new_page()
        page.goto(URL)
        app = page.frame_locator("#root iframe")
        app.locator(".navbar").wait_for(timeout=240000)
        page.wait_for_selector("#mg-loader", state="detached", timeout=60000)
        frame = [f for f in page.frames if f != page.main_frame][0]

        def idle(extra=1.5):
            time.sleep(0.5)
            quiet, end = 0.0, time.time() + 40
            while time.time() < end and quiet < 1.2:
                busy = frame.evaluate("document.documentElement.classList.contains('shiny-busy')")
                quiet = 0.0 if busy else quiet + 0.2
                time.sleep(0.2)
            time.sleep(extra)

        def top():
            frame.evaluate("window.scrollTo(0, 0)")
            time.sleep(0.4)

        def shot(name):
            page.screenshot(path=str(OUT / name))
            print("wrote", name)

        def tab(value):
            app.locator(f".navbar a[data-value='{value}']").click()
            idle(1)
            top()

        def sub(navset, value, extra=2.5):
            app.locator(f"#{navset} a[data-value='{value}']").click()
            idle(extra)

        # 01 Home
        idle(1)
        shot("01_home.png")

        # 02 Gene search with the curation columns and two ticked rows
        tab("search")
        app.locator("#s_query").fill("PMK1, MST12\nmelanin")
        app.locator("#s_go").click()
        idle(2)
        boxes = app.locator("#s_table .rt-tbody input[type=checkbox]")
        boxes.nth(0).click(); boxes.nth(1).click()
        idle(1)
        app.locator(".bslib-sidebar-layout .collapse-toggle").first.click()
        time.sleep(1.5)
        top()
        shot("02_gene_search.png")
        app.locator(".bslib-sidebar-layout .collapse-toggle").first.click()
        time.sleep(1)

        # 03 Gene profile (PMK1, first result)
        sub("s_view", "profile", 3)
        top()
        shot("03_gene_profile.png")

        # 04-06 Transcriptomics example
        tab("transcriptomics")
        app.locator("#tx_example").click(); idle(1)
        app.locator("#tx_go").click(); idle(4)
        top()
        shot("04_tx_lineplots.png")
        sub("tx_view", "heatmap", 3); top()
        shot("05_tx_heatmap.png")
        sub("tx_view", "corr", 4); top()
        shot("06_tx_correlation.png")

        # 07-08 Phosphoproteomics example (MST7)
        tab("phospho")
        app.locator("#ph_example").click(); idle(5)
        top()
        shot("07_phospho_peptides.png")
        sub("ph_view", "map", 4); top()
        shot("08_phospho_sitemap.png")

        # 09 Phylogenomics example
        tab("phylo")
        app.locator("#og_example").click(); idle(1)
        app.locator("#og_go").click(); idle(5)
        top()
        shot("09_phylogenomics.png")

        # 10-11 GO enrichment, WGCNA module 7
        tab("go")
        # Module 7 is the default example module (EX_GO_MODULE)
        app.locator("#go_example").click(); idle(1)
        app.locator("#go_run").click(); idle(5)
        top()
        shot("10_go_dotplot.png")
        sub("go_view", "matrix", 3); top()
        shot("11_go_matrix.png")
        ctx.close()


if __name__ == "__main__":
    main()
