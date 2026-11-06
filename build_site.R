# =============================================================================
# Build the static MagnaGenes site with Shinylive.
#
#   Rscript build_site.R            # writes site/
#
# Adds a loading screen (M. oryzae infection cycle, progress bar and a short
# note) that stays up while webR and the R packages download, then fades out
# once the app has drawn its navigation bar.
# =============================================================================

out_dir <- "site"
unlink(out_dir, recursive = TRUE)

loader_head <- '
<style>
  #mg-loader { position: fixed; inset: 0; z-index: 9999; background: #ffffff;
    display: flex; align-items: center; justify-content: center;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
    transition: opacity .6s ease; }
  #mg-loader.done { opacity: 0; pointer-events: none; }
  #mg-loader .mg-cycle { position: absolute; inset: 0; background: url("loading_cycle.png") center / contain no-repeat;
    opacity: .9; margin: 4vh 4vw; animation: mg-fade 1.2s ease both; }
  #mg-loader .mg-card { position: relative; width: min(340px, 78vw); text-align: center;
    background: rgba(255,255,255,.92); border-radius: 14px; padding: 1.1rem 1.3rem;
    box-shadow: 0 4px 18px rgba(0,0,0,.08); }
  #mg-loader h1 { margin: 0 0 .15rem; font-size: 1.55rem; font-weight: 800; color: #1d2b2a; letter-spacing: -.01em; }
  #mg-loader h1 span { font-weight: 400; color: #8a8984; }
  #mg-loader .mg-sub { font-size: .86rem; color: #52514e; margin-bottom: .9rem; }
  #mg-loader .mg-track { height: 8px; border-radius: 999px; background: #e3efe9; overflow: hidden; }
  #mg-loader .mg-bar { height: 100%; width: 0%; border-radius: 999px; background: #1B9E77; transition: width .4s ease; }
  #mg-loader .mg-status { margin-top: .55rem; font-size: .8rem; color: #137a5c; font-weight: 600; min-height: 1.1em; }
  #mg-loader .mg-note { margin-top: .5rem; font-size: .76rem; color: #8a8984; line-height: 1.4; }
  @keyframes mg-fade { from { opacity: 0; transform: scale(.97); } to { opacity: .9; transform: none; } }
  @media (prefers-reduced-motion: reduce) { #mg-loader .mg-cycle { animation: none; } #mg-loader .mg-bar { transition: none; } }
</style>'

loader_body <- '
<div id="mg-loader" role="status" aria-live="polite">
  <div class="mg-cycle" aria-hidden="true"></div>
  <div class="mg-card">
    <h1>MagnaGenes <span>v2.0</span></h1>
    <div class="mg-sub">Multi-omics explorer for <em>Magnaporthe oryzae</em></div>
    <div class="mg-track"><div class="mg-bar" id="mg-bar"></div></div>
    <div class="mg-status" id="mg-status">Starting R in your browser…</div>
    <div class="mg-note">The first visit takes less than a minute while the R packages load.
      After that MagnaGenes opens quickly.</div>
  </div>
</div>
<script>
(function () {
  var start = Date.now(), loader = document.getElementById("mg-loader"),
      bar = document.getElementById("mg-bar"), status = document.getElementById("mg-status");
  var steps = [[0, "Starting R in your browser…"], [6, "Downloading R packages…"],
               [22, "Loading MagnaGenes data…"], [38, "Almost there…"]];
  function appReady() {
    var f = document.querySelector("#root iframe");
    try { var d = f && f.contentDocument; return !!(d && d.querySelector(".navbar")); } catch (e) { return false; }
  }
  var timer = setInterval(function () {
    var t = (Date.now() - start) / 1000;
    bar.style.width = (95 * (1 - Math.exp(-t / 20))).toFixed(1) + "%";
    for (var i = steps.length - 1; i >= 0; i--) { if (t >= steps[i][0]) { status.textContent = steps[i][1]; break; } }
    if (appReady()) {
      clearInterval(timer);
      bar.style.width = "100%"; status.textContent = "Ready";
      setTimeout(function () { loader.classList.add("done"); setTimeout(function () { loader.remove(); }, 700); }, 350);
    }
  }, 250);
})();
</script>'

# The loader text contains … escapes meant for JavaScript; keep them literal
shinylive::export(
  "app", out_dir,
  template_params = list(
    title = "MagnaGenes v2.0",
    include_in_head = loader_head,
    include_before_body = loader_body
  )
)

file.copy("site-assets/loading_cycle.png", file.path(out_dir, "loading_cycle.png"), overwrite = TRUE)
message("Site written to ", out_dir, "/")
