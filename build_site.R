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
  #mg-loader .mg-stage { position: relative; width: min(520px, 86vw); aspect-ratio: 1400 / 991;
    display: flex; align-items: center; justify-content: center; }
  /* The infection cycle: faint outline, drawn in clockwise as loading progresses */
  #mg-loader .mg-cycle { position: absolute; inset: 0; background: url("loading_cycle.png") center / contain no-repeat; }
  #mg-loader .mg-cycle.base { opacity: .08; }
  #mg-loader .mg-cycle.reveal { opacity: .42;
    -webkit-mask-image: conic-gradient(from 0deg at 50% 50%, #000 var(--p, 0%), transparent var(--p, 0%));
            mask-image: conic-gradient(from 0deg at 50% 50%, #000 var(--p, 0%), transparent var(--p, 0%)); }
  #mg-loader .mg-card { position: relative; width: min(250px, 60vw); text-align: center;
    background: rgba(255,255,255,.78); border-radius: 12px; padding: .8rem 1rem; }
  #mg-loader h1 { margin: 0 0 .1rem; font-size: 1.35rem; font-weight: 800; color: #1d2b2a; letter-spacing: -.01em; }
  #mg-loader h1 span { font-weight: 400; color: #8a8984; }
  #mg-loader .mg-sub { font-size: .78rem; color: #52514e; margin-bottom: .7rem; }
  #mg-loader .mg-track { height: 7px; border-radius: 999px; background: #e3efe9; overflow: hidden; }
  #mg-loader .mg-bar { height: 100%; width: 0%; border-radius: 999px; background: #1B9E77; transition: width .4s ease; }
  #mg-loader .mg-status { margin-top: .45rem; font-size: .76rem; color: #137a5c; font-weight: 600; min-height: 1.1em; }
  #mg-loader .mg-note { margin-top: .4rem; font-size: .7rem; color: #8a8984; line-height: 1.35; }
  @media (prefers-color-scheme: dark) {
    #mg-loader { background: #121716; }
    #mg-loader .mg-card { background: rgba(27,33,32,.8); }
    #mg-loader h1 { color: #e6ebe9; } #mg-loader .mg-sub { color: #b4bdb9; } #mg-loader .mg-track { background: #23332d; }
    #mg-loader .mg-status { color: #5fd3ad; } #mg-loader .mg-cycle.reveal { opacity: .55; filter: invert(1) hue-rotate(180deg); }
    #mg-loader .mg-cycle.base { filter: invert(1) hue-rotate(180deg); }
  }
</style>'

loader_body <- '
<div id="mg-loader" role="status" aria-live="polite">
  <div class="mg-stage">
    <div class="mg-cycle base" aria-hidden="true"></div>
    <div class="mg-cycle reveal" id="mg-cycle" aria-hidden="true"></div>
    <div class="mg-card">
      <h1>MagnaGenes <span>v2.0</span></h1>
      <div class="mg-sub">Multi-omics explorer for <em>Magnaporthe oryzae</em></div>
      <div class="mg-track"><div class="mg-bar" id="mg-bar"></div></div>
      <div class="mg-status" id="mg-status">Starting R in your browser…</div>
      <div class="mg-note">The first visit takes less than a minute while the R packages load.
        After that MagnaGenes opens quickly.</div>
    </div>
  </div>
</div>
<script>
(function () {
  var start = Date.now(), loader = document.getElementById("mg-loader"),
      bar = document.getElementById("mg-bar"), status = document.getElementById("mg-status"),
      cycle = document.getElementById("mg-cycle");
  var steps = [[0, "Starting R in your browser…"], [6, "Downloading R packages…"],
               [22, "Loading MagnaGenes data…"], [38, "Almost there…"]];
  function appReady() {
    var f = document.querySelector("#root iframe");
    try { var d = f && f.contentDocument; return !!(d && d.querySelector(".navbar")); } catch (e) { return false; }
  }
  function show(p) { bar.style.width = p + "%"; cycle.style.setProperty("--p", p + "%"); }
  var timer = setInterval(function () {
    var t = (Date.now() - start) / 1000;
    show((95 * (1 - Math.exp(-t / 20))).toFixed(1));
    for (var i = steps.length - 1; i >= 0; i--) { if (t >= steps[i][0]) { status.textContent = steps[i][1]; break; } }
    if (appReady()) {
      clearInterval(timer);
      show(100); status.textContent = "Ready";
      setTimeout(function () { loader.classList.add("done"); setTimeout(function () { loader.remove(); }, 700); }, 450);
    }
  }, 250);
})();
</script>'

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
