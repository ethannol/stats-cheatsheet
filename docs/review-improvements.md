# Review improvements, 2026-10-12

The LOO example formerly refitted the regression but reused variograms fitted using
all station temperatures. Each fold now refits the trend and both raw-temperature and
residual variograms using only training stations. The small held-out-target perturbation
fixture verifies that held-out responses cannot change those models. LOO remains an
interpolation assessment; it does not establish transfer to an unsampled region. The
full numerical rerun uses the Mac mini, with source/version/input hashes recorded in
`reports/results.json` and `reports/data-manifest.json`. Old headline RMSEs were removed
rather than retained before that rerun.

Data loaders now clip time-series dates at the configured cutoff, record request metadata
and hashes, and invalidate parsed caches after loader/cutoff changes. Known cache bytes
are integrity-checked before use. Legacy raw snapshots are retained where possible,
with `retrieved_at: null` and `origin: legacy-cache`; the loader's request is recorded but
is not claimed to reconstruct its original acquisition. Provider revisions to historical observations still require retaining the exact raw
bytes for reproducibility.

`make results` executes every topic, preserves failed transcripts, and publishes only
a complete successful set. Numerical excerpts and full transcripts replace the hand-edited
README talking points. The method chooser links to code, and the portfolio shows selected
figures copied from the successful run, with source/input/output hashes in the manifest.
