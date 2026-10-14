# Review improvements

## earth-stats-cheatsheet (sections 01-11), 2026-10-12

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

## mining-stats-cheatsheet (sections 12-21), 2026-10-12

The block estimate remains the average of four point estimates, so its predictions
and comparisons with exhaustive block truth are unchanged by this correction. The
uncertainty panel formerly averaged point variances and labelled their square root as
kriging SD without distinguishing support. It now says **RMS point kriging SD (not block
SE)**. This is an interim descriptive panel, not uncertainty of the block estimate:
a block SE requires the cross-point prediction-error covariances or a direct block
kriging solution. The full numerical rerun uses the Mac mini, with source/version/input
hashes recorded in `reports/results.json` and `reports/data-manifest.json`.

Data loaders now clip time-series dates at the configured cutoff, record request metadata
and hashes, and invalidate parsed caches after loader/cutoff changes. Known cache bytes
are integrity-checked before use. Legacy raw snapshots are retained where possible,
with `retrieved_at: null` and `origin: legacy-cache`; the loader's request is recorded but
is not claimed to reconstruct its original acquisition. GitHub teaching data are newly
fetched from a pinned upstream revision. Provider revisions to historical observations
still require retaining the exact raw bytes for reproducibility.

`make results` executes every topic, preserves failed transcripts, and publishes only
a complete successful set. Numerical excerpts and full transcripts replace the hand-edited
README talking points. The method chooser links to code, and the portfolio shows selected
figures copied from the successful run, with source/input/output hashes in the manifest.
