# TaK Explorer

**TaK Explorer** is an R Shiny app for diagnosing the taxonomic quality of
a biodiversity dataset — how deep the identifications go (species? genus?
only phylum?) and how much of the data that depth actually covers. It
implements the **TaK (Taxonomic Knowledge) index**, described in
*"TaK – a novel index for assessing taxonomic knowledge in biodiversity
datasets"* (Fonseca et al., submitted to *Ecological Indicators*).

No installation beyond R itself is required — the app installs its own
dependencies on first run (see [Quick start](#quick-start)).

## The idea in one paragraph

A biodiversity dataset is rarely identified uniformly: some records reach
species level, others stop at phylum. TaK Explorer turns that unevenness
into two numbers per group (a dataset, a phylum, whatever column you
choose to compare by):

- **TR — Taxonomic Resolution**: *of the taxa you've named, how deep did
  you go?* It's the average identification depth across **unique
  lineages** (one lineage = one distinct path through the hierarchy, e.g.
  `Arthropoda > Copepoda > Calanoida`), regardless of how many individuals
  each one represents.
- **TC — Taxonomic Completeness**: *of your actual specimens/records, how
  many are covered by a good identification?* It's the same
  identification-depth score, but weighted by abundance — so a single
  coarsely-identified but hugely abundant lineage can pull TC down a lot
  even when TR looks fine.

Because TR and TC answer different questions, a dataset can score high on
one and low on the other. Plotting TR against TC in a biplot (the app's
main output) sorts groups into four diagnostic quadrants:

| | TC < 0.5 | TC ≥ 0.5 |
|---|---|---|
| **TR ≥ 0.5** | Under-resolved (deep IDs, but only for rare taxa) | Well-resolved |
| **TR < 0.5** | Data-deficient | Abundance-biased (a few dominant taxa carry the completeness score) |

## Quick start

```r
install.packages('shiny')
shiny::runGitHub('Tak_Explorer', 'gfonseca-unifesp', ref = 'main')
```

The app opens with a small built-in benchmarking dataset (five toy
sub-datasets designed to land in each of the four quadrants above), so you
can explore the interface before uploading your own data.

## How the app is structured

TaK Explorer is a single `page_navbar` app with four tabs:

1. **Data Editor** — where your data lives. Upload a CSV, edit any cell
   directly in the table (double-click), add extra taxonomic-rank columns,
   set the weight vector (see below), and export what you have at any
   point. Everything downstream reacts to this table.
2. **Visualization** — the Biplot (TR vs. TC, one point per group, sized
   by abundance) and a stacked bar chart showing what fraction of each
   Dataset falls in each of the four quadrants. Both are downloadable as
   publication-ready PNGs.
3. **Summary** — one row per Dataset: mean TR/TC, total individuals, and
   the count of unique taxa — a quick numeric readout to accompany the
   plots, downloadable as CSV.
4. **Sampling Sensitivity** — a rarefaction/bootstrap check of *how much
   TR and TC depend on how much data you have*. See below — this is worth
   reading before you trust a TC value from a small or patchy dataset.

### Why there's a "Sampling Sensitivity" tab

TR is an **unweighted** mean over unique lineages, so — mathematically —
it shouldn't matter much how many individuals or records back each
identification. TC, being **abundance-weighted**, can behave very
differently: if a dataset's abundance is dominated by one or two lineages
that only got identified to a coarse rank, TC can swing hard depending on
whether your sample happened to capture those dominant records.

This isn't hypothetical — it's the reason the tab exists. Testing it
against real ISA DeepData occurrence records, one single record (out of
2,155) accounted for 72% of the Clarion-Clipperton Zone's total measured
abundance and was identified only to Phylum; whether that one record was
included in a random subsample moved the region's TC from **0.17 to
0.50** — a three-fold difference from a single line of data — while TR for
the same region barely moved at all (0.752 → 0.752) across the same
subsampling range.

**How to use it:** upload your data, go to the Sampling Sensitivity tab,
set how many bootstrap draws and subsampling levels you want (defaults are
reasonable for exploration), and click "Run rarefaction analysis" (it's a
deliberate button, not a live-updating chart — a full bootstrap sweep
takes a few seconds and shouldn't re-run on every keystroke). You get a
TR/TC-vs-sampling-depth curve per Dataset, plus the underlying numbers as
a downloadable CSV. A flat TR curve and a flat TC curve mean your indices
are trustworthy regardless of how much more data you might collect; a TC
curve that's still moving as you approach 100% of your data is a sign
that your completeness score is sensitive to a small number of records,
and should be reported with that caveat.

## Weight vector

Each taxonomic rank column (Phylum, Class, Order, ...) gets a weight, and
those weights set how much finer identifications are rewarded when
computing TR/TC. The default (`1, 2, 3, 4, 5, 6`) is a **linear** scheme —
each additional rank of depth is worth one more point. You can type any
comma-separated vector of the same length as your rank columns:

- **Uniform** (e.g. `1, 1, 1, 1, 1, 1`): every rank counts equally: only
  *how many* ranks got identified matters, not which ones.
- **Linear** (the default): deeper ranks are worth proportionally more.
- **Custom / exponential**: type any increasing sequence (e.g.
  `1, 2, 4, 8, 16, 32`) to weight the deepest ranks much more heavily than
  the shallow ones — useful when species-level identification is
  disproportionately harder to obtain than the ranks above it.

There's no dropdown for named presets — the weight field is free text on
purpose, so it can match datasets with any number of taxonomic-rank
columns (the app adapts to however many rank columns your data has, via
"Add taxonomic column").

## Data format

Upload a CSV (comma- or semicolon-delimited, auto-detected) with:

- One column per taxonomic rank, in order from coarsest to finest (e.g.
  `Phylum, Class, Order, Family, Genus, Species`). Leave a cell blank (or
  `NA`) if that record wasn't identified that deep — don't guess or
  backfill.
- A `Dataset` column, whatever grouping you want to compare by (a region,
  a station, a year, a survey method). If missing, everything is treated
  as one Dataset called `"Uploaded"`.
- An abundance/count column, auto-detected under any of these names:
  `individualCount`, `n`, `Abundance`, `abundance`, `count`. If none of
  these is found, every record is assumed to represent 1 individual.

`example_data.csv` and `ISA_DeepData_2026.csv` (a real 4-region occurrence
dataset from the International Seabed Authority's DeepData/OBIS node,
used in the accompanying manuscript) are both in this repo as worked
examples you can upload directly.

## Repository structure

```
app.R              The Shiny app (UI + server). Sources TaK_fun_1.R for
                    all the actual math -- it does not redefine it.
TaK_fun_1.R         Core R functions: calculate_TaK_shiny() (the TR/TC
                    engine), rarefy_tak() (the bootstrap behind the
                    Sampling Sensitivity tab), and draw_rarefaction_plot()
                    (shared by the in-app plot and its PNG download).
TaK_fun_2.py        Standalone Python re-implementation of the same TR/TC
                    math (function TaK()), for users who work in Python
                    instead of R. Not used by the Shiny app.
Vignette_app.Rmd    A runnable R Markdown walkthrough of the same app
                    (useful as a teaching/demo document).
example_data.csv    The benchmarking dataset (five scenarios engineered
                    to land in each biplot quadrant).
ISA_DeepData_2026.csv
                    Real occurrence data (CCZ, Indian Ocean, NW Pacific,
                    Mid-Atlantic Ridge) used as the empirical case study
                    in the manuscript.
```

## Citing TaK

If you use TaK or TaK Explorer, please cite:

> Fonseca, G., Paula, F.S., Ávila, A.K.F., Brustolin, M.C., Watanabe, W.B.,
> Pape, E., Vandepitte, L., Genio, L. TaK – a novel index for assessing
> taxonomic knowledge in biodiversity datasets. *Ecological Indicators*
> (submitted).

## License

MIT — see [LICENSE](LICENSE).
