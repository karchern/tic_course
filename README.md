# TIC screen analysis

Analysis of a Toxin Inhibition Conjugation (TIC) screen: *E. coli* colonies each
carry the RcaT toxin (arabinose-inducible) and a designed blocker from the
library (IPTG-inducible). A blocker works if its colonies still grow when both
are switched on.

Two screens are in here, with the same plate layout:

| | Iris measurements | plate photographs | colony table |
|---|---|---|---|
| original | `data/iris` | `data/photos` | `data/colonies.tsv` |
| predoc course | `data/iris_predocs` | `data/plate_images_predocs` | `data/colonies_predocs.tsv` |

## Running it

The environment is in `envs/r_tic.yml` (conda-forge and bioconda only):

```bash
conda env create -f envs/r_tic.yml
conda activate r_tic
```

Then, from the repository root:

```bash
Rscript preprocess.R     # reads data/iris, writes data/colonies.tsv
Rscript tic_analysis.R   # sanity checks, effect sizes, t-tests, volcano plot
```

Output lands in `output/`, which is not tracked.

## The scripts

- **`preprocess.R`** — reads every `.iris` file, attaches the plate map so each
  colony knows its strain, decodes the plate name into the experimental
  conditions, and writes one tidy table of colonies.
- **`tic_analysis.R`** — sanity checks on the control plates, plate-edge
  correction, normalisation to the GFP controls, effect size per blocker,
  t-tests, volcano plot.
- **`inspect_gene.R`** — `inspect_gene("50aa_1_55")` shows one gene's colonies:
  the Iris numbers on the left, the colonies cut out of the plate photographs
  on the right. Pass `"iris_predocs"` as the second argument for the other
  screen.
- **`plate_photos.R`** — finds the colony grid in a plate photograph so the
  colonies can be cut out. Every grid is checked against the Iris measurements
  of the same plate and refused below r = 0.85; the fits are cached in
  `data/grid_fits`.

## The data

One 96-well source plate is pinned into 384 format, so every well becomes a
2x2 quartet: two clones, each pinned twice. Each plate is photographed under
four conditions — no inducer, IPTG (blocker only), Ara (toxin only), and
Ara + IPTG — and Iris turns each photograph into one opacity value per colony.

`data/AlphaDesign_384wp.identifiers.csv` says which strain sits at each of the
384 positions. `data/HitList_Paper.txt` lists the hits reported in the paper.
