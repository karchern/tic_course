# Cut the raw colony images out of the plate photographs.
#
# Iris turned each plate photo into one opacity number per colony and threw the
# picture away, so to look at the colonies themselves we have to find the
# colony grid in the photo again. Iris's own grid coordinates are no help: it
# worked on a cropped and resized copy of the photo.
#
# So the grid is found from scratch, in four steps:
#   1. the plate is located as the bright area in the photo;
#   2. the spacing between colonies comes from the periodicity of the plate's
#      own brightness profile;
#   3. the position of the grid comes from lining the profile up with the
#      opacities Iris measured, column by column and row by row;
#   4. the whole grid is then scored against all 384 Iris measurements and
#      polished. Below r = 0.85 we refuse to cut anything out, so a grid that
#      is off by a colony cannot pass unnoticed.
# The answer is cached in data/grid_fits as <source>__<plate>.grid.csv.
#
# There are two screens, each with its own Iris output and its own plate
# photographs, named slightly differently. plate_source() says where a screen
# keeps its files; everything else here takes one.
library(jpeg)

grid_dir <- "data/grid_fits"

plate_sources <- list(
  iris = list(
    colonies = "data/colonies.tsv",
    photos   = "data/photos",
    iris     = "data/iris",
    # data/photos/SpecGentAraIPTG-10050201-3A_1.jpg
    photo_of = function(plate) paste0(plate, ".jpg")),
  iris_predocs = list(
    colonies = "data/colonies_predocs.tsv",
    photos   = "data/plate_images_predocs",
    iris     = "data/iris_predocs",
    photo_of = function(plate) paste0(plate, ".jpg"))
)

plate_source <- function(source = "iris") {
  if (is.list(source)) return(source)
  if (!source %in% names(plate_sources))
    stop(source, " is not a screen - use one of: ", paste(names(plate_sources), collapse = ", "))
  c(plate_sources[[source]], name = source)
}

plates_with_photos <- function(src = plate_source()) {
  photos <- list.files(src$photos, "\\.jpe?g$", ignore.case = TRUE)
  photos <- photos[!str_detect(photos, "\\.grid\\.")]   # Iris's own overlays
  plates <- str_remove(basename(list.files(src$iris, "\\.iris$")), "\\.jpg\\.iris$")
  plates[map_chr(plates, src$photo_of) %in% photos]
}

plate_photo_path <- function(plate, src = plate_source()) file.path(src$photos, src$photo_of(plate))
plate_iris_path  <- function(plate, src = plate_source()) file.path(src$iris, paste0(plate, ".jpg.iris"))
grid_cache_path  <- function(plate, src = plate_source()) {
  file.path(grid_dir, paste0(src$name, "__", plate, ".grid.csv"))
}

# Searching the grid reads thousands of little patches, so it is done on a
# thinned copy of the photo. Colonies are over a hundred pixels across, so
# nothing that matters is lost; the answer is scaled back up at the end.
read_grey <- function(path, target_px = 2.2e6) {
  full <- readJPEG(path)
  if (length(dim(full)) == 3) full <- (full[, , 1] + full[, , 2] + full[, , 3]) / 3
  step <- max(1, floor(sqrt(length(full) / target_px)))
  list(grey = full[seq(1, nrow(full), step), seq(1, ncol(full), step)], step = step)
}

# the bright rectangle in the middle of the photo: the plate
plate_box <- function(profile, shrink = 0.02) {
  edges <- range(which(profile > (min(profile) + max(profile)) / 2))
  round(c(edges[1] + diff(edges) * shrink, edges[2] - diff(edges) * shrink))
}

# flatten the slow brightness gradient, leaving the colonies
detrend <- function(profile, window) {
  smooth <- stats::filter(profile, rep(1 / (2 * window + 1), 2 * window + 1), sides = 2)
  profile - ifelse(is.na(smooth), median(profile), smooth)
}


# Grids already fitted for other plates, in this photo's own coordinates: the
# median of the good ones, and each of them on its own.
rig_grids <- function(plate, step, src = plate_source(), min_agreement = 0.95) {
  files <- setdiff(list.files(grid_dir, paste0("^", src$name, "__.*\\.grid\\.csv$"),
                              full.names = TRUE),
                   grid_cache_path(plate, src))
  if (length(files) == 0) return(list())
  fitted <- map(files, read.csv) %>% list_rbind() %>% filter(agreement >= min_agreement)
  if (nrow(fitted) == 0) return(list())
  shrink <- function(g) list(x0 = (g$x0 - 1) / step + 1, px = g$px / step,
                             y0 = (g$y0 - 1) / step + 1, py = g$py / step)
  c(list(shrink(summarise(fitted, across(c(x0, px, y0, py), median)))),
    map(seq_len(nrow(fitted)), \(i) shrink(fitted[i, ])))
}


# Fitting the colony grid ----------------------------------------------------

fit_colony_grid <- function(plate, src = plate_source(), quiet = FALSE) {
  src   <- plate_source(src)
  cache <- grid_cache_path(plate, src)
  if (file.exists(cache)) return(as.list(read.csv(cache)))

  img  <- read_grey(plate_photo_path(plate, src))
  grey <- img$grey
  iris <- read.delim(plate_iris_path(plate, src), skip = 6)

  cols <- plate_box(colMeans(grey))
  rows <- plate_box(rowMeans(grey))
  inside <- grey[rows[1]:rows[2], cols[1]:cols[2]]
  col_profile <- detrend(colMeans(inside), 60)
  row_profile <- detrend(rowMeans(inside), 60)

  # Iris's mean opacity per column (and per row) is a fingerprint of the plate:
  # dim columns are dim in the photo too. So we try every plausible spacing and
  # starting point, and keep the few that reproduce that fingerprint best.
  col_opacity <- tapply(iris$opacity, iris$column, mean)
  row_opacity <- tapply(iris$opacity, iris$row, mean)
  # On a plate where the toxin killed most colonies the fingerprint is faint,
  # so we keep several candidates rather than one, and make them genuinely
  # different grids instead of neighbours of the same peak.
  candidates <- function(profile, n, target, keep = 8, pitches = NULL) {
    span <- length(profile)
    if (is.null(pitches)) pitches <- seq(span / (n + 6), (span - 1) / (n - 1), by = 0.25)
    found <- NULL
    for (pitch in pitches) {
      starts <- seq(1, span - pitch * (n - 1), by = 0.5)
      fit <- vapply(starts, function(s) cor(profile[round(s + pitch * (0:(n - 1)))], target),
                    numeric(1))
      found <- rbind(found, data.frame(pitch = pitch, start = starts, fit = fit))
    }
    found <- found[order(found$fit, decreasing = TRUE), ]
    out <- found[0, ]
    for (k in seq_len(nrow(found))) {
      if (nrow(out) == keep) break
      far <- abs(out$start - found$start[k]) > 3 | abs(out$pitch - found$pitch[k]) > 1
      if (all(far)) out <- rbind(out, found[k, ])
    }
    out
  }
  # 24 columns give a far sharper fingerprint than 16 rows, so the columns
  # decide the spacing; the 384 format is square, so the rows reuse it and only
  # have to find where they start.
  col_fit <- candidates(col_profile, 24, col_opacity)
  row_fit <- candidates(row_profile, 16, row_opacity, keep = 3,
                        pitches = seq(min(col_fit$pitch) - 1, max(col_fit$pitch) + 1, by = 0.25))

  # Brightness of a predicted colony above its local agar background. This
  # mirrors what Iris calls opacity, so the two are directly comparable.
  patch <- function(x, y, d) {
    x <- round(x); y <- round(y)
    if (x - 2 * d < 1 || y - 2 * d < 1 ||
        x + 2 * d > ncol(grey) || y + 2 * d > nrow(grey)) return(NA_real_)
    mean(grey[(y - d):(y + d), (x - d):(x + d)]) -
      median(grey[(y - 2 * d):(y + 2 * d), (x - 2 * d):(x + 2 * d)])
  }
  score_grid <- function(g) {
    d <- round(g$px * 0.25)
    b <- mapply(function(r, cc) patch(g$x0 + g$px * (cc - 1), g$y0 + g$py * (r - 1), d),
                iris$row, iris$column)
    if (anyNA(b)) return(-1)
    cor(b, iris$opacity)
  }

  # the best candidate grid, judged on all 384 colonies at once
  best <- list(g = NULL, s = -Inf)
  for (i in seq_len(nrow(col_fit))) {
    for (j in seq_len(nrow(row_fit))) {
      if (abs(col_fit$pitch[i] - row_fit$pitch[j]) > 1.5) next
      g <- list(x0 = col_fit$start[i] + cols[1] - 1, px = col_fit$pitch[i],
                y0 = row_fit$start[j] + rows[1] - 1, py = row_fit$pitch[j])
      s <- score_grid(g)
      if (s > best$s) best <- list(g = g, s = s)
    }
  }
  # On a plate where the toxin suppressed nearly every colony there may be too
  # little left to recognise. The plates were photographed in the same rig, so
  # a grid that already fitted other plates is a fair guess; it still has to
  # win on the Iris agreement like any other candidate.
  for (g in rig_grids(plate, img$step, src)) {
    s <- score_grid(g)
    if (s > best$s) best <- list(g = g, s = s)
  }

  grid <- best$g
  if (is.null(grid)) stop(plate, ": could not place a colony grid in the photo")

  # polish each parameter in turn
  refine <- function(g, field, values) {
    best <- g; best_score <- score_grid(g)
    for (v in values) {
      cand <- g; cand[[field]] <- v
      s <- score_grid(cand)
      if (s > best_score) { best <- cand; best_score <- s }
    }
    best
  }
  # coarse first, in case the guess came from another plate's grid
  grid <- refine(grid, "x0", seq(grid$x0 - 25, grid$x0 + 25, by = 1))
  grid <- refine(grid, "y0", seq(grid$y0 - 25, grid$y0 + 25, by = 1))
  for (pass in 1:2) {
    grid <- refine(grid, "x0", seq(grid$x0 - 4, grid$x0 + 4, by = 0.25))
    grid <- refine(grid, "px", seq(grid$px - 0.3, grid$px + 0.3, by = 0.02))
    grid <- refine(grid, "y0", seq(grid$y0 - 4, grid$y0 + 4, by = 0.25))
    grid <- refine(grid, "py", seq(grid$py - 0.3, grid$py + 0.3, by = 0.02))
  }

  agreement <- score_grid(grid)
  if (!quiet) cat(sprintf("%s: colony grid found, r = %.3f with Iris opacity\n",
                          plate, agreement))
  if (agreement < 0.85)
    stop(plate, ": the grid does not match the Iris measurements (r = ",
         round(agreement, 2), ") - refusing to cut out colonies")

  # back to the coordinates of the full-size photo
  step <- img$step
  grid <- list(x0 = (grid$x0 - 1) * step + 1, px = grid$px * step,
               y0 = (grid$y0 - 1) * step + 1, py = grid$py * step,
               agreement = agreement)
  dir.create(grid_dir, showWarnings = FALSE)
  write.csv(as.data.frame(grid), cache, row.names = FALSE)
  grid
}


# Cutting out the colonies ---------------------------------------------------
# One 96-well plate is pinned into 384 format, so each well is a 2x2 quartet of
# colonies. quartet_pixels() returns the photo pixels of one quartet as a long
# table, sampled onto a fixed square grid with x and y running 0-1, so that
# quartets from differently sized photos can be drawn side by side.

quartet_pixels <- function(plate, well_row, well_col, src = plate_source(),
                           pad = 0.6, out_px = 220) {
  src  <- plate_source(src)
  grid <- fit_colony_grid(plate, src, quiet = TRUE)
  grey <- readJPEG(plate_photo_path(plate, src))
  if (length(dim(grey)) == 3) grey <- (grey[, , 1] + grey[, , 2] + grey[, , 3]) / 3

  centre <- function(i, first, pitch) first + pitch * (i - 1)
  x_lo <- centre(2 * well_col - 1, grid$x0, grid$px) - grid$px * pad
  x_hi <- centre(2 * well_col,     grid$x0, grid$px) + grid$px * pad
  y_lo <- centre(2 * well_row - 1, grid$y0, grid$py) - grid$py * pad
  y_hi <- centre(2 * well_row,     grid$y0, grid$py) + grid$py * pad

  xs <- pmin(ncol(grey), pmax(1, round(seq(x_lo, x_hi, length.out = out_px))))
  ys <- pmin(nrow(grey), pmax(1, round(seq(y_lo, y_hi, length.out = out_px))))

  expand_grid(iy = seq_along(ys), ix = seq_along(xs)) %>%
    mutate(grey = grey[cbind(ys[iy], xs[ix])],
           x = (ix - 1) / (out_px - 1),
           y = (iy - 1) / (out_px - 1),
           .keep = "unused")
}
