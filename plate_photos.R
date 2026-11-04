# Cut the raw colony images out of the plate photographs.
#
# Iris turned each plate photo into one opacity number per colony. To look at
# the colonies themselves we need to find the colony grid in the photo again.
# That is done by lining the photo up with the Iris measurements of the same
# plate: a grid that is off by even one colony correlates badly, so a misplaced
# grid cannot pass unnoticed. The fit is cached next to the photo.
#
# Photos live in data/photos, named after the plate exactly like the .iris
# files, e.g. data/photos/SpecGentAraIPTG-10050201-3A_1.png. Drop a plate's
# photo in there and it becomes available to inspect_gene(); plates without a
# photo are simply left out.
library(png)

photo_dir <- "data/photos"
iris_dir  <- "data/iris"

# Which plates we have a photo for.
plates_with_photos <- function() {
  str_remove(basename(list.files(photo_dir, "\\.png$")), "\\.png$")
}

plate_photo_path <- function(plate) file.path(photo_dir, paste0(plate, ".png"))


# Fitting the colony grid ----------------------------------------------------
# A grid is four numbers: the centre of the first colony and the spacing
# between colonies, in pixels, for each axis.

fit_colony_grid <- function(plate, quiet = FALSE, iris_file = NULL) {
  cache <- file.path(photo_dir, paste0(plate, ".grid.csv"))
  if (file.exists(cache)) return(as.list(read.csv(cache)))

  if (is.null(iris_file)) iris_file <- file.path(iris_dir, paste0(plate, ".jpg.iris"))
  iris <- read.delim(iris_file, skip = 6)
  full <- readPNG(plate_photo_path(plate))[, , 1]

  # Searching the grid means reading thousands of little patches, so a large
  # photo is fitted on a thinned-out copy and the answer scaled back up. The
  # colonies are tens of pixels across, so nothing is lost.
  step <- max(1, floor(sqrt(length(full) / 2.2e6)))
  grey <- full[seq(1, nrow(full), step), seq(1, ncol(full), step)]

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

  # starting guess as fractions of the image, so it survives a different size
  grid <- list(x0 = 0.1286 * ncol(grey), px = 0.03293 * ncol(grey),
               y0 = 0.1265 * nrow(grey), py = 0.0493  * nrow(grey))

  refine <- function(g, field, values) {
    best <- g; best_score <- score_grid(g)
    for (v in values) {
      cand <- g; cand[[field]] <- v
      s <- score_grid(cand)
      if (s > best_score) { best <- cand; best_score <- s }
    }
    best
  }
  for (pass in 1:2) {
    grid <- refine(grid, "x0", seq(grid$x0 - 25, grid$x0 + 25, by = 0.5))
    grid <- refine(grid, "px", seq(grid$px - 1.5, grid$px + 1.5, by = 0.02))
    grid <- refine(grid, "y0", seq(grid$y0 - 25, grid$y0 + 25, by = 0.5))
    grid <- refine(grid, "py", seq(grid$py - 1.5, grid$py + 1.5, by = 0.02))
  }

  agreement <- score_grid(grid)
  if (!quiet) cat(sprintf("%s: colony grid fitted, r = %.3f with Iris opacity\n",
                          plate, agreement))
  if (agreement < 0.85)
    stop(plate, ": the fitted grid does not match the Iris measurements (r = ",
         round(agreement, 2), ") - refusing to cut out colonies")

  # back to the coordinates of the full-size photo
  grid <- list(x0 = (grid$x0 - 1) * step + 1, px = grid$px * step,
               y0 = (grid$y0 - 1) * step + 1, py = grid$py * step,
               agreement = agreement)
  write.csv(as.data.frame(grid), cache, row.names = FALSE)
  grid
}


# Cutting out the colonies ---------------------------------------------------
# One 96-well plate is pinned into 384 format, so each well is a 2x2 quartet of
# colonies. quartet_pixels() returns the photo pixels of one quartet as a long
# table, with x and y running 0-1 inside the cut-out so that quartets from
# different plates can be drawn side by side.

quartet_pixels <- function(plate, well_row, well_col, pad = 0.6, out_px = 220) {
  grid <- fit_colony_grid(plate)
  grey <- readPNG(plate_photo_path(plate))[, , 1]

  centre <- function(i, first, pitch) first + pitch * (i - 1)
  x_lo <- centre(2 * well_col - 1, grid$x0, grid$px) - grid$px * pad
  x_hi <- centre(2 * well_col,     grid$x0, grid$px) + grid$px * pad
  y_lo <- centre(2 * well_row - 1, grid$y0, grid$py) - grid$py * pad
  y_hi <- centre(2 * well_row,     grid$y0, grid$py) + grid$py * pad

  # Sampled onto a fixed out_px x out_px grid, so cut-outs from plates
  # photographed at different sizes can be drawn next to each other.
  xs <- pmin(ncol(grey), pmax(1, round(seq(x_lo, x_hi, length.out = out_px))))
  ys <- pmin(nrow(grey), pmax(1, round(seq(y_lo, y_hi, length.out = out_px))))

  expand_grid(iy = seq_along(ys), ix = seq_along(xs)) %>%
    mutate(grey = grey[cbind(ys[iy], xs[ix])],
           x = (ix - 1) / (out_px - 1),
           y = (iy - 1) / (out_px - 1),
           .keep = "unused")
}
