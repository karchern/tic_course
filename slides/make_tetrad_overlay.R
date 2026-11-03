# Draw a box around every 2x2 tetrad of colonies on the plate photo:
# red for library designs, green for the GFP controls.
# The colony grid is fitted by lining the photo up with the Iris measurements of
# the same plate, so a misplaced grid cannot pass unnoticed.
# Rscript slides/make_tetrad_overlay.R
library(png)

photo     <- "slides/figures/plate_photo_AraIPTG_Sen2_nobox.png"
iris_file <- "data/iris/SpecGentAraIPTG-10050201-3A_1.jpg.iris"

img  <- readPNG(photo)
grey <- img[, , 1]
iris <- read.delim(iris_file, skip = 6)

# Brightness of a predicted colony above its local agar background. This mirrors
# what Iris calls opacity, so the two are directly comparable.
patch <- function(x, y, d) {
  x <- round(x); y <- round(y)
  if (x - 2 * d < 1 || y - 2 * d < 1 || x + 2 * d > ncol(grey) || y + 2 * d > nrow(grey)) return(NA_real_)
  mean(grey[(y - d):(y + d), (x - d):(x + d)]) -
    median(grey[(y - 2 * d):(y + 2 * d), (x - 2 * d):(x + 2 * d)])
}

# how well a candidate grid explains the plate: brightness at each predicted
# colony against the opacity Iris measured there
score_grid <- function(g) {
  d <- round(g$px * 0.25)
  b <- mapply(function(r, cc) patch(g$x0 + g$px * (cc - 1), g$y0 + g$py * (r - 1), d), iris$row, iris$column)
  if (anyNA(b)) return(-1)
  cor(b, iris$opacity)
}

# starting guess, as fractions of the image so it survives a different export size
grid <- list(x0 = 0.1286 * ncol(grey), px = 0.03293 * ncol(grey),
             y0 = 0.1265 * nrow(grey), py = 0.0493  * nrow(grey))

# coordinate descent: refine one parameter at a time, twice over
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
cat(sprintf("columns: pitch %.2f px, first %.1f, last %.1f\nrows:    pitch %.2f px, first %.1f, last %.1f\nagreement with Iris opacity: r = %.3f\n",
            grid$px, grid$x0, grid$x0 + grid$px * 23,
            grid$py, grid$y0, grid$y0 + grid$py * 15, agreement))
if (agreement < 0.85) stop("grid fit does not match the Iris measurements - not drawing boxes")

plate_map <- read.csv("data/AlphaDesign_384wp.identifiers.csv")
gfp <- unique(data.frame(tx = ceiling(plate_map$column.384[plate_map$genename == "gfp"] / 2),
                         ty = ceiling(plate_map$row.384[plate_map$genename == "gfp"] / 2)))
cat("GFP tetrads:", nrow(gfp), "of 96\n")

draw_box <- function(im, x0, x1, y0, y1, rgb, w) {
  x0 <- max(1, x0); y0 <- max(1, y0); x1 <- min(ncol(im), x1); y1 <- min(nrow(im), y1)
  for (k in 1:3) {
    im[y0:(y0 + w), x0:x1, k] <- rgb[k]
    im[(y1 - w):y1, x0:x1, k] <- rgb[k]
    im[y0:y1, x0:(x0 + w), k] <- rgb[k]
    im[y0:y1, (x1 - w):x1, k] <- rgb[k]
  }
  im
}

colony_x <- function(c) grid$x0 + grid$px * (c - 1)
colony_y <- function(r) grid$y0 + grid$py * (r - 1)

out <- img[, , 1:3]
for (tx in 1:12) for (ty in 1:8) {
  is_gfp <- any(gfp$tx == tx & gfp$ty == ty)
  out <- draw_box(out,
                  round(colony_x(2 * tx - 1) - grid$px * 0.42), round(colony_x(2 * tx) + grid$px * 0.42),
                  round(colony_y(2 * ty - 1) - grid$py * 0.42), round(colony_y(2 * ty) + grid$py * 0.42),
                  if (is_gfp) c(0.10, 0.75, 0.25) else c(0.90, 0.10, 0.10), 4)
}

writePNG(out, "slides/figures/plate_photo_tetrads.png")
cat("wrote slides/figures/plate_photo_tetrads.png\n")
