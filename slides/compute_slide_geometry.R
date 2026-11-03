# Work out where to place each image on its slide so that the colony grid lands
# in the same rectangle on every slide: the plate photo and all the heatmaps.
#
# The images are not interchangeable - the photo has the plate's plastic rim
# around the colonies, a heatmap has axis labels and a legend - so the same
# on-screen box would NOT line the grids up. For each image we find the grid
# inside it and then scale and offset the image so the grid lands on target.
# Rscript slides/compute_slide_geometry.R

# the on-screen rectangle the colony grid should occupy, on the 1920x1080 slide
TARGET <- list(left = 246, top = 310, width = 820, height = 546.7)   # 24:16 = 1.5

boxes <- rbind(read.csv("slides/figures/photo_grid.csv"),
               read.csv("slides/figures/heatmap_grid.csv"))

boxes$grid_w <- boxes$right - boxes$left
boxes$grid_h <- boxes$bottom - boxes$top
boxes$scale  <- TARGET$width / boxes$grid_w
boxes$css_width  <- round(boxes$width  * boxes$scale, 1)
boxes$css_height <- round(boxes$height * boxes$scale, 1)
boxes$css_left   <- round(TARGET$left - boxes$left * boxes$scale, 1)
boxes$css_top    <- round(TARGET$top  - boxes$top  * boxes$scale, 1)
# a check: with the width matched, does the height land on target too?
boxes$height_err <- round(boxes$grid_h * boxes$scale - TARGET$height, 1)

print(boxes[, c("image", "grid_w", "grid_h", "css_left", "css_top", "css_width", "css_height", "height_err")],
      row.names = FALSE)
cat("\nUse: position:absolute; left:<css_left>px; top:<css_top>px; width:<css_width>px; height:<css_height>px\n")
cat("object-fit must be 'fill' or omitted, not 'contain': the box already has the image's aspect ratio.\n")
