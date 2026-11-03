# Paint out the red box the paper drew on its Fig. 7A plate photo: every red pixel is
# replaced by the mean grey of its non-red neighbours (repeated until none are left).
# Rscript slides/remove_red_box.R <in.png> <out.png>
library(png)
args <- commandArgs(trailingOnly = TRUE)
img <- readPNG(args[1])
r <- img[, , 1]; g <- img[, , 2]; b <- img[, , 3]
red <- r > 0.5 & r - g > 0.25 & r - b > 0.25
red <- red | (cbind(red[, -1], FALSE) | cbind(FALSE, red[, -ncol(red)]) |
              rbind(red[-1, ], FALSE) | rbind(FALSE, red[-nrow(red), ]))   # include anti-aliased edge
grey <- (r + g + b) / 3
grey[red] <- NA
while (any(is.na(grey))) {
  idx <- which(is.na(grey), arr.ind = TRUE)
  vals <- apply(idx, 1, function(p) {
    rr <- max(1, p[1] - 2):min(nrow(grey), p[1] + 2); cc <- max(1, p[2] - 2):min(ncol(grey), p[2] + 2)
    mean(grey[rr, cc], na.rm = TRUE)
  })
  grey[idx] <- vals   # NaN where no neighbour known yet -> stays NA, filled next round
  grey[is.nan(grey)] <- NA
}
out <- img[, , 1:3]
for (k in 1:3) out[, , k] <- ifelse(red, grey, img[, , k])
writePNG(out, args[2])
cat("painted out", sum(red), "pixels\n")
