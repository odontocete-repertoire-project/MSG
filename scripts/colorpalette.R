library(RColorBrewer)
library(colorspace)

pinks <- c("#ffe8ee", "#ffd6e0", "#feb4c7", "#fe7799","#fe0849", "#710a25")
blues   <- brewer.pal(6, "Blues")
purples <- brewer.pal(6, "Purples")
greens  <- brewer.pal(6, "Greens")

all_colors <- c(pinks, blues, purples, greens)

sg_standard_color <- pinks[4]
tt_standard_color <- blues[4]
msg_standard_color <- purples[4]

main_colors <- c(sg_standard_color,tt_standard_color, msg_standard_color)

swatchplot(
  "Pinks" = pinks,
  "Blues" = blues,
  "Purples" = purples,
  "Greens" = greens
)

swatchplot(main_colors)
