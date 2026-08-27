#!/usr/bin/env Rscript
# C9 RNA-seq expression figure: copy-number-coloured tree | per-copy read counts
# (log10-per-copy stacked) | fractional locus coverage (gradient + colourbar below) | tissue.
suppressMessages(library(ape))

## ---- data ----
tree <- read.tree("figures/rna_seq_sps.nwk")
tree$tip.label[tree$tip.label == "Bassiana_duperreyi"] <- "Acritoscincus_duperreyi"
tree <- ladderize(tree)
Ntip <- length(tree$tip.label)

d <- read.csv("figures/rnaseq_species_summary.csv", stringsAsFactors = FALSE)
d <- d[match(tree$tip.label, d$species), ]
copies   <- c("C9A","C9B","C9C")
col_copy <- c(C9A = "steelblue3", C9B = "darkgoldenrod1", C9C = "mediumseagreen")
col_state <- c("1" = "#6699cc", "2" = "#1b9e77", "3" = "#e7298a")

rmat <- as.matrix(d[, c("reads_C9A","reads_C9B","reads_C9C")]); rmat[is.na(rmat)] <- 0
cmat <- as.matrix(d[, c("cov_C9A","cov_C9B","cov_C9C")])
lg   <- log10(rmat + 1)
ylim <- c(0.5, Ntip + 0.5)
h  <- 0.336; hc <- 0.204; wc <- 0.265        # cells +2%
cex_head <- 0.82; cex_sub <- 0.62
y_head <- Ntip + 2.3; y_sub <- Ntip + 1.0

anch_v <- c(0, 0.2, 0.6, 1.0); anch <- col2rgb(c("white","aquamarine","turquoise4","purple4"))
covcol <- function(v){ v<-pmin(pmax(v,0),1)
  rgb(approx(anch_v,anch[1,],v)$y, approx(anch_v,anch[2,],v)$y, approx(anch_v,anch[3,],v)$y, maxColorValue=255) }

state <- setNames(as.integer(d$n_copies), d$species)
getDesc <- function(n) if (n <= Ntip) n else unlist(lapply(tree$edge[tree$edge[,1]==n,2], getDesc))
edge_state <- sapply(tree$edge[,2], function(ch) min(state[tree$tip.label[getDesc(ch)]]))
edge_col <- col_state[as.character(edge_state)]

## ---- draw (called for each device) ----
draw <- function() {
  layout(matrix(1:3, nrow = 1), widths = c(1.5, 1.95, 1.5))
  par(oma = c(3.5, 0, 2.4, 0))

  # panel 1: tree
  par(mar = c(1, 0.3, 3.4, 0))
  dep <- max(node.depth.edgelength(tree))
  plot(tree, cex = 0.9, font = 3, label.offset = dep*0.015,
       x.lim = c(0, dep*2.05), y.lim = ylim, edge.color = edge_col, edge.width = 1)
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tipy <- setNames(lastPP$yy[1:Ntip], tree$tip.label)
  legend("bottomleft", legend = c("single copy","2 copies","3 copies"),
         col = col_state, lwd = 3, bty = "n", cex = 0.78, title = "copy number", inset = c(0, 0.01))

  # panel 2: per-copy read counts (stacked); enlarged spaced copy key above
  par(mar = c(1, 0.1, 3.4, 0.6))
  rowmax <- apply(lg, 1, sum); xmax <- max(rowmax)
  plot.new(); plot.window(xlim = c(0, xmax*1.02), ylim = ylim)
  gl <- seq(0, ceiling(xmax), by = 2)
  abline(v = gl, col = "grey93", lwd = 0.6); axis(1, at = gl, labels = gl, cex.axis = 0.72)
  mtext("per copy expression (log10)", side = 1, line = 2.2, cex = 0.7)
  for (i in seq_len(Ntip)) {
    y <- tipy[d$species[i]]
    if (rowmax[i] == 0) { text(0.1, y, "no expression detected", adj = 0, cex = 0.6, font = 3, col = "grey45")
    } else { x0 <- 0
      for (j in 1:3) { s <- lg[i,j]; if (s>0){ rect(x0,y-h,x0+s,y+h,col=col_copy[j],border=NA); x0<-x0+s } } }
  }
  text(xmax/2, y_head, "Read counts (uniquely mapped)", cex = cex_head, font = 2, xpd = NA)
  keyx <- xmax*c(0.30, 0.47, 0.64)
  for (j in 1:3) { rect(keyx[j]-xmax*0.012, y_sub-0.45, keyx[j]+xmax*0.012, y_sub+0.45,
                        col = col_copy[j], border = NA, xpd = NA)
    text(keyx[j]+xmax*0.02, y_sub, copies[j], adj = 0, cex = cex_sub, xpd = NA) }

  # panel 3: coverage cells (gradient) + colourbar below + tissue
  par(mar = c(1, 0.4, 3.4, 0.4))
  plot.new(); plot.window(xlim = c(0.3, 6), ylim = ylim)
  xc <- c(0.7, 1.26, 1.82); xt <- c(4.3, 5.1)
  for (i in seq_len(Ntip)) {
    y <- tipy[d$species[i]]; present <- strsplit(d$copies_present[i], ";")[[1]]
    for (j in 1:3) {
      has <- copies[j] %in% present; v <- cmat[i, j]
      rect(xc[j]-wc, y-hc, xc[j]+wc, y+hc,
           col = if (!has || is.na(v)) "grey96" else covcol(v), border = "grey85", lwd = 0.4)
    }
    tis <- d$tissue[i]
    if (grepl("liver", tis)) points(xt[1], y, pch = 21, bg = "brown3", col = "black", lwd = 0.4, cex = 1.35)
    if (grepl("lung",  tis)) points(xt[2], y, pch = 21, bg = "plum",   col = "black", lwd = 0.4, cex = 1.35)
  }
  # colourbar: beside the cells; ticks 0,0.2,...,1
  xbar <- 2.85; yb0 <- Ntip*0.30; yb1 <- Ntip*0.72; nb <- 60
  for (k in 1:nb) rect(xbar-0.14, yb0+(k-1)/nb*(yb1-yb0), xbar+0.14, yb0+k/nb*(yb1-yb0),
                       col = covcol((k-0.5)/nb), border = NA)
  rect(xbar-0.14, yb0, xbar+0.14, yb1, border = "grey60")
  for (val in seq(0,1,0.2)) { yy <- yb0+val*(yb1-yb0)
    segments(xbar+0.14, yy, xbar+0.22, yy, col="grey40"); text(xbar+0.30, yy, sprintf("%.1f",val), cex=0.5, adj=0) }
  # headings & sub-labels (matched sizes/lines)
  text(mean(xc), y_head, "Fractional locus", cex = cex_head, font = 2, xpd = NA)
  text(mean(xt), y_head, "Tissue",           cex = cex_head, font = 2, xpd = NA)
  text(xc, y_sub, copies, cex = cex_sub, xpd = NA)
  text(xt, y_sub, c("Liver","Lung"), cex = cex_sub, xpd = NA)

  mtext("C9 transcriptional expression across squamates (RNA-seq)",
        outer = TRUE, cex = 1.05, font = 2, line = 0.4)
}

png("figures/c9_rnaseq_expression_figure.png", width = 1230, height = 1350, res = 150); draw(); invisible(dev.off())
pdf("figures/c9_rnaseq_expression_figure.pdf", width = 8.2, height = 9);                draw(); invisible(dev.off())
cat("wrote figures/c9_rnaseq_expression_figure.{png,pdf}\n")
