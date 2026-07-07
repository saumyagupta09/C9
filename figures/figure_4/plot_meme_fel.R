#!/usr/bin/env Rscript
# Selection across C9 protein: domain architecture + FEL (pervasive) & MEME (episodic) site tracks.
# Input: a per-residue CSV built by build_selection_sites.py.
# Args (optional): 1=input csv  2=output basename  3=dataset label
args    <- commandArgs(trailingOnly = TRUE)
input   <- if (length(args) >= 1) args[1] else "figure_4/c9_selection_sites.csv"
outbase <- if (length(args) >= 2) args[2] else "figure_4/c9_domain_meme_fel_figure"
dslab   <- if (length(args) >= 3) args[3] else "snake C9A"
d <- read.csv(input, stringsAsFactors = FALSE)
L <- max(d$ref_pos)

## domains (InterProScan on the Pituophis C9A reference) — CCP/FIMAC are C6/C7, not present in C9
dom <- data.frame(
  name  = c("signal","TSP1","LDLRA","MACPF","EGF","TSP1"),
  start = c(1, 38, 96, 135, 496, 537),
  end   = c(20, 91, 133, 495, 536, 583),
  col   = c("#BFBFBF","cornflowerblue","#F0C000","sandybrown","mediumseagreen","cornflowerblue"),
  stringsAsFactors = FALSE)

col_fel  <- c(positive="#C0392B", purifying="#2CA02C", ns="#C4C4C4")
col_meme <- c(episodic="#7B3FA0", ns="#C4C4C4")

FCAP_LO <- -6; FCAP_HI <- 5                       # clip extreme FEL dN-dS for display
d$fel_plot <- pmax(pmin(d$fel_dNdS, FCAP_HI), FCAP_LO)
d$meme_lp  <- -log10(pmax(d$meme_p, 1e-4))

xlim <- c(-2, L+2)
band <- function() for (i in seq_len(nrow(dom)))
  rect(dom$start[i], par("usr")[3], dom$end[i], par("usr")[4],
       col = adjustcolor(dom$col[i], 0.10), border = NA)
MARL <- 12                                         # shared left margin (labels live here)

draw <- function() {
  layout(matrix(1:4, ncol = 1), heights = c(1.7, 3, 3, 0.9))
  par(oma = c(1, 3, 3, 1))

  ## main title
  ## panel 1: domain architecture, with its label placed BESIDE the bar (left margin)
  par(mar = c(0.5, MARL, 3.0, 2), xaxs = "i")
  plot.new(); plot.window(xlim = xlim, ylim = c(0, 1))
  rect(1, 0.30, L, 0.70, col = "grey93", border = "grey70")
  for (i in seq_len(nrow(dom))) {
    rect(dom$start[i], 0.12, dom$end[i], 0.88, col = dom$col[i], border = "grey30")
    if (dom$name[i] != "signal")
      text(mean(c(dom$start[i], dom$end[i])), 0.5, dom$name[i], cex = 0.8, font = 2, col = "grey10")
  }
  text(mean(c(135,495)), -0.05, "(PF01823) pore-forming domain", cex = 0.6, col = "grey30", xpd = NA)
  ## position ruler spanning the FULL protein length (N-terminus .. end of last TSP1)
  axis(3, at = seq(0, 550, 50), cex.axis = 0.6, tcl = -0.3, mgp = c(2,0.3,0))
  axis(3, at = L, labels = L, cex.axis = 0.6, tcl = -0.3, mgp = c(2,0.3,0))
  axis(3, at = c(0, L), labels = FALSE, tcl = 0)                 # extend ruler line to residue 584
  text(-3, 0.5, "N", cex = 0.85, font = 2, adj = 1, xpd = NA)
  text(L+4, 0.5, "C", cex = 0.85, font = 2, adj = 0, xpd = NA)
  ## the (former) heading, now a label beside the domain figure; species name italicised
  mtext("C9A protein",        side = 2, las = 1, at = 0.86, line = 1.2, adj = 1, cex = 0.6, font = 2)
  mtext("domain architecture", side = 2, las = 1, at = 0.66, line = 1.2, adj = 1, cex = 0.6, font = 2)
  mtext(expression("("*italic("Pituophis catenifer")), side = 2, las = 1, at = 0.40, line = 1.2, adj = 1, cex = 0.58)
  mtext("C9A reference)",       side = 2, las = 1, at = 0.20, line = 1.2, adj = 1, cex = 0.58)

  ## main title, shifted right so "C9" sits directly above residue 300 of the bar
  tcex <- 1.45
  text(300 - strwidth("Selection across ", cex = tcex, font = 2), 1.82,
       "Selection across C9 protein", adj = 0, cex = tcex, font = 2, xpd = NA)

  ## panel 2: FEL — heading lifted above the panel with a gap
  par(mar = c(0.5, MARL, 2.3, 2))
  plot.new(); plot.window(xlim = xlim, ylim = c(FCAP_LO*1.03, FCAP_HI*1.05)); band()
  abline(h = 0, col = "grey55", lwd = 0.7)
  for (i in seq_len(nrow(d)))
    rect(d$ref_pos[i]-0.45, 0, d$ref_pos[i]+0.45, d$fel_plot[i], col = col_fel[d$fel_class[i]], border = NA)
  axis(2, cex.axis = 0.62, las = 1, mgp = c(2,0.5,0)); box(col="grey70")
  mtext(expression(d[N]-d[S]), side = 2, line = 2.6, cex = 0.72)
  mtext("FEL - pervasive selection", side = 3, line = 0.7, adj = 0, cex = 0.82, font = 2)
  legend("bottomright", bg = "white", box.col = "grey80", cex = 0.6, fill = col_fel,
         legend = c("Pervasive positive (p<0.05)","Pervasive negative (p<0.05)","Not significant"))

  ## panel 3: MEME
  par(mar = c(2.4, MARL, 2.3, 2))
  plot.new(); plot.window(xlim = xlim, ylim = c(0, max(d$meme_lp)*1.10)); band()
  for (i in seq_len(nrow(d)))
    rect(d$ref_pos[i]-0.45, 0, d$ref_pos[i]+0.45, d$meme_lp[i], col = col_meme[d$meme_class[i]], border = NA)
  abline(h = -log10(0.05), col = "grey55", lty = 2, lwd = 0.7)
  axis(2, cex.axis = 0.62, las = 1, mgp = c(2,0.5,0)); box(col="grey70")
  axis(1, at = seq(0, L, 50), cex.axis = 0.62, mgp = c(2,0.4,0))
  mtext("-log10(p)", side = 2, line = 2.6, cex = 0.72)
  mtext("residue position (Pituophis catenifer C9A)", side = 1, line = 1.4, cex = 0.66)
  mtext("MEME - episodic diversifying selection", side = 3, line = 0.7, adj = 0, cex = 0.82, font = 2)
  legend("topright", bg = "white", box.col = "grey80", cex = 0.6, fill = col_meme,
         legend = c("Episodic diversifying (p<0.05)","Not significant"))

  ## panel 4: domain legend + note
  par(mar = c(0.2, MARL, 0.2, 2))
  plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
  legend("left", horiz = TRUE, bty = "n", cex = 0.62,
         fill = c("cornflowerblue","#F0C000","sandybrown","mediumseagreen"),
         legend = c("TSP1 (PF00090)","LDLRA (PF00057)","MACPF (PF01823)","EGF (PF07645)"))
  mtext(paste0("selection dataset: ", dslab,
               "   |   residue numbering per the Pituophis catenifer C9A InterProScan reference; extreme FEL dN-dS clipped for display."),
        side = 1, line = -0.7, adj = 1, cex = 0.52, col = "grey35")
}

png(paste0(outbase, ".png"), width = 1500, height = 1000, res = 150); draw(); invisible(dev.off())
pdf(paste0(outbase, ".pdf"), width = 10, height = 6.7);               draw(); invisible(dev.off())
cat("wrote ", outbase, ".{png,pdf}\n", sep = "")
