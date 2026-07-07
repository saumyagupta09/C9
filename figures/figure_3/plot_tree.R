#!/usr/bin/env Rscript
# Figure 3: C8/C9 nucleotide ML tree (Fig-5 style), rooted on C8, tips coloured by category,
# ultrafast-bootstrap support shown on well-supported nodes.
suppressMessages(library(ape))

t <- read.tree("figure_3/c8_c9_tree.contree")
c8 <- t$tip.label[grepl("_C8[AB]$", t$tip.label)]
t <- root(t, outgroup = c8, resolve.root = TRUE)
t <- ladderize(t)
Ntip <- length(t$tip.label)

cat_of <- function(l) if(grepl("_C8A$",l))"C8a" else if(grepl("_C8B$",l))"C8b" else
  if(grepl("_C9A$",l))"C9A" else if(grepl("_C9B$",l))"C9B" else if(grepl("_C9C$",l))"C9C" else "C9out"
catcol <- c(C8a="#8c6d31", C8b="#c2a83e", C9A="cadetblue4", C9B="mediumseagreen",
            C9C="steelblue3", C9out="grey35")
tc <- catcol[sapply(t$tip.label, cat_of)]
lab <- gsub("_"," ", t$tip.label)

bs <- suppressWarnings(as.numeric(t$node.label))   # UFBoot support

draw <- function(){
  par(mar=c(1,0.5,2,0.5))
  plot(t, tip.color=tc, cex=0.42, label.offset=0.01, edge.width=0.7, no.margin=FALSE, font=3)
  # bootstrap on internal nodes with support >= 70
  nl <- rep("", length(bs)); nl[!is.na(bs) & bs>=70] <- bs[!is.na(bs) & bs>=70]
  nodelabels(nl, cex=0.32, frame="none", col="grey20", adj=c(1.2,-0.35))
  add.scale.bar(cex=0.6, lwd=2)
  mtext("C8 / C9 maximum-likelihood tree (nucleotide, GTR+F+R5; UFBoot x1000)  — rooted on C8",
        side=3, line=0.5, cex=0.95, font=2, adj=0)
  legend("bottomleft", inset=c(0.02,0.02), bty="n", cex=0.62, pch=15, pt.cex=1.3,
    col=catcol, legend=c("C8α","C8β","C9A (DAB2-prox)","C9B","C9C","C9 outgroup (frog/bird/turtle/croc)"))
}
png("figure_3/figure3_c8c9_tree.png", width=1400, height=2000, res=170); draw(); invisible(dev.off())
pdf("figure_3/figure3_c8c9_tree.pdf", width=9, height=13);              draw(); invisible(dev.off())
cat("wrote figure_3/figure3_c8c9_tree.{png,pdf}\n")
