#!/usr/bin/env Rscript
# Figure 1 (updated / professional): soft pastel clade-highlight blocks (Chapter-2 style, no names),
# thicker branches, larger labels/boxes/borders, bigger canvas. Separate output: Figure_1_updated.
suppressMessages(library(ape))

tree <- read.tree("figure_1/rep_species.nwk"); tree <- ladderize(tree, right = FALSE)
Ntip <- length(tree$tip.label)
d <- read.csv("figure_1/c9_status_data.csv", stringsAsFactors = FALSE)
d <- d[match(tree$tip.label, d$species), ]

gid <- c(DAB2="gray70", C9A="cadetblue4", C9B="mediumseagreen", C9C="steelblue3", FYB1="gray70")
PSEUDO <- "firebrick2"; selcol <- c(relaxed="darkorange2", ns="grey55")
colubrid <- c("Hemorrhois_hippocrepis","Arizona_elegans","Pantherophis_obsoletus","Pituophis_catenifer","Thamnophis_elegans")
elapid   <- c("Pseudonaja_textilis","Notechis_scutatus","Naja_naja")
snakes   <- c(colubrid, elapid, "Crotalus_tigris","Crotalus_adamanteus","Protobothrops_mucrosquamatus")
sel_of <- function(g,sp,o) if(g=="C9B" && o=="functional" && sp %in% snakes) "relaxed" else "ns"

tip_col  <- setNames(ifelse(d$absrel=="yes","red2","gray10"), d$species)
edge_col <- sapply(tree$edge[,2], function(ch) if(ch<=Ntip) tip_col[tree$tip.label[ch]] else "gray10")
edge_w   <- ifelse(edge_col=="red2", 3.8, 2.8)               # aBSREL red +1, black +2

## clades as contiguous tip blocks (Colubridae split: nested block + basal Hemorrhois) + pastel colours
clades <- list(
  Gekkota=list(c("Asaccus_caudivolvulus","Correlophus_ciliatus")),
  Scincoidea=list(c("Carinascincus_ocellatus","Cryptoblepharus_egeriae")),
  "Boidae/Pythonidae"=list(c("Python_bivittatus","Candoia_aspera")),
  Colubridae=list(c("Arizona_elegans","Pantherophis_obsoletus","Pituophis_catenifer","Thamnophis_elegans"),
                  c("Hemorrhois_hippocrepis")),
  Elapidae=list(c("Pseudonaja_textilis","Notechis_scutatus","Naja_naja")),
  Viperidae=list(c("Crotalus_tigris","Crotalus_adamanteus","Protobothrops_mucrosquamatus")),
  Anguimorpha=list(c("Dopasia_gracilis","Elgaria_multicarinata")),
  Iguania=list(c("Anolis_carolinensis","Pogona_vitticeps")),
  Lacertoidea=list(c("Zootoca_vivipara","Rhineura_floridana")))
clade_col <- c(Gekkota="lavenderblush2", Scincoidea="#E2EFDA", "Boidae/Pythonidae"="#EDE7F6",
  Colubridae="#DFF3F5", Elapidae="#FBE1DA", Viperidae="snow3", Anguimorpha="#FCE1EC",
  Iguania="#FFF3CC", Lacertoidea="#DEEBF7")

roundrect <- function(xc,yc,w,h,rx,ry,col,border="grey25",lwd=1,lty=1){
  BR<-cbind((xc+w/2-rx)+rx*cos(seq(-pi/2,0,len=8)),(yc-h/2+ry)+ry*sin(seq(-pi/2,0,len=8)))
  TR<-cbind((xc+w/2-rx)+rx*cos(seq(0,pi/2,len=8)),(yc+h/2-ry)+ry*sin(seq(0,pi/2,len=8)))
  TL<-cbind((xc-w/2+rx)+rx*cos(seq(pi/2,pi,len=8)),(yc+h/2-ry)+ry*sin(seq(pi/2,pi,len=8)))
  BL<-cbind((xc-w/2+rx)+rx*cos(seq(pi,3*pi/2,len=8)),(yc-h/2+ry)+ry*sin(seq(pi,3*pi/2,len=8)))
  polygon(rbind(BR,TR,TL,BL),col=col,border=border,lwd=lwd,lty=lty) }

dep <- max(node.depth.edgelength(tree)); xlab <- dep*1.04
bw <- dep*0.10; bh <- 0.80; rx <- bw*0.14; ry <- bh*0.14
LCEX <- 1.0                                                   # tree label size (bigger)
BB   <- 2.0                                                   # box border
BLREL<- 3.4                                                   # relaxed-outline border
GLINE<- 2.4                                                   # gene backbone line
XR   <- dep*2.55

draw <- function(){
  par(mar=c(5.5, 1, 3.6, 0.8))
  plot(tree, show.tip.label=FALSE, edge.color=edge_col, edge.width=edge_w, x.lim=c(0, XR+bw), y.lim=c(0.4, Ntip+1.3))
  pp <- get("last_plot.phylo", envir=.PlotPhyloEnv); tipy<-pp$yy[1:Ntip]; xx<-pp$xx; yy<-pp$yy

  ## clade highlight blocks (behind), hugging each clade's branches (short); then redraw the tree on top
  for(cl in names(clades)) for(blk in clades[[cl]]){
    idx<-match(blk,tree$tip.label); yr<-range(tipy[idx])
    x0 <- if(length(blk)>1) xx[getMRCA(tree,blk)] else xx[tree$edge[tree$edge[,2]==idx,1]]
    bcol <- if("Hemorrhois_hippocrepis" %in% blk) "white" else clade_col[cl]   # Hemorrhois left white
    rect(x0, yr[1]-0.45, dep*1.02, yr[2]+0.45, col=bcol, border=NA)             # square clade blocks
  }
  ee<-tree$edge
  for(e in seq_len(nrow(ee))) segments(xx[ee[e,1]], yy[ee[e,2]], xx[ee[e,2]], yy[ee[e,2]], col=edge_col[e], lwd=edge_w[e], lend=1)
  for(nd in unique(ee[,1])){ ch<-ee[ee[,1]==nd,2]; segments(xx[nd], min(yy[ch]), xx[nd], max(yy[ch]), col="gray10", lwd=2.8, lend=1) }

  for(i in seq_len(Ntip)) text(xlab, tipy[i], gsub("_"," ",tree$tip.label[i]), adj=0, cex=LCEX, font=3, xpd=NA)
  maxnw <- max(strwidth(gsub("_"," ",tree$tip.label), cex=LCEX, font=3))
  base_x <- xlab + maxnw + dep*0.15
  xf <- base_x + dep*0.66; SPAN <- xf - base_x
  hx <- base_x + (0:4)*SPAN/4
  text(hx, Ntip+1.0, c("DAB2","C9A","C9B","C9C","FYB1"), font=4, cex=0.82, col="black", xpd=NA)
  for(i in seq_len(Ntip)){ y<-tipy[i]; sp<-tree$tip.label[i]
    cps<-c(); sts<-c()
    for(cc in c("C9A","C9B","C9C")) if(d[[tolower(cc)]][i]!="absent"){ cps<-c(cps,cc); sts<-c(sts,d[[tolower(cc)]][i]) }
    n<-length(cps)
    segments(base_x, y, xf, y, col="grey55", lwd=GLINE)
    roundrect(base_x,y,bw,bh,rx,ry,col=gid["DAB2"],border="grey35",lwd=BB)
    if(n>0) for(k in 1:n){ x<-base_x + k*SPAN/(n+1); cc<-cps[k]; st<-sts[k]
      if(st=="lost") text(x,y,"locus deleted",cex=0.52,font=3,col="firebrick2")
      else if(st=="remnant") roundrect(x,y,bw,bh,rx,ry,col=adjustcolor(gid[cc],0.40),border="grey45",lwd=1.6,lty=3)
      else if(st=="stops") roundrect(x,y,bw,bh,rx,ry,col=PSEUDO,border="grey25",lwd=BB)
      else roundrect(x,y,bw,bh,rx,ry,col=gid[cc],border=if(sel_of(cc,sp,st)=="relaxed")selcol["relaxed"] else "grey35",
                     lwd=if(sel_of(cc,sp,st)=="relaxed")BLREL else BB) }
    roundrect(xf,y,bw,bh,rx,ry,col=gid["FYB1"],border="grey35",lwd=BB) }

  mark <- function(node,txt){ e<-which(tree$edge[,2]==node); if(!length(e))return()
    x<-(xx[tree$edge[e,1]]+xx[node])/2; yv<-yy[node]
    points(x,yv,pch=23,bg="#333333",col="black",cex=1.9,lwd=1.4); text(x,yv+0.95,txt,cex=0.82,font=2) }
  mark(getMRCA(tree,snakes),"+C9B"); mark(getMRCA(tree,c("Arizona_elegans","Pituophis_catenifer")),"+C9C")
  for(sp in c("Asaccus_caudivolvulus","Dopasia_gracilis","Protobothrops_mucrosquamatus")){ i<-match(sp,tree$tip.label); points(dep*0.99,tipy[i],pch=4,col="#B2182B",lwd=3,cex=1.3) }

  mtext("Evolution of C9 copy number and genomic organization across Squamata",
        side=3, line=0.8, adj=0.5, cex=1.55, font=2)
  try(add.scale.bar(0,-0.7,cex=0.8,lwd=2.5),silent=TRUE)
  yl <- -1.6
  legend(x=0, y=yl, title="Branch (aBSREL)", bty="n", cex=0.82, lwd=4, seg.len=1.2, xpd=NA,
    col=c("red2","gray10"), legend=c("aBSREL positive","not significant"))
  legend(x=dep*0.66, y=yl, title="Gene (box fill)", bty="n", cex=0.82, xpd=NA, pch=22, pt.cex=2.2,
    pt.bg=c("gray70","cadetblue4","mediumseagreen","steelblue3","firebrick2"), col="grey40",
    legend=c("DAB2 / FYB1 (flanking)","C9A","C9B","C9C","pseudogene (stops)"))
  legend(x=dep*1.42, y=yl, title="ORF / outline", bty="n", cex=0.82, xpd=NA, pch=22, pt.cex=2.2, pt.lwd=c(1,1.4,3),
    pt.bg=c("grey85",adjustcolor("grey70",0.4),"white"), col=c("grey40","grey55","darkorange2"),
    legend=c("functional","remnant (dotted)","C9B relaxed (orange)"))
  legend(x=dep*1.95, y=yl, title="Events", bty="n", cex=0.82, xpd=NA, pch=c(23,4), pt.bg=c("#333333",NA),
    col=c("black","#B2182B"), pt.cex=c(1.7,1.4), legend=c("duplication (+C9B/+C9C)","loss / pseudogenisation"))
}
png("figure_1/Figure_1_updated.png", width=2300, height=1900, res=180); draw(); invisible(dev.off())
pdf("figure_1/Figure_1_updated.pdf", width=14.8, height=12.2);         draw(); invisible(dev.off())
cat("wrote figure_1/Figure_1_updated.{png,pdf}\n")
