#!/usr/bin/env Rscript
# Figure 1 (variant): DAB2 & FYB1 fixed flanking columns, C9 copies in tight fixed columns between
# them (mimics another_example.png). New palette; whitespace trimmed. Separate output file.
suppressMessages(library(ape))

tree <- read.tree("figure_1/rep_species.nwk"); tree <- ladderize(tree, right = FALSE)
Ntip <- length(tree$tip.label)
d <- read.csv("figure_1/c9_status_data.csv", stringsAsFactors = FALSE)
d <- d[match(tree$tip.label, d$species), ]

gid <- c(DAB2="gray70", C9A="cadetblue4", C9B="mediumseagreen", C9C="steelblue3", FYB1="gray70")
PSEUDO <- "firebrick2"
selcol <- c(relaxed="darkorange2", ns="grey55")
colubrid <- c("Hemorrhois_hippocrepis","Arizona_elegans","Pantherophis_obsoletus","Pituophis_catenifer","Thamnophis_elegans")
elapid   <- c("Pseudonaja_textilis","Notechis_scutatus","Naja_naja")
snakes   <- c(colubrid, elapid, "Crotalus_tigris","Crotalus_adamanteus","Protobothrops_mucrosquamatus")
# orange box outline = C9B relaxation only
sel_of <- function(g,sp,o) if(g=="C9B" && o=="functional" && sp %in% snakes) "relaxed" else "ns"

# branches: aBSREL only (blue if positive, else black); internal branches black
tip_col  <- setNames(ifelse(d$absrel=="yes","seagreen4","gray10"), d$species)
edge_col <- sapply(tree$edge[,2], function(ch) if(ch<=Ntip) tip_col[tree$tip.label[ch]] else "gray10")
edge_w  <- ifelse(edge_col=="seagreen4", 2.8, 0.8)

roundrect <- function(xc,yc,w,h,rx,ry,col,border="grey25",lwd=1,lty=1){
  BR<-cbind((xc+w/2-rx)+rx*cos(seq(-pi/2,0,len=6)),(yc-h/2+ry)+ry*sin(seq(-pi/2,0,len=6)))
  TR<-cbind((xc+w/2-rx)+rx*cos(seq(0,pi/2,len=6)),(yc+h/2-ry)+ry*sin(seq(0,pi/2,len=6)))
  TL<-cbind((xc-w/2+rx)+rx*cos(seq(pi/2,pi,len=6)),(yc+h/2-ry)+ry*sin(seq(pi/2,pi,len=6)))
  BL<-cbind((xc-w/2+rx)+rx*cos(seq(pi,3*pi/2,len=6)),(yc-h/2+ry)+ry*sin(seq(pi,3*pi/2,len=6)))
  polygon(rbind(BR,TR,TL,BL),col=col,border=border,lwd=lwd,lty=lty) }

dep <- max(node.depth.edgelength(tree)); xlab <- dep*1.03
bw <- dep*0.077; bh <- 0.62; rx <- bw*0.14; ry <- bh*0.16
XR <- dep*2.13                                   # right edge of content (trim whitespace)

draw <- function(){
  par(mar=c(5, 1, 3.2, 0.6))
  plot(tree, show.tip.label=FALSE, edge.color=edge_col, edge.width=edge_w, x.lim=c(0, XR+bw), y.lim=c(0.4, Ntip+1.3))
  lastPP <- get("last_plot.phylo", envir=.PlotPhyloEnv); tipy <- lastPP$yy[1:Ntip]
  for(i in seq_len(Ntip)) text(xlab, tipy[i], gsub("_"," ",tree$tip.label[i]), adj=0, cex=0.72, font=3, xpd=NA)

  maxnw <- max(strwidth(gsub("_"," ",tree$tip.label), cex=0.72, font=3))
  base_x <- xlab + maxnw + dep*0.15                # DAB2 fixed here
  xf <- base_x + dep*0.40                          # FYB1 fixed here (compact span)
  SPAN <- xf - base_x
  hx <- base_x + (0:4)*SPAN/4                       # header reference (triplicate positions)
  text(hx, Ntip+0.9, c("DAB2","C9A","C9B","C9C","FYB1"), font=4, cex=0.6, col="black", xpd=NA)
  for(i in seq_len(Ntip)){ y<-tipy[i]; sp<-tree$tip.label[i]
    cps<-c(); sts<-c()
    for(cc in c("C9A","C9B","C9C")) if(d[[tolower(cc)]][i]!="absent"){ cps<-c(cps,cc); sts<-c(sts,d[[tolower(cc)]][i]) }
    n<-length(cps)
    segments(base_x, y, xf, y, col="grey70", lwd=0.9)
    roundrect(base_x,y,bw,bh,rx,ry,col=gid["DAB2"],border="grey45",lwd=0.7)
    if(n>0) for(k in 1:n){ x<-base_x + k*SPAN/(n+1); cc<-cps[k]; st<-sts[k]   # genes evenly distributed
      if(st=="lost") text(x,y,"locus deleted",cex=0.40,font=3,col="firebrick2")
      else if(st=="remnant") roundrect(x,y,bw,bh,rx,ry,col=adjustcolor(gid[cc],0.40),border="grey55",lwd=0.9,lty=3)
      else if(st=="stops") roundrect(x,y,bw,bh,rx,ry,col=PSEUDO,border="grey30",lwd=0.7)
      else roundrect(x,y,bw,bh,rx,ry,col=gid[cc],border=selcol[sel_of(cc,sp,st)],lwd=if(sel_of(cc,sp,st)=="ns")0.7 else 1.8) }
    roundrect(xf,y,bw,bh,rx,ry,col=gid["FYB1"],border="grey45",lwd=0.7) }

  mark <- function(node,txt){ e<-which(tree$edge[,2]==node); if(!length(e))return()
    xx<-(lastPP$xx[tree$edge[e,1]]+lastPP$xx[node])/2; yy<-lastPP$yy[node]
    points(xx,yy,pch=23,bg="#4C4C4C",col="black",cex=1.4); text(xx,yy+0.7,txt,cex=0.6,font=2) }
  mark(getMRCA(tree,snakes),"+C9B"); mark(getMRCA(tree,c("Arizona_elegans","Pituophis_catenifer")),"+C9C")
  for(sp in c("Asaccus_caudivolvulus","Dopasia_gracilis","Protobothrops_mucrosquamatus")){ i<-match(sp,tree$tip.label); points(dep*0.99,tipy[i],pch=4,col="#B2182B",lwd=2,cex=0.9) }

  mtext("Evolution of C9 copy number and genomic organization across Squamata", side=3, line=0.6, adj=0.5, cex=1.3, font=2)
  try(add.scale.bar(0,-0.6,cex=0.55,lwd=2),silent=TRUE)

  yl <- -1.2
  legend(x=0,        y=yl, title="Branch (aBSREL)", bty="n", cex=0.58, lwd=3, seg.len=1, xpd=NA,
    col=c("seagreen4","gray10"),
    legend=c("aBSREL positive","not significant"))
  legend(x=dep*0.62, y=yl, title="Gene (box fill)", bty="n", cex=0.58, xpd=NA, pch=22, pt.cex=1.6,
    pt.bg=c("gray70","cadetblue4","mediumseagreen","steelblue3","firebrick2"), col="grey40",
    legend=c("DAB2 / FYB1 (flanking)","C9A","C9B","C9C","pseudogene (stops)"))
  legend(x=dep*1.25, y=yl, title="ORF / outline", bty="n", cex=0.58, xpd=NA, pch=22, pt.cex=1.6, pt.lwd=c(0.7,0.9,1.8),
    pt.bg=c("grey85",adjustcolor("grey70",0.4),"white"), col=c("grey40","grey55","darkorange2"),
    legend=c("functional","remnant (dotted)","C9B relaxed (orange)"))
  legend(x=dep*1.68, y=yl, title="Events", bty="n", cex=0.58, xpd=NA, pch=c(23,4), pt.bg=c("#4C4C4C",NA),
    col=c("black","#B2182B"), pt.cex=c(1.2,1), legend=c("duplication (+C9B/+C9C)","loss / pseudogenisation"))
}

png("figure_1/figure1_c9_status_scaled.png", width=1420, height=1400, res=150); draw(); invisible(dev.off())
pdf("figure_1/figure1_c9_status_scaled.pdf", width=8.8, height=9.3);           draw(); invisible(dev.off())
cat("wrote figure_1/figure1_c9_status_scaled.{png,pdf}\n")
