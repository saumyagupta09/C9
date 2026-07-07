#!/usr/bin/env Rscript
# Figure 1: Status & selection of C9 across squamate representatives.
# Tree branches: RELAX x aBSREL ; faint clade bands behind the tree.
# Synteny boxes DAB2-C9A-C9B-C9C-FYB1: FILL = gene identity colour; STYLE = ORF status
#   (solid=functional, red=multiple stop codons, faded+dotted=remnant, "No remnant"=lost);
#   OUTLINE = per-copy RELAX (mustard relaxed / red intensified).
suppressMessages(library(ape))

tree <- read.tree("figure_1/rep_species.nwk"); tree <- ladderize(tree, right = FALSE)
Ntip <- length(tree$tip.label)
d <- read.csv("figure_1/c9_status_data.csv", stringsAsFactors = FALSE)
d <- d[match(tree$tip.label, d$species), ]

## gene identity colours
gid <- c(DAB2="slategray2", C9A="slateblue3", C9B="orchid3", C9C="plum2", FYB1="slategray2")
selcol <- c(relaxed="#D4A017", intensified="#D7191C", ns="grey55")

## selection groups
elapid   <- c("Pseudonaja_textilis","Notechis_scutatus","Naja_naja")
colubrid <- c("Hemorrhois_hippocrepis","Arizona_elegans","Pantherophis_obsoletus","Pituophis_catenifer","Thamnophis_elegans")
snakes   <- c(colubrid, elapid, "Crotalus_tigris","Crotalus_adamanteus","Protobothrops_mucrosquamatus")
sel_a <- function(sp,o) if(o!="functional") "ns" else if(sp%in%elapid)"intensified" else if(sp%in%colubrid)"relaxed" else "ns"
sel_b <- function(sp,o) if(o=="functional" && sp%in%snakes) "relaxed" else "ns"
sel_c <- function(sp,o) "ns"

branchcol <- function(relax, absrel) if(absrel=="yes") switch(relax, relaxed="#1A9850", intensified="#5E3C99", "#2C7FB8") else switch(relax, relaxed="#D4A017", intensified="#D7191C", "black")
tip_status <- setNames(d$relax, d$species)
tip_col    <- setNames(mapply(branchcol, d$relax, d$absrel), d$species)
getDesc <- function(n) if(n<=Ntip) n else unlist(lapply(tree$edge[tree$edge[,1]==n,2], getDesc))
edge_col <- sapply(tree$edge[,2], function(ch){ if(ch<=Ntip) return(tip_col[tree$tip.label[ch]])
  st<-unique(tip_status[tree$tip.label[getDesc(ch)]])
  if(length(st)==1 && st=="relaxed")"#D4A017" else if(length(st)==1 && st=="intensified")"#D7191C" else "grey35" })

## clades (Hemorrhois deliberately excluded from Colubridae -> stays white)
clades <- list(Gekkota=c("Asaccus_caudivolvulus","Correlophus_ciliatus"),
  Scincoidea=c("Carinascincus_ocellatus","Cryptoblepharus_egeriae"),
  "Boidae/Pythonidae"=c("Python_bivittatus","Candoia_aspera"),
  Colubridae=c("Arizona_elegans","Pantherophis_obsoletus","Pituophis_catenifer","Thamnophis_elegans"),
  Elapidae=c("Pseudonaja_textilis","Notechis_scutatus","Naja_naja"),
  Viperidae=c("Crotalus_tigris","Crotalus_adamanteus","Protobothrops_mucrosquamatus"),
  Anguimorpha=c("Dopasia_gracilis","Elgaria_multicarinata"),
  Iguania=c("Anolis_carolinensis","Pogona_vitticeps"),
  Lacertoidea=c("Zootoca_vivipara","Rhineura_floridana"))
clade_col <- c(Gekkota="#FFF3CC", Scincoidea="#EFEFEF", "Boidae/Pythonidae"="#EAE3F2",
  Colubridae="#DCEAF7", Elapidae="#FCE3D6", Viperidae="#E2F0DA", Anguimorpha="#FBE2EF",
  Iguania="#E8F3E0", Lacertoidea="#DFF1F6")

roundrect <- function(xc,yc,w,h,rx,ry,col,border="grey25",lwd=1,lty=1){
  BR<-cbind((xc+w/2-rx)+rx*cos(seq(-pi/2,0,len=6)),(yc-h/2+ry)+ry*sin(seq(-pi/2,0,len=6)))
  TR<-cbind((xc+w/2-rx)+rx*cos(seq(0,pi/2,len=6)),(yc+h/2-ry)+ry*sin(seq(0,pi/2,len=6)))
  TL<-cbind((xc-w/2+rx)+rx*cos(seq(pi/2,pi,len=6)),(yc+h/2-ry)+ry*sin(seq(pi/2,pi,len=6)))
  BL<-cbind((xc-w/2+rx)+rx*cos(seq(pi,3*pi/2,len=6)),(yc-h/2+ry)+ry*sin(seq(pi,3*pi/2,len=6)))
  polygon(rbind(BR,TR,TL,BL),col=col,border=border,lwd=lwd,lty=lty) }

dep <- max(node.depth.edgelength(tree)); xlab <- dep*1.03
bw <- dep*0.09; bh <- 0.60; rx <- bw*0.13; ry <- bh*0.16

draw <- function(){
  par(mar=c(8, 3.2, 3.5, 0.6))
  plot(tree, show.tip.label=FALSE, edge.color=edge_col, edge.width=2.8, x.lim=c(0, dep*2.42), y.lim=c(0.4, Ntip+1.3))
  lastPP <- get("last_plot.phylo", envir=.PlotPhyloEnv); tipy <- lastPP$yy[1:Ntip]

  ## clade labels only (no coloured bands)
  for(cl in names(clades)){ idx<-match(clades[[cl]],tree$tip.label)
    mtext(cl, side=2, at=mean(range(tipy[idx])), line=0.4, las=0, cex=0.55, col="grey35") }

  ## species names
  for(i in seq_len(Ntip)) text(xlab, tipy[i], gsub("_"," ",tree$tip.label[i]), adj=0, cex=0.72, font=3, xpd=NA)
  maxnw <- max(strwidth(gsub("_"," ",tree$tip.label), cex=0.72, font=3))
  cx <- (xlab + maxnw + dep*0.10) + (0:4)*dep*0.135

  draw_gene <- function(x,y,colid,status,sel){
    if(status=="absent") return(invisible())
    if(status=="lost"){ text(x,y,"No remnant",cex=0.44,font=3,col="grey45"); return() }
    if(status=="remnant"){ roundrect(x,y,bw,bh,rx,ry,col=adjustcolor(colid,0.40),border="grey55",lwd=0.9,lty=3); return() }
    if(status=="stops"){ roundrect(x,y,bw,bh,rx,ry,col="#E03531",border="grey30",lwd=0.7); return() }
    roundrect(x,y,bw,bh,rx,ry,col=colid,border=selcol[sel],lwd=if(sel=="ns")0.7 else 1.8) }
  for(i in seq_len(Ntip)){ y<-tipy[i]; sp<-tree$tip.label[i]
    present<-c(TRUE, d$c9a[i]!="absent", d$c9b[i]!="absent", d$c9c[i]!="absent", TRUE)
    segments(min(cx[present]),y,max(cx[present]),y,col="grey70",lwd=0.9)
    roundrect(cx[1],y,bw,bh,rx,ry,col=gid["DAB2"],border="grey45",lwd=0.7)
    draw_gene(cx[2],y,gid["C9A"],d$c9a[i],sel_a(sp,d$c9a[i]))
    draw_gene(cx[3],y,gid["C9B"],d$c9b[i],sel_b(sp,d$c9b[i]))
    draw_gene(cx[4],y,gid["C9C"],d$c9c[i],sel_c(sp,d$c9c[i]))
    roundrect(cx[5],y,bw,bh,rx,ry,col=gid["FYB1"],border="grey45",lwd=0.7) }
  text(cx, Ntip+0.9, c("DAB2","C9A","C9B","C9C","FYB1"), font=4, cex=0.8, col="black", xpd=NA)

  ## events
  mark <- function(node,txt){ e<-which(tree$edge[,2]==node); if(!length(e))return()
    xx<-(lastPP$xx[tree$edge[e,1]]+lastPP$xx[node])/2; yy<-lastPP$yy[node]
    points(xx,yy,pch=23,bg="#4C4C4C",col="black",cex=1.5); text(xx,yy+0.7,txt,cex=0.62,font=2) }
  mark(getMRCA(tree,snakes),"+C9B"); mark(getMRCA(tree,c("Arizona_elegans","Pituophis_catenifer")),"+C9C")
  for(sp in c("Asaccus_caudivolvulus","Dopasia_gracilis","Protobothrops_mucrosquamatus")){ i<-match(sp,tree$tip.label); points(dep*0.99,tipy[i],pch=4,col="#B2182B",lwd=2,cex=0.9) }

  mtext("Status and selection of C9 across squamate representatives", side=3, line=1.5, adj=0, cex=1.15, font=2, at=0)
  try(add.scale.bar(0,0.0,cex=0.55,lwd=2),silent=TRUE)

  ## legends (bottom margin)
  yl <- -2.2
  legend(x=0, y=yl, title="Branch (RELAX x aBSREL)", bty="n", cex=0.6, lwd=3, seg.len=1, xpd=NA,
    col=c("black","#D4A017","#D7191C","#2C7FB8","#1A9850","#5E3C99"),
    legend=c("not significant","relaxed","intensified","aBSREL +","aBSREL+relaxed","aBSREL+intensified"))
  legend(x=dep*0.78, y=yl, title="Gene (box fill = identity)", bty="n", cex=0.6, xpd=NA, pch=22, pt.cex=1.7,
    pt.bg=c("slategray2","slateblue3","orchid3","plum2"), col="grey40",
    legend=c("DAB2 / FYB1 (flanking)","C9A","C9B","C9C"))
  legend(x=dep*1.5, y=yl, title="ORF status", bty="n", cex=0.6, xpd=NA, pch=c(22,22,22,NA), pt.cex=1.7, pt.lwd=c(1,1,1,1),
    pt.bg=c("grey85","#E03531",adjustcolor("grey70",0.4),NA), col=c("grey40","grey30","grey55","white"),
    legend=c("functional","multiple stops","remnant (dotted)","'No remnant' = lost"))
  legend(x=dep*1.5, y=yl-Ntip*0.16, title="Outline = per-copy RELAX", bty="n", cex=0.6, xpd=NA, pch=22, pt.cex=1.7, pt.lwd=2.4,
    pt.bg="white", col=c("#D4A017","#D7191C"), legend=c("relaxed","intensified"))
  legend(x=dep*2.05, y=yl, title="Events", bty="n", cex=0.6, xpd=NA, pch=c(23,4), pt.bg=c("#4C4C4C",NA),
    col=c("black","#B2182B"), pt.cex=c(1.3,1), legend=c("duplication (+C9B/+C9C)","loss / pseudogenisation"))
}

png("figure_1/figure1_c9_status.png", width=1600, height=1450, res=150); draw(); invisible(dev.off())
pdf("figure_1/figure1_c9_status.pdf", width=10.6, height=9.6);          draw(); invisible(dev.off())
cat("wrote figure_1/figure1_c9_status.{png,pdf}\n")
