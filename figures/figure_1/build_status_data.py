#!/usr/bin/env python3
# Build figure_1/c9_status_data.csv for the C9 status figure.
# Per copy ORF status: functional / stops (multiple stop codons) / remnant (partial) / lost / absent.
# RELAX (snakes_C9A per-clade) + aBSREL branch status per species.
import csv
cat={}
for r in csv.DictReader(open("Deliverables/tables/c9_all_copies_status_table.csv")):
    c=r['category'].split(' (')[0]
    if   c=="functional": st="functional"
    elif c=="lost":       st="lost"          # whole gene gone -> "No remnant"
    elif "remnant" in c:  st="remnant"       # partial remnant -> faded box + dotted line
    elif "degraded" in c or "pseudo" in c: st="stops"   # multiple stop codons -> red box
    else: st="fragment"
    copy=r['copy'].split()[0]
    if copy=="(none)": copy="C9A"            # whole-locus loss recorded on the C9A slot
    cat.setdefault(r['species'],{})[copy]=st

sp=["Asaccus_caudivolvulus","Correlophus_ciliatus","Carinascincus_ocellatus","Cryptoblepharus_egeriae",
"Python_bivittatus","Candoia_aspera","Hemorrhois_hippocrepis","Pseudonaja_textilis","Notechis_scutatus",
"Naja_naja","Arizona_elegans","Pantherophis_obsoletus","Pituophis_catenifer","Thamnophis_elegans",
"Crotalus_tigris","Crotalus_adamanteus","Protobothrops_mucrosquamatus","Dopasia_gracilis",
"Elgaria_multicarinata","Anolis_carolinensis","Pogona_vitticeps","Zootoca_vivipara","Rhineura_floridana"]
intensified={"Pseudonaja_textilis","Notechis_scutatus","Naja_naja"}                 # elapid C9A K=1.73
relaxed    ={"Hemorrhois_hippocrepis","Arizona_elegans","Pantherophis_obsoletus","Pituophis_catenifer","Thamnophis_elegans"}
absrel_yes ={"Naja_naja","Hemorrhois_hippocrepis"}
def stat(s,c):
    v=cat.get(s,{}).get(c,"absent")
    return v if v in ("functional","stops","remnant","lost") else "absent"
rows=[]
for s in sp:
    relax="intensified" if s in intensified else ("relaxed" if s in relaxed else "ns")
    rows.append(dict(species=s, c9a=stat(s,"C9A"), c9b=stat(s,"C9B"), c9c=stat(s,"C9C"),
                     relax=relax, absrel=("yes" if s in absrel_yes else "no")))
with open("figure_1/c9_status_data.csv","w",newline="") as f:
    w=csv.DictWriter(f, fieldnames=["species","c9a","c9b","c9c","relax","absrel"]); w.writeheader(); w.writerows(rows)
for r in rows: print(f"  {r['species']:30} A={r['c9a']:11} B={r['c9b']:11} C={r['c9c']:9} {r['relax']:11} absrel={r['absrel']}")
