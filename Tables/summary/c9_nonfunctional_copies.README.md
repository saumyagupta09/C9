# Non-functional C9 copies (`c9_nonfunctional_copies.csv`)

The 11 non-functional C9 copies (subset of `c9_all_copies_status_table.csv`), by clade.
Categories: **lost** (1), **degraded/pseudogene** (4), **fragment** (6). All other copies are functional.

- **lost** — no C9 anywhere in the genome (Asaccus; chromosome-level assembly => genuine, not a gap).
- **degraded/pseudogene** — a synteny-placed copy with a disrupted ORF.
    - Pantherophis C9B & Protobothrops C9B: **SRA-read confirmed genuine** (2026-07-02, see CORRECTION_LOG).
    - Dopasia C9A: frameshift(%3=1)+39 stops — signature resembles the Naja/Notechis assembly artifacts;
      **NOT yet read-validated** (candidate for re-check).
    - Anolis_sagrei C9A: C-terminal truncation (no stop), parallels A. carolinensis C9A.
- **fragment** — only a partial C9 hit on a secondary scaffold (assembly fragment), not a full copy.

Corrected AWAY from this list after read validation: Naja_naja C9B and Notechis_scutatus C9B
(were "pseudogene/remnant", now functional; assembly-indel/gap artifacts).

See `c9_alt_fragments_assembly_context.csv` (+README): all 6 fragments are on unplaced/non-chromosomal scaffolds; only Naja & Notechis have chromosome-level assemblies and both put the fragment on an unplaced scaffold => assembly artifacts, not independent copies.
