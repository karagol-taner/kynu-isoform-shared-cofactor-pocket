# Splice truncation in human kynureninase separates homodimer assembly from completion of the pyridoxal 5'-phosphate pocket

Taner Karagöl1,*,‡, Alper Karagöl2,3,*,‡

1McGovern Institute for Brain Research, Massachusetts Institute of Technology, Cambridge, MA, USA

2Broad Institute of MIT and Harvard, Cambridge, MA, USA

3Department of Mechanical Engineering, Massachusetts Institute of Technology, Cambridge, MA, USA

## Abstract

Human kynureninase (KYNU) hydrolyses 3-hydroxykynurenine on the route from tryptophan to NAD, and it works as an obligate homodimer in which pyridoxal 5'-phosphate lies at a pocket completed across the subunit boundary. The UniProt entry for KYNU annotates a second isoform, Q16719-2, in which the last 164 residues of the 465 residue canonical chain are replaced by a six residue sequence, giving a 307 residue product. We asked what that truncation does to the dimer, using the two wild type human crystal structures and values taken from coordinates. The isoform keeps the residues that carry most of the homodimer interface. Of the per subunit buried surface area, 71.3 per cent in 2HZP and 71.2 per cent in 3E9K comes from the retained block 1 to 301. The six residues the isoform replaces carry 8.1 per cent and the deleted span 308 to 465 carries 20.6 per cent. The partition holds across 90 combinations of surface area implementation, probe radius and burial cutoff, with the removed share between 27.8 and 29.4 per cent. The retained interface is a single contiguous patch of 58 residues spanning 47 to 301. The cofactor behaves in the opposite way. In the assembled dimer the partner subunit buries 58.59 square angstroms of the cofactor in 2HZP and 54.04 square angstroms in 3E9K. All of that burial comes from residues 302 to 465, with 0.00 square angstroms from the retained block. Truncating the partner at residue 301 returns the cofactor to the exposure it has in the isolated monomer. The same analysis on selenocysteine lyase, a family member with the same fold and an equivalent annotated truncation, buries 41.47 square angstroms entirely through its retained block, so the KYNU result is not a property of the family. AlphaFold3 predicts the canonical with isoform heterodimer at an ipTM of 0.89 compared with 0.91 for the canonical homodimer, over an interface the isoform forms from its retained block alone. Four tryptic peptides spanning the novel junction are unique to Q16719-2 across 42,545 reviewed human entries. Re-searching a liver dataset whose own search database contained the isoform entry returned no detectable diagnostic peptide in 443,786 MS/MS rows. Molecular dynamics of both dimers over 100 nanoseconds leaves the mixed dimer holding 0.85 of the canonical buried interface and 0.82 of its contacts. The isoform chain fluctuates 3.04 times as much as the same residues of a full length partner, at 296 of 296 positions. These values describe an isoform that can enter the dimer and bring no competent active site to it.

Keywords: Kynurenine pathway, Truncated isoforms, Homodimer interface, Pyridoxal 5'-phosphate, Complex structure prediction

## Repository layout

```
README.md        this file
codes/           the analysis, in the order below (30 files)
codes/colab/     the simulation notebooks (8 files)
data/            constructs, pairings, transitions and the peptide panel (10 files)
data/af3/        per prediction, the ranked model, its confidence file and the request (32 files)
results/         the derived tables behind every value in the paper (37 files)
results/dcs/     per residue identity, fluctuation and the two conserved scores (4 files)
results/md/      what the simulations measured, by series, frame and residue (33 files)
figures/         the figures, PNG at 400 dpi and editable SVG (18 files)
supplementary/   the supplementary information and the tables workbook (2 files)
```

## Dataset summary

| File | Rows | Description |
| --- | --- | --- |
| `data/constructs.fasta` | 7 | The canonical chain, the annotated isoform and every control construct, as the predictions were given them. |
| `data/construct_provenance.tsv` | 7 | Where each construct comes from, residue by residue, and the seed the permuted tails were drawn on. |
| `data/pairings.tsv` | 8 | The complexes predicted, one row per job. |
| `data/transitions.tsv` | 264 | Tryptic peptides of the canonical and isoform sequences with the transitions a targeted run would use. |
| `data/peptide_panel.tsv` | 20 | The junction peptides unique to the isoform. |
| `results/af3_scores.tsv` | 8 | Interface confidence of every predicted complex, with the interface residue count and both chain lengths. |
| `results/interface_residue_bsa.tsv` | 352 | Buried surface area per residue in each deposited entry, with the block each residue falls in. |
| `results/interface_parameter_sweep.tsv` | 90 | The interface partition under every combination of surface area implementation, probe radius and burial cutoff. |
| `results/contact_conservation.tsv` | 10 | Orthologue identity at the cofactor contact positions, with the residue each orthologue carries. |
| `results/per_position_conservation.tsv` | 465 | Orthologue identity at every position of the human chain, over the same orthologue set. |
| `results/variants_by_block.tsv` | 3 | Reported human variants, assigned to blocks. |
| `results/peptide_uniqueness.tsv` | 4 | Each junction peptide against the reviewed human proteome. |
| `results/wang2019_kynu_peptides.tsv` | 82 | What the re-search of the public liver and kidney datasets returned. |
| `results/dcs/canonical_sub1.tsv` | 465 | Per residue identity, fluctuation and the two conserved scores, one file per chain and build. |
| `results/md/series/canonical/sasa_dimer.xvg` | 2001 | The measured series of each run, one folder per build. |

## Reproducing the analysis

```bash
pip install biopython freesasa gemmi matplotlib numpy pillow pyteomics scipy

python3 codes/score_af3.py                 # interface confidence of every complex
python3 codes/interface_partition.py       # fetches each entry and partitions the interface
python3 codes/interface_parameter_sweep.py # the partition across 90 parameter combinations
python3 codes/cofactor_pocket.py           # the partner contribution to the cofactor pocket
python3 codes/own_subunit_pocket.py        # which site of a mixed dimer survives
python3 codes/phosphate_hbonds.py          # hydrogen bonds to the cofactor phosphate
python3 codes/variants_and_conservation.py # orthologue alignment and human variants
python3 codes/per_position_conservation.py # identity at every position
python3 codes/extract_peptides.py
python3 codes/peptide_uniqueness.py        # the junction peptides against the proteome
python3 codes/transcript_support.py
python3 codes/existence_evidence.py
python3 codes/public_spectra_search.py     # candidate projects in the public archive
python3 codes/triage_datasets.py           # the public proteomics search
python3 codes/dataset_searchability.py
python3 codes/sweep_search_spaces.py
python3 codes/search_public_results.py
python3 codes/search_maxquant_archives.py
python3 codes/md_setup.py                  # the simulation setup, read back out of the builds
python3 codes/compute_environment.py       # the hardware of the runs, from the mdrun output
python3 codes/dcs_analysis.py              # the dynamic and rigid conserved scores
python3 codes/render_core_struct.py        # the structural panels, needs PyMOL
python3 codes/render_md_struct.py          # the panels of the simulated systems, needs PyMOL
python3 codes/figure_core.py               # Figures 1 to 4
python3 codes/figure_md.py                 # Figure 5
python3 codes/figure_dcs.py                # Figure 6
python3 codes/figure_supp.py               # Supplementary Figures S1 to S3
```

Every script exits non-zero on a fatal condition instead of continuing.
Two folders are not committed here: `md/`, the CHARMM-GUI builds, written by the notebooks in `codes/colab/`, and `structures/`, the deposited entries, fetched from the Protein Data Bank. The values read from them are committed under `results/`. The
constructs, the pairings, the prediction inputs and the tryptic transitions are
committed under `data/`, so the analysis starts from them. The
simulations themselves are run from the notebooks in `codes/colab/`, in the order
`kynu_MD_1_prepare`, `kynu_MD_2_canonical` and `kynu_MD_2_mixed`, then
`kynu_MD_3_compare`, `kynu_MD_3b_cofactor_drift`, `kynu_MD_3c_frames` and
`kynu_MD_4_rmsf`. Each reads and writes a Google Drive folder and resumes from the
last chunk it wrote. Their measured output is committed under `results/md/`, so every
number and figure here can be reproduced without repeating the simulations. The card,
the thread counts, the compiled vector width and the throughput that gmx mdrun reported are
extracted from the executed notebooks into `data/md_run_log_canonical.tsv` and
`data/md_run_log_mixed.tsv`, which `codes/compute_environment.py` reads.

## Software dependencies

| Domain | Packages |
| --- | --- |
| Arrays and statistics | `numpy`, `scipy` |
| Plotting | `matplotlib`, `pillow` |
| Structures | `biopython`, `freesasa`, `gemmi` |
| Proteomics | `pyteomics` |
| Molecular dynamics | GROMACS 2026.3 and CHARMM-GUI, through the notebooks in `codes/colab/` |
| Rendering | PyMOL, for `codes/render_core_struct.py` and `codes/render_md_struct.py` |

Tested on Python 3.10.12.

## Data availability

* The crystal structures 2HZP, 3E9K and 7S3V for kynureninase and 3GZC and 3GZD for
  selenocysteine lyase are public depositions and are not redistributed here.
  `codes/interface_partition.py` downloads each one from the Protein Data Bank and
  checks it against the chain count and residue span the analyses were run against.
* Sequences and annotations are from UniProt, entry Q16719 for human kynureninase,
  and the orthologue set is the reviewed membership of InterPro family IPR010111.
  The sequences the alignment used are committed as
  `results/orthologue_sequences.json` so the alignment repeats without the network.
* Human variants are from the EBI Proteins API.
* The public proteomics datasets re-searched here are third-party depositions in
  PRIDE; what the search returned is committed under `results/`.
* Of each AlphaFold3 job, `data/af3/` carries the model its own ranking placed first,
  the confidence file that model is scored from, and the request the server returned.
  The remaining samples and the alignments are not needed to reproduce a value here.
  The job folders sit directly under `data/af3/` rather than under the dated folder
  the server delivered, so `models_root` in the committed `results/af3_scores.json`
  names that original folder while a rescore here writes `data/af3`. Every score the
  rescore produces is the same.

## Citation

If you use this code or dataset, please cite:

> Karagöl T., Karagöl A. *Splice truncation in human kynureninase separates homodimer assembly from completion of the pyridoxal 5'-phosphate pocket* (2026).

A BibTeX entry will be added upon publication of the journal record.

## Contact

Questions and issues are welcome, please open a GitHub Issue or contact the
corresponding authors.

## License

The code in this repository is released under the MIT License. Curated datasets and
figures are released under CC-BY 4.0. Protein Data Bank entries, UniProt records and
the third-party proteomics depositions are subject to the licenses of those upstream
providers.
