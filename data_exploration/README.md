# Data Exploration — Evo-Hex

Exploratory analysis of alpha-helical secondary structure across the CATH S40
non-redundant subset. The debug run covers **100 structures**, yielding **1,461
helices** (1,146 alpha / 281 3-10 / 34 pi) and **28,970 residues** with zero
DSSP failures. Full-dataset results (≈8,288 structures) are pending.

---

## 1. Helix-type distribution

| Type | Count | % |
|------|------:|----:|
| Alpha (H) | 1,146 | 78.4% |
| 3-10 (G) | 281 | 19.2% |
| Pi (I) | 34 | 2.3% |

**Notes:**
- The dominance of alpha helices is expected for a mainly-alpha CATH class.
- Pi helices are rare but present; their distinct amino acid profile (I and V
  enriched) warrants attention.
- 3-10 helices appear predominantly at helix termini — consistent with their
  role as helix caps or short linkers.

**Open questions:**
- Does the 3-10/alpha ratio change with protein length or structural family?
- Are pi helices systematically associated with particular functional sites?

---

## 2. Residue propensity

Top helix formers (observed propensity > 1.0, ranked):

| AA | Obs. propensity | Theo. (Chou-Fasman) |
|----|----------------:|--------------------:|
| A  | 1.311 | 1.42 |
| L  | 1.263 | 1.21 |
| Q  | 1.183 | 1.11 |
| M  | 1.182 | 1.45 |
| E  | 1.181 | 1.51 |
| I  | 1.152 | 1.08 |
| W  | 1.135 | 1.08 |
| R  | 1.130 | 0.98 |

Strong helix breakers: **P** (0.420) and **G** (0.512).

**Notes:**
- A and L match the classical expectation. The relatively high observed
  propensity of R (vs. theoretical 0.98) may reflect salt-bridge stabilization
  at solvent-exposed positions.
- M observed propensity (1.18) is substantially below its Chou-Fasman value
  (1.45) — worth investigating whether this is a dataset bias or a real
  discrepancy in structural contexts.
- W appears more helix-compatible than classically assumed (1.135 vs. 1.08).

**Open questions:**
- Does the propensity ranking hold across helix-length bins (short / medium /
  long)?
- Is the R over-representation driven by coiled-coil or membrane-spanning
  contexts?

---

## 3. Proteome comparison (enrichment)

Residues significantly enriched in helices vs. the human proteome reference:

| AA | Enrichment |
|----|----------:|
| A  | 1.407 |
| L  | 1.396 |
| I  | 1.344 |
| Q  | 1.230 |
| E  | 1.246 |

Depleted in helices: **G** (0.472), **P** (0.350), **S** (0.647), **C**
(0.603).

**Notes:**
- G and P depletion is mechanistically clear (backbone rigidity / lack of
  side-chain contacts).
- S depletion is noteworthy — serine is abundant in loops and turns and
  disfavors the helical backbone geometry.
- C depletion likely reflects its preferential role in disulfide bridges and
  beta-sheet contexts in this dataset.

**Open questions:**
- How does cysteine depletion correlate with the redox environment of the
  structures?

---

## 4. Alpha vs. 3-10 helix composition

Largest differences in residue frequency between helix types:

| AA | Alpha (%) | 3-10 (%) | Difference |
|----|----------:|---------:|-----------:|
| I  | 7.38 | 2.60 | +4.78 (alpha) |
| L  | 13.52 | 9.54 | +3.97 (alpha) |
| A  | 9.81 | 6.18 | +3.62 (alpha) |
| P  | 1.76 | 8.24 | +6.49 (3-10) |
| D  | 4.47 | 8.03 | +3.55 (3-10) |
| S  | 4.62 | 7.81 | +3.19 (3-10) |

**Notes:**
- The enrichment of P, D, and S in 3-10 helices is consistent with their
  known role as helix-terminal or transition structures — P can initiate a
  helix cap, and D/S are common N-cap residues.
- I, L, A enrichment in alpha helices reflects the hydrophobic core packing
  typical of longer helices.

**Open questions:**
- Are the P-rich 3-10 helices exclusively at N-termini, or also C-terminal?
- Does D enrichment in 3-10 helices co-locate with N-cap hydrogen bond
  donors?

---

## 5. N-cap / C-cap preferences

The heatmap (cell-p11) shows residue preference at positions N1-N3 and C1-C3
after z-score normalization per residue.

**Notes (preliminary, 100 structures):**
- D and N show z-score enrichment at N1-N2 positions, consistent with the
  classical "N-cap motif" involving side-chain H-bonds to backbone NH groups.
- Hydrophobic residues (L, I, V) are enriched at C-cap positions, supporting
  the idea that the C-terminus of a helix is stabilized by hydrophobic
  contacts rather than H-bond donors.

**Open questions:**
- Is the D/N enrichment at N1 strong enough to be statistically significant
  in the full dataset?
- Does the C-cap hydrophobic preference vary between alpha and 3-10 helices?

---

## 6. Heptad repeat and hydrophobic periodicity

Shannon entropy is nearly uniform across all seven heptad positions
(H ≈ 4.15-4.18 bits), with L as the top residue at every position (~10-11%).

**Notes:**
- The lack of entropy differentiation between positions a/d (hydrophobic
  core) and b/c/f (solvent-exposed) is unexpected for coiled-coil sequences.
  It may indicate that the dataset is not coiled-coil-enriched, or that the
  heptad assignment in this subset is noisy.
- Mean hydrophobic moment (μH, Eisenberg) should be computed once the full
  dataset is available to test amphipathicity at scale.

**Open questions:**
- Are the heptad registers being assigned correctly for non-coiled-coil helices?
- Does entropy differentiation emerge in the full S40 dataset?

---

## 7. Amino acid composition by helix length

| Group | Short (4-9) | Medium (10-19) | Long (≥20) |
|-------|------------:|---------------:|-----------:|
| Hydrophobic | 47.5% | 46.8% | 41.3% |
| Polar | 22.1% | 22.5% | 25.6% |
| Charged+ | 13.7% | 14.7% | 12.3% |
| Charged- | 13.6% | 13.2% | 16.3% |
| Special (G/P) | 6.0% | 4.6% | 10.2% |

**Notes:**
- Hydrophobic frequency decreases as helices get longer, while Polar and
  Charged- increase. This is consistent with longer helices being more likely
  to span membrane interfaces or to be solvent-exposed over their length.
- The spike in Special (G/P) in long helices is surprising and may indicate
  kinks or proline-induced bends in longer helical segments.

**Open questions:**
- Do long helices with high G/P frequency correspond to kinked or curved
  helices as assessed by geometric descriptors?
- Is the Charged- increase in long helices driven by glutamate (common in
  surface-exposed, solvent-accessible helices)?

---

## 8. UMAP — amino acid composition space

The UMAP projection (cell-p12) embeds each helix as a 20-dimensional
amino acid composition vector (normalized, StandardScaler).

**Hypotheses to test:**
- Helices dominated by A or L should cluster tightly (high-frequency
  residues, low diversity).
- Pi helices (rare, I/V-enriched) may appear as an isolated island.
- 3-10 helices (P/D/S-enriched) should be separable from alpha helices if
  the composition signal is strong enough.

**Open questions:**
- Once the full dataset is available, do UMAP clusters correlate with CATH
  superfamily labels?
- Is the density relief showing genuine compositional archetypes or just
  the marginal distribution of a continuous composition space?

---

## 9. Codon degeneracy vs. helix propensity

Residues with high codon degeneracy (L: 6, R: 6, S: 6) show a range of
propensities, suggesting codon availability does not directly drive helix
composition.

**Notes:**
- A (4 codons, propensity 1.31) is the clearest case of a highly degenerate
  residue that is also a strong helix former.
- M (1 codon, propensity 1.18) is a notable outlier: a rare codon residue
  with above-average helix propensity.

**Open questions:**
- Is there a phylogenetic signal in codon usage that correlates with helix
  content across taxa represented in CATH?

---

## Status

| Stage | Status |
|-------|--------|
| Debug run (100 structures) | Done |
| Full S40 pipeline run | Running |
| Full dataset plots | Pending |
| Statistical tests (chi-square, KS) | Planned |
| Comparison across CATH topologies | Planned |
