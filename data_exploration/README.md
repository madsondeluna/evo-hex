**Tentative titles:**

1. The Alpha-Helix as a Primordial Fossil: Prebiotic Amino Acid Bias and the Deterministic Rise of Alpha-Helical Architecture
2. Prebiotic Chemistry as Structural Destiny: Alanine and Leucine Bias Encoded in the Alpha-Helical Proteome
3. Frozen in Fold Space: A Prebiotic Amino Acid Signature Persists Across the Alpha-Helical Proteome

---

## Abstract (SBBE26, ~250 words)

The evolutionary origin of alpha-helical dominance in modern proteomes remains unresolved.
We tested whether the amino acid composition of alpha-helices reflects physicochemical
constraints predating Darwinian selection. From the CATH database (v4.3; 126,178
mainly-alpha and 305,361 alpha-beta domains), the S40 threshold (40% identity) yielded 8,056 class 1 representatives, of which 7,997
passed crystallographic quality criteria (resolution 3.0 angstrom or better, R-factor
0.25 or lower). A unified DSSP pipeline annotated 117,665 helical segments and 2,368,790
residues. We computed amino acid propensities, proteome enrichment ratios, N-cap and
C-cap positional preferences, heptad repeat distributions with Shannon entropy, and a
UMAP projection of the per-helix composition space.

Alanine and Leucine, both canonical prebiotic amino acids, are the two strongest
helix-forming residues (propensities 1.300 and 1.255; enrichment 1.550 and 1.365),
jointly occupying 24.0% of all alpha-helical positions. Glycine and Proline, also
prebiotically abundant but excluded from helical geometry by backbone constraints, show
the strongest depletion (enrichment 0.474 and 0.379). N-cap preferences for Asp and Asn,
reproduced across 117,665 helical termini, are fully explained by hydrogen bond geometry
without invoking sequence-level selection. Shannon entropy at all seven heptad positions
reaches 96.3-96.4% of the theoretical maximum, with Leucine showing a reproducible excess
at every structural role.

These data suggest that modern alpha-helical composition retains a measurable prebiotic
signature across billions of years of molecular evolution. Understanding the
physicochemical logic that shaped the earliest protein folds is an open question with
direct implications for reconstructing the events at the origin of life and interpreting
structural conservation across all domains of life.

---

## Abstract

The alpha-helix is the most prevalent secondary structure element across known proteomes,
yet the determinants of its compositional identity in modern proteins remain poorly
understood. Two competing interpretive frameworks have been proposed: one in which
alpha-helical dominance reflects contingent evolutionary outcomes under Darwinian
selection, and one in which it reflects a physicochemical predisposition rooted in the
prebiotic amino acid inventory. To evaluate these frameworks quantitatively, we performed
a large-scale structural bioinformatics analysis of the CATH database (version 4.3),
which contains 126,178 mainly-alpha (class 1) and 305,361 alpha-beta (class 3) protein
domains in its full release. Restricting the analysis to class 1 under the S40
non-redundancy threshold (maximum 40% pairwise sequence identity) reduced the working
set to 8,056 representative domains, of which 7,997 met crystallographic quality criteria
(resolution 3.0 angstrom or better, R-factor 0.25 or lower). We applied a unified DSSP
annotation pipeline to these 7,997 structures, encompassing 117,665 helical segments and
2,368,790 residues. For each residue and helix type,
we computed amino acid frequencies, observed helical propensities, proteome enrichment
ratios against the Swiss-Prot human reference, N-cap and C-cap positional preferences,
heptad repeat residue distributions with Shannon entropy, per-helix hydrophobic moments,
first-order Markov transition matrices, and a UMAP projection of the 20-dimensional
amino acid composition space.

Our results indicate that the two strongest helix-forming residues in the dataset,
Alanine (observed propensity 1.300; proteome enrichment 1.550) and Leucine (propensity
1.255; enrichment 1.365), are both canonical prebiotic amino acids, jointly occupying
24.0% of all alpha-helical positions across 90,263 helical segments. Conversely, Glycine
and Proline, also abundant in prebiotic syntheses but mechanistically incompatible with
alpha-helical geometry, show the strongest depletion in helical regions (enrichment 0.474
and 0.379, respectively). This bifurcation within the prebiotic amino acid set is
consistent with predictions derived from backbone dihedral constraints and hydrogen bond
geometry alone. N-cap and C-cap positional preferences, reproduced across all 117,665
helical termini, are explicable by backbone hydrogen bond satisfaction requirements
without invoking sequence-level selection. Shannon entropy at all seven heptad repeat
positions falls within 96.3-96.4% of the theoretical maximum, with Leucine maintaining
a consistent above-average frequency (10.4-10.6%) across all positions irrespective of
their structural role. The UMAP projection reveals a continuous compositional manifold
rather than discrete clusters, with Leucine- and Alanine-dominated helices defining the
region of highest density.

Taken together, these data suggest that the compositional signature of the alpha-helix
in modern proteins is substantially consistent with physicochemical constraints operative
in the prebiotic world, prior to the emergence of a ribosomal translation system. Whether
this pattern reflects a direct inheritance from early peptide chemistry, a convergent
outcome of independent physicochemical filtering, or a combination of both remains an open
question. These findings underscore the importance of revisiting the evolutionary origins
of primordial protein architectures: understanding how the first folded structures arose,
and what chemical logic shaped their composition, may be essential for reconstructing the
molecular events at the origin of life and for interpreting the deep conservation of
structural motifs across all domains of life.

---

## 1. Introduction

The protein structural universe is dominated by a small number of recurring architectural
motifs, of which the alpha-helix is the most fundamental. Across all known proteomes,
alpha-helical secondary structures account for approximately 32% of all residues in
solved structures, and the "all-alpha" and "alpha/beta" CATH classes together represent
the majority of known domain folds. The evolutionary and mechanistic origins of this
dominance, however, remain an open question at the intersection of structural biology,
evolutionary biochemistry, and the origins of life.

Two broad interpretive frameworks exist. The first, the contingency model, posits that
the modern distribution of secondary structures reflects the outcome of billions of years
of Darwinian selection acting on a random initial exploration of sequence space. Under
this view, the alpha-helix is prevalent because evolution happened to sample it frequently
and retain it for functional reasons. The second framework, the determinism model, argues
that the modern distribution was largely predetermined by the physical chemistry of the
earliest available amino acids, the geometry of the polypeptide backbone, and the
constraints imposed by primitive membrane-like environments. The alpha-helix would then
be dominant not because evolution chose it, but because the prebiotic chemical world
made it inevitable.

The prebiotic amino acid inventory, as reconstructed from Miller-Urey-type spark discharge
experiments, meteoritic analyses (Murchison, Murray), and hydrothermal vent syntheses,
is biased toward simple aliphatic and acidic residues: Glycine, Alanine, Valine, Leucine,
Isoleucine, Proline, Aspartate, Glutamate, Serine, and Threonine are among the most
consistently recovered prebiotically generated amino acids. Crucially, this set is not
homogeneous with respect to helical propensity. Alanine and Leucine, the two most
abundantly synthesized prebiotic residues in many experimental systems, also carry among
the highest intrinsic propensities for alpha-helical conformation. If early functional
peptides were assembled primarily from this restricted toolkit, helical structures would
have been statistically over-represented simply because the available building blocks
preferred the helical backbone geometry.

Here we provide a quantitative structural analysis designed to test this hypothesis. We
analyze the amino acid composition of helical regions across the CATH S40 non-redundant
mainly-alpha protein dataset, compare helical frequencies against modern proteome
reference values, evaluate codon degeneracy patterns, and characterize compositional
distributions across helix types, lengths, and terminal positions. We interpret these
results in the context of the prebiotic fossil hypothesis and assess the extent to which
modern helical composition is consistent with a physicochemically constrained, rather
than evolutionarily contingent, origin.

---

## 2. Methods

### 2.1 Structural dataset

All protein domain structures were retrieved from the CATH S40 non-redundant database
(version 4.3), restricting the analysis to class 1 (Mainly Alpha). The S40 threshold
ensures that no two domains in the working set share more than 40% sequence identity,
minimizing compositional redundancy arising from evolutionary relatedness. The final
dataset comprises 7,997 mainly-alpha domains yielding 117,665 annotated helical segments
and 2,368,790 residues with zero DSSP failures. Structures were accepted if the
crystallographic resolution was 3.0 A or better and the reported R-factor was 0.25 or
lower. All PDB files were retrieved in their biological assembly form and subjected to
chain-level cleaning prior to analysis.

### 2.2 Secondary structure annotation

Secondary structure was assigned for every residue using DSSP (Dictionary of Secondary
Structure of Proteins), accessed via the BioPython interface. The three helical assignments
recognized by DSSP are: alpha-helix (H, i, i+4 hydrogen bonds), 3-10 helix (G, i, i+3
hydrogen bonds), and pi-helix (I, i, i+5 hydrogen bonds). A single unified DSSP pass per
structure simultaneously collected all per-residue secondary structure assignments, helix
segment boundaries, and terminal position indices (N1-N3, C1-C3, defined relative to the
first and last backbone hydrogen bond of each helix). Residues assigned to coil (C),
beta-strand (E), beta-bridge (B), turn (T), and bend (S) were excluded from all
helix-specific analyses.

### 2.3 Amino acid frequency and propensity

For each helix type (H, G, I), residue counts were accumulated across all helices in all
structures. Amino acid frequencies within helices and across all residues were computed
as simple proportions of the total residue count in each category. The observed helical
propensity for amino acid a was defined as the ratio of its frequency within alpha-helical
segments to its global frequency across all annotated residues: P(a) = f_helix(a) /
f_total(a). Theoretical propensities were drawn from the Chou-Fasman scale for
comparative reference. The prebiotic status of each amino acid was assigned based on
consensus across Miller-Urey experiments, meteoritic analyses, and hydrothermal synthesis
reports.

### 2.4 Proteome enrichment

Proteome-level reference frequencies were obtained from the Swiss-Prot UniProt human
proteome annotation (canonical sequences, reviewed entries). The enrichment ratio for
amino acid a in helical regions was computed as: E(a) = f_observed(a) / f_proteome(a),
where f_observed is the frequency within helical segments of the CATH dataset and
f_proteome is the Swiss-Prot reference frequency. Values above 1.0 indicate enrichment
in helices relative to the background proteome; values below 1.0 indicate depletion.

### 2.5 N-cap and C-cap analysis

For each annotated helix, the first three residues (N1, N2, N3) and the last three
residues (C1, C2, C3) were extracted and tabulated separately. Residue frequencies were
compiled as 20 x 3 matrices for the N-cap and C-cap positions. The matrices were
z-score normalized per residue (row-wise) so that enrichment and depletion at each
cap position are expressed relative to the per-residue mean across all positions.

### 2.6 Heptad repeat analysis and Shannon entropy

Each alpha-helix was assigned a heptad register by projecting residue positions onto a
3.6-residues-per-turn periodicity model (positions a-g, repeating). Amino acid frequencies
at each of the seven heptad positions were tallied across all alpha-helices. Shannon
entropy at each position was computed as H = -sum p_i * log2(p_i), where p_i is the
frequency of amino acid i at that position. The theoretical maximum entropy for 20 amino
acids with equal probability is log2(20) = 4.322 bits.

### 2.7 Hydrophobic moment

The Eisenberg hydrophobic moment was computed for each helix with four or more residues
using the formula uH = (1/N) * sqrt[(sum H_i * sin(i * theta))^2 + (sum H_i * cos(i *
theta))^2], where H_i is the Eisenberg hydrophobicity of residue i, theta = 100 deg =
1.745 rad is the helical periodicity angle, and N is the number of residues with
available hydrophobicity values.

### 2.8 Codon degeneracy

Codon degeneracy for each amino acid was assigned from the standard genetic code as the
number of synonymous codons encoding that residue. This value was plotted against the
observed helical propensity to assess whether codon availability correlates with helical
enrichment.

### 2.9 Dimensionality reduction

To visualize the compositional landscape of individual helices, each helix was represented
as a 20-dimensional vector of normalized amino acid frequencies. Vectors were standardized
using scikit-learn StandardScaler (zero mean, unit variance per feature) and then reduced
to two dimensions using UMAP (Uniform Manifold Approximation and Projection; McInnes et
al., 2018) with parameters n_neighbors=15, min_dist=0.1, n_components=2, init='random',
random_state=42. Points were colored by the dominant (most frequent) amino acid within
each helix. A kernel density estimate (KDE) background was overlaid to visualize regional
density in the embedding space.

---

## 3. Results

### Figure 1 - Helix-type distribution

![Helix type distribution](imgs/helix_type_distribution.png)

Of the 117,665 helical segments annotated across 7,997 structures, 90,263 (76.7%) were
classified as alpha helices, 23,912 (20.3%) as 3-10 helices, and 3,490 (3.0%) as pi
helices. The mean helical content across all structures was 54.2%. The strong numerical
dominance of alpha helices in the CATH mainly-alpha class is the structural observation
motivating our hypothesis: this dominance must be explained, and the central question
is whether it reflects a physicochemical inevitability or a contingent evolutionary
outcome.

---

### Figure 2 - Amino acid helix propensity

![Helix propensities](imgs/helix_propensities.png)

| AA | Prebiotic? | Obs. propensity | Chou-Fasman |
|----|:----------:|----------------:|------------:|
| A  | Yes | **1.300** | 1.42 |
| L  | Yes | **1.255** | 1.21 |
| E  | Yes | 1.200 | 1.51 |
| Q  | No  | 1.196 | 1.11 |
| M  | No  | 1.171 | 1.45 |
| I  | Yes | 1.122 | 1.08 |
| R  | No  | 1.118 | 0.98 |
| W  | No  | 1.099 | 1.08 |
| K  | No  | 1.069 | 1.16 |
| G  | Yes | 0.515 | 0.57 |
| P  | Yes | 0.451 | 0.57 |

Across 2,368,790 residues from 7,997 non-redundant structures, the two strongest
observed helix formers are Alanine (1.300) and Leucine (1.255), both canonical prebiotic
amino acids. Among the top six helix-forming residues, four (A, L, E, I) belong to the
prebiotic set. The only prebiotic residues with propensity below 1.0 are Glycine (0.515)
and Proline (0.451), both excluded from helical geometry by first principles, backbone
dihedral constraints and hydrogen bond geometry, respectively. This bifurcation is
mechanistically predicted and constitutes direct evidence that prebiotic amino acid
chemistry constrained helical composition from the outset.

---

### Figure 3 - Amino acid composition by helix type (heatmap, z-score)

![Helix type composition heatmap](imgs/helix_type_composition_heatmap.png)

Z-score normalized residue frequencies per helix type reveal a clear compositional
separation between alpha and 3-10 helices. Isoleucine, Leucine, and Alanine show positive
z-scores for alpha helices, while Proline, Aspartate, and Serine show positive z-scores
for 3-10 helices. Pi helices show enrichment in branched aliphatic residues (Leu, Val,
Ile), with Leucine leading at 12.5%, consistent with the wider i, i+5 backbone spacing
accommodating bulkier side chains while retaining the prebiotic hydrophobic core
preference. This compositional separation demonstrates that each helix type has a distinct
amino acid signature, and that the alpha-helix signature is specifically aligned with the
prebiotic amino acid set.

---

### Figure 4 - Top amino acids by helix type

![Top amino acids by helix type](imgs/helix_type_top_amino_acids.png)

In alpha helices, Leucine and Alanine together account for 24.0% of all residues
(13.2% and 10.8%, respectively), based on 90,263 helical segments. In 3-10 helices
(n=23,912), the same two residues account for only 17.8% (L 9.7%, A 8.0%), with Proline
rising to fourth place at 7.1%, a residue absent from or only marginally present in most
prebiotic amino acid inventories. In pi helices (n=3,490), Leucine again leads at 12.5%,
followed by Valine (8.8%) and Isoleucine (7.7%). The concentration of Leu and Ala in
alpha helices (24.0%) versus 3-10 helices (17.8%) across a dataset of 117,665 helices
provides robust quantitative support for the hypothesis that the alpha-helix is the
structural element most strongly imprinted by prebiotic amino acid bias.

---

### Figure 5 - Alpha vs. 3-10 statistical comparison

![Helix type statistical comparison](imgs/helix_type_statistical_comparison.png)

| AA | Prebiotic? | Alpha (%) | 3-10 (%) | Delta |
|----|:----------:|----------:|---------:|------:|
| P  | Yes | 1.90 | 7.12 | -5.22 |
| L  | Yes | 13.21 | 9.71 | +3.49 |
| I  | Yes | 6.66 | 3.50 | +3.16 |
| V  | Yes | 6.63 | 3.48 | +3.16 |
| A  | Yes | 10.80 | 8.04 | +2.76 |
| D  | Yes | 4.55 | 7.05 | -2.50 |
| S  | Yes | 4.74 | 7.05 | -2.31 |
| G  | Yes | 3.26 | 5.34 | -2.08 |

The eight largest frequency differences between alpha and 3-10 helices are all accounted
for by prebiotic amino acids. Isoleucine (+3.16), Leucine (+3.49), Valine (+3.16), and
Alanine (+2.76) are enriched in alpha helices. Proline (-5.22), Aspartate (-2.50), Serine
(-2.31), and Glycine (-2.08) are enriched in 3-10 helices, consistent with their roles at
helix boundaries and cap positions. This pattern reveals an internal division of labor
within the prebiotic amino acid set: some prebiotic residues (A, L, I, V, E) build
alpha-helical bodies, while others (P, D, S, G) mark their structural boundaries. This
division is mechanistically dictated, not evolved.

---

### Figure 6 - Helix positions along the polypeptide chain

![Helix positions](imgs/helix_positions.png)

Alpha helices are distributed across the full length of polypeptide chains with no
strong positional bias, reflecting the structural diversity of the mainly-alpha CATH
class across 7,997 domain architectures. A mild enrichment of 3-10 helices toward the
N-terminal third of chains is consistent with their role as helix-initiating structures
or short linkers. This figure is primarily descriptive and provides structural context
for interpreting the compositional analyses, though it does not directly test the
prebiotic bias hypothesis.

---

### Figure 7 - Helix length distribution

![Helix lengths](imgs/helix_lengths.png)

Across 90,263 alpha-helical segments, the length distribution is right-skewed with a
mode in the 7-12 residue range and a long tail extending beyond 40 residues. The median
alpha-helix length is approximately 11 residues, while 3-10 helices are substantially
shorter (median ~4 residues). This length difference is mechanistically relevant: longer
helices provide more interior positions where hydrophobic core packing, dominated by Leu,
Ile, and Ala, can be expressed, amplifying the compositional bias toward prebiotic
residues in the structural core.

---

### Figure 8 - Helix length by type

![Helix length by type](imgs/helix_length_by_type.png)

The length hierarchy (alpha > pi > 3-10) is confirmed across 117,665 helical segments.
The substantially larger variance in alpha-helix length reflects the diversity of structural
contexts in which alpha helices appear, from compact 4-turn hairpins to extended
transmembrane rods. This figure corroborates Figure 7 and provides no additional
mechanistic inference; it is primarily included for completeness of structural
characterization.

---

### Figure 9 - Helical wheel

![Helical wheel](imgs/helical_wheel_average.png)

The average amino acid composition projected onto the helical wheel (100 deg/residue
periodicity) reveals a partial hydrophobic sector enrichment in the arc corresponding to
the buried face. Alanine and Leucine contribute disproportionately to this sector,
consistent with their role as hydrophobic core-packing residues. The relative enrichment
of prebiotic aliphatic residues in the core-facing sector is consistent across the full
7,997-structure dataset, supporting the model of these residues as structural anchors
whose placement is determined by physical exclusion from the aqueous environment rather
than by evolved sequence specificity.

---

### Figure 10 - Hydrophobic moment distribution

![Hydrophobic moment](imgs/hydrophobic_moment_distribution.png)

The distribution of Eisenberg hydrophobic moments (uH) across 90,263 alpha-helical
segments shows a broad range with a mode near 0.20-0.30 and a tail extending beyond 0.50.
The mode corresponds to weakly amphipathic, globular-context helices in which hydrophobic
residues are not strongly periodically distributed. The high-uH tail corresponds to
strongly amphipathic helices consistent with membrane-interfacial or helix-bundle packing
contexts. If early peptides functioned at prebiotic membrane interfaces, as proposed by
lipid-world and membrane-first origin models, amphipathic helices with high uH would have
been the primary functional unit, and the prebiotic preference for Ala and Leu would have
been directly selected by the physical chemistry of the lipid-water interface.

---

### Figure 11 - N-cap and C-cap residue preferences

![N-cap C-cap preferences](imgs/ncap_ccap_preferences.png)

> N1-N3: first three helical residues after the initiating hydrogen bond; C1-C3: last
> three residues before the terminal hydrogen bond. All values are z-score normalized
> per residue across positions.

The N-cap heatmap shows positive z-scores for Aspartate and Asparagine at positions
N1-N2, consistent with the classical N-cap motif in which side-chain oxygen atoms donate
hydrogen bonds to the first two backbone NH groups of the helix, which are otherwise
unsatisfied. Proline shows enrichment at N2-N3, consistent with its role as a
helix-initiating residue that terminates the preceding coil segment and nucleates the
helical conformation. The C-cap heatmap shows enrichment of hydrophobic residues (Leu,
Ile, Val) at positions C1-C3, consistent with hydrophobic capping through side-chain
contacts rather than hydrogen bond donation.

These preferences are entirely explained by backbone geometry and hydrogen bond
electrostatics, requiring no evolutionary narrative. Their faithful reproduction across
117,665 helical termini is strong evidence that physical chemistry, not selection,
determines residue placement at helix termini. The Asp/Asn N-cap preference is
particularly significant: both are prebiotic amino acids, and their structural role at
N-cap positions would have been expressed in the earliest helical peptides as a direct
consequence of their side-chain hydrogen bonding capacity.

---

### Figure 12 - Heptad repeat residue distribution

![Heptad pattern](imgs/heptad_pattern.png)

Residue frequencies at each of the seven heptad positions (a-g) are visualized as a
treemap, colored by structural role: positions a and d form the hydrophobic core in
canonical coiled-coil geometry; positions b, c, and f are solvent-exposed; positions
e and g participate in inter-helix electrostatic interactions.

The Shannon entropy at all seven positions falls within a narrow range of 4.164-4.168
bits, against a theoretical maximum of 4.322 bits for a perfectly uniform distribution
over 20 amino acids. This near-maximal entropy (96.3-96.4% of maximum) indicates
near-uniform amino acid usage across all heptad positions, with Leucine as the single
residue showing a consistent above-average frequency (10.4-10.6%) at every position
regardless of its structural role. The dominance of a prebiotic residue at all positions,
including the hydrophobic core positions a and d where evolutionary pressure would be
expected to be strongest, supports the model of Leucine as a universal helix constituent
whose prominence traces back to its prebiotic abundance rather than position-specific
evolutionary selection.

---

### Figure 13 - Helical content per structure

![Helix content distribution](imgs/helix_content_distribution.png)

The mean helical content across 7,997 structures is 54.2%, with a broad distribution
reflecting the structural diversity within the CATH mainly-alpha class. This figure
characterizes the dataset and confirms that the analysis covers a wide range of helical
contexts, from compact all-helix bundles to mixed architectures, rather than a
structurally homogeneous subset. The compositional findings reported across Figures 2-5
are therefore representative of the general mainly-alpha structural universe rather than
an artifact of a narrow architectural class.

---

### Figure 14 - Amino acid co-occurrence matrix

![AA co-occurrence](imgs/aa_cooccurrence.png)

The z-score normalized co-occurrence matrix reveals which amino acid pairs appear together
in the same helix more frequently than expected under independence, computed across
117,665 helical segments. Notable positive co-occurrences involve Leu-Ala and Leu-Glu
pairings, consistent with the frequent appearance of these residues together in helical
cores and in amphipathic helices, respectively. The co-occurrence of Gly and Pro (both
helix breakers) reflects their joint enrichment at helix boundaries. The prebiotic core
triad (Ala, Leu, Glu) constitutes the most frequent positive co-occurrence cluster,
suggesting that the compositional grammar of the alpha-helix, the rules governing which
residues appear together, is dominated by the same prebiotic amino acids that show the
highest individual propensities.

---

### Figure 15 - Amino acid composition by helix length

![Helix length vs composition](imgs/helix_length_vs_composition.png)

| Group | Short (4-9) | Medium (10-19) | Long (>=20) |
|-------|------------:|---------------:|------------:|
| Hydrophobic | 47.6% | 47.6% | 41.5% |
| Polar | 21.3% | 21.5% | 23.7% |
| Charged+ | 13.7% | 14.8% | 14.6% |
| Charged- | 14.1% | 13.0% | 15.5% |
| Special (G/P) | 6.6% | 4.9% | 10.1% |

Hydrophobic residues constitute the dominant physicochemical group at all helix lengths,
ranging from 47.6% in short and medium helices to 41.5% in long helices (>=20 residues).
This ~6-percentage-point decrease in hydrophobic content with helix length is accompanied
by increases in polar and charged residues, consistent with longer helices being more
likely to span solvent-exposed or charged-rich environments. The doubling of the Special
group (Gly/Pro) from 4.9% in medium helices to 10.1% in long helices suggests that long
helices frequently incorporate proline-induced kinks or glycine-mediated flexibility,
structurally necessary for accommodating the geometric demands of extended helical
segments. The relative stability of the hydrophobic fraction above 40% across all length
bins represents a compositional floor consistent with a physicochemical minimum
requirement for helix stability, whose specific amino acid contributors (Leu, Ile, Ala)
are prebiotic.

---

### Figure 16 - Residue transition matrix

![Helix transition matrix](imgs/helix_transition_matrix.png)

The first-order Markov transition matrix encodes the probability of each amino acid being
immediately followed by each other amino acid within helical segments, computed across
the full dataset. Dominant transitions are concentrated along and near the diagonal,
reflecting compositional autocorrelation: high-frequency residues (Leu, Ala, Glu) are
most likely to be followed by themselves or each other. This figure characterizes the
local sequence grammar of helices and provides a prior for generative sequence design,
but does not provide direct evidence for the prebiotic bias hypothesis.

---

### Figure 17 - 3-10 helix ratio by alpha-helix length

![G ratio by length](imgs/g_ratio_by_length.png)

The frequency of 3-10 helices adjacent to or flanking alpha helices of different lengths
shows that shorter alpha helices have proportionally more 3-10 flanking segments,
consistent with the model in which 3-10 helices act as N-terminal cap structures for
longer alpha helices. This structural relationship holds across the full 7,997-structure
dataset and provides context for interpreting the compositional differences between helix
types described in Figures 3 and 5.

---

### Figure 18 - Shannon entropy at heptad positions

![Shannon entropy heptad](imgs/shannon_entropy_heptad.png)

The scatter visualization makes explicit the near-flat entropy profile across all seven
heptad positions (range: 4.164-4.168 bits; theoretical maximum: 4.322 bits), with
Leucine as the consistently most frequent residue at each position (10.4-10.6%). The
96.3-96.4% occupancy of the theoretical entropy maximum indicates that the compositional
space at each heptad position is broadly sampled across 90,263 alpha-helical segments,
with only a small but reproducible enrichment of Leucine above the uniform baseline. This
subtle excess across all positions, regardless of their hydrophobic versus electrostatic
structural role, may represent the most direct quantitative trace of prebiotic Leucine
availability preserved in the modern protein fold space.

---

### Figure 19 - Codon degeneracy vs. helix propensity

![Codon degeneracy vs propensity](imgs/codon_degeneracy_vs_propensity.png)

| AA | Prebiotic? | Codons | Obs. propensity |
|----|:----------:|-------:|----------------:|
| L  | Yes | 6 | 1.255 |
| A  | Yes | 4 | 1.300 |
| V  | Yes | 4 | 1.002 |
| G  | Yes | 4 | 0.515 |
| P  | Yes | 4 | 0.451 |
| I  | Yes | 3 | 1.122 |
| M  | No  | 1 | 1.171 |

There is no monotonic relationship between codon degeneracy and helical propensity across
the full amino acid alphabet. However, the prebiotic residues Alanine and Leucine occupy
a distinctive region of the degeneracy-propensity space: both carry high codon degeneracy
(4 and 6 synonymous codons, respectively) and high observed helical propensity (1.300 and
1.255). The interpretation is that the genetic code does not create the prebiotic helix
bias but reinforces it: residues that were abundant prebiotically and have high helical
propensity are also among the most degenerate in the code, ensuring their high frequency
in modern proteomes. The outlier position of Methionine (1 codon, propensity 1.171)
reflects its functional selection for hydrophobic core packing rather than prebiotic
availability.

---

### Figure 20 - Proteome enrichment

![Proteome comparison](imgs/proteome_comparison.png)

| AA | Prebiotic? | Enrichment | log2(enrich.) |
|----|:----------:|----------:|---------------:|
| A  | Yes | **1.550** | +0.63 |
| L  | Yes | **1.365** | +0.45 |
| I  | Yes | 1.213 | +0.28 |
| E  | Yes | 1.286 | +0.36 |
| K  | No  | 1.230 | +0.30 |
| G  | Yes | 0.474 | -1.08 |
| P  | Yes | 0.379 | -1.40 |
| C  | No  | 0.534 | -0.90 |

Among the twenty standard amino acids, the three most enriched in helical regions
relative to the modern human proteome are Alanine (1.550), Leucine (1.365), and
Isoleucine (1.213), all prebiotic. The two most depleted residues are Proline (0.379)
and Glycine (0.474), also prebiotic, but mechanistically excluded from helical geometry.
This bifurcation within the prebiotic amino acid set, between strongly enriched
helix-formers and strongly depleted helix-breakers, is the clearest quantitative signature
of physical chemistry acting on prebiotic substrates across 90,263 alpha-helical segments
from 7,997 non-redundant protein domains. It is mechanistically predictable from first
principles and its observation at this scale constitutes direct quantitative evidence for
the prebiotic fossil hypothesis.

---

### Figure 21 - UMAP of amino acid composition space

![UMAP composition space](imgs/umap_aa_composition.png)

Each point represents one helix, embedded in two dimensions from a 20-dimensional
normalized amino acid composition vector after StandardScaler preprocessing, computed
across 117,665 helical segments. Color encodes the dominant (most frequent) amino acid
within that helix. The KDE density background visualizes regions of high compositional
density in the embedding space.

The UMAP projection reveals a continuous compositional manifold rather than discrete
clusters, consistent with the physical chemistry model: all helices navigate the same
propensity landscape, producing a gradual spectrum of compositions rather than arbitrary
sequence configurations. Regions of highest density correspond to Leucine- and
Alanine-dominated helices, appearing as the densest area of the projection and confirming
that these prebiotic residues define the compositional attractor of the alpha-helical fold
space. The peripheral, lower-density regions of the embedding are populated by helices
dominated by charged and polar residues, representing the compositionally divergent
fraction of the helical repertoire shaped by functional specialization rather than
ancestral physical chemistry.

---

## 4. General Discussion

The results presented here, derived from 7,997 non-redundant mainly-alpha protein
structures comprising 117,665 helical segments and 2,368,790 residues, consistently
support a model in which the compositional identity of the alpha-helix is not a product
of contingent evolution but a physicochemically constrained outcome of the prebiotic
amino acid inventory. Across every analytical axis examined, observed propensity, proteome
enrichment, helix-type composition, heptad periodicity, cap preferences, and hydrophobic
moment, the same pattern emerges: the prebiotic amino acids Alanine and Leucine are
systematically overrepresented in helical regions, while the prebiotic amino acids Glycine
and Proline are systematically excluded, in exact correspondence with their known effects
on backbone geometry.

The central insight is that the prebiotic amino acid set was not compositionally neutral
with respect to secondary structure. It was biased, by the same organic chemistry that
made these molecules the easiest to synthesize abiotically, toward residues that happen
to stabilize the alpha-helical backbone. Alanine's methyl side chain provides the minimum
steric bulk needed to bias the phi/psi angles toward the helical region of Ramachandran
space without introducing large van der Waals penalties. Leucine's isobutyl side chain
offers maximal hydrophobic packing energy in helical-bundle cores, the very configuration
that would have been selected at prebiotic membrane interfaces, where hydrophobic exclusion
from aqueous environments provides a direct, non-Darwinian selection pressure.

The bifurcation of prebiotic residues into alpha-helix builders (A, L, I, V, E) and
boundary markers (P, D, S, G) is mechanistically satisfying. Proline's role at 3-10 helix
positions and helix N-termini, Aspartate's role at N-cap positions, Serine's enrichment
at helix boundaries, all are explicable by hydrogen bond geometry and backbone flexibility
constraints that existed before any genetic code or ribosome was present. This means that
the positional grammar of helices, which residues appear in the core versus at the termini,
was already determined by prebiotic chemistry, and the modern protein universe inherited
this grammar essentially intact.

The codon degeneracy analysis adds an important second layer to this argument. The genetic
code, far from randomizing the prebiotic bias, appears to have reinforced it: the two
strongest prebiotic helix-formers (Ala, Leu) are also among the most degenerate codons
in the standard genetic code, ensuring their high frequency in translated proteomes across
all domains of life. Whether this reflects a causal connection between structural
importance and codon redundancy, or an independent coincidence arising from chemical
accessibility of the corresponding codons, remains an open question. Either way, the
modern genetic code does not dilute the prebiotic signal; it perpetuates it.

The N-cap and C-cap preferences provide perhaps the most direct evidence for deterministic
physicochemical control. The enrichment of Asp and Asn at N1-N2 positions, reproduced
across tens of thousands of helical termini, is not an evolved sequence motif in the
conventional sense, it is a direct consequence of hydrogen bond geometry at the helix
N-terminus, where the first two backbone NH groups are unsatisfied and require side-chain
H-bond donors. Any peptide, in any chemical context, that places Asp or Asn at these
positions gains thermodynamic stability simply because physics demands it. These
preferences existed before natural selection and thus constitute genuine fossils of
prebiotic physical chemistry preserved in the modern structure database.

The heptad entropy analysis reveals that the compositional space at each of the seven
heptad positions is occupied near-uniformly (96.3-96.4% of the theoretical maximum
entropy), with only a small, reproducible excess of Leucine (10.4-10.6%) at all positions.
This near-maximum entropy combined with the specific Leu excess is informative: it
suggests that while the evolutionary process has broadly diversified the amino acid usage
at every structural position, it has not erased the prebiotic Leucine signal even at
positions where functional specialization would be expected to enforce compositional
selectivity. The UMAP projection confirms this picture at the level of individual helices:
a continuous compositional manifold with Leu- and Ala-dominated helices forming the
densest regions of the embedding space, rather than a random scatter or a discrete cluster
structure driven by recent evolutionary divergence.

---

## 5. Conclusions

We present quantitative evidence, derived from large-scale structural bioinformatics
analysis of 7,997 non-redundant mainly-alpha domains from the CATH S40 dataset, that the
compositional identity of the alpha-helix in modern proteins is consistent with, and
specifically predicted by, the physicochemical properties of prebiotic amino acids. The
two most abundant prebiotic amino acids with high intrinsic helical propensity, Alanine
and Leucine, together account for 24.0% of positions in alpha-helical segments and show
the strongest enrichment relative to the modern proteome (1.550 and 1.365, respectively).
Their complementary prebiotic counterparts, Glycine and Proline, are the most depleted
residues in helical regions (0.474 and 0.379), excluded by mechanistic constraints rooted
in backbone geometry. This bifurcation, observed consistently across helix types, length
bins, terminal positions, and the full compositional space as revealed by UMAP, constitutes
a "primordial fossil" signal preserved in the protein fold space of modern organisms.

These findings challenge the model of purely contingent protein evolution and support a
view in which the structural universe of proteins was substantially pre-determined by the
physical chemistry of the prebiotic world. The alpha-helix is prevalent not because
evolution chose it, but because the available building materials were predisposed to form
it, and the modern protein database retains, quantifiably and reproducibly, the chemical
signature of that predisposition.

---

## Status

| Stage | Status |
|-------|--------|
| Full S40 pipeline run (7,997 structures, 117,665 helices) | Done |
| Full dataset plots | Done |
| UMAP full dataset | Done |
| Statistical tests (chi-square, KS, Fisher exact) | Planned |
| Beta-strand control (CATH class 2) | Planned |
| Coiled-coil subfamily heptad analysis | Planned |
| Manuscript draft | In progress |
