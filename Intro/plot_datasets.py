"""
Scatter plot of GRN datasets: genes perturbed × cells/samples,
with pie-chart markers (perturbation type makeup) and species-coded shapes.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Wedge
from matplotlib.transforms import Affine2D

# ── Curated dataset table (from TSV + primary sources) ──
# genes_measured convention:
#   - bulk microarray/RNA-seq: size of the profiled feature space (array probes / annotated genes)
#   - droplet scRNA-seq (Perturb-seq etc.): MEDIAN genes detected per cell, the honest
#     comparable metric (droplet capture detects only a few thousand of ~20k genes per cell),
#     not the reference-annotation feature space
#   - targeted panels: the panel size actually assayed
datasets = pd.DataFrame([
    # name (Author Year), genes_perturbed, genes_measured, total_cells, species, pert_makeup, year
    # ── Early bulk (model organisms) ──
    ("Gardner 2003", 9, 9, 99, "E. coli", {"KD": 1.0}, 2003),
    ("Gardner 2007", 50, 4297, 121, "E. coli", {"KO": 0.5, "OE": 0.5}, 2007),
    ("Parker 2019", 71, 4400, 71, "E. coli", {"KO": 1.0}, 2019),
    ("Kemmeren 2014", 1484, 6170, 1484, "S. cerevisiae", {"KO": 1.0}, 2014),
    ("Hackett 2020", 201, 6170, 1608, "S. cerevisiae", {"OE": 1.0}, 2020),
    ("Nishiyama 2009", 50, 25030, 340, "M. musculus", {"KD": 1.0}, 2009),
    ("De Cegli 2010", 32, 45101, 114, "M. musculus", {"OE": 1.0}, 2010),
    ("Lee 2016", 46, 15000, 138, "D. melanogaster", {"KD": 1.0}, 2016),
    # ── ENCODE bulk shRNA (whole-transcriptome bulk RNA-seq) ──
    ("ENCODE-fly", 23, 15000, 400, "D. melanogaster", {"KD": 1.0}, 2014),
    ("Graveley 2020", 232, 20000, 464, "H. sapiens", {"KD": 1.0}, 2020),
    ("ENCODE-HepG2 2020", 259, 20000, 518, "H. sapiens", {"KD": 1.0}, 2020),
    # ── Connectivity maps ──
    ("Subramanian 2017", 978, 978, 1000000, "H. sapiens", {"KD": 0.4, "OE": 0.3, "Compound": 0.3}, 2017),
    ("Feng 2020", 308, 20000, 570, "H. sapiens", {"KD": 0.5, "KO": 0.5}, 2020),
    # ── Single-cell perturbation screens (median genes/cell) ──
    ("Dixit 2016", 10, 3000, 200000, "H. sapiens", {"CRISPRi": 1.0}, 2016),
    ("Datlinger 2017", 33, 2500, 5905, "H. sapiens", {"CRISPRi": 1.0}, 2017),
    ("Jaitin 2016", 22, 2000, 5000, "M. musculus", {"CRISPRi": 1.0}, 2016),
    ("Adamson 2016", 82, 4000, 65000, "H. sapiens", {"CRISPRi": 1.0}, 2016),
    ("Rubin 2018", 12, 3000, 4300, "H. sapiens", {"CRISPRi": 1.0}, 2018),
    ("Tian 2019", 92, 2325, 100000, "H. sapiens", {"CRISPRi": 1.0}, 2019),
    ("Sunshine 2022", 183, 4000, 90380, "H. sapiens", {"CRISPRi": 1.0}, 2022),
    ("Replogle-Essential 2022", 2057, 5000, 310385, "H. sapiens", {"CRISPRi": 1.0}, 2022),
    ("Replogle-GWPS 2022", 11056, 5000, 1989578, "H. sapiens", {"CRISPRi": 1.0}, 2022),
    ("Joung 2023", 1836, 4000, 1000000, "H. sapiens", {"OE": 1.0}, 2023),
    ("Huang-HCT116 2025", 18903, 5387, 3409169, "H. sapiens", {"CRISPRi": 1.0}, 2025),
    ("Huang-HEK293T 2025", 18903, 5871, 4534299, "H. sapiens", {"CRISPRi": 1.0}, 2025),
    ("Wang 2025", 19000, 5000, 25600000, "H. sapiens", {"CRISPRi": 1.0}, 2025),
    ("10x-K562 2023", 567, 3200, 83943, "H. sapiens", {"CRISPRi": 1.0}, 2023),
    ("Jackson 2020", 11, 299, 11000, "S. cerevisiae", {"Barcode": 1.0}, 2020),
    # ── Large observational (non-perturbation) atlases ──
    ("GTEx 2020", 0, 17000, 17382, "H. sapiens", {"Observational": 1.0}, 2020),
    ("FANTOM5 2014", 0, 21000, 1829, "H. sapiens", {"Observational": 1.0}, 2014),
    ("Tabula Sapiens 2022", 0, 4000, 483152, "H. sapiens", {"Observational": 1.0}, 2022),
    ("Pijuan-Sala 2019", 0, 3500, 116312, "M. musculus", {"Observational": 1.0}, 2019),
    ("CELLxGENE 2023", 0, 4000, 50000000, "H. sapiens", {"Observational": 1.0}, 2023),
], columns=["name", "genes_perturbed", "genes_measured", "total_cells", "species", "pert_makeup", "year"])

# ── Visual encoding ──
species_markers = {
    "H. sapiens": "o",       # circle
    "M. musculus": "s",      # square
    "S. cerevisiae": "D",    # diamond
    "E. coli": "^",          # triangle
    "D. melanogaster": "P",  # plus (filled)
}

pert_colors = {
    "KD": "#4C72B0",        # blue - knockdown
    "KO": "#DD8452",        # orange - knockout
    "OE": "#55A868",        # green - overexpression
    "CRISPRi": "#C44E52",   # red - CRISPRi
    "Compound": "#8172B3",  # purple - chemical
    "Barcode": "#937860",   # brown - barcoding
    "Observational": "#64B5CD",  # light blue - no perturbation
}

# ── Size scaling ──
min_r, max_r = 6, 28
log_measured = np.log10(datasets["genes_measured"].values)
radii = min_r + (max_r - min_r) * (log_measured - log_measured.min()) / (log_measured.max() - log_measured.min())

# ── Label offsets to avoid overlap (manually tuned for known clusters) ──
label_offsets = {
    "Gardner 2003": (-14, -20),
    "Gardner 2007": (8, 4),
    "Parker 2019": (8, -14),
    "De Cegli 2010": (-20, -20),
    "Nishiyama 2009": (6, 6),
    "Lee 2016": (-52, -6),
    "ENCODE-fly": (8, -14),
    "Graveley 2020": (8, -12),
    "ENCODE-HepG2 2020": (8, 6),
    "Datlinger 2017": (8, -12),
    "Jaitin 2016": (-58, 6),
    "Subramanian 2017": (-90, 8),
    "Feng 2020": (8, 6),
    "Dixit 2016": (-48, -12),
    "Rubin 2018": (-46, 6),
    "Jackson 2020": (6, -12),
    "Huang-HCT116 2025": (-120, -16),
    "Huang-HEK293T 2025": (8, -4),
    "Wang 2025": (8, 8),
    "Joung 2023": (8, -14),
    "Sunshine 2022": (8, 6),
    "Adamson 2016": (8, 6),
    "Replogle-Essential 2022": (8, -12),
    "Replogle-GWPS 2022": (8, 6),
    "10x-K562 2023": (8, -12),
    "GTEx 2020": (5, -12),
    "FANTOM5 2014": (5, 6),
    "Tabula Sapiens 2022": (5, 6),
    "Pijuan-Sala 2019": (5, -12),
    "CELLxGENE 2023": (5, 6),
}

# ── Plot ──
fig, ax = plt.subplots(figsize=(16, 9))
ax.set_xscale("log")
ax.set_yscale("log")

for idx, row in datasets.iterrows():
    x = row["genes_perturbed"] if row["genes_perturbed"] > 0 else 0.8  # place at left edge for observational
    y = row["total_cells"]
    r = radii[idx]
    makeup = row["pert_makeup"]
    species = row["species"]
    marker = species_markers[species]

    # Draw species-shaped scatter marker (white background, neutral edge)
    ax.scatter(x, y, s=r**2, marker=marker, facecolors="white",
               edgecolors="#444444", linewidths=1.5, zorder=2)

    # Fill with perturbation color(s)
    if len(makeup) == 1:
        color = pert_colors[list(makeup.keys())[0]]
        ax.scatter(x, y, s=(r * 0.82)**2, marker=marker,
                   facecolors=color, edgecolors="none", alpha=0.7, zorder=3)
    else:
        # Multi-type: concentric rings from outside in, largest fraction outermost
        sorted_types = sorted(makeup.items(), key=lambda kv: kv[1], reverse=True)
        n = len(sorted_types)
        for i, (pert_type, frac) in enumerate(sorted_types):
            # Scale radius from full down to smallest ring
            ring_r = r * 0.82 * (1.0 - i * 0.28)
            ax.scatter(x, y, s=ring_r**2, marker=marker,
                       facecolors=pert_colors[pert_type],
                       edgecolors="white", linewidths=0.4,
                       alpha=0.8, zorder=3 + i)

    # Label with custom offset
    offset = label_offsets.get(row["name"], (5, 5))
    ax.annotate(row["name"], (x, y), fontsize=11,
                xytext=offset, textcoords="offset points",
                ha="left", va="bottom", alpha=0.9,
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.4)
                if abs(offset[0]) > 30 else None)


# ── Axis labels ──
ax.set_xlabel("Genes Perturbed", fontsize=16)
ax.set_ylabel("Total Cells / Samples", fontsize=16)
ax.tick_params(axis="both", labelsize=13)
ax.set_title("Datasets for GRN Inference\n"
             "(size = genes measured; shape = species; fill = perturbation type)",
             fontsize=16, pad=15)

# ── Legends (outside plot area) ──
# Species
species_handles = [plt.Line2D([0], [0], marker=m, color="gray", linestyle="",
                              markersize=10, markerfacecolor="lightgray", label=s)
                   for s, m in species_markers.items()]
leg1 = ax.legend(handles=species_handles, title="Species",
                 loc="upper left", bbox_to_anchor=(1.02, 1.0),
                 fontsize=12, title_fontsize=13, framealpha=0.9)
ax.add_artist(leg1)

# Perturbation type
pert_handles = [mpatches.Patch(facecolor=c, edgecolor="gray", label=k)
                for k, c in pert_colors.items()]
leg2 = ax.legend(handles=pert_handles, title="Perturbation Type",
                 loc="center left", bbox_to_anchor=(1.02, 0.55),
                 fontsize=12, title_fontsize=13, framealpha=0.9)
ax.add_artist(leg2)

# Size (genes measured)
gm_min, gm_max = log_measured.min(), log_measured.max()
for sz_val in [10, 100, 1000, 10000]:
    r_leg = min_r + (max_r - min_r) * (np.log10(sz_val) - gm_min) / (gm_max - gm_min)
    ax.scatter([], [], s=r_leg**2, c="lightgray", edgecolors="gray",
               label=f"{sz_val:,}")
leg3 = ax.legend(title="Genes Measured",
                 loc="lower left", bbox_to_anchor=(1.02, 0.15),
                 fontsize=12, title_fontsize=13, framealpha=0.9)
ax.add_artist(leg3)

ax.grid(True, alpha=0.3, which="both")
ax.set_xlim(0.5, 30000)
ax.set_ylim(50, 100000000)

plt.tight_layout()
plt.subplots_adjust(right=0.80)
plt.savefig("datasets_scatter.png", dpi=150)
plt.savefig("datasets_scatter.pdf")
print("Saved: datasets_scatter.png, datasets_scatter.pdf")
plt.show()
