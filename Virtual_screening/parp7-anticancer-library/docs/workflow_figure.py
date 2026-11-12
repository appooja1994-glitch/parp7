import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(9.6, 11.2))
fig.patch.set_facecolor("white")
ink, body = "#1b1b1b", "#3a3a3a"
SRC_F, SRC_E = "#dce8f3", "#4a7fad"
CUR_F, CUR_E = "#e2eee2", "#5a8a5a"
OUT_F, OUT_E = "#4a6b8a", "#2e4a66"

def box(x, y, w, h, title, sub, fill, edge, tcol=None, fs=10.0, sfs=8.2):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h,
        boxstyle="round,pad=0.035,rounding_size=0.09",
        facecolor=fill, edgecolor=edge, linewidth=1.2, zorder=3))
    tc = tcol or ink
    if sub:
        ax.text(x, y + h*0.19, title, ha="center", va="center",
                fontsize=fs, fontweight="bold", color=tc, zorder=4)
        ax.text(x, y - h*0.20, sub, ha="center", va="center",
                fontsize=sfs, color=tc, linespacing=1.38, zorder=4)
    else:
        ax.text(x, y, title, ha="center", va="center",
                fontsize=fs, fontweight="bold", color=tc, zorder=4)

def arrow(y1, y2, x=0.0, color="#5a6b7a"):
    ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle="-|>",
        mutation_scale=14, linewidth=1.35, color=color, zorder=2))

W = 7.4
steps = [
 (10.30, 1.08, "NPASS 3.0", "204,023 natural products  •  1,048,756 activity records\n"
                            "8,764 targets, of which 1,658 are cancer cell lines", SRC_F, SRC_E),
 (8.85, 1.08, "Step 1.  Select cancer cell-line actives",
              "target is a cancer cell line  •  IC50, GI50, EC50, CC50, LC50 or ED50\n"
              "values normalised to µM  •  threshold 10 µM  •  most potent record kept",
              SRC_F, SRC_E),
 (7.40, 1.00, "Step 2.  Join activities to structures",
              "surviving compound identifiers matched to SMILES and InChI", SRC_F, SRC_E),
 (6.02, 1.00, "Step 3.  Merge any further sources",
              "deduplicated on full InChIKey  •  provenance recorded per compound",
              SRC_F, SRC_E),
 (4.58, 1.14, "Step 4.  Chemical curation",
              "sanitise  •  desalt  •  neutralise  •  organic elements only  •  deglycosylate\n"
              "MW 150–600 Da  •  at most one Lipinski violation\n"
              "PAINS flagged, not removed  •  InChIKey deduplication", CUR_F, CUR_E),
 (3.12, 1.00, "Step 5.  Three dimensional preparation",
              "ETKDGv3 embedding  •  MMFF94 optimisation\n"
              "Meeko PDBQT with explicit torsion tree", CUR_F, CUR_E),
]
for y, h, t, s, f, e in steps:
    box(0, y, W, h, t, s, f, e)
for i in range(len(steps) - 1):
    y1 = steps[i][0] - steps[i][1] / 2
    y2 = steps[i+1][0] + steps[i+1][1] / 2
    arrow(y1, y2)

arrow(steps[-1][0] - steps[-1][1]/2, 1.72 + 0.50)
box(0, 1.72, W, 1.00, "SCREENING LIBRARY",
    "final_molecules.smi  •  anticancer_final.sdf  •  pdbqt/\n"
    "annotations.tsv  •  provenance and attrition logs",
    OUT_F, OUT_E, tcol="white", fs=10.6, sfs=8.4)

arrow(1.72 - 0.50, 0.46 + 0.44)
box(0, 0.46, W, 0.88, "Molecular docking against the PARP7 catalytic domain model",
    "AutoDock Vina  •  NAD⁺ site defined by the RBN-2397 pose\n"
    "hits counter docked against PARP1, PARP2 and PARP12",
    "#eef1f4", "#8a9aa8", fs=9.8, sfs=8.0)

# benchmark actives, deliberately outside the flow
ax.add_patch(FancyBboxPatch((4.02, -0.12), 3.20, 1.26,
    boxstyle="round,pad=0.03,rounding_size=0.07",
    facecolor="#faf0f0", edgecolor="#b07a7a", linewidth=1.0,
    linestyle=(0, (4, 2)), zorder=3))
ax.text(5.62, 0.86, "Benchmark actives", ha="center", fontsize=8.6,
        fontweight="bold", color="#8a3b3b", zorder=4)
ax.text(5.62, 0.33, "built separately from ChEMBL\n"
                    "and BindingDB for Q7Z3E1,\nwith property matched decoys",
        ha="center", fontsize=7.4, color="#8a3b3b", linespacing=1.45, zorder=4)
ax.add_patch(FancyArrowPatch((4.02, 0.46), (W/2 + 0.04, 0.46),
    arrowstyle="-|>", mutation_scale=12, linewidth=1.0, color="#b07a7a", zorder=2))

ax.text(0, 11.18, "Construction of the anticancer natural product screening library",
        ha="center", fontsize=12.6, fontweight="bold", color=ink)
ax.text(0, 10.92, "for structure based virtual screening against PARP7",
        ha="center", fontsize=9.4, color=body)

ax.set_xlim(-4.6, 7.5)
ax.set_ylim(-0.28, 11.45)
ax.axis("off")
plt.tight_layout()
plt.savefig("workflow.png", dpi=220, bbox_inches="tight", facecolor="white")
print("wrote workflow.png")
