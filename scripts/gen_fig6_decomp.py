import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np

# --- Data from Table 4 (GB model, DELTAS section) ---
residues = ['ARG68', 'VAL110', 'ILE126', 'GLN129', 'HIS130', 'LYS133']

compounds = ['CID 101685610', 'AKBA (AF7)', 'CID 101663382', u'\u03b1-Amyrone', 'EPIC']

# Total free energy per residue per compound
total = {
    'ARG68':   [-0.57, -1.12, -1.08, -0.65, -1.36],
    'VAL110':  [-2.60, -1.85, -1.94, -1.50, -1.76],
    'ILE126':  [-1.18, -1.23, -1.30, -1.51, -1.66],
    'GLN129':  [-0.99, -1.18, -0.90, -0.90, -0.77],
    'HIS130':  [-1.56, -1.42, -1.43, -1.51, -1.23],
    'LYS133':  [-1.61, -2.21, -1.35, -1.27, -1.12],
}

# Energy components
vdw = {
    'ARG68':   [-1.19, -1.79, -1.79, -1.22, -1.57],
    'VAL110':  [-2.52, -1.87, -1.90, -1.47, -1.70],
    'ILE126':  [-1.37, -1.53, -1.50, -1.65, -1.83],
    'GLN129':  [-1.70, -1.72, -1.96, -1.58, -1.55],
    'HIS130':  [-2.59, -2.10, -2.22, -2.21, -1.93],
    'LYS133':  [-1.81, -1.90, -1.77, -1.54, -1.58],
}
ele = {
    'ARG68':   [0.19,  0.85,  1.12,  0.37,  0.35],
    'VAL110':  [-0.28, 0.22, -0.05,  0.03,  0.00],
    'ILE126':  [-0.08,-0.16, -0.15,  0.00, -0.06],
    'GLN129':  [-0.04,-0.10, -0.04,  0.07,  0.00],
    'HIS130':  [-0.26, 0.19, -0.02, -0.16,  0.02],
    'LYS133':  [0.87, -2.06,  0.43,  0.27,  0.66],
}
polar_solv = {
    'ARG68':   [0.64,  0.09, -0.13,  0.41,  0.07],
    'VAL110':  [0.53,  0.06,  0.26,  0.24,  0.24],
    'ILE126':  [0.44,  0.65,  0.53,  0.36,  0.48],
    'GLN129':  [0.98,  0.85,  1.33,  0.82,  0.95],
    'HIS130':  [1.45,  0.65,  0.95,  1.04,  0.83],
    'LYS133':  [-0.41, 2.06,  0.24,  0.24,  0.06],
}
nonpolar_solv = {
    'ARG68':   [-0.21, -0.28, -0.28, -0.20, -0.21],
    'VAL110':  [-0.33, -0.26, -0.26, -0.29, -0.29],
    'ILE126':  [-0.17, -0.19, -0.19, -0.22, -0.25],
    'GLN129':  [-0.23, -0.21, -0.23, -0.22, -0.17],
    'HIS130':  [-0.16, -0.15, -0.14, -0.18, -0.16],
    'LYS133':  [-0.25, -0.31, -0.24, -0.23, -0.26],
}

n_residues = len(residues)
n_compounds = len(compounds)
x = np.arange(n_residues)

# ============================================================
# PANEL A: Grouped bar chart — Total DG per residue
# ============================================================
fig, ax = plt.subplots(figsize=(10, 5.5), dpi=300)

bar_width = 0.15
offsets = np.arange(n_compounds) - (n_compounds - 1) / 2
colors = ['#2166AC', '#67A9CF', '#B2182B', '#EF8A62', '#F7C94C']

for i, comp in enumerate(compounds):
    vals = [total[r][i] for r in residues]
    ax.bar(x + offsets[i] * bar_width, vals, bar_width * 0.9,
           label=comp, color=colors[i], edgecolor='black', linewidth=0.4)

ax.set_xlabel('Residue', fontsize=12, fontweight='bold')
ax.set_ylabel(u'\u0394G (kcal/mol)', fontsize=12, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(residues, fontsize=10, fontweight='bold')
ax.tick_params(axis='y', labelsize=10)
ax.axhline(y=0, color='black', linewidth=0.5, linestyle='-')
ax.set_ylim(-3.2, 0.3)
ax.legend(fontsize=8.5, loc='lower right', framealpha=0.9, edgecolor='gray')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_title('(A) Per-Residue Total Binding Free Energy', fontsize=13, fontweight='bold', pad=12)

plt.tight_layout()
out_a = r'e:\Documents\DEGOO\amirona\dinamica-molecular\Figure_6A_decomp_total.png'
plt.savefig(out_a, dpi=300, bbox_inches='tight')
print(f'Saved: {out_a}')
plt.close()

# ============================================================
# PANEL B: Energy components for CID 101685610 vs AF7
# ============================================================
fig2, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300, sharey=True)

focus_idx = [0, 1]
focus_labels = ['(B) CID 101685610', '(C) AKBA (AF7)']
comp_labels = [u'\u0394E$_{vdW}$', u'\u0394E$_{ele}$', u'\u0394G$_{polar}$', u'\u0394G$_{nonpolar}$']
comp_colors = ['#4393C3', '#D6604D', '#92C5DE', '#F4A582']

for panel, ci in enumerate(focus_idx):
    ax2 = axes[panel]
    data_sets = [
        [vdw[r][ci] for r in residues],
        [ele[r][ci] for r in residues],
        [polar_solv[r][ci] for r in residues],
        [nonpolar_solv[r][ci] for r in residues],
    ]
    bar_w = 0.18
    off2 = np.arange(4) - 1.5

    for j, (vals, lab, col) in enumerate(zip(data_sets, comp_labels, comp_colors)):
        ax2.bar(x + off2[j] * bar_w, vals, bar_w * 0.9,
                label=lab, color=col, edgecolor='black', linewidth=0.4)

    ax2.set_title(focus_labels[panel], fontsize=11, fontweight='bold')
    ax2.set_xlabel('Residue', fontsize=10, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(residues, fontsize=9, fontweight='bold', rotation=30, ha='right')
    ax2.axhline(y=0, color='black', linewidth=0.5)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.tick_params(axis='y', labelsize=9)

axes[0].set_ylabel('Energy (kcal/mol)', fontsize=10, fontweight='bold')
axes[1].legend(fontsize=8, loc='upper right', framealpha=0.9, edgecolor='gray')

plt.tight_layout()
out_b = r'e:\Documents\DEGOO\amirona\dinamica-molecular\Figure_6BC_decomp_components.png'
plt.savefig(out_b, dpi=300, bbox_inches='tight')
print(f'Saved: {out_b}')
plt.close()

print('Done. Both panels generated.')
