import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import json
from cycler import cycler
import os 

# Set up directories
cwd = os.getcwd()
file_dir =  os.path.dirname(cwd)
plots_dir = os.path.join(cwd, "plots")

with open(os.path.join(cwd, 'results.json'), 'r') as file:
        data = json.load(file)

clifford_results = data["clifford_shadows"]
pauli_results = data["pauli_shadows"]
big_clifford_results = data["big_clifford_shadows"]
big_pauli_results = data["big_pauli_shadows"]
vqe_1000_results = data["vqe_1000"]
vqe_10000_results = data["vqe_10000"]

molecule_list = pauli_results.keys()
bk_molecules = [m for m in molecule_list if m[-7:-5] == "BK"]
jw_molecules = [m for m in molecule_list if m[-7:-5] == "JW"]
minimal_basis = [m for m in molecule_list if m.split("_")[1] == "STO-3G"]
not_minimal_basis = [m for m in molecule_list if m not in minimal_basis] 

big_molecule_list = big_pauli_results.keys()
big_bk_molecules = [m for m in big_molecule_list if m[-7:-5] == "BK"]
big_jw_molecules = [m for m in big_molecule_list if m[-7:-5] == "JW"]
big_minimal_basis = [m for m in big_molecule_list if m.split("_")[1] == "STO-3G"]
big_not_minimal_basis = [m for m in big_molecule_list if m not in big_minimal_basis] 

clifford_cst_rel_error = []
clifford_cst_up_rel_error = []
for molecule in clifford_results.keys():
        exact_energy = complex(clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst_up"])
        clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        clifford_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

pauli_cst_rel_error = []
pauli_cst_up_rel_error = []
for molecule in pauli_results.keys():
        exact_energy = complex(pauli_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst_up"])
        pauli_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        pauli_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

big_pauli_cst_rel_error = []
big_pauli_cst_up_rel_error = []
big_pauli_cst_gc_rel_error = []
for molecule in big_pauli_results.keys():
        exact_energy = complex(big_pauli_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(big_pauli_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(big_pauli_results[molecule]["gs_nrg_tap_cst_up"])
        cst_gc_energy = complex(big_pauli_results[molecule]["gs_nrg_tap_cst_gc"])
        big_pauli_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        big_pauli_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))
        big_pauli_cst_gc_rel_error.append(np.abs(1 - (cst_gc_energy / cst_energy)))


big_clifford_cst_rel_error = []
big_clifford_cst_up_rel_error = []
big_clifford_cst_gc_rel_error = []
for molecule in big_clifford_results.keys():
        exact_energy = complex(big_clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(big_clifford_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(big_clifford_results[molecule]["gs_nrg_tap_cst_up"])
        cst_gc_energy = complex(big_clifford_results[molecule]["gs_nrg_tap_cst_gc"])
        big_clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        big_clifford_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))
        big_clifford_cst_gc_rel_error.append(np.abs(1 - (cst_gc_energy / cst_energy)))

################################################

vqe_1000_rel_error = []
for molecule in vqe_1000_results:
        exact_energy = complex(vqe_1000_results[molecule]["gs_nrg_tap"])
        vqe_energy = complex(vqe_1000_results[molecule]["gs_nrg_basic_vqe"])
        vqe_1000_rel_error.append(np.abs(1 - (vqe_energy / exact_energy)))

vqe_10000_rel_error = []
for molecule in vqe_10000_results:
        exact_energy = complex(vqe_10000_results[molecule]["gs_nrg_tap"])
        vqe_energy = complex(vqe_10000_results[molecule]["gs_nrg_basic_vqe"])
        vqe_10000_rel_error.append(np.abs(1 - (vqe_energy / exact_energy)))

################################################

jw_clifford_cst_rel_error = []
jw_clifford_cst_up_rel_error = []
for molecule in jw_molecules:
        exact_energy = complex(clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst_up"])
        jw_clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        jw_clifford_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

bk_clifford_cst_rel_error = []
bk_clifford_cst_up_rel_error = []
for molecule in bk_molecules:
        exact_energy = complex(clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst_up"])
        bk_clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        bk_clifford_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

jw_pauli_cst_rel_error = []
jw_pauli_cst_up_rel_error = []
for molecule in jw_molecules:
        exact_energy = complex(pauli_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst_up"])
        jw_pauli_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        jw_pauli_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

bk_pauli_cst_rel_error = []
bk_pauli_cst_up_rel_error = []
for molecule in bk_molecules:
        exact_energy = complex(pauli_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst_up"])
        bk_pauli_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        bk_pauli_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

jw_big_clifford_cst_rel_error = []
jw_big_clifford_cst_up_rel_error = []
for molecule in big_jw_molecules:
        exact_energy = complex(big_clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(big_clifford_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(big_clifford_results[molecule]["gs_nrg_tap_cst_up"])
        jw_big_clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        jw_big_clifford_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

bk_big_clifford_cst_rel_error = []
bk_big_clifford_cst_up_rel_error = []
for molecule in big_bk_molecules:
        exact_energy = complex(big_clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(big_clifford_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(big_clifford_results[molecule]["gs_nrg_tap_cst_up"])
        bk_big_clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        bk_big_clifford_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

jw_big_pauli_cst_rel_error = []
jw_big_pauli_cst_up_rel_error = []
for molecule in big_jw_molecules:
        exact_energy = complex(big_pauli_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(big_pauli_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(big_pauli_results[molecule]["gs_nrg_tap_cst_up"])
        jw_big_pauli_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        jw_big_pauli_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

bk_big_pauli_cst_rel_error = []
bk_big_pauli_cst_up_rel_error = []
for molecule in big_bk_molecules:
        exact_energy = complex(big_pauli_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(big_pauli_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(big_pauli_results[molecule]["gs_nrg_tap_cst_up"])
        bk_big_pauli_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        bk_big_pauli_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

#################################

mb_pauli_cst_rel_error = []
mb_pauli_cst_up_rel_error = []
for molecule in minimal_basis:
        exact_energy = complex(pauli_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst_up"])
        mb_pauli_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        mb_pauli_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

mb_clifford_cst_rel_error = []
mb_clifford_cst_up_rel_error = []
for molecule in minimal_basis:
        exact_energy = complex(clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst_up"])
        mb_clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        mb_clifford_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

nmb_pauli_cst_rel_error = []
nmb_pauli_cst_up_rel_error = []
for molecule in not_minimal_basis:
        exact_energy = complex(pauli_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(pauli_results[molecule]["gs_nrg_tap_cst_up"])
        nmb_pauli_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        nmb_pauli_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

nmb_clifford_cst_rel_error = []
nmb_clifford_cst_up_rel_error = []
for molecule in not_minimal_basis:
        exact_energy = complex(clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(clifford_results[molecule]["gs_nrg_tap_cst_up"])
        nmb_clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        nmb_clifford_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

mb_big_pauli_cst_rel_error = []
mb_big_pauli_cst_up_rel_error = []
for molecule in big_minimal_basis:
        exact_energy = complex(big_pauli_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(big_pauli_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(big_pauli_results[molecule]["gs_nrg_tap_cst_up"])
        mb_big_pauli_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        mb_big_pauli_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

mb_big_clifford_cst_rel_error = []
mb_big_clifford_cst_up_rel_error = []
for molecule in big_minimal_basis:
        exact_energy = complex(big_clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(big_clifford_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(big_clifford_results[molecule]["gs_nrg_tap_cst_up"])
        mb_big_clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        mb_big_clifford_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

nmb_big_pauli_cst_rel_error = []
nmb_big_pauli_cst_up_rel_error = []
for molecule in big_not_minimal_basis:
        exact_energy = complex(big_pauli_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(big_pauli_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(big_pauli_results[molecule]["gs_nrg_tap_cst_up"])
        nmb_big_pauli_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        nmb_big_pauli_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))

nmb_big_clifford_cst_rel_error = []
nmb_big_clifford_cst_up_rel_error = []
for molecule in big_not_minimal_basis:
        exact_energy = complex(big_clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(big_clifford_results[molecule]["gs_nrg_tap_cst"])
        cst_up_energy = complex(big_clifford_results[molecule]["gs_nrg_tap_cst_up"])
        nmb_big_clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))
        nmb_big_clifford_cst_up_rel_error.append(np.abs(1 - (cst_up_energy / cst_energy)))
        

#################################

# small shadows - minimal basis vs. not minimal basis

plt.rcParams.update({'axes.labelsize': 18, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 16})
mpl.rcParams['axes.prop_cycle'] = cycler(color=['b', 'g', 'r'])

xx_list = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',
           'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
xx = {}
for mol_idx in range(len(molecule_list)):
        xx[list(molecule_list)[mol_idx]] = xx_list[mol_idx]

fig, ax = plt.subplots(2,1)
mb = minimal_basis
nmb = not_minimal_basis
ax[0].scatter(mb, mb_pauli_cst_rel_error, c='b', marker='x')
ax[0].scatter(mb, mb_clifford_cst_rel_error, c='g', marker='o')
ax[1].scatter(nmb, nmb_pauli_cst_rel_error, c='b', marker='x')
ax[1].scatter(nmb, nmb_clifford_cst_rel_error, c='g', marker='o')
ax[0].set_xlabel("Molecule (Minimal Basis Set)")
ax[1].set_xlabel("Molecule (Non-Minimal Basis Set)")
xlabels0 = [xx[m] for m in mb]
xlabels1 = [xx[m] for m in nmb]
ax[0].set_xticklabels(xlabels0, rotation=0)
ax[1].set_xticklabels(xlabels1, rotation=0)
ax[0].set_ylabel("Relative Error")
ax[1].set_ylabel("Relative Error")
ax[0].grid(True)
ax[1].grid(True)
fig.set_figwidth(10.0)
fig.set_figheight(8.0)
fig.tight_layout(w_pad=0.0)
ax[0].legend(['Pauli CST', 'Clifford CST'])
ax[1].legend(["Pauli CST", "Clifford CST"])
plt.savefig(os.path.join(plots_dir, "small_shadows_mb_vs_nmb"))

# small shadows - jw vs. bk

# fig, ax = plt.subplots(2,1)
# jw = jw_molecules
# bk = bk_molecules
# ax[0].plot(jw, jw_pauli_cst_rel_error, jw, jw_clifford_cst_rel_error, marker='x', linestyle="None")
# ax[1].plot(bk, bk_pauli_cst_rel_error, bk, bk_clifford_cst_rel_error, marker='x', linestyle="None")
# ax[0].set_xlabel("Molecule (JW Encoding)")
# ax[1].set_xlabel("Molecule (BK Encoding)")
# xlabels0 = [xx[m] for m in jw]
# xlabels1 = [xx[m] for m in bk]
# ax[0].set_xticklabels(xlabels0, rotation=0)
# ax[1].set_xticklabels(xlabels1, rotation=0)
# ax[0].set_ylabel("Relative Error")
# ax[1].set_ylabel("Relative Error")
# ax[0].grid(True)
# ax[1].grid(True)
# fig.set_figwidth(10.0)
# fig.set_figheight(8.0)
# fig.tight_layout(w_pad=0.0)
# ax[0].legend(['Pauli CST', 'Clifford CST'])
# ax[1].legend(["Pauli CST", "Clifford CST"])
# plt.savefig(os.path.join(plots_dir, "small_shadows_jw_vs_bk"))

# small shadows vs. VQE

fig, ax = plt.subplots(2,1)
m = molecule_list
ax[0].scatter(m, vqe_1000_rel_error[:50], c='b', marker='x')
ax[0].scatter(m, pauli_cst_rel_error, c='g', marker='o')
ax[1].scatter(m, vqe_1000_rel_error[:50], c='b', marker='x')
ax[1].scatter(m, clifford_cst_rel_error, c='g', marker='o')
ax[0].set_xlabel("Molecule")
ax[1].set_xlabel("Molecule")
xlabels0 = [xx[s] for s in m]
ax[0].set_xticklabels(xlabels0, rotation=0)
ax[1].set_xticklabels(xlabels0, rotation=0)
ax[0].set_ylabel("Relative Error")
ax[1].set_ylabel("Relative Error")
ax[0].grid(True)
ax[1].grid(True)
fig.set_figwidth(10.0)
fig.set_figheight(8.0)
fig.tight_layout(w_pad=0.0)
ax[0].legend(['VQE', 'Pauli CST'])
ax[1].legend(["VQE", "Clifford CST"])
plt.savefig(os.path.join(plots_dir, "small_shadows_vs_vqe"))

# small shadows + up - minimal basis vs. not minimal basis

# fig, ax = plt.subplots(2,1)
# mb = minimal_basis
# nmb = not_minimal_basis
# ax[0].plot(mb, mb_pauli_cst_up_rel_error, mb, mb_clifford_cst_up_rel_error, marker='x', linestyle="None")
# ax[1].plot(nmb, nmb_pauli_cst_up_rel_error, nmb, nmb_clifford_cst_up_rel_error, marker='x', linestyle="None")
# ax[0].set_xlabel("Molecule (Minimal Basis Set)")
# ax[1].set_xlabel("Molecule (Non-Minimal Basis Set)")
# xlabels0 = [xx[m] for m in mb]
# xlabels1 = [xx[m] for m in nmb]
# ax[0].set_xticklabels(xlabels0, rotation=0)
# ax[1].set_xticklabels(xlabels1, rotation=0)
# ax[0].set_ylabel("Relative Error")
# ax[1].set_ylabel("Relative Error")
# ax[0].grid(True)
# ax[1].grid(True)
# fig.set_figwidth(10.0)
# fig.set_figheight(8.0)
# fig.tight_layout(w_pad=0.0)
# ax[0].legend(['Pauli CST + UP', 'Clifford CST + UP'])
# ax[1].legend(["Pauli CST + UP", "Clifford CST + UP"])
# plt.savefig(os.path.join(plots_dir, "small_shadows_up_mb_vs_nmb"))

# small shadows + up - jw vs. bk

# fig, ax = plt.subplots(2,1)
# jw = jw_molecules
# bk = bk_molecules
# ax[0].plot(jw, jw_pauli_cst_up_rel_error, jw, jw_clifford_cst_up_rel_error, marker='x', linestyle="None")
# ax[1].plot(bk, bk_pauli_cst_up_rel_error, bk, bk_clifford_cst_up_rel_error, marker='x', linestyle="None")
# ax[0].set_xlabel("Molecule (JW Encoding)")
# ax[1].set_xlabel("Molecule (BK Encoding)")
# xlabels0 = [xx[m] for m in jw]
# xlabels1 = [xx[m] for m in bk]
# ax[0].set_xticklabels(xlabels0, rotation=0)
# ax[1].set_xticklabels(xlabels1, rotation=0)
# ax[0].set_ylabel("Relative Error")
# ax[1].set_ylabel("Relative Error")
# ax[0].grid(True)
# ax[1].grid(True)
# fig.set_figwidth(10.0)
# fig.set_figheight(8.0)
# fig.tight_layout(w_pad=0.0)
# ax[0].legend(['Pauli CST + UP', 'Clifford CST + UP'])
# ax[1].legend(["Pauli CST + UP", "Clifford CST + UP"])
# plt.savefig(os.path.join(plots_dir, "small_shadows_up_jw_vs_bk"))

# small shadows - pauli vs. clifford

# fig, ax = plt.subplots(2,1)
# m = molecule_list
# ax[0].plot(m, pauli_cst_rel_error, m, clifford_cst_rel_error, marker='x', linestyle="None")
# ax[1].plot(m, pauli_cst_up_rel_error, m, clifford_cst_up_rel_error, marker='x', linestyle="None")
# ax[0].set_xlabel("Molecule")
# ax[1].set_xlabel("Molecule")
# xlabels0 = [xx[s] for s in m]
# ax[0].set_xticklabels(xlabels0, rotation=0)
# ax[1].set_xticklabels(xlabels0, rotation=0)
# ax[0].set_ylabel("Relative Error")
# ax[1].set_ylabel("Relative Error")
# ax[0].grid(True)
# ax[1].grid(True)
# fig.set_figwidth(10.0)
# fig.set_figheight(8.0)
# fig.tight_layout(w_pad=0.0)
# ax[0].legend(['Pauli CST', 'Clifford CST'])
# ax[1].legend(["Pauli CST + UP", "Clifford CST + UP"])
# plt.savefig(os.path.join(plots_dir, "small_shadows_pauli_vs_clifford"))

# big shadows - minimal basis vs. not minimal basis

fig, ax = plt.subplots(2,1)
mb = big_minimal_basis
nmb = big_not_minimal_basis
ax[0].scatter(mb, mb_big_pauli_cst_rel_error, c='b', marker='x')
ax[0].scatter(mb, mb_big_clifford_cst_rel_error, c='g', marker='o')
ax[1].scatter(nmb, nmb_big_pauli_cst_rel_error, c='b', marker='x')
ax[1].scatter(nmb, nmb_big_clifford_cst_rel_error, c='g', marker='o')
ax[0].set_xlabel("Molecule (Minimal Basis Set)")
ax[1].set_xlabel("Molecule (Non-Minimal Basis Set)")
xlabels0 = [xx[m] for m in mb]
xlabels1 = [xx[m] for m in nmb]
ax[0].set_xticklabels(xlabels0, rotation=0)
ax[1].set_xticklabels(xlabels1, rotation=0)
ax[0].set_ylabel("Relative Error")
ax[1].set_ylabel("Relative Error")
ax[0].grid(True)
ax[1].grid(True)
fig.set_figwidth(10.0)
fig.set_figheight(8.0)
fig.tight_layout(w_pad=0.0)
ax[0].legend(['Pauli CST', 'Clifford CST'])
ax[1].legend(["Pauli CST", "Clifford CST"])
plt.savefig(os.path.join(plots_dir, "big_shadows_mb_vs_nmb"))

# big shadows - jw vs. bk

fig, ax = plt.subplots(2,1)
jw = big_jw_molecules
bk = big_bk_molecules
ax[0].scatter(jw, jw_big_pauli_cst_rel_error, c='b', marker='x')
ax[0].scatter(jw, jw_big_clifford_cst_rel_error, c='g', marker='o')
ax[1].scatter(bk, bk_big_pauli_cst_rel_error, c='b', marker='x')
ax[1].scatter(bk, bk_big_clifford_cst_rel_error, c='g', marker='o')
ax[0].set_xlabel("Molecule (JW Encoding)")
ax[1].set_xlabel("Molecule (BK Encoding)")
xlabels0 = [xx[m] for m in jw]
xlabels1 = [xx[m] for m in bk]
ax[0].set_xticklabels(xlabels0, rotation=0)
ax[1].set_xticklabels(xlabels1, rotation=0)
ax[0].set_ylabel("Relative Error")
ax[1].set_ylabel("Relative Error")
ax[0].grid(True)
ax[1].grid(True)
fig.set_figwidth(10.0)
fig.set_figheight(8.0)
fig.tight_layout(w_pad=0.0)
ax[0].legend(['Pauli CST', 'Clifford CST'])
ax[1].legend(["Pauli CST", "Clifford CST"])
plt.savefig(os.path.join(plots_dir, "big_shadows_jw_vs_bk"))

# big shadows vs. VQE

fig, ax = plt.subplots(2,1)
m = big_molecule_list
ax[0].scatter(m, vqe_10000_rel_error[:27], c='b', marker='x')
ax[0].scatter(m, big_pauli_cst_rel_error, c='g', marker='o')
ax[1].scatter(m, vqe_10000_rel_error[:27], c='b', marker='x')
ax[1].scatter(m, big_clifford_cst_rel_error, c='g', marker='o')
ax[0].set_xlabel("Molecule")
ax[1].set_xlabel("Molecule")
xlabels0 = [xx[s] for s in m]
ax[0].set_xticklabels(xlabels0, rotation=0)
ax[1].set_xticklabels(xlabels0, rotation=0)
ax[0].set_ylabel("Relative Error")
ax[1].set_ylabel("Relative Error")
ax[0].grid(True)
ax[1].grid(True)
fig.set_figwidth(10.0)
fig.set_figheight(8.0)
fig.tight_layout(w_pad=0.0)
ax[0].legend(['VQE', 'Pauli CST'])
ax[1].legend(["VQE", "Clifford CST"])
plt.savefig(os.path.join(plots_dir, "big_shadows_vs_vqe"))

# big shadows + up - minimal basis vs. not minimal basis

# fig, ax = plt.subplots(2,1)
# mb = big_minimal_basis
# nmb = big_not_minimal_basis
# ax[0].plot(mb, mb_big_pauli_cst_up_rel_error, mb, mb_big_clifford_cst_up_rel_error, marker='x', linestyle="None")
# ax[1].plot(nmb, nmb_big_pauli_cst_up_rel_error, nmb, nmb_big_clifford_cst_up_rel_error, marker='x', linestyle="None")
# ax[0].set_xlabel("Molecule (Minimal Basis Set)")
# ax[1].set_xlabel("Molecule (Non-Minimal Basis Set)")
# xlabels0 = [xx[m] for m in mb]
# xlabels1 = [xx[m] for m in nmb]
# ax[0].set_xticklabels(xlabels0, rotation=0)
# ax[1].set_xticklabels(xlabels1, rotation=0)
# ax[0].set_ylabel("Relative Error")
# ax[1].set_ylabel("Relative Error")
# ax[0].grid(True)
# ax[1].grid(True)
# fig.set_figwidth(10.0)
# fig.set_figheight(8.0)
# fig.tight_layout(w_pad=0.0)
# ax[0].legend(['Pauli CST + UP', 'Clifford CST + UP'])
# ax[1].legend(["Pauli CST + UP", "Clifford CST + UP"])
# plt.savefig(os.path.join(plots_dir, "big_shadows_up_mb_vs_nmb"))

# big shadows + up - jw vs. bk

# fig, ax = plt.subplots(2,1)
# jw = big_jw_molecules
# bk = big_bk_molecules
# ax[0].plot(jw, jw_big_pauli_cst_up_rel_error, jw, jw_big_clifford_cst_up_rel_error, marker='x', linestyle="None")
# ax[1].plot(bk, bk_big_pauli_cst_up_rel_error, bk, bk_big_clifford_cst_up_rel_error, marker='x', linestyle="None")
# ax[0].set_xlabel("Molecule (JW Encoding)")
# ax[1].set_xlabel("Molecule (BK Encoding)")
# xlabels0 = [xx[m] for m in jw]
# xlabels1 = [xx[m] for m in bk]
# ax[0].set_xticklabels(xlabels0, rotation=0)
# ax[1].set_xticklabels(xlabels1, rotation=0)
# ax[0].set_ylabel("Relative Error")
# ax[1].set_ylabel("Relative Error")
# ax[0].grid(True)
# ax[1].grid(True)
# fig.set_figwidth(10.0)
# fig.set_figheight(8.0)
# fig.tight_layout(w_pad=0.0)
# ax[0].legend(['Pauli CST + UP', 'Clifford CST + UP'])
# ax[1].legend(["Pauli CST + UP", "Clifford CST + UP"])
# plt.savefig(os.path.join(plots_dir, "big_shadows_up_jw_vs_bk"))

# big shadows - pauli vs. clifford

# fig, ax = plt.subplots(2,1)
# m = big_molecule_list
# ax[0].scatter(m, big_pauli_cst_rel_error, c='b', marker='x')
# ax[0].scatter(m, big_clifford_cst_rel_error, c='g', marker='o')
# ax[0].scatter(m, big_pauli_cst_rel_error, c='b', marker='x')
# ax[0].scatter(m, big_clifford_cst_rel_error, c='g', marker='o')
# ax[1].plot(m, big_pauli_cst_up_rel_error, m, big_clifford_cst_up_rel_error, marker='x', linestyle="None")
# ax[0].set_xlabel("Molecule")
# ax[1].set_xlabel("Molecule")
# xlabels0 = [xx[s] for s in m]
# ax[0].set_xticklabels(xlabels0, rotation=0)
# ax[1].set_xticklabels(xlabels0, rotation=0)
# ax[0].set_ylabel("Relative Error")
# ax[1].set_ylabel("Relative Error")
# ax[0].grid(True)
# ax[1].grid(True)
# fig.set_figwidth(10.0)
# fig.set_figheight(8.0)
# fig.tight_layout(w_pad=0.0)
# ax[0].legend(['Pauli CST', 'Clifford CST'])
# ax[1].legend(["Pauli CST + UP", "Clifford CST + UP"])
# plt.savefig(os.path.join(plots_dir, "big_shadows_pauli_vs_clifford"))

# small shadows - cst vs. cst + up

# fig, ax = plt.subplots(2,1)
# m = molecule_list
# ax[0].plot(m, pauli_cst_rel_error, m, pauli_cst_up_rel_error, marker='x', linestyle="None")
# ax[1].plot(m, clifford_cst_rel_error, m, clifford_cst_up_rel_error, marker='x', linestyle="None")
# ax[0].set_xlabel("Molecule")
# ax[1].set_xlabel("Molecule")
# xlabels0 = [xx[s] for s in m]
# ax[0].set_xticklabels(xlabels0, rotation=0)
# ax[1].set_xticklabels(xlabels0, rotation=0)
# ax[0].set_ylabel("Relative Error")
# ax[1].set_ylabel("Relative Error")
# ax[0].grid(True)
# ax[1].grid(True)
# fig.set_figwidth(10.0)
# fig.set_figheight(8.0)
# fig.tight_layout(w_pad=0.0)
# ax[0].legend(['Pauli CST', 'Pauli CST + UP'])
# ax[1].legend(["Clifford CST", "Clifford CST + UP"])
# plt.savefig(os.path.join(plots_dir, "small_shadows_cst_vs_cst_up"))

# big shadows - cst vs. cst + up vs. cst + gc

fig, ax = plt.subplots(2,1)
m = big_molecule_list
ax[0].scatter(m, big_pauli_cst_rel_error, c='b', marker='x')
ax[0].scatter(m, big_pauli_cst_up_rel_error, c='g', marker='o')
ax[0].scatter(m, big_pauli_cst_gc_rel_error, c='r', marker='+')
ax[1].scatter(m, big_clifford_cst_rel_error, c='b', marker='x')
ax[1].scatter(m, big_clifford_cst_up_rel_error, c='g', marker='o')
ax[1].scatter(m, big_clifford_cst_gc_rel_error, c='r', marker='+')
ax[0].set_xlabel("Molecule")
ax[1].set_xlabel("Molecule")
xlabels0 = [xx[s] for s in m]
ax[0].set_xticklabels(xlabels0, rotation=0)
ax[1].set_xticklabels(xlabels0, rotation=0)
ax[0].set_ylabel("Relative Error")
ax[1].set_ylabel("Relative Error")
ax[0].grid(True)
ax[1].grid(True)
fig.set_figwidth(10.0)
fig.set_figheight(8.0)
fig.tight_layout(w_pad=0.0)
ax[0].legend(['Pauli CST', 'Pauli CST + UP', 'Pauli CST + GC'])
ax[1].legend(["Clifford CST", "Clifford CST + UP", "Clifford CST + GC"])
plt.savefig(os.path.join(plots_dir, "big_shadows_cst_vs_cst_up_vs_cst_gc"))

# # big shadows + up  vs. VQE

# fig, ax = plt.subplots(2,1)
# m = big_molecule_list
# ax[0].plot(m, vqe_10000_rel_error[:27], m, big_pauli_cst_up_rel_error, marker='x', linestyle="None")
# ax[1].plot(m, vqe_10000_rel_error[:27], m, big_clifford_cst_up_rel_error, marker='x', linestyle="None")
# ax[0].set_xlabel("Molecule")
# ax[1].set_xlabel("Molecule")
# xlabels0 = [xx[s] for s in m]
# ax[0].set_xticklabels(xlabels0, rotation=0)
# ax[1].set_xticklabels(xlabels0, rotation=0)
# ax[0].set_ylabel("Relative Error")
# ax[1].set_ylabel("Relative Error")
# ax[0].grid(True)
# ax[1].grid(True)
# fig.set_figwidth(10.0)
# fig.set_figheight(8.0)
# fig.tight_layout(w_pad=0.0)
# ax[0].legend(['VQE', 'Pauli CST + UP'])
# ax[1].legend(["VQE", "Clifford CST + UP"])
# plt.savefig(os.path.join(plots_dir, "big_shadows_up_vs_vqe"))

# small shadows + up vs. VQE

# fig, ax = plt.subplots(2,1)
# m = molecule_list
# ax[0].plot(m, vqe_1000_rel_error[:50], m, pauli_cst_up_rel_error, marker='x', linestyle="None")
# ax[1].plot(m, vqe_1000_rel_error[:50], m, clifford_cst_up_rel_error, marker='x', linestyle="None")
# ax[0].set_xlabel("Molecule")
# ax[1].set_xlabel("Molecule")
# xlabels0 = [xx[s] for s in m]
# ax[0].set_xticklabels(xlabels0, rotation=0)
# ax[1].set_xticklabels(xlabels0, rotation=0)
# ax[0].set_ylabel("Relative Error")
# ax[1].set_ylabel("Relative Error")
# ax[0].grid(True)
# ax[1].grid(True)
# fig.set_figwidth(10.0)
# fig.set_figheight(8.0)
# fig.tight_layout(w_pad=0.0)
# ax[0].legend(['VQE', 'Pauli CST + UP'])
# ax[1].legend(["VQE", "Clifford CST + UP"])
# plt.savefig(os.path.join(plots_dir, "small_shadows_up_vs_vqe"))

# small and big shadows - Pauli vs. Clifford

fig, ax = plt.subplots(2,1)
m = molecule_list
m2 = big_molecule_list
ax[0].scatter(m, pauli_cst_rel_error, c='b', marker='x')
ax[0].scatter(m, clifford_cst_rel_error, c='g', marker='o')
ax[1].scatter(m2, big_pauli_cst_rel_error, c='b', marker='x')
ax[1].scatter(m2, big_clifford_cst_rel_error, c='g', marker='o')
ax[0].set_xlabel("Molecule")
ax[1].set_xlabel("Molecule")
xlabels0 = [xx[s] for s in m]
xlabels1 = [xx[s] for s in m2]
ax[0].set_xticklabels(xlabels0, rotation=0)
ax[1].set_xticklabels(xlabels1, rotation=0)
ax[0].set_ylabel("Relative Error")
ax[1].set_ylabel("Relative Error")
ax[0].grid(True)
ax[1].grid(True)
fig.set_figwidth(10.0)
fig.set_figheight(8.0)
fig.tight_layout(w_pad=0.0)
ax[0].legend(['Pauli CST (1000)', 'Clifford CST (1000)'])
ax[1].legend(["Pauli CST (10000)", "Clifford CST (10000)"])
plt.savefig(os.path.join(plots_dir, "shadows_pauli_vs_clifford"))
