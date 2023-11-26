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

with open(os.path.join(cwd, 'biased_results.json'), 'r') as file:
        data = json.load(file)

big_clifford_results = data["big_clifford_shadows"]
big_biased_clifford_results = data["big_biased_clifford_shadows"]

molecule_list = big_clifford_results.keys()

random_clifford_cst_rel_error = []
biased_clifford_cst_rel_error = []
for molecule in molecule_list:
        exact_energy = complex(big_clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(big_clifford_results[molecule]["gs_nrg_tap_cst"])
        random_clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))

        exact_energy = complex(big_biased_clifford_results[molecule]["gs_nrg_tap"])
        cst_energy = complex(big_biased_clifford_results[molecule]["gs_nrg_tap_cst"])
        biased_clifford_cst_rel_error.append(np.abs(1 - (cst_energy / exact_energy)))


# small shadows - minimal basis vs. not minimal basis

plt.rcParams.update({'axes.labelsize': 18, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 16})
mpl.rcParams['axes.prop_cycle'] = cycler(color=['b', 'g', 'r'])

xx = {"B_STO-3G_DOUBLET_JW.json" : "b", "Be_STO-3G_SINGLET_BK.json" : "g", "H2_6-31G_SINGLET_BK.json" : "i", "H2_3-21G_SINGLET_BK.json" : "k", 
"H2_6-31G_SINGLET_JW.json" : "m", "H2_3-21G_SINGLET_JW.json" : "p", "Be_STO-3G_SINGLET_JW.json" : "q", "B_STO-3G_DOUBLET_BK.json" : "v", 
"HeH+_3-21G_SINGLET_BK.json" : "x", "H3+_STO-3G_SINGLET_JW.json" : "y", "H4_STO-3G_SINGLET_BK.json" : "z", "C_STO-3G_TRIPLET_JW.json" : "D", 
"N_STO-3G_QUARTET_JW.json" : "F", "B+_STO-3G_SINGLET_JW.json" : "H", "Li_STO-3G_DOUBLET_BK.json" : "I", "B+_STO-3G_SINGLET_BK.json" : "M",
"O_STO-3G_TRIPLET_JW.json" : "N", "H3+_STO-3G_SINGLET_BK.json" : "V", "H4_STO-3G_SINGLET_JW.json" : "W", "HeH+_3-21G_SINGLET_JW.json" : "X"}

fig, ax = plt.subplots(1,1)
x = molecule_list
ax.scatter(x, random_clifford_cst_rel_error, c='b', marker='x')
ax.scatter(x, biased_clifford_cst_rel_error, c='g', marker='o')
ax.set_xlabel("Molecule")
xlabels0 = [xx[m] for m in x]
ax.set_xticklabels(xlabels0, rotation=0)
ax.set_ylabel("Relative Error")
ax.grid(True)
fig.set_figwidth(10.0)
fig.set_figheight(4.0)
fig.tight_layout(w_pad=0.0)
ax.legend(['Random Clifford CST', 'Globally Biased CST'])
plt.savefig(os.path.join(plots_dir, "random_vs_biased_clifford_shadows"))