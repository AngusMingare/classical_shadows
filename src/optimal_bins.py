import os
import json
from symmer import PauliwordOp, QubitTapering, QuantumState
from tomography.cst import ClassicalShadow
import numpy as np
from qiskit import QuantumCircuit
from symmer.utils import exact_gs_energy
import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler

# Set up directories
cwd = os.getcwd()
plots_dir = os.path.join(cwd, "plots")
file_dir =  os.path.dirname(cwd)
symmer_dir = os.path.join(cwd, 'symmer')
test_dir = os.path.join(symmer_dir, 'tests')
ham_data_dir = os.path.join(test_dir, 'hamiltonian_data')
shadow_data_dir = os.path.join(file_dir, 'shadow_data')
random_clifford_shadows_dir = os.path.join(shadow_data_dir, 'random_clifford')
big_random_clifford_shadow_data_dir = os.path.join(shadow_data_dir, 'clifford_10000')

random_clifford_shadow_list = os.listdir(big_random_clifford_shadow_data_dir)
big_clifford_shadow_list = os.listdir(big_random_clifford_shadow_data_dir)

results = {"random clifford" : {}, "big random clifford" : {}}

# print("Clifford Shadows")
# cliff = {}
# for cliff_file in random_clifford_shadow_list:
#     print(cliff_file)

#     with open(os.path.join(random_clifford_shadows_dir, cliff_file), 'rb') as infile:
#         data_dict = json.load(infile)

#     with open(os.path.join(ham_data_dir, cliff_file), 'r') as infile:
#         ham_file = json.load(infile)
#     H = PauliwordOp.from_dictionary(ham_file['hamiltonian'])
#     QT = QubitTapering(H)
#     hf_state   = QuantumState(np.asarray(ham_file['data']['hf_array'])) # Hartree-Fock state
#     hf_energy  = ham_file['data']['calculated_properties']['HF']['energy']
#     H_taper   = QT.taper_it(ref_state=hf_state) 
#     gs_nrg_tap, gs_psi_tap = exact_gs_energy(H_taper.to_sparse_matrix)
#     print("gs_nrg_tap = ", str(gs_nrg_tap))

#     classical_shadow = ClassicalShadow(QuantumCircuit(gs_psi_tap.n_qubits), [])
#     classical_shadow.unitary_ensemble = "random clifford"
#     classical_shadow.shadows = []
#     for s, b in data_dict["shadow"].values():
#         classical_shadow.shadows.append((s,b))

#     H_taper_dict = H_taper.to_dictionary
#     classical_shadow.observables = [PauliwordOp.from_dictionary({o : w}).to_sparse_matrix for o,w in zip(list(H_taper_dict.keys()), list(H_taper_dict.values()))]

#     classical_shadow.num_shadows = len(classical_shadow.shadows)

#     temp = []
#     for i in [1, 4, 7, 10, 13]:
#         print("num bins = ", str(i))
        
#         _, res = classical_shadow.linearPredictions(i)
#         gs_nrg_tap_cst = 0
#         for w, exp in zip(list(H_taper_dict.values()), res):
#             gs_nrg_tap_cst += exp 
        
#         rel_error = np.abs(1 - (gs_nrg_tap_cst / gs_nrg_tap))
#         temp.append((str(i), str(rel_error)))
    
#     cliff[cliff_file] = temp

# results["random clifford"] = cliff
# print(results)

##########################################

# results = {'random clifford': {'OH-_STO-3G_SINGLET_JW.json': [('1', '0.09453796101088674'), ('4', '0.10031804682671641'), ('7', '0.20907253294169248'), ('10', '0.11104239115387582'), ('13', '0.4179250060605516')], 'B_STO-3G_DOUBLET_JW.json': [('1', '0.028338186769521045'), ('4', '0.017196834161514962'), ('7', '0.033559509130874354'), ('10', '0.07270420759484564'), ('13', '0.06009855180437229')], 'CH+_STO-3G_SINGLET_BK.json': [('1', '0.05324553105717511'), ('4', '0.04221745736860072'), ('7', '0.18429456153625412'), ('10', '0.14584024989462074'), ('13', '0.33125699897280614')], 'NH_STO-3G_SINGLET_JW.json': [('1', '0.10949174604500134'), ('4', '0.14810299746165823'), ('7', '0.35365668803774886'), ('10', '0.3760600934176622'), ('13', '0.3382827184288064')], 'BH_STO-3G_SINGLET_BK.json': [('1', '0.12272305498909342'), ('4', '0.13688792324523236'), ('7', '0.20695417074793765'), ('10', '0.3247894963760919'), ('13', '0.3863673954277139')], 'Be_STO-3G_SINGLET_BK.json': [('1', '0.02083098351034085'), ('4', '0.01116446316786357'), ('7', '0.03020115585551386'), ('10', '0.004781235860807476'), ('13', '0.024611237375028572')], 'H2_6-31G_SINGLET_BK.json': [('1', '0.5020456791973829'), ('4', '0.5292001832060589'), ('7', '0.6170474775269565'), ('10', '0.6620137181830104'), ('13', '0.3137853621343669')], 'H2_3-21G_SINGLET_BK.json': [('1', '0.3124717286324351'), ('4', '0.5912077183357325'), ('7', '0.3851074484135978'), ('10', '0.6314010497554139'), ('13', '0.4922573651509786')], 'H2_6-31G_SINGLET_JW.json': [('1', '0.14682177541485542'), ('4', '0.19636461927849258'), ('7', '0.05697575732428306'), ('10', '0.19319258605278922'), ('13', '0.7124592027212753')], 'H2_3-21G_SINGLET_JW.json': [('1', '0.2107688172688429'), ('4', '0.25739843997916023'), ('7', '0.04060226676916545'), ('10', '0.3152035614754478'), ('13', '0.06881686134573761')], 'Be_STO-3G_SINGLET_JW.json': [('1', '0.09549971201826435'), ('4', '0.07535821258583564'), ('7', '0.06632405687720411'), ('10', '0.008518992425219585'), ('13', '0.03369157264800804')], 'B_STO-3G_DOUBLET_BK.json': [('1', '0.08813499592888185'), ('4', '0.08290543835989883'), ('7', '0.09248514423707721'), ('10', '0.07464387751596069'), ('13', '0.12238866809956173')], 'HeH+_3-21G_SINGLET_BK.json': [('1', '0.03116431526014063'), ('4', '0.04093434312055522'), ('7', '0.060106349594464925'), ('10', '0.001029013637585563'), ('13', '0.359036295165167')], 'H3+_STO-3G_SINGLET_JW.json': [('1', '0.021768466761377336'), ('4', '0.030925017989923864'), ('7', '0.02713892703777765'), ('10', '0.01277113577550315'), ('13', '0.001639640960494404')], 'H4_STO-3G_SINGLET_BK.json': [('1', '0.006765676095095929'), ('4', '0.02240180056394636'), ('7', '0.01802209221779627'), ('10', '0.026145247625469437'), ('13', '0.02715281304915773')], 'C_STO-3G_TRIPLET_JW.json': [('1', '0.047852532149464944'), ('4', '0.07140913553108397'), ('7', '0.03341710537362452'), ('10', '0.06980600757683209'), ('13', '0.04014524883822279')], 'N_STO-3G_QUARTET_JW.json': [('1', '0.0318017800856909'), ('4', '0.015405296221345788'), ('7', '0.02516229713909457'), ('10', '0.0032643778971714887'), ('13', '0.024871727991802484')], 'B+_STO-3G_SINGLET_JW.json': [('1', '0.04408256432362556'), ('4', '0.03437037590970604'), ('7', '0.055183967540906265'), ('10', '0.0057285663624975225'), ('13', '0.017561719841874046')], 'Li_STO-3G_DOUBLET_BK.json': [('1', '0.03796366150064934'), ('4', '0.04152403778092473'), ('7', '0.04118775115499762'), ('10', '0.057870136571073316'), ('13', '0.037572407127897')], 'B+_STO-3G_SINGLET_BK.json': [('1', '0.03741937648123095'), ('4', '0.033425336000818984'), ('7', '0.05144761651138685'), ('10', '0.032365752822335336'), ('13', '0.009909394771822178')], 'O_STO-3G_TRIPLET_JW.json': [('1', '0.005401668710456997'), ('4', '0.0042217220065545735'), ('7', '0.006503156381420627'), ('10', '0.04961220846605541'), ('13', '0.01712349622622167')], 'Li_STO-3G_DOUBLET_JW.json': [('1', '0.038909772965092015'), ('4', '0.07369150983714512'), ('7', '0.028849146110014656'), ('10', '0.05635643533946855'), ('13', '0.061722938620721424')], 'N_STO-3G_QUARTET_BK.json': [('1', '0.041261003506414706'), ('4', '0.0775931571791244'), ('7', '0.11128420858879373'), ('10', '0.058689731946488966'), ('13', '0.032606232996316264')], 'C_STO-3G_TRIPLET_BK.json': [('1', '4.39625097312657e-05'), ('4', '0.011150651758129015'), ('7', '0.004718678473607074'), ('10', '0.0803478063638684'), ('13', '0.012570096535717057')], 'H3+_STO-3G_SINGLET_BK.json': [('1', '0.07915133563254861'), ('4', '0.08113258488359798'), ('7', '0.12231220614148375'), ('10', '0.08345865753234838'), ('13', '0.08500430932436376')], 'H4_STO-3G_SINGLET_JW.json': [('1', '0.027431520279767252'), ('4', '0.07838873956926273'), ('7', '0.05412293392343992'), ('10', '0.08172469699876472'), ('13', '0.0923676564627568')], 'HeH+_3-21G_SINGLET_JW.json': [('1', '0.025129578359811777'), ('4', '0.1420573630956048'), ('7', '0.1088393105962987'), ('10', '0.23639366435963605'), ('13', '0.18425530734476414')]}, 'big random clifford': {}}

# print("Big Clifford Shadows")
# big_cliff = {}
# for cliff_file in big_clifford_shadow_list:
#     print(cliff_file)

#     with open(os.path.join(big_random_clifford_shadow_data_dir, cliff_file), 'rb') as infile:
#         data_dict = json.load(infile)

#     with open(os.path.join(ham_data_dir, cliff_file), 'r') as infile:
#         ham_file = json.load(infile)
#     H = PauliwordOp.from_dictionary(ham_file['hamiltonian'])
#     QT = QubitTapering(H)
#     hf_state   = QuantumState(np.asarray(ham_file['data']['hf_array'])) # Hartree-Fock state
#     hf_energy  = ham_file['data']['calculated_properties']['HF']['energy']
#     H_taper   = QT.taper_it(ref_state=hf_state) 
#     gs_nrg_tap, gs_psi_tap = exact_gs_energy(H_taper.to_sparse_matrix)

#     classical_shadow = ClassicalShadow(QuantumCircuit(gs_psi_tap.n_qubits), [])
#     classical_shadow.unitary_ensemble = "random clifford"
#     classical_shadow.shadows = []
#     for s, b in data_dict["shadow"].values():
#         classical_shadow.shadows.append((s,b))

#     H_taper_dict = H_taper.to_dictionary
#     classical_shadow.observables = [PauliwordOp.from_dictionary({o : w}).to_sparse_matrix for o,w in zip(list(H_taper_dict.keys()), list(H_taper_dict.values()))]

#     classical_shadow.num_shadows = len(classical_shadow.shadows)

#     temp = []
#     for i in [1, 4, 7, 10, 13]:
#         print("num bins = ", str(i))
        
#         _, res = classical_shadow.linearPredictions(i)
#         gs_nrg_tap_cst = 0
#         for w, exp in zip(list(H_taper_dict.values()), res):
#             gs_nrg_tap_cst += exp 

#         rel_error = np.abs(1 - (gs_nrg_tap_cst / gs_nrg_tap))
#         temp.append((str(i), str(rel_error)))
    
#     big_cliff[cliff_file] = temp

# results["big random clifford"] = big_cliff
# print(results)

# ###########################

# with open(os.path.join(cwd, 'optimal_bins_results.json'), 'w') as file:
#         json.dump(results, file)

# ############################

# with open(os.path.join(cwd, 'optimal_bins_results.json'), 'w') as file:
#         data_dict = json.load(file)

data_dict = {'random clifford': {'OH-_STO-3G_SINGLET_JW.json': [('1', '0.09453796101088674'), ('4', '0.10031804682671641'), ('7', '0.20907253294169248'), ('10', '0.11104239115387582'), ('13', '0.4179250060605516')], 'B_STO-3G_DOUBLET_JW.json': [('1', '0.028338186769521045'), ('4', '0.017196834161514962'), ('7', '0.033559509130874354'), ('10', '0.07270420759484564'), ('13', '0.06009855180437229')], 'CH+_STO-3G_SINGLET_BK.json': [('1', '0.05324553105717511'), ('4', '0.04221745736860072'), ('7', '0.18429456153625412'), ('10', '0.14584024989462074'), ('13', '0.33125699897280614')], 'NH_STO-3G_SINGLET_JW.json': [('1', '0.10949174604500134'), ('4', '0.14810299746165823'), ('7', '0.35365668803774886'), ('10', '0.3760600934176622'), ('13', '0.3382827184288064')], 'BH_STO-3G_SINGLET_BK.json': [('1', '0.12272305498909342'), ('4', '0.13688792324523236'), ('7', '0.20695417074793765'), ('10', '0.3247894963760919'), ('13', '0.3863673954277139')], 'Be_STO-3G_SINGLET_BK.json': [('1', '0.02083098351034085'), ('4', '0.01116446316786357'), ('7', '0.03020115585551386'), ('10', '0.004781235860807476'), ('13', '0.024611237375028572')], 'H2_6-31G_SINGLET_BK.json': [('1', '0.5020456791973829'), ('4', '0.5292001832060589'), ('7', '0.6170474775269565'), ('10', '0.6620137181830104'), ('13', '0.3137853621343669')], 'H2_3-21G_SINGLET_BK.json': [('1', '0.3124717286324351'), ('4', '0.5912077183357325'), ('7', '0.3851074484135978'), ('10', '0.6314010497554139'), ('13', '0.4922573651509786')], 'H2_6-31G_SINGLET_JW.json': [('1', '0.14682177541485542'), ('4', '0.19636461927849258'), ('7', '0.05697575732428306'), ('10', '0.19319258605278922'), ('13', '0.7124592027212753')], 'H2_3-21G_SINGLET_JW.json': [('1', '0.2107688172688429'), ('4', '0.25739843997916023'), ('7', '0.04060226676916545'), ('10', '0.3152035614754478'), ('13', '0.06881686134573761')], 'Be_STO-3G_SINGLET_JW.json': [('1', '0.09549971201826435'), ('4', '0.07535821258583564'), ('7', '0.06632405687720411'), ('10', '0.008518992425219585'), ('13', '0.03369157264800804')], 'B_STO-3G_DOUBLET_BK.json': [('1', '0.08813499592888185'), ('4', '0.08290543835989883'), ('7', '0.09248514423707721'), ('10', '0.07464387751596069'), ('13', '0.12238866809956173')], 'HeH+_3-21G_SINGLET_BK.json': [('1', '0.03116431526014063'), ('4', '0.04093434312055522'), ('7', '0.060106349594464925'), ('10', '0.001029013637585563'), ('13', '0.359036295165167')], 'H3+_STO-3G_SINGLET_JW.json': [('1', '0.021768466761377336'), ('4', '0.030925017989923864'), ('7', '0.02713892703777765'), ('10', '0.01277113577550315'), ('13', '0.001639640960494404')], 'H4_STO-3G_SINGLET_BK.json': [('1', '0.006765676095095929'), ('4', '0.02240180056394636'), ('7', '0.01802209221779627'), ('10', '0.026145247625469437'), ('13', '0.02715281304915773')], 'C_STO-3G_TRIPLET_JW.json': [('1', '0.047852532149464944'), ('4', '0.07140913553108397'), ('7', '0.03341710537362452'), ('10', '0.06980600757683209'), ('13', '0.04014524883822279')], 'N_STO-3G_QUARTET_JW.json': [('1', '0.0318017800856909'), ('4', '0.015405296221345788'), ('7', '0.02516229713909457'), ('10', '0.0032643778971714887'), ('13', '0.024871727991802484')], 'B+_STO-3G_SINGLET_JW.json': [('1', '0.04408256432362556'), ('4', '0.03437037590970604'), ('7', '0.055183967540906265'), ('10', '0.0057285663624975225'), ('13', '0.017561719841874046')], 'Li_STO-3G_DOUBLET_BK.json': [('1', '0.03796366150064934'), ('4', '0.04152403778092473'), ('7', '0.04118775115499762'), ('10', '0.057870136571073316'), ('13', '0.037572407127897')], 'B+_STO-3G_SINGLET_BK.json': [('1', '0.03741937648123095'), ('4', '0.033425336000818984'), ('7', '0.05144761651138685'), ('10', '0.032365752822335336'), ('13', '0.009909394771822178')], 'O_STO-3G_TRIPLET_JW.json': [('1', '0.005401668710456997'), ('4', '0.0042217220065545735'), ('7', '0.006503156381420627'), ('10', '0.04961220846605541'), ('13', '0.01712349622622167')], 'Li_STO-3G_DOUBLET_JW.json': [('1', '0.038909772965092015'), ('4', '0.07369150983714512'), ('7', '0.028849146110014656'), ('10', '0.05635643533946855'), ('13', '0.061722938620721424')], 'N_STO-3G_QUARTET_BK.json': [('1', '0.041261003506414706'), ('4', '0.0775931571791244'), ('7', '0.11128420858879373'), ('10', '0.058689731946488966'), ('13', '0.032606232996316264')], 'C_STO-3G_TRIPLET_BK.json': [('1', '4.39625097312657e-05'), ('4', '0.011150651758129015'), ('7', '0.004718678473607074'), ('10', '0.0803478063638684'), ('13', '0.012570096535717057')], 'H3+_STO-3G_SINGLET_BK.json': [('1', '0.07915133563254861'), ('4', '0.08113258488359798'), ('7', '0.12231220614148375'), ('10', '0.08345865753234838'), ('13', '0.08500430932436376')], 'H4_STO-3G_SINGLET_JW.json': [('1', '0.027431520279767252'), ('4', '0.07838873956926273'), ('7', '0.05412293392343992'), ('10', '0.08172469699876472'), ('13', '0.0923676564627568')], 'HeH+_3-21G_SINGLET_JW.json': [('1', '0.025129578359811777'), ('4', '0.1420573630956048'), ('7', '0.1088393105962987'), ('10', '0.23639366435963605'), ('13', '0.18425530734476414')]}, 'big random clifford': {'OH-_STO-3G_SINGLET_JW.json': [('1', '0.03118809424105473'), ('4', '0.045296999123040305'), ('7', '0.027102697321306213'), ('10', '0.04990129019816547'), ('13', '0.04269878505346747')], 'B_STO-3G_DOUBLET_JW.json': [('1', '0.0026327546433718707'), ('4', '0.002331162460346814'), ('7', '0.009730798088450343'), ('10', '0.006932329317425934'), ('13', '0.022804911353216717')], 'CH+_STO-3G_SINGLET_BK.json': [('1', '0.008738030824773091'), ('4', '0.013157898712032101'), ('7', '0.029526510629561686'), ('10', '0.03869455766280028'), ('13', '0.0040742608002943825')], 'NH_STO-3G_SINGLET_JW.json': [('1', '0.02721145873394004'), ('4', '0.008883535830833011'), ('7', '0.04125672305968653'), ('10', '0.04850301848342564'), ('13', '0.00608805029775028')], 'BH_STO-3G_SINGLET_BK.json': [('1', '0.021008196616478347'), ('4', '0.025118767740017645'), ('7', '0.013023568914816819'), ('10', '0.018984940812080353'), ('13', '0.0006125963459389094')], 'Be_STO-3G_SINGLET_BK.json': [('1', '0.00026588284420103747'), ('4', '0.00037617352478314103'), ('7', '0.002266156244255413'), ('10', '0.006550673150366282'), ('13', '0.008555378006829706')], 'H2_6-31G_SINGLET_BK.json': [('1', '0.11208742721622178'), ('4', '0.09428814760105164'), ('7', '0.1235383867879909'), ('10', '0.06515677382008722'), ('13', '0.07943201406694766')], 'H2_3-21G_SINGLET_BK.json': [('1', '0.12968744718909386'), ('4', '0.1773857296704452'), ('7', '0.1620405499273183'), ('10', '0.19155971342326916'), ('13', '0.10812187420375663')], 'H2_6-31G_SINGLET_JW.json': [('1', '0.052687988824701826'), ('4', '0.10776874631762867'), ('7', '0.06823856860130606'), ('10', '0.029515656687977065'), ('13', '0.07546298003799068')], 'H2_3-21G_SINGLET_JW.json': [('1', '0.025062393243552528'), ('4', '0.04997543660639758'), ('7', '0.04157627137954212'), ('10', '0.050517235543423245'), ('13', '0.013040559582321198')], 'Be_STO-3G_SINGLET_JW.json': [('1', '0.007476411746774492'), ('4', '0.009466417240956204'), ('7', '0.013997170514713941'), ('10', '0.008910376181828106'), ('13', '0.012704836378344408')], 'B_STO-3G_DOUBLET_BK.json': [('1', '0.022297582876187327'), ('4', '0.02441397343984142'), ('7', '0.01976708796438087'), ('10', '0.02277918519511135'), ('13', '0.031091130888479923')], 'HeH+_3-21G_SINGLET_BK.json': [('1', '0.0009051222578488893'), ('4', '0.005667844887454532'), ('7', '0.014350670241874974'), ('10', '0.01997010999686122'), ('13', '0.037819232910833334')], 'H3+_STO-3G_SINGLET_JW.json': [('1', '0.014133140106396658'), ('4', '0.0022646962404828663'), ('7', '0.018839254709377085'), ('10', '0.023609479571473635'), ('13', '0.017457840600729746')], 'H4_STO-3G_SINGLET_BK.json': [('1', '0.03129603313656815'), ('4', '0.03139351631127907'), ('7', '0.03335906034635405'), ('10', '0.025885698232720022'), ('13', '0.036092596654884845')], 'C_STO-3G_TRIPLET_JW.json': [('1', '0.022941877918807307'), ('4', '0.02056005930478888'), ('7', '0.03620086198449346'), ('10', '0.02515865573015641'), ('13', '0.017121877988214407')], 'N_STO-3G_QUARTET_JW.json': [('1', '0.012154567323625298'), ('4', '0.008823064663485147'), ('7', '0.009295934208395917'), ('10', '0.0038521138122262766'), ('13', '0.009491249105721411')], 'B+_STO-3G_SINGLET_JW.json': [('1', '0.017311058383037148'), ('4', '0.02051436829605735'), ('7', '0.01826540284945133'), ('10', '0.012670284305829949'), ('13', '0.009434520224316012')], 'Li_STO-3G_DOUBLET_BK.json': [('1', '0.03297541726445985'), ('4', '0.029310287724593653'), ('7', '0.029730528929118538'), ('10', '0.036138612486086075'), ('13', '0.028384942742296193')], 'B+_STO-3G_SINGLET_BK.json': [('1', '0.015209603376713954'), ('4', '0.01703433040655422'), ('7', '0.019840222813713893'), ('10', '0.02120705108517118'), ('13', '0.01709942509454654')], 'O_STO-3G_TRIPLET_JW.json': [('1', '0.01577680395975123'), ('4', '0.021582252354753306'), ('7', '0.015527230224887223'), ('10', '0.004593050312061697'), ('13', '0.00963050635047269')], 'Li_STO-3G_DOUBLET_JW.json': [('1', '0.020458006083394653'), ('4', '0.023116560298926836'), ('7', '0.03079289731051249'), ('10', '0.025617703804822645'), ('13', '0.026985477025894378')], 'N_STO-3G_QUARTET_BK.json': [('1', '0.008478670970336166'), ('4', '0.0029870163663830995'), ('7', '0.003152402086959105'), ('10', '0.004518207705232369'), ('13', '0.013703099890796988')], 'C_STO-3G_TRIPLET_BK.json': [('1', '0.002908760537478905'), ('4', '0.004934573436164991'), ('7', '0.004417062937615235'), ('10', '0.009090792776158096'), ('13', '0.018077750590604014')], 'H3+_STO-3G_SINGLET_BK.json': [('1', '0.006704362965173072'), ('4', '0.009200689792460803'), ('7', '0.012079916133890256'), ('10', '0.014744798770880863'), ('13', '0.0027839747821649663')], 'H4_STO-3G_SINGLET_JW.json': [('1', '0.0002876620349346881'), ('4', '0.01309884219247226'), ('7', '0.023399358573962026'), ('10', '0.012490143794825048'), ('13', '0.022956602827971073')], 'HeH+_3-21G_SINGLET_JW.json': [('1', '0.02768731858814688'), ('4', '0.03565809590209379'), ('7', '0.05020268179932186'), ('10', '0.0734093526444054'), ('13', '0.034871738932427454')]}}

plt.rcParams.update({'axes.labelsize': 18, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 16})
mpl.rcParams['axes.prop_cycle'] = cycler(color=['b', 'g', 'r', 'c', 'm', 'y'])

cliff_dict = data_dict["random clifford"]
big_cliff_dict = data_dict["big random clifford"]

molecule_list = cliff_dict.keys()
big_molecule_list = big_cliff_dict.keys()

xx_list = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',
           'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
xx = {}
for mol_idx in range(len(molecule_list)):
        xx[list(molecule_list)[mol_idx]] = xx_list[mol_idx]

rel_error_1_bins = []
rel_error_4_bins = []
rel_error_7_bins = []
rel_error_10_bins = []
rel_error_13_bins = []

for vals in cliff_dict.values():
    for bins, err in vals:
        if bins == "1":
            rel_error_1_bins.append(float(err))
        elif bins == "4":
            rel_error_4_bins.append(float(err))
        elif bins == "7":
            rel_error_7_bins.append(float(err))
        elif bins == "10":
            rel_error_10_bins.append(float(err))
        elif bins == "13":
            rel_error_13_bins.append(float(err))
     

big_rel_error_1_bins = []
big_rel_error_4_bins = []
big_rel_error_7_bins = []
big_rel_error_10_bins = []
big_rel_error_13_bins = []

for vals in big_cliff_dict.values():
    for bins, err in vals:
        if bins == "1":
            big_rel_error_1_bins.append(float(err))
        elif bins == "4":
            big_rel_error_4_bins.append(float(err))
        elif bins == "7":
            big_rel_error_7_bins.append(float(err))
        elif bins == "10":
            big_rel_error_10_bins.append(float(err))
        elif bins == "13":
            big_rel_error_13_bins.append(float(err))

# fig, ax = plt.subplots(2,1)
# x0 = molecule_list
# x1 = big_molecule_list
# ax[0].plot(x0, rel_error_1_bins, x0, rel_error_4_bins, x0, rel_error_7_bins, x0, rel_error_10_bins, x0, rel_error_13_bins, marker='x', linestyle="None")
# ax[0].set_xlabel("Molecule")
# ax[1].plot(x1, big_rel_error_1_bins, x1, big_rel_error_4_bins, x1, big_rel_error_7_bins, x1, big_rel_error_10_bins, x1, big_rel_error_13_bins, marker='x', linestyle="None")
# ax[1].set_xlabel("Molecule")

# ax[0].axhline(y=np.mean(rel_error_1_bins),c="b",linewidth=0.5,zorder=0)
# ax[0].axhline(y=np.mean(rel_error_4_bins),c="g",linewidth=0.5,zorder=0)
# ax[0].axhline(y=np.mean(rel_error_7_bins),c="r",linewidth=0.5,zorder=0)
# ax[0].axhline(y=np.mean(rel_error_10_bins),c="c",linewidth=0.5,zorder=0)
# ax[0].axhline(y=np.mean(rel_error_13_bins),c="m",linewidth=0.5,zorder=0)

# ax[1].axhline(y=np.mean(big_rel_error_1_bins),c="b",linewidth=0.5,zorder=0)
# ax[1].axhline(y=np.mean(big_rel_error_4_bins),c="g",linewidth=0.5,zorder=0)
# ax[1].axhline(y=np.mean(big_rel_error_7_bins),c="r",linewidth=0.5,zorder=0)
# ax[1].axhline(y=np.mean(big_rel_error_10_bins),c="c",linewidth=0.5,zorder=0)
# ax[1].axhline(y=np.mean(big_rel_error_13_bins),c="m",linewidth=0.5,zorder=0)

# xlabels0 = [xx[m] for m in x0]
# xlabels1 = [xx[m] for m in x1]
# ax[0].set_xticklabels(xlabels0, rotation=0)
# ax[0].set_ylabel("Relative Error")
# ax[0].grid(True)
# ax[1].set_xticklabels(xlabels1, rotation=0)
# ax[1].set_ylabel("Relative Error")
# ax[1].grid(True)
# fig.set_figwidth(10.0)
# fig.set_figheight(8.0)
# fig.tight_layout(w_pad=0.0)
# ax[0].legend(['1 bin', '4 bins', '7 bins', '10 bins', '13 bins'])
# ax[1].legend(['1 bin', '4 bins', '7 bins', '10 bins', '13 bins'])
# plt.savefig(os.path.join(plots_dir, "optimal_bins_both.png"))

fig, ax = plt.subplots(2,1)
x0 = molecule_list
x1 = big_molecule_list

ax[0].axhline(y=np.mean(rel_error_1_bins),c="b",linewidth=0.5,zorder=0)
ax[0].axhline(y=np.mean(rel_error_4_bins),c="g",linewidth=0.5,zorder=0)
ax[0].axhline(y=np.mean(rel_error_7_bins),c="r",linewidth=0.5,zorder=0)
ax[0].axhline(y=np.mean(rel_error_10_bins),c="gold",linewidth=0.5,zorder=0)
ax[0].axhline(y=np.mean(rel_error_13_bins),c="indigo",linewidth=0.5,zorder=0)

ax[1].axhline(y=np.mean(big_rel_error_1_bins),c="b",linewidth=0.5,zorder=0)
ax[1].axhline(y=np.mean(big_rel_error_4_bins),c="g",linewidth=0.5,zorder=0)
ax[1].axhline(y=np.mean(big_rel_error_7_bins),c="r",linewidth=0.5,zorder=0)
ax[1].axhline(y=np.mean(big_rel_error_10_bins),c="c",linewidth=0.5,zorder=0)
ax[1].axhline(y=np.mean(big_rel_error_13_bins),c="m",linewidth=0.5,zorder=0)

ax[0].set_xlabel("Molecule")
ax[1].set_xlabel("Molecule")
ax[0].set_ylim(0.0, 0.4)
ax[1].set_ylim(0.0, 0.1)
xlabels0 = [xx[m] for m in x0]
xlabels1 = [xx[m] for m in x1]
ax[0].set_xticklabels(xlabels0, rotation=0)
ax[0].set_ylabel("Relative Error")
ax[0].grid(True)
ax[1].set_xticklabels(xlabels1, rotation=0)
ax[1].set_ylabel("Relative Error")
ax[1].grid(True)
fig.set_figwidth(10.0)
fig.set_figheight(8.0)
fig.tight_layout(w_pad=0.0)
ax[0].legend(['1 bin', '4 bins', '7 bins', '10 bins', '13 bins'])
ax[1].legend(['1 bin', '4 bins', '7 bins', '10 bins', '13 bins'])
plt.savefig(os.path.join(plots_dir, "optimal_bins_means.png"))