import math

import pandas as pd
import mdtraj as md
import numpy as np

import warnings
import glob

warnings.simplefilter('ignore', DeprecationWarning)

########################################################################################################################
########################################################################################################################

path = "input_path"
files = sorted(glob.glob(path + "*pdb"))

out_path = "path_to_save"
dic = {}
def main():

    for comp in files:
        target_name = comp.split("/")[-1].split("_")[0]
        output_name = out_path + target_name + "_PL_Hbonds.csv"

        traj = md.load(comp)
        top = traj.topology
        lig = "UNK"

        results = run(top, traj, lig)
        int_res = summary(results, output_name)
        dic[target_name] = int_res


########################################################################################################################
########################################################################################################################

class find_PL_int:
    def __init__(self, topology, trajectory, ligand=None):

        self.topo = topology
        self.traj = trajectory

        if ligand is not None:
            self.ligand = self.topo.select("resname " + ligand)
        else:
            self.ligand = self.topo.select("resname UNK")

        self.lig_don_h = {}      # {H-index: don_N_index}
        self.lig_acc = []        # {index}
        self.lig_hv_at = {}      # {index: element}
        self.prot_don_h = {}     # {H-index: don_N_index}
        self.prot_acc = []       # {index}
        self.prot_hv_at = {}     # {index: element}

        self.to_outputs = []     # interactions pairs for the outputs
        self.PL_hbond = []       # Pairs of [P, L] within the criteria of distance and angle.
        self.PL_no_angle = []    # Pairs of [P, L] within the distance criteria.
        self.prot_res_list = []


    # Function to return pairs within distance of cutoff.
    def find_pair(self, pair, criteria=1.2):
        pair_list = np.array([[pair[0], pair[1]]])
        dist = md.compute_distances(self.traj, pair_list)
        roundup = round(dist[0][0] * 10, 1)

        # if dist * 10 <= criteria:
        if roundup <= criteria:
            return pair


    def ligand_fcn(self):

        lig_hydrogen = {}

        # From complex file, extract ligand N, O, H index.
        for at in self.ligand:
            element = self.topo.atom(at).element.symbol
            if element in ["N", "O", "S"]:
                self.lig_hv_at[at] = element
            if element == "H":
                lig_hydrogen[at] = element

        # Possible complex of ligand N, O and H.
        lig_at_h = [[x, y] for x in lig_hydrogen.keys() for y in self.lig_hv_at.keys()]

        # Use find_pair function to find bonded N,O,S-H pairs. Fill the dictionary, lig_don_h.
        for pair in lig_at_h:
            bond_pair = self.find_pair(pair)
            if bond_pair:
                self.lig_don_h[pair[0]] = pair[1]

        # Among heavy atoms in lig_hv_at dictionary, if the atom does not have a bonded hydrogen, add to the list, lig_acc.
        for at in self.lig_hv_at:
            if at not in self.lig_don_h.values():
                self.lig_acc.append(at)


        return self.lig_don_h, self.lig_acc


    def protein_fcn(self):

        prot_hydrogen = {}

        solutes = md.compute_neighbors(self.traj, 0.4, query_indices=self.lig_hv_at)  # List of proteins within 4 A of ligand heavy atoms.
        solutes = [idx for idx in solutes[0] if idx not in self.ligand]  # Exclude indices of ligand.

        # From complex file, extract protein N, O, H index.
        for at in solutes:
            element = self.topo.atom(at).element.symbol

            if element in ["N", "O", "S"]:
                self.prot_hv_at[at] = element
            if element == "H":
                prot_hydrogen[at] = element

        # Possible complex of ligand N, O and H.
        prot_at_h = [[x, y] for x in prot_hydrogen.keys() for y in self.prot_hv_at.keys()]

        # Use find_pair function to find bonded N-H or O-H pairs. Fill the dictionary, lig_don_h.
        for pair in prot_at_h:
            bond_pair = self.find_pair(pair, 1.2)
            if bond_pair:
                self.prot_don_h[pair[0]] = pair[1]

        # Among heavy atoms in lig_hv_at dictionary, if the atom does not have a bonded hydrogen, add to the list, lig_acc.
        for at in self.prot_hv_at:
            if at not in self.prot_don_h.values():
                self.prot_acc.append(at)

        return self.prot_don_h, self.prot_acc


    def prot_lig_interaction(self, distance=3.6, angle=30):

        don_acc_element = ["oxygen", "sulfur"]
        prot_resn_list = []
        PL_within_d = []        # Including pairs of {distance : [protein heavy atom, ligand heavy atom]} within cutoff. Always [P, L].

        # Possible complex of protein heavy atoms and ligand heavy atoms. Always [P, L]
        prot_lig_hv_at = [[x, y] for x in self.prot_hv_at for y in self.lig_hv_at]

        # Use find_pair function to find protein heavy atom and ligand heavy atom pairs within input distance (default=3.6 A). Add to the list, PL_within_d.
        for pair in prot_lig_hv_at:
            hb_dist_pair = self.find_pair(pair, distance)
            if hb_dist_pair:
                PL_within_d.append(hb_dist_pair)

        # From possible P-distance-L pairs, if either of a pair is in *_don_h dictionary, find the angle A-D-H.
        # If donor atom is found, the other is always counted as an acceptor atom.
        for pair in PL_within_d:
            if pair[0] in self.prot_don_h.values():  # Check protein (donor) atom.
                if pair[1] in self.lig_acc or (str(self.topo.atom(pair[1]).element) in don_acc_element):
                    hydrogens = [k for k, v in self.prot_don_h.items() if v == pair[0]]

                    if len(hydrogens) >= 1:
                        for n in range(len(hydrogens)):
                            h1 = hydrogens[n]
                            triplet_pd = np.array([[pair[1], pair[0], h1]])
                            tri_angle = md.compute_angles(self.traj, triplet_pd)
                            deg = round(math.degrees(tri_angle), 3)
                            self.PL_no_angle.append([deg, pair[0], pair[1], "pd"])
                            if deg <= angle:
                                self.PL_hbond.append([deg, pair[0], pair[1], h1, "pd"])
                                prot_res = int(self.topo.atom(pair[0]).residue.index)
                                prot_resn_list.append(str(self.topo.atom(pair[0]).residue.name) + str(prot_res))
                                self.prot_res_list.append(prot_res+1)
                                int_data = ["prot-to-lig", self.topo.atom(pair[0]), self.topo.atom(pair[1])]
                                self.to_outputs.append(int_data)


            if pair[1] in self.lig_don_h.values():  # Check Ligand (donor) atom.
                if pair[0] in self.prot_acc or (str(self.topo.atom(pair[0]).element) in don_acc_element):
                    hydrogens = [k for k, v in self.lig_don_h.items() if v == pair[1]]

                    if len(hydrogens) >= 1:
                        for n in range(len(hydrogens)):
                            h1 = hydrogens[n]
                            triplet_ld = np.array([[pair[0], pair[1], h1]])
                            deg = md.compute_angles(self.traj, triplet_ld)
                            deg = round(math.degrees(deg), 3)
                            self.PL_no_angle.append([deg, pair[0], pair[1], "ld"])
                            if deg <= angle:
                                self.PL_hbond.append([deg, pair[0], pair[1], h1, "ld"])
                                prot_res = int(self.topo.atom(pair[0]).residue.index)
                                prot_resn_list.append(str(self.topo.atom(pair[0]).residue.name) + str(prot_res))
                                self.prot_res_list.append(prot_res+1)
                                int_data = ["lig-to-prot", self.topo.atom(pair[0]), self.topo.atom(pair[1])]
                                self.to_outputs.append(int_data)

        self.prot_res_list = list(set(self.prot_res_list))

        return self.to_outputs, self.PL_hbond, self.PL_no_angle, self.prot_res_list

########################################################################################################################
########################################################################################################################

def run(top, traj, lig):

    PL_INT = find_PL_int(top, traj, lig)
    PL_INT.ligand_fcn()
    PL_INT.protein_fcn()

    results = PL_INT.prot_lig_interaction()

    return results

def summary(results, output_name):

    data, hbonds, within_d, int_res = results[0], results[1], results[2], results[3]

    df = pd.DataFrame(columns=["H-bonds", "Protein atom", "Ligand atom"])
    for n, hb in enumerate(data):
        df.loc[n, "H-bonds"] = hb[0]
        df.loc[n, "Protein atom"] = hb[1]
        df.loc[n, "Ligand atom"] = hb[2]
    # df.to_csv(output_name)

    return int_res


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    main()

print(dic)
