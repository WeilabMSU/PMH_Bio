import gudhi
import cmath
import numpy as np
from Bio import PDB
from biopandas.mol2 import PandasMol2
import matplotlib.pyplot as plt
import time
import os
import argparse
import warnings

warnings.filterwarnings("ignore")

el_p = ["C", "N", "O", "S"]  # elements in protein to be considered
el_l = [
    "C",
    "N",
    "O",
    "S",
    "P",
    "F",
    "Cl",
    "Br",
    "I",
    "H",
]  # elements in ligand to be considered


####################################################################Alpha complex
class alpha_complex:
    def __init__(self, pointcloud):
        self.simplex_list = []
        self.filtration_list = []
        comp = gudhi.AlphaComplex(points=pointcloud)
        simplex_tree = comp.create_simplex_tree()
        for filtered_value in simplex_tree.get_filtration():
            self.simplex_list.append(filtered_value[0])
            self.filtration_list.append(filtered_value[1])

    def get_complex(self, filtration_distance):
        fil_list = np.array(self.filtration_list)
        index = np.array(fil_list) <= (filtration_distance / 2) ** 2
        simplexes = np.array(self.simplex_list, dtype=object)[index]
        return simplexes.tolist()


################################################################### Mayer homology


def croot(k, n):
    if n <= 0:
        return None
    return cmath.exp((2 * cmath.pi * 1j * k) / n)


def boundary(complex, p):

    t1 = time.time()

    maxdim = len(max(complex, key=len))

    print("maxdim", maxdim)

    simplices = [
        sorted([spx for spx in complex if len(spx) == i]) for i in range(1, maxdim + 1)
    ]

    bnd = [np.zeros((0, len(simplices[0])))]

    for spx_k, spx_kp1 in zip(simplices, simplices[1:]):

        # print(spx_k,spx_kp1)
        mtx = []
        for sigma in spx_kp1:
            faces = get_faces(sigma)
            mtx.append([get_coeff(spx, faces, p) for spx in spx_k])
            # print(mtx)
        bnd.append(np.array(mtx).T)

    print(f"boundary time: {time.time() - t1}")

    return bnd


def get_faces(lst):
    return [lst[:i] + lst[i + 1 :] for i in range(len(lst))]


def get_coeff(simplex, faces, p):
    if simplex in faces:
        idxs = [i for i in range(len(faces)) if faces[i] == simplex]
        # print(idx,simplex)
        return sum([croot(idx, p) for idx in idxs])
    else:
        return 0


def matrix_rank(A):
    if A.size == 0:
        return 0
    else:
        return np.linalg.matrix_rank(A, tol=1e-10)


def betti_number(f, g):  # kerf/img
    rk_f = matrix_rank(f)

    dim_kerf = f.shape[1] - rk_f
    # print("f",f.size,f.shape,dim_kerf)
    rk_g = matrix_rank(g)
    # print("g",g.size,rk_g)
    return dim_kerf - rk_g


def mat_multiply(bnd_op, dim, q):
    f = bnd_op.get(dim - q + 1)
    # print(dim-q+1)
    for i in range(1, q):
        # print(dim+i-q+1)
        f = f @ bnd_op.get(dim + i - q + 1)

    return f


def calculate_betti(bnd_op, n, q, max_dim):  # ker d^q/imd^(n-q) for q<n

    betti = []
    for i in range(max_dim + 1):  # only consider b0,b1,l0,l1 if max_dim = 1

        f = mat_multiply(bnd_op, i, q)  # f = d^q
        g = mat_multiply(bnd_op, i + n - q, n - q)  # g = d^(n-q)

        betti.append(betti_number(f, g))

    return betti


class boundary_operators:
    def __init__(self, complex, p):
        self.boundaries = {}
        bnd = boundary(complex, p)
        self.max_dim = len(bnd)
        for i in range(self.max_dim):
            self.boundaries[i] = bnd[i]

    def get(self, index):
        if index in self.boundaries:
            return self.boundaries[index]
        elif index == self.max_dim:
            return np.zeros((self.boundaries[self.max_dim - 1].shape[1], 0))
        else:
            return np.zeros((0, 0))


def get_betti(X, N, max_dim):
    # returns Betti numberat given stage q and certain dimension
    # row number is N-2; column number is max_dim+1.

    Bettis = []
    bnd_op = boundary_operators(X, N)
    for q in range(1, N):
        betti = calculate_betti(bnd_op, N, q, max_dim=max_dim)
        Bettis.append(betti)

    return Bettis


####################################################################persistence
def calculate(pointcloud, distances, N, max_dim, feature_path, xyzfile_name, suffix):

    # max_dim=1 indicates b0 and b1
    if pointcloud.shape[0] > 0:
        # B is a dictionary indexed by filtration. Correspondingly, each value is a list of betti value indexed by q stage from 1 to N-1
        B = {}
        pre_complex_num = 0
        alpha = alpha_complex(pointcloud)
        for i, distance in enumerate(distances):
            cur_complex = alpha.get_complex(distance)
            if len(cur_complex) == pre_complex_num:
                B[distance] = B[distances[i - 1]]
            else:
                print(distance)
                bettis = get_betti(cur_complex, N, max_dim)
                print(np.shape(bettis))
                B[distance] = bettis
            pre_complex_num = len(cur_complex)

        for q_stage in range(1, N):
            for i_dim in range(max_dim + 1):
                bi = []
                for i, distance in enumerate(distances):
                    bi.append(B[distance][q_stage - 1][i_dim])
                filename = f"{feature_path}/{xyzfile_name}_n{N}_q{q_stage}_b{i_dim}_{suffix}.npy"
                np.save(filename, np.array(bi))

    else:

        for q_stage in range(1, N):
            for i_dim in range(max_dim + 1):
                filename = f"{feature_path}/{xyzfile_name}_n{N}_q{q_stage}_b{i_dim}_{suffix}.npy"
                bi = [0] * len(distances)
                np.save(filename, np.array(bi))


########################################################## element specific
class get_cloudpoint:
    def __init__(self, pdb_path, cut=12):
        self.cut = cut
        self.pdb_path = pdb_path

    def get_protein_atom_coordinate(self, pdbid):
        parser = PDB.PDBParser()
        pdbfile_path = f"{self.pdb_path}/{pdbid}/{pdbid}_protein.pdb"
        struct = parser.get_structure(pdbid, pdbfile_path)

        for model in struct:
            coor = []
            for chain in model:
                for residue in chain:
                    if residue.id[0] == " ":
                        for atom in residue:
                            XYZ = atom.get_coord()

                            coor.append(np.append(XYZ, atom.element))
            break

        return np.array(coor)

    def get_ligand_atom_coordinate(self, pdbid):
        mol2file = f"{self.pdb_path}/{pdbid}/{pdbid}_ligand.mol2"
        pmol = PandasMol2().read_mol2(mol2file)
        x = pmol.df["x"].values.reshape(-1, 1)
        y = pmol.df["y"].values.reshape(-1, 1)
        z = pmol.df["z"].values.reshape(-1, 1)
        atom_type = np.array(
            [x.split(".")[0] for x in pmol.df["atom_type"].values]
        ).reshape(-1, 1)
        return np.concatenate((x, y, z, atom_type), axis=1)

    def fit_cutoff(self, pdbid):
        # eliminate points in protein whose distance with points in ligands > 12A.
        p_coor = self.get_protein_atom_coordinate(pdbid)
        l_coor = self.get_ligand_atom_coordinate(pdbid)

        remainder = []
        for elm1 in p_coor:
            for elm2 in l_coor:
                dis = np.linalg.norm(
                    elm1[:3].astype(np.float32) - elm2[:3].astype(np.float32)
                )

                if dis <= self.cut:
                    remainder.append(elm1)
                    break

        PRO = np.array(remainder)
        LIG = l_coor

        return PRO, LIG


def generate_xyz(PRO, LIG, e_p, e_l):
    PRO_xyz = PRO[PRO[:, 3] == e_p][:, :3]
    LIG_xyz = LIG[LIG[:, 3] == e_l][:, :3]
    xyz_coor = np.concatenate([PRO_xyz, LIG_xyz])
    return xyz_coor


################## main program
def process(args, pdbid, cutoff, el_p, el_l, radius):

    N = args.N
    max_dim = args.max_dim

    feature_path = f"{args.feature_path}/{pdbid}"
    if not os.path.exists(feature_path):
        os.makedirs(feature_path)

    cloudpoint_generator = get_cloudpoint(args.pdb_path, cut=cutoff)
    PRO, LIG = cloudpoint_generator.fit_cutoff(pdbid)

    suffix = f"dmax{args.filtration_upper}-dr{args.dr:.2f}"
    for e_p in el_p:
        for e_l in el_l:
            xyzfile_name = f"{pdbid}_p{e_p}_l{e_l}"
            print(e_p, e_l)
            xyz_coor = generate_xyz(PRO, LIG, e_p, e_l)
            distances = radius * 2
            calculate(
                xyz_coor.astype(np.float32),
                distances,
                N,
                max_dim,
                feature_path,
                xyzfile_name,
                suffix,
            )


def main(args):

    cutoff = 12
    filtration_radius = np.arange(0, args.filtration_upper + args.dr, args.dr)
    pdbid = "2p7z"
    process(
        args,
        pdbid,
        cutoff,
        el_p,
        el_l,
        filtration_radius,
    )


def parse_args():

    parser = argparse.ArgumentParser(
        description="Get PMH features for protein-ligand complex"
    )
    parser.add_argument("--feature_path", type=str, default="features")
    parser.add_argument("--pdb_path", type=str, default="PDB")
    parser.add_argument("--max_dim", type=int, default=1)
    parser.add_argument("--N", type=int, default=3)
    parser.add_argument("--filtration_upper", type=int, default=10)
    parser.add_argument("--dr", type=float, default=0.25)
    args = parser.parse_args()

    print(args)
    main(args)


if __name__ == "__main__":
    parse_args()
    print("End!")
