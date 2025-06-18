import numpy as np
import os
import sys

# Add parent directory to sys.path to allow imports from one level up
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from count_cusp_off import read_off_file, find_cusps, get_angle, get_individual_cusps

import subprocess, tempfile, os, json, sys, shutil

import argparse

from sklearn.neighbors import NearestNeighbors
import math


def parse_arguments():
    parser = argparse.ArgumentParser(description='Create inputs for tooth evolution simulation')
    parser.add_argument('filename', nargs='?', default='', help='Input filename to process (leave empty to generate random parameters)')
    parser.add_argument('--exe', default="./newrunt.e", help='Path to the Fortran binary (default: ./newrunt.e)')
    parser.add_argument('--steps', type=int, default=9000, help='Number of simulation steps (default: 9000)')
    parser.add_argument('--rand_output', default="random_tooth.txt", help='Output filename (default: random_tooth.txt)')
    parser.add_argument('--jitter', type=bool, default=True, help='Use jitter method (default: True)')
    parser.add_argument('r_jitter', type=int, default=5, help='Number of replicas per evaluation (default: 5)')
    parser.add_argument('sigma_rel', type=float, default=5e-3, help='Sigma relative (default: 5e-3)')
    return parser.parse_args()

# Parse command line arguments
args = parse_arguments()


# ---------------------------------------------------------------------
# 1.  Configuration you may need to adjust
# ---------------------------------------------------------------------
EXE          = args.exe        # path to the Fortran binary
STEPS        = args.steps                   # 2nd CLI arg in your call
WORKDIR_FLAG = "."                  # leave outputs in cwd

PARAM_NAMES = [
    "Egr_Epithelial_proliferation_rate",
    "Mgr_Mesenchymal_proliferation_rate",
    "Rep_Repulsion",
    "Nothing",
    "ADH_Traction_between_neighbors",
    "ACT_BMP4_auto-activation",
    "INH_Inhibtion_of_SHH_over_BMP4",
    "Nothing",
    "Sec_FGF4_secretion_rate",
    "Nothing",
    "Da_BMP4_diffusion_rate",
    "Di_SHH_diffusion_rate",
    "Ds_FGF4_diffusion_rate",
    "Nothing",
    "Int_Initial(SHH)_threshold",
    "Set_secondary(FGF4)_threshold",
    "Boy_Boyancy_Mesenchyme_mechanic_resistance",
    "Dff_Differentiation_rate",
    "Bgr_Border_growth_Amount_of_Mesenchyme_in_AP",
    "Abi_Anterior_bias",
    "Pbi_Posterior_bias",
    "Bbi_Buccal_bias_(invariable_in_seals)",
    "Lbi_Lingual_bias_(invariable)",
    "Rad_Radius_of_initial_conditions(2_for_having_7_cells)_(invariable_in_seals)",
    "Deg_Protein_degradation_rate",
    "Dgr_Downward_vector_of_growth",
    "Ntr_Nuclear_atraction_Mechanical_traction_from_the_borders_to_the_nucleus",
    "Bwi_Width_of_border"
]

# Define the parameters and their ranges
parameter_ranges = {
    'Egr_Epithelial_proliferation_rate': (0.004, 0.031),
    'Mgr_Mesenchymal_proliferation_rate': (0, 3352.557),
    'Rep_Repulsion': (0,112.481),
    'ADH_Traction_between_neighbors': (0, 1.026),
    'ACT_BMP4_auto-activation': (0.121, 0.258),
    'INH_Inhibtion_of_SHH_over_BMP4': (1.046, 21.99),
    'Sec_FGF4_secretion_rate': (0.022,0.107),
    'Da_BMP4_diffusion_rate': (0.001, 1.122),
    'Di_SHH_diffusion_rate': (0.001, 0.299),
    'Ds_FGF4_diffusion_rate': (0, 0.333),
    'Int_Initial(SHH)_threshold': (0, 0.195),
    'Set_secondary(FGF4)_threshold': (0.288, 0.999),
    'Boy_Boyancy_Mesenchyme_mechanic_resistance': (0, 0.201),
    'Dff_Differentiation_rate': (0.0004, 0.033),
    'Deg_Protein_degradation_rate': (0.01, 0.199),
    'Ntr_Nuclear_atraction_Mechanical_traction_from_the_borders_to_the_nucleus': (0, 1.599),
    'Bgr_Border_growth_Amount_of_Mesenchyme_in_AP': (0.125, 2.303),  # New sampling range
    'Pbi_Posterior_bias': (0, 32.88),  # New sampling range
    'Dgr_Downward_vector_of_growth': (733.901, 12500)  # New sampling range
}

# Define the unchanged parameters
unchanged_parameters = {
    'Nothing': 0.0,
    'Abi_Anterior_bias': 12.0,
    'Bbi_Buccal_bias_(invariable_in_seals)': 1.0,
    'Lbi_Lingual_bias_(invariable)': 1.0,
    'Rad_Radius_of_initial_conditions(2_for_having_7_cells)_(invariable_in_seals)': 2.0,
    'Bwi_Width_of_border': 0.8
}

# indices (0-based) of **free** parameters you want to vary
FREE = np.array([0,1,2,4,5,6,8,10,11,12,14,15,16,17,18,19,20,24,25,26])

def generate_random_theta():
    print("Generating random parameters")
    theta0 = np.zeros(len(PARAM_NAMES))
    for i, param_name in enumerate(PARAM_NAMES):
        if param_name in parameter_ranges:
            low, high = parameter_ranges[param_name]
            theta0[i] = rng.uniform(low, high)
        elif param_name in unchanged_parameters:
            theta0[i] = unchanged_parameters[param_name]
        else:
            # Fallback for any parameters not in either dictionary
            theta0[i] = rng.uniform(0.05, 1.0)
    return theta0

def parse_theta_from_file(input_file="1.txt"):
    if os.path.exists(args.filename): ##change this to use argparse
        theta0 = np.array([float(line.split()[0]) for line in open("1.txt") if not line.startswith("1,2,3") and line.strip()])
        return theta0
    else:
        raise ValueError(f"File {args.filename} does not exist")
    

# ---------------------------------------------------------------------
# 2.  Wrapper to run ToothMaker once
# ---------------------------------------------------------------------
def generate_toothmaker_file(theta: np.ndarray, input_file):
    """Generate a ToothMaker input file with the given parameter vector."""
    with open(input_file, "w") as f:
        # Write parameter lines
        for val, name in zip(theta, PARAM_NAMES):
            f.write(f"{val:20.13f}      {name}\n")
        
        # Write blank line and footer lines (exactly as ToothMaker expects)
        f.write("\n1,2,3,NA,4,5,6,NA,7,NA,8,9,10,NA,11,12,13,14,15,16,17,18,19,NA,20,21,NA\n")
        f.write("1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27\n")


def run_once(input_file) -> np.ndarray:
    """Run ToothMaker from a given input file, return observable y."""
    # Command line: newrunt.e  ./1.txt  .   1  9000  1
    input_file_stem = os.path.splitext(input_file)[0]
    cmd = [EXE, f"./{input_file}", WORKDIR_FLAG, f"{input_file_stem}", str(STEPS), "1"]
    print(cmd)  
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if res.returncode != 0:
        print(res.stderr.decode(), file=sys.stderr);  sys.exit(1)

    ##now, we need to parse the output file
    # Generate filename with exactly enough underscores to make total length 34
    base_name = f"{str(STEPS)}_{input_file_stem}"
    underscores_needed = 34 - len(base_name) - 4  # -4 for ".off"
    filename = f"{base_name}{'_' * underscores_needed}.off"
    vertices, faces = read_off_file(filename)
    return vertices, faces


# ---------------------------------------------------------------------
# 3.  Light-weight geometry features
# ---------------------------------------------------------------------

# --- Smoothed cusp & face complexity ---


##cusps are not particularly great as a measure, let's switch to the number of faces
##also not convinced that this makes sense
def _soft_cusp_count(cusp_heights: np.ndarray) -> float:
    """Return a **continuous** surrogate for the number of cusps.

    Each cusp contributes a weight between 0 and 1 using a logistic
    transfer-function centred at the median height.  This makes the
    value differentiable with respect to vertex heights (and therefore
    model parameters).
    """
    if cusp_heights.size == 0:
        return 0.0

    # Threshold (tau) is the median cusp height; width (kappa) is 10 % of the range
    tau   = float(np.median(cusp_heights))
    kappa = 0.1 * float(np.ptp(cusp_heights) + 1e-6)

    weights = 1.0 / (1.0 + np.exp(-(cusp_heights - tau) / kappa))
    return float(np.sum(weights))


##this doesn't work yet because everything is discrete
def extract_features(vertices, faces):
    # ---- cusps ----
    cusps = find_cusps(vertices, faces)
    cusp_heights = -np.asarray(cusps)[:, 2] if len(cusps) else np.empty(0)
    soft_cusps   = _soft_cusp_count(cusp_heights)

    # ---- faces ----
    face_count = float(len(faces))
    return np.array([face_count], dtype=float)

# --- parameters ---------------------------------
R_JITTER   = 5       # replicas per evaluation
SIGMA_REL  = 5e-6       # 0.5 % relative noise

def run_with_features(theta_vec: np.ndarray, jitter: bool = True) -> np.ndarray:
    """Return expectation of the feature vector under Gaussian jitter."""
    acc = 0.0
    if jitter:
        for _ in range(R_JITTER):
            theta_pert = theta_vec + rng.normal(scale=SIGMA_REL * np.maximum(1.0, np.abs(theta_vec)))
            generate_toothmaker_file(theta_pert, "temp_theta.txt")
            verts, faces = run_once("temp_theta.txt")
            acc += len(faces)
        acc = acc / R_JITTER
    else:
        generate_toothmaker_file(theta_vec, "temp_theta.txt")
        verts, faces = run_once("temp_theta.txt")
        acc = len(faces)

    return np.array([acc], dtype=float)


# ---------------------------------------------------------------------
# 4.  Finite-difference Jacobian
# ---------------------------------------------------------------------
EPS_REL = 1e-6
EPS_ABS = 1e-8
rng     = np.random.default_rng(42)

def calculate_jacobian(theta0, y0):
    M, P = y0.size, FREE.size
    J = np.empty((M, P))

    for j,k in enumerate(FREE):
        eps = max(abs(theta0[k])*EPS_REL, EPS_ABS)
        th_p, th_m = theta0.copy(), theta0.copy()
        th_p[k] += eps;  th_m[k] -= eps

        y_p = run_with_features(th_p)
        y_m = run_with_features(th_m)

        ##alternative where we don't use the jitter method

        J[:,j] = (y_p - y_m)/(2*eps)
        print(f"param {PARAM_NAMES[k]:35s}  done")

    H = J.T @ J
    w = np.linalg.eigvalsh(H)[::-1]          # descending
    return J, H, w



def main():
    # ---- sample one genotype ------------------------- 
    if args.filename:
        theta0 = parse_theta_from_file(args.filename)
        y0 = extract_features(*run_once(args.filename))
        yN = run_with_features(theta0)

    else:
        theta0 = generate_random_theta()
        generate_toothmaker_file(theta0, args.rand_output)
        y0 = extract_features(*run_once(args.rand_output))
        yN = run_with_features(theta0)

    
    J, H, w = calculate_jacobian(theta0, y0)

    print("True Features:")
    print(y0)
    print("Noised Features:")
    print(yN)
    print("Jacobian:")
    print(J)
    print("Hessian:")
    print(H)
    print("Eigenvalues:")
    print(w)

if __name__ == "__main__":
    main()
