import matplotlib.pyplot as plt
import numpy as np
import os
import re
import csv
import glob
def extract_soc_matrix(text, mode='Real'):
    lines = text.splitlines()

    # locate SOC block
    start = None
    for i, line in enumerate(lines):
        if "SOC MATRIX" in line:
            start = i
        if start is not None and "Image part:" in line:
            start = i + 2
            break

    if start is None:
        return None

    matrix = []

    for line in lines[start:]:
        if line.strip() == "":
            break

        # match numeric rows like: 0  -169.4  0.0 ...
        parts = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)
        # skip header lines
        if len(parts) < 2:
            continue

        # first entry is row index, rest is matrix values
        matrix.append([float(x) for x in parts[1:]])
    for i, line in enumerate(lines):
        if "SOC MATRIX" in line:
            start = i
        if start is not None and "Real part:" in line:
            start = i + 2
            break

    if start is None:
        return None

    matrix_real = []

    for line in lines[start:]:
        if line.strip() == '' or "Image Part" in line:
            break

        # match numeric rows like: 0  -169.4  0.0 ...
        parts = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)
        # skip header lines
        if len(parts) < 2:
            continue

        # first entry is row index, rest is matrix values
        matrix_real.append([float(x) for x in parts[1:]])
    soc_energy = float(re.search(r"Energy stabilization:\s*([-\d.]+)", text).group(1))
    return np.array(matrix), np.array(matrix_real), soc_energy

# Folder containing ORCA output files
folder = "./"   # change if needed

# Output CSV file
output_csv = "orca_s1_energies.csv"

results = []

# Grab all ORCA output files
files = glob.glob(os.path.join(folder, "*_deg.log"))

for file in files:
    with open(file, "r") as f:
        text = f.read()
        soc_imag, soc_real, energy_soc = extract_soc_matrix(text)
        inside_block=False
        for line in text.splitlines():
            if "NEVPT2 TOTAL ENERGIES" in line:
                inside_block = True
                continue

            if "NEVPT2 TRANSITION ENERGIES" in line:
                break

            if inside_block:
                fields = line.split()

        # Data lines look like:
        # ['0:', '0', '2', '-169.968801', 'EDIAG[0]', '-169.968800935']
                if len(fields) >= 6 and fields[0].endswith(":"):

                    state = int(fields[0].replace(":", ""))
                    energy = float(fields[5])   # full precision EDIAG value

                    if state == 0:
                        nevpt2_s0 = energy
                    elif state == 1:
                        nevpt2_s1 = energy
    casscf_matches = re.findall(r"ROOT\s+(\d+):\s+E=\s+(-?\d+\.\d+)\s+Eh",text)

    casscf_energies = {int(root): float(E) for root, E in casscf_matches}

    casscf_s0 = casscf_energies.get(0, np.nan)
    casscf_s1 = casscf_energies.get(1, np.nan)
    #nevpt2_s1 = nevpt2_energies.get(1, np.nan)
    results.append([
        os.path.basename(file),
        casscf_s0,
        casscf_s1,
        nevpt2_s0,
        nevpt2_s1,
        abs(soc_real[0][3]*219474.3), 
        abs(soc_imag[0][3]*219474.3),
        abs(energy_soc)
        ])

# -----------------------------
# Write CSV
# -----------------------------

energy_s0_cas = []
energy_s1_cas = []
energy_s0_nevpt2 = []
energy_s1_nevpt2 = []
angle = []
soc_imag_l = []
soc_real_l = []
energy_soc_l = []
for row in results:
    energy_s0_cas.append(row[1])
    energy_s1_cas.append(row[2])
    energy_s0_nevpt2.append(row[3])
    energy_s1_nevpt2.append(row[4])

    name = row[0].split('_')
    angle.append(float(name[0][4:]))
    soc_imag_l.append(row[6])
    soc_real_l.append(row[5])
    energy_soc_l.append(row[7])
#exit()
energy_s0_cas = np.array(energy_s0_cas)*219474.63
energy_s1_cas = np.array(energy_s1_cas)*219474.63

energy_s1_cas = -min(energy_s0_cas) + energy_s1_cas
energy_s0_cas = -min(energy_s0_cas) + energy_s0_cas


energy_s0_nev = np.array(energy_s0_nevpt2)*219474.63
energy_s1_nev = np.array(energy_s1_nevpt2)*219474.63
energy_s1_nev = -min(energy_s0_nev) + energy_s1_nev
energy_s0_nev = -min(energy_s0_nev) + energy_s0_nev


plt.figure()
plt.scatter(angle, soc_imag_l, label='SOC-Off-Diagonal Imaginary')
plt.scatter(angle, soc_real_l, label='SOC-Off-Diagonal Real')
#plt.scatter(angle, energy_soc_l, label='SOC Stabillization')
plt.xlabel('NO Rotation Angle (Degrees)')
plt.ylabel('Relative Energy (cm^-1)')
plt.legend()
plt.savefig('spin_orbit')
plt.show()

plt.figure()
plt.scatter(angle, energy_s0_cas, label='S0')
plt.scatter(angle, energy_s1_cas, label='S1')
plt.xlabel('NO Rotation Angle (Degrees)')
plt.ylabel('Relative Energy CASSCF (cm^-1)')
plt.legend()
plt.show()
plt.savefig('CASSCF_Surf.png')

plt.figure()
plt.scatter(angle, energy_s0_nev, label='S0')
plt.scatter(angle, energy_s1_nev, label='S1')
plt.xlabel('NO Rotation Angle (Degrees)')
plt.ylabel('Relative Energy NEVPT2 (cm^-1)')
plt.legend()
plt.show()
plt.savefig('NEVPT2_Surf.png')

print(f"Saved results to {output_csv}")
