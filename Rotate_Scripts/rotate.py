import numpy as np
import os
import quaternion as quat
import math
class Rotater():
    def __init__(self, input_file, target_atom_index, output_directory):
        self.filename = input_file
        self.target_atom_index=target_atom_index
        self.od = output_directory
    def kernel(self):
        os.makedirs(self.od,  exist_ok=True)
        os.makedirs(self.od+'/perp',  exist_ok=True)
        os.makedirs(self.od+'/para',  exist_ok=True)
        m1, m2, c1, c2, c_l, a_l = self.read_coords()
        CM = self.calc_center_of_mass_coordinate(m1, m2, c1, c2)
        new_coords = self.center_to_CM(c_l, CM)
        saved_coords = new_coords.copy()
        for theta in np.linspace(0, 2*math.pi, 50):
            theta_deg = theta * 180/math.pi
            v_N, v_O = self.calculate_rotation_vector(saved_coords, theta)
            new_coords[0] = v_N
            new_coords[1] = v_O
            self.write(new_coords, a_l, filename=self.od+'/perp/perp'+str(round(theta_deg))+'_deg.inp', mode='Perp' )
        for theta in np.linspace(0, 2*math.pi, 50):
            theta_deg = theta * 180/math.pi
            v_N, v_O = self.calculate_rotation_vector(saved_coords, theta, mode='Ortho')
            new_coords[0] = v_N
            new_coords[1] = v_O
            self.write(new_coords, a_l, filename=self.od+'/para/para'+str(round(theta_deg))+'_deg.inp', mode='Para')
        self.write_workup_submit( )
    def read_coords(self):
        with open(self.filename) as f:
            full_list = []
            for line in f:
                full_list.append(line[0:-1])
        mass_dict = {'H':1, 'C':12, 'N':14, 'O':16}
        atom_list = []
        coord_list = []
        coord_dict = {}
        for line in full_list:
            z = line.split()
            atom_list.append(z[0])
            coord_list.append([float(i) for i in z[1:]])
            coord_dict[z[0]] = [float(i) for i in z[1:]]
        
        m1 = mass_dict[atom_list[1]]
        m2 = mass_dict[atom_list[0]]
        c1 = coord_list[1]
        c2 = coord_list[0]
        self.atom_count = len(full_list)
        return m1, m2, c1, c2, coord_list, atom_list
    def calc_center_of_mass_coordinate(self, m1, m2, c1, c2):
        CM_x = (m1*c1[0] + m2*c2[0])/(m1+m2)
        CM_y = (m1*c1[1] + m2*c2[1])/(m1+m2)
        CM_z = (m1*c1[2] + m2*c2[2])/(m1+m2)
        CM = [CM_x, CM_y, CM_z]
        return CM
    def center_to_CM(self, all_coords,CM):
        coords = np.array(all_coords)
        CM = np.array(CM)
        coords = coords - CM
        new_CM = self.calc_center_of_mass_coordinate(16, 14, coords[1], coords[0])
        return coords
    def calculate_rotation_vector(self, all_coords, theta, mode=''):
        O_v = np.array(all_coords[1])
        N_v = np.array(all_coords[0])
        if mode == 'Ortho':
            axis = np.cross(np.array(O_v), np.array(all_coords[self.target_atom_index]))
        else:
            axis = all_coords[self.target_atom_index]
        axis = np.array(axis)
        axis = axis/np.linalg.norm(axis)

        q = quat.quaternion(np.cos(theta/2), *(np.sin(theta/2))*axis)
        vecO = quat.quaternion(0, *O_v)
        vecN = quat.quaternion(0, *N_v)
        v_prime_O = q * vecO * q.conjugate()
        v_prime_N = q * vecN * q.conjugate()
        v_prime_N = v_prime_N.imag
        v_prime_O = v_prime_O.imag
        return v_prime_N, v_prime_O
    def write(self, coords, atoms, filename, mode=''):
        with open(filename, "w") as w:
            w.write('!aug-cc-PVTZ\n')
            w.write('%PAL NPROCS 10 END\n')
            w.write('%CASSCF\n')
            w.write('nel 1\n')
            w.write('norb 2\n')
            w.write('nroots 2\n')
            w.write('PTMethod FIC_NEVPT2\n')
            w.write('weights 1 1\n')
            w.write('rel\n')
            w.write('  dosoc true\n')
            w.write('  PrintLevel 2\n')
            w.write('end\n')
            w.write('maxiter 500\n')
            w.write('end\n')
            w.write('* xyz 0 2\n')
            for line, atom in zip(coords, atoms):
                str_line = [atom]
                for val in line:
                    str_line.append(str(val))
                w.write(" ".join(str_line))
                w.write('\n')
            w.write('*')
        with open(self.od+'/trajectory'+mode+'.xyz', 'a') as traj:
            traj.write(str(self.atom_count)+'\n')
            traj.write(filename+'\n')
            for line, atom in zip(coords, atoms):
                str_line = [atom]
                for val in line:
                    str_line.append(str(val))
                traj.write(" ".join(str_line))
                traj.write('\n')
    def write_workup_submit(self):
        with open(self.od+"/para/submit.sub", "w") as w:
            w.write('#!/bin/bash\nfor f in *.inp; do\njobname=$(basename "$f" .inp)\nsbatch <<EOF\n#!/bin/bash\n')
            w.write('#SBATCH --job-name=$jobname\n')
            w.write('#SBATCH --output=${jobname}.out\n')
            w.write('#SBATCH --error=${jobname}.err\n')
            w.write('#SBATCH --time=24:00:00\n')
            w.write('#SBATCH --nodes=1\n')
            w.write('#SBATCH --ntasks=10\n')
            w.write('#SBATCH --mem=64G\n')
            w.write('module load orca \n')
            w.write('module load openmpi-ib/intel-2024.0/4.1.6\n')
            w.write('/sciclone/apps/orca/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca $f > ${jobname}.log\n')
            w.write('EOF\n')
            w.write('done\n')
        with open(self.od+"/perp/submit.sub", "w") as w:
            w.write('#!/bin/bash\nfor f in *.inp; do\njobname=$(basename "$f" .inp)\nsbatch <<EOF\n#!/bin/bash\n')
            w.write('#SBATCH --job-name=$jobname\n')
            w.write('#SBATCH --output=${jobname}.out\n')
            w.write('#SBATCH --error=${jobname}.err\n')
            w.write('#SBATCH --time=24:00:00\n')
            w.write('#SBATCH --nodes=1\n')
            w.write('#SBATCH --ntasks=10\n')
            w.write('#SBATCH --mem=64G\n')
            w.write('module load orca \n')
            w.write('module load openmpi-ib/intel-2024.0/4.1.6\n')
            w.write('/sciclone/apps/orca/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca $f > ${jobname}.log\n')
            w.write('EOF\n')
            w.write('done\n')
        with open(self.od+'/para/1d_conv.py', 'w') as w:
            script_text = r'''
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
'''
            w.write(script_text)
        with open(self.od+'/perp/1d_conv.py', 'w') as w:
            w.write(script_text)
        
k = Rotater(input_file='input_geom_file', target_atom_index=2, output_directory='output_directory_name')
k.kernel()