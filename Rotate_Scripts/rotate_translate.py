import numpy as np
import quaternion as quat
import math
import os
class Rotater():
    def __init__(self, input_file, target_atom_index, output_directory):
        self.filename = input_file
        self.target_atom_index=target_atom_index
        self.od = output_directory
    def kernel(self, r_step, t_step):
        os.makedirs(self.od,  exist_ok=True)
        os.makedirs(self.od+'/perp',  exist_ok=True)
        os.makedirs(self.od+'/para',  exist_ok=True)
        m1, m2, c1, c2, c_l, a_l = self.read_coords()
        CM = self.calc_center_of_mass_coordinate(m1, m2, c1, c2)
        new_coords = self.center_to_CM(c_l, CM)
        for t in np.linspace(0, -1, t_step):
            c_new = self.translate_NO(new_coords, t)
            CM = self.calc_center_of_mass_coordinate(m1, m2, c_new[1], c_new[0])
            fin_coords = self.center_to_CM(c_new, CM)
            dist = round(np.linalg.norm(fin_coords[self.target_atom_index]),2)
            saved_coords = fin_coords.copy()
            os.makedirs(self.od+'/para/'+str(dist)+'_A', exist_ok=True)
            os.makedirs(self.od+'/perp/'+str(dist)+'_A', exist_ok=True)
            for theta in np.linspace(0, 2*math.pi, r_step):
                theta_deg = theta * 180/math.pi
                v_N, v_O = self.calculate_rotation_vector(saved_coords, theta)
                fin_coords[0] = v_N
                fin_coords[1] = v_O
                self.write(fin_coords, a_l, filename=self.od+'/perp/'+str(dist)+'_A/perp'+str(round(theta_deg))+'_deg.inp', mode='Perp_'+str(dist))
            for theta in np.linspace(0, 2*math.pi, 50):
                theta_deg = theta * 180/math.pi
                v_N, v_O = self.calculate_rotation_vector(saved_coords, theta, mode='Ortho')
                fin_coords[0] = v_N
                fin_coords[1] = v_O
                self.write(fin_coords, a_l, filename=self.od+'/para/'+str(dist)+'_A/para'+str(round(theta_deg))+'_deg.inp', mode='Para_'+str(dist))
        self.write_workup_submit()
    def translate_NO(self, coords, translation):
        new_coords = coords.copy()
        dir_vector = coords[self.target_atom_index]
        #print(np.linalg.norm(dir_vector))
        N_coord = coords[0] - (dir_vector/np.linalg.norm(dir_vector))*translation
        O_coord = coords[1] - (dir_vector/np.linalg.norm(dir_vector))*translation
        new_coords[0] = N_coord
        new_coords[1]= O_coord
        return new_coords
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
            w.write('!aug-cc-PVDZ\n')
            w.write('%PAL NPROCS 20 END\n')
            w.write('%CASSCF\n')
            w.write('nel 1\n')
            w.write('norb 2\n')
            w.write('nroots 2\n')
            w.write('weights 1 1\n')
            w.write('PTMethod FIC_NEVPT2\n')
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
        para_sub = r'''
#!/bin/bash
for dir in para/*/; do
    echo "Processing $dir"

    for inp in "$dir"/*.inp; do
        [ -f "$inp" ] || continue

        echo "Submitting $inp"

        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=$(basename "${inp%.inp}")
#SBATCH --output=${inp%.inp}.out
#SBATCH --error=${inp%.inp}.err
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --ntasks=20
module load orca  # adjust for your cluster
module load openmpi-ib/intel-2024.0/4.1.6
/sciclone/apps/orca/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca "$inp" > "${inp%.inp}.log"
EOF
    done
done
        '''
        perp_sub = r'''
#!/bin/bash
for dir in perp/*/; do
    echo "Processing $dir"

    for inp in "$dir"/*.inp; do
        [ -f "$inp" ] || continue

        echo "Submitting $inp"

        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=$(basename "${inp%.inp}")
#SBATCH --output=${inp%.inp}.out
#SBATCH --error=${inp%.inp}.err
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --ntasks=20
module load orca  # adjust for your cluster
module load openmpi-ib/intel-2024.0/4.1.6
/sciclone/apps/orca/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca "$inp" > "${inp%.inp}.log"
EOF
    done
done
        '''
        with open(self.od+'/submit_para.sub', 'w') as w:
            w.write(para_sub)
        with open(self.od+'/submit_perp.sub', 'w') as w:
            w.write(perp_sub)
        
        workup_script = r'''
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from itertools import zip_longest
import seaborn as sns
import pandas as pd
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
       parts = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)
       if len(parts) < 2:
           continue
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
   return abs(np.array(matrix)[0][3])*219474.63, abs(np.array(matrix_real)[0][3])*219474.63, soc_energy

# Folder containing ORCA output files
folders = [d for d in glob.glob("*_A/") if os.path.isdir(d)]

def read(directory):
   # Output CSV file
   results = []
   # Grab all ORCA output files
   files = glob.glob(os.path.join(directory, "*deg.log"))
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
                           #print(nevpt2_s0)
                       elif state == 1:
                           nevpt2_s1 = energy
                          # print(nevpt2_s1)
       casscf_matches = re.findall(r"ROOT\s+(\d+):\s+E=\s+(-?\d+\.\d+)\s+Eh",text)
       casscf_energies = {int(root): float(E) for root, E in casscf_matches}
       casscf_s0 = casscf_energies.get(0, np.nan)
       casscf_s1 = casscf_energies.get(1, np.nan)
       results.append([os.path.basename(file),casscf_s0,casscf_s1,soc_imag,soc_real, nevpt2_s0, nevpt2_s1])
      
   energy_s0_cas = []
   energy_s1_cas = []
   energy_s0_nev = []
   energy_s1_nev = []
   angle = []
   soc_imag = []
   soc_real = []
   for row in results:
       energy_s0_cas.append(row[1])
       energy_s1_cas.append(row[2])
       soc_imag.append(row[3])
       soc_real.append(row[4])
       energy_s0_nev.append(row[5])
       energy_s1_nev.append(row[6])  
       name = row[0].split('_')
       angle.append(float(name[0][4:]))
  
   energy_s0_cas = np.array(energy_s0_cas)*219474.63
   energy_s1_cas = np.array(energy_s1_cas)*219474.63
   energy_s0_nev = np.array(energy_s0_nev)*219474.63
   energy_s1_nev = np.array(energy_s1_nev)*219474.63
   energy_s1_nev_cop = energy_s0_nev
   energy_s0_nev_cop = energy_s1_nev
   energy_s1_nev = -min(energy_s0_nev) + energy_s1_nev
   energy_s0_nev = -min(energy_s0_nev) + energy_s0_nev
   energy_s1_cop = energy_s1_cas
   energy_s0_cop = energy_s0_cas
   energy_s1_cas = -min(energy_s0_cas) + energy_s1_cas
   energy_s0_cas = -min(energy_s0_cas) + energy_s0_cas

   plt.figure()
   plt.scatter(angle, energy_s0_cas, label='S0')
   plt.scatter(angle, energy_s1_cas, label='S1')
   plt.xlabel('NO Rotation Angle (Degrees)')
   plt.ylabel('Relative Energy CASSCF (cm^-1)')
   plt.savefig(directory[0:3]+'CAS.png')
   plt.legend()
   plt.figure()
   plt.scatter(angle, energy_s0_nev, label='S0')
   plt.scatter(angle, energy_s1_nev, label='S1')
   plt.xlabel('NO Rotation Angle (Degrees)')
   plt.ylabel('Relative Energy NEVPT2 (cm^-1)')
   plt.savefig(directory[0:3]+'NEV.png')
   plt.legend()

   return energy_s0_cop, energy_s1_cop, angle, soc_imag, soc_real, energy_s0_nev_cop, energy_s1_nev_cop
def make_surface(folders):
   row_s0 = []
   row_s1 = []
   row_s0nev = []
   row_s1nev = []
   row_imag = []
   row_real = []
   row_diff_nev = []
   row_diff_cas = []

   for i in folders:
       e_s0, e_s1, angle, imag, real, e_s0_nev, e_s1_nev = read(i)
       dist = float(i[0:3])
       for es0, es1, a, i, r, e0nev, e1nev in zip(e_s0, e_s1, angle, imag, real, e_s0_nev, e_s1_nev):
           row_s0.append({
                   "Distance": dist,
                   "Angle": a,
                   "Energy_S0": es0})
           row_s1.append({
                   "Distance": dist,
                   "Angle": a,
                   "Energy_S0": es1})
           row_imag.append(({
                   "Distance": dist,
                   "Angle": a,
                   "Energy_S0": i}))
           row_real.append(({
                   "Distance": dist,
                   "Angle": a,
                   "Energy_S0": r}))
           row_s1nev.append(({
                   "Distance": dist,
                   "Angle": a,
                   "Energy_S0": e1nev}))
           row_s0nev.append(({
                   "Distance": dist,
                   "Angle": a,
                   "Energy_S0": e0nev}))
           row_diff_nev.append(({
                   "Distance": dist,
                   "Angle": a,
                   "Energy_S0": e1nev-e0nev}))
           row_diff_cas.append(({
                   "Distance": dist,
                   "Angle": a,
                   "Energy_S0": es1-es0}))
   plotsurf(row_s0, 'S0_surf_cas')
   plotsurf(row_s1, 'S1_surf_cas')
   plotsurf(row_diff_nev, 'S1S0_diff_nev')
   plotsurf(row_diff_cas, 'S1S0_diff_cas')
   plotsurf(row_s1nev, 'S1_surf_nev')
   plotsurf(row_s0nev, 'S0_surf_nev')
   plotsurf(row_imag,'Imag_soc_surf')
   plotsurf(row_real, 'Real_soc_surf')

def plotsurf(rows, name):
   df = pd.DataFrame(rows)
   k= df["Energy_S0"].min()
   df["Energy_S0"] = df["Energy_S0"] - k 
   df = df.drop_duplicates()
   heatmap_data = df.pivot_table(
       index="Distance",
       columns="Angle",
       values="Energy_S0"
   )

   angles = heatmap_data.columns.astype(float)
   distances = heatmap_data.index.astype(float)

   X, Y = np.meshgrid(angles, distances)
   Z = heatmap_data.values
   plt.figure()
   cf = plt.contourf(
       X, Y, Z,
       levels=30,
       cmap='viridis',
       #norm=PowerNorm(gamma=0.5),
   )

   plt.contour(
       X, Y, Z,
       levels=20,
       colors='black',
       linewidths=0.5
   )

   plt.colorbar(cf, label='Relative Energy (cm^-1)')
   #sns.heatmap(heatmap_data,cmap='viridis',vmin=10,vmax=100.5)
   plt.xlabel('Angle (Degrees)')
   plt.ylabel('Distance (Å)')
   plt.savefig(name)

make_surface(folders)'''
        with open(self.od+'/para/workup_para.py', 'w') as w:
            w.write(workup_script)
        with open(self.od+'/perp/workup_perp.py', 'w') as w:
            w.write(workup_script)

k = Rotater(input_file='tz_geom.xyz', target_atom_index=2, output_directory='dz_casscf_nevpt2_surf')
k.kernel(r_step=50, t_step=10)