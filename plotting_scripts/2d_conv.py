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
folders = ['3.42_A', '3.53_A', '3.64_A', '3.75_A', '3.86_A','3.97_A']

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
   row_diff = []
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
           row_diff.append(({
                   "Distance": dist,
                   "Angle": a,
                   "Energy_S0": abs(e0nev-e1nev)}))
   plotsurf(row_s0, 'S0_surf_cas')
   plotsurf(row_s1, 'S1_surf_cas')
   plotsurf(row_diff, 'S1S0_diff')
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

make_surface(folders)


