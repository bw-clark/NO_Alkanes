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
        with open(self.od+'para/1d_conv.py', 'w') as w:
            script_text = r'''
            '''
            w.write(script_text)
        with open(self.od+'perp/1d_conv.py', 'w') as w:
            w.write(script_text)
        
k = Rotater(input_file='input_geom_file', target_atom_index=2, output_directory='output_directory_name')
k.kernel()