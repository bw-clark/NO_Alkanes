# ORCA ROTATION SIMULATION

1) In order to generate rotational dynamics videos for your molecules, you must first obtain an energy minimized geometry. Better geometries include **NEVPT2** and **MP2**, while **CASSCF** may be used. Once you have obtained the gemoetry, find the XYZ coordinates that represents the system in its energy minimized state. This can be accomplished by navigating to the directory in which you ran the optimization and locating a file ending in `.xyz`. Open the `.xyz` file. **If there are multiple trajectories**, copy the coordinates for the atoms from the bottom most trajectory. **If there is one trajectory** copy the coordinates for the geometry. Be sure to only copy the coordinates and atom label like this: 

        O -3.48765083192002     -0.42445509478239	  0.17021825494316

        N  -2.72622306718642     -0.22518484770389     -0.62579766694305

    Once you have the trajectory copied, type

        nano ComplexName_xyz

    and paste your coordinates into this file and save it. 
> **Note:** Click **Control + X**, then **Y**, and then **Enter** to stop viewing the file. Ensure that the `xyz` file is formatted with the N and O as the first two atoms as shown below for NO-CH4:

        N -0.0007482666666666637 0.6172133333333334 -0.0351232

        O 0.0006547333333333516 -0.5400616666666667 0.0307328

        C 3.6404999672145255 -0.007335493931674476 -0.003875630790721231

        H 3.2526469672145257 -0.48252149393167454 0.8926253692092788

        H 3.2531509672145256 -0.5174994939316744 -0.8811636307907212

        H 4.7252149672145265 -0.06333849393167447 -0.002448630790721231

        H 3.3308939672145255 1.0339795060683257 -0.024477630790721232

Finally, open this link: **https://liwt31.github.io/2022/01/02/online_viewer/** and copy the contents of the `_xyz` file into the text box. Locate the atom in the render that appears that you intend to rotate around, and make note of the number associated with its position (The label will be displayed on the atom).

2) For this next section, you will need to have Conda installed with some important packages. Now that you have your coordinates, place your `_xyz` file in the same directory as the `rotate.py` file which is located in the `Rotate_Scripts` folder within the GitHub repository. Open the `rotate.py` file in VSCode and locate the line at the bottom which should look something like this

        k = Rotater(input_file='input_geom_file', target_atom_index=2, output_directory='output_directory_name')

    Change `'input_geom_file'` to the name of your `_xyz` file that you chose earlier. Next, change the value that `target_atom_index` equals to the number (minus 1 b/c of 0-indexing) that you noted earlier in step 1. Finally, name your output directory appropriately where it says `'output_directory_name'`.

3) In your terminal, navigate to the directory with both the `rotate.py` file and your `_xyz` file. Create a new conda environment, activate your conda environment,and install quaternion:

    `conda activate EnvironmentName`

    `conda install quaternion`

You are now ready to run `rotate.py`:

    python rotate.py

You should now see new files created in your directory under the name that you chose for your `'output_directory_name'`. Open that directory and you should see four files:

    para
    perp
    trajectoryPara.xyz
    trajectoryPerp.xyz

You can drag the xyz files into the following website to visualize the rotation to ensure that it is what you are looking for: **https://www.prs.wiki/utilities/xyz-trajectory**.

4) Navigate to the HPC Cluster of your choice (works on Femto and Bora) and type `module load orca` & `module load openmpi-ib/intel-2024.0/4.1.6` into the command line. 

5) Now navigate into either the `para` or `perp` folders. You should see a submit file labeled something like `paraSubmit.sub`. Make sure that at this point, your `para` folders are in a scratch directory within the HPC (scr10, scr20, or scr30). If your files are on your personal machine, you can transfer them over using Globus (If you haven't used Globus, ask someone for help). On the HPC, navigate to your `para` folder and type

        sbatch paraSubmit.sub
(or generally, `sbatch` + whatever your submit file is named in the folder). If everything goes correctly, your jobs should be succesfully running.

6) Once your jobs have finished, you should see several `.log` files in your directory. View one of the log files and make sure that it ran correctly by checking out the bottom, which should have a line that says:

                                     ****ORCA TERMINATED NORMALLY****                             
