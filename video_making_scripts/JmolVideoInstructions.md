# Once you have gathered the perpendicular and parallel rotations for your molecule, you can start.

1) Begin by opening the github repository containing the video creation files as well as your folder containg the para or perp rotations. Without loss of generality, the following instructions will refer to just the para, but can be applied just the same to the perp files. 

2) Then, you must create .molden files that are to be used by the Jmol software to generate your orbitals. To do so copy the file in the git repository called conv_molden into the para rotations file.

3) In the terminal, navigate to the directory that contains the para files. Type:
    bash conv_molden
which will create molden files for each of your "para###_deg.log" files. Then, type
    mkdir MoldenFiles
and then
    mv *.molden.output MoldenFiles
then move to the MoldenFiles directory and add from the git repository the following files:
    make_frames
    conv_png
    make_gif.py
Furthermore, create two directories in the current directory called 'homo' and 'lumo.' Finally, add your copy of Jmol (16.3.27 is a known, working version) into the MoldenFiles directory.

4) Before your run any programs, view the 'make_frames' file and find the line that says "mo HOMO-1" or "mo HOMO." When running this script, writing mo HOMO-1 tells Jmol to view the LUMO and thus leaving it like this will generate a video of the LUMO orbitals. To get the HOMO orbitals, remove the -1 so the line says: mo HOMO.


5) Now, in the terminal we can run the script to generate the images that will be turned into a video. Type:
    bash make_frames
(Expect this to take a few minutes. If this process stops midway through, sometimes clicking control + c fixes it).
Once this process finishes, we have to fix the .png formatting. Type:
bash conv_png
You should now see your directory populated with a bunch of png files of your molecule(s) and orbitals.

6) Now that all of your files are in the correct .png formatting, you can greate the gif video. In your terminal, initiate conda and type:
    python make_gif.py
An error may occur telling you that you don't have some conda packages installed. In that event, type:
    conda install NameOfMissingPackage
and rerun the python make_gif.py line.
You should now see a file labeled something like "animated_output.gif"

Repeat steps 5-6 and change the '-1' in the make_frames file to get the .gif for the missing of the HOMO or LUMO.

Repeat steps 2-6 for whichever of the perp/para directories that you did not do originally.