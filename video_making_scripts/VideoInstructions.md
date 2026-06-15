# Generating HOMO/LUMO Orbital GIFs

Once you have gathered the perpendicular and parallel rotations for your molecule, you can start.



## 1. Open the Required Directories

Begin by opening:

- The GitHub repository containing the video creation files.
- Your folder containing the parallel (`para`) or perpendicular (`perp`) rotations.

Without loss of generality, the following instructions will refer only to the `para` files, but the same procedure applies to the `perp` files.



## 2. Generate Molden Files

You must create `.molden` files that will be used by Jmol to generate your orbitals using the `.gbw` files.

Copy the file named `conv_molden` from the Git repository into the `para` rotations directory.



## 3. Convert ORCA Outputs to Molden Format

In the terminal, navigate to the directory containing the `para` files and run:

```bash
bash conv_molden
```

This will create Molden files for each of your log files.

Next, create a directory to store the Molden files:

```bash
mkdir MoldenFiles
```

Move the generated files into it:

```bash
mv *.molden.output MoldenFiles
```

Now navigate into the `MoldenFiles` directory and add the following files from the Git repository:

```text
make_frames
conv_png
make_gif.py
```

Create two directories by running:

```bash
mkdir homo
mkdir lumo
```

Finally, place your copy of **Jmol** into the `MoldenFiles` directory.

> **Note:** Jmol version **16.3.27** is a known working version; we've had issues with newer versions of Jmol.



## 4. Select HOMO or LUMO Generation

Before running any programs, open the `make_frames` file and find the line containing:

```text
mo HOMO-1
```

or

```text
mo HOMO
```

### LUMO Generation

If the script contains:

```text
mo HOMO-1
```

Jmol will generate images for the **LUMO**.

### HOMO Generation

To generate images for the **HOMO**, remove the `-1` so the line becomes:

```text
mo HOMO
```



## 5. Generate Orbital Frames

Run:

```bash
bash make_frames
```

> **Note:** This process may take several minutes.

If the process stops midway through, sometimes pressing:

```text
Ctrl + C
```

resolves the issue.

Once the process finishes, fix the PNG formatting:

```bash
bash conv_png
```

You should now see the directory populated with PNG images of your molecule(s) and orbitals.



## 6. Create the GIF Animation

Now that all PNG files are correctly formatted, generate the GIF.

First, initialize your Conda environment, then run:

```bash
python make_gif.py
```

If an error occurs indicating that a required package is missing, install it:

```bash
conda install NameOfMissingPackage
```

Then rerun:

```bash
python make_gif.py
```

You should now see a file named something like ` animated_output.gif`



## 7. Generate the Other Orbital

Repeat **Steps 5–6** after changing the orbital selection in `make_frames`:

- `mo HOMO-1` → generates the **LUMO**
- `mo HOMO` → generates the **HOMO**

This will create the GIF for whichever orbital you have not yet generated.



## 8. Generate the Other Rotation Set

Repeat **Steps 2–6** for whichever rotation directory (`perp` or `para`) you did not process originally.

