# FEC
## Description
This directory primarily contains the original script RNAfec used to calculate force, extension, and energy data using the modified ViennaRNA package. The force, extension, and energy data calculated with this file can be used to construct force-extension curves and perform related analyses. Output obtained with the RNAfec file enclosed in this directory was used in the work "Computational modeling of RNA-protein binding interactions under an external force" by Danielle Wampler and Ralf Bundschuh, which is under review and currently accessible at https://arxiv.org/pdf/2603.22269.

## Contents
The following folders and files are included in this directory :
  - RNAfec.c
  - example_sequences
  - fec_outfiles
  - protein_files

RNAfec.c: c file that takes an RNA sequence as input and uses the modified ViennaRNA package to model folding the specified RNA under an incrementing external force [pN]. The code calculates at each force (between a minimum and maximum, with discrete steps) the ensemble free energy [pNnm] and molecule end-to-end extension [nm]. These forces, extensions, and ensemble free energies are respectively written to columns in a txt file in fec_outfiles. RNAfec.c can also take a protein and free protein conentration [nM] as additional input to output force, extension, and ensemble free energy data in the presence of the specified protein at the specified concentration.

example_sequences: Example RNA sequences that can be used as sequence input with RNAfec.c. p5ab and PolyU were both run as examples with no protein and NM_005857 was an example sequence that was run both with 0 nM HuR and 5 nM HuR.

fec_outfiles: This is the directory to which RNAfec.c directs the force [pN], extension [nm], and ensemble free energy [pNnm] output. The files contained in the folder correspond to the sequences contained in 'example_sequnces' and can be used as a check that RNAfec.c and the modified ViennaRNA are running as expected.

protein_files: This directory contains inverse Kd data for the proteins HuR, RBFOX1, U2AF2, and KHDRBS3 obtained from relative RNAcompete (https://doi.org/10.1038/nbt.1550) data which have been converted to absolute binding affintiy data (for more details, please see Wampler and Bundschuh cited above). This directory also conatins the txt file RBDls which contain the physiological binding domain lengths [10⁻² nm] for HuR, RBFOX1, U2AF2, and KHDRBS3 as determined from crystal structures (see Wampler and Bundschuh cited above for more details). Binding affinity data and protein binding domain length are necessary to run RNAfec.c with protein at a nonzero concentration.

## Requirements and Usage
The modified ViennaRNA package should be installed, configured, and compiled as described in the appropriate ViennaRNA documentation. 

After ViennaRNA is compiled, RNAfold can be compiled with:

"g++ RNAfec.c -I ../src -L ../src/ViennaRNA -o RNAfec.x -lm -lgomp -lRNA -fopenmp"

After compiling, RNAfec can be ran by calling "./RNAfold.x" in the terminal followed by the following arguments:
  1. RNA sequence 
  2. Protein name (Optional)
  3. Protein concentation [nM] (Optional)