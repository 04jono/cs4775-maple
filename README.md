# Data
Tree generated from 500 samples from December 2019 to February 4 2020 is in the Feb42020/ directory. Aligned FASTA files may be found in the files/ directory.

Tree generated from 7000 samples from May 1 2021 is in the May12021/ directory. The aligned FASTA files are too large to upload to GitHub.

# Data pipeline

1. Download raw FASTA file from GISAID

2. Perform multiple sequence alignment with MAFFT to make sure nucleotides are in the correct positions and sequences are the same length:
   `mafft --thread 4 --auto feb4.fasta > aligned_feb4.fasta`

3. Run `createMapleFile.py` using an aligned reference (Wuhan-Hu-1) to create a MAPLE alignment file:

`python createMapleFile.py --input <input fasta> --output <output file>`

4. Run `MAPLEv0.6.11.py` on the MAPLE alignment file:

` pypy3 MAPLEv0.6.11.py --input <input file> --output <output directory>`

5. Use `rf_compare_phylo_trees.py` or `tree_visualizer.py` for analysis.
