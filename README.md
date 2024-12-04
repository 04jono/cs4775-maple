# Data pipeline

1. Download raw FASTA file from GISAID

2. Perform multiple sequence alignment with MAFFT to make sure nucleotides are in the correct positions and sequences are the same length:
   `mafft --thread -4 --auto feb4.fasta > aligned_feb4.fasta`

3. Run `createMapleFile.py` using an aligned reference (Wuhan-Hu-1) to create a MAPLE alignment file.

4. Run `MAPLEv0.6.11.py` on the MAPLE alignment file.
