ASSESSMe.py
===============================================
Assess Statistical Estimates from Shotgun Proteomics Scoring Methods
-----------------------------------------------
The program helps to validate the accuracy, or calibration, of statistical scores estimated for peptide-spectrum matches (PSMs) or peptides in Shotgun proteomics. The procedure, called the semi-labelled calibration test, is described in:  "On Using Samples of Known Protein Content to Assess the Statistical Calibration of Scores Assigned to Peptide-Spectrum Matches in Shotgun Proteomics" by **Granholm et al.**, *Journal of Proteome Research*, 2011. The program, ASSESSMe.py, is a python program ran from the command line.

Workflow:
- Get MS/MS spectra from a sample of known protein content. One example of such data is the ISB18 mix, available on https://regis-web.systemsbiology.net/PublicDatasets/

- Get a fasta file containing only those proteins known to be in the purified sample, the *sample* database.

- Produce a target *biparite* database by appending *entrapment* sequences to the *sample* database. This is done with ASSESSMe.py in the 'database' mode, which also produces a reversed decoy database.

- Search the spectra against the target *bipartite* database and the decoy database, and estimate *p* values for each PSM or peptide.

- Run ASSESSMe.py in 'calibration' mode, and give it a file with a *p* value and a protein id, on each line.
