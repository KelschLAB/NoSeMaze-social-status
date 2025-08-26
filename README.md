"Individual differences drive social hierarchies in mouse societies" (Reinwald et al., 2025)

OVERVIEW

This repository contains the scripts used for the manuscript "Individual differences drive social hierarchies in mouse societies" (Reinwald et al., 2025).
It includes MATLAB and R code for preprocessing, analysis, and visualization.

Note: Data (raw and processed) and results are not included (to be shared after publication).
The repository preserves the folder structure so that scripts can be directly reused when data become available.

REPOSITORY STRUCTURE

NoSeMaze_Experiment/
├── analysis/            # analysis scripts (MATLAB and R)
├── config/              # basic config files (meta information on cohorts)       
├── data/                # (empty) raw and processed data will be placed here
├── results/             # (empty) figures, tables, and output files will be stored here
└── src/                 # source scripts (MATLAB and R)

USAGE

Clone the repository:
git clone https://github.com/USERNAME/NoSeMaze_Experiment.git
Set the root path inside the scripts to match your local setup:
myRootPath = 'C:/your/local/path/NoSeMaze_Experiment';
or in R:
myRootPath <- "C:/your/local/path/NoSeMaze_Experiment"
Place your data in the data/ folder (once available).
Run the scripts from within MATLAB or R.

License

Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)
This work is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
You are free to:
  - Share — copy and redistribute the material in any medium or format
Under the following terms:
  - Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made. 
    You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
  - NonCommercial — You may not use the material for commercial purposes.
  - NoDerivatives — If you remix, transform, or build upon the material, you may not distribute the modified material.
No additional restrictions — You may not apply legal terms or technological measures that legally restrict others 
from doing anything the license permits.
Full license text available at:
https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode


Citation

If you use these scripts in your work, please cite:
Title
Individual differences drive social hierarchies in mouse societies 
Authors 
Jonathan R. Reinwald1,2,3,#, Sarah Ghanayem1,2,#, David Wolf1,2, Julia Lebedeva1,2, Philipp Lebhardt4, Oliver Gölz4, Corentin Nelias1,2, Wolfgang Kelsch1,2,*
Affiliations
1 Dept. of Psychiatry and Psychotherapy, University Medical Center Mainz, Johannes-Gutenberg University, Untere Zahlbacher Strasse 8, 55131 Mainz, Germany
2 Dept. of Psychiatry and Psychotherapy, Central Institute of Mental Health, Medical Faculty Mannheim, Heidelberg University, Square J5, 68159 Mannheim, Germany
3 Dept. of Neuroimaging, Central Institute of Mental Health, Medical Faculty Mannheim, Heidelberg University, Square J5, 68159 Mannheim, Germany
4 Dept. of Clinical Health Technologies, Institute for Manufacturing Engineering and Automation, Fraunhofer Society, Theodor-Kutzer-Ufer 1-3, 68167 Mannheim, Germany
* corresponding author: Wolfgang Kelsch, wokelsch@uni-mainz.de 
# these authors shared last authorship 
(DOI will be added after publication.)

Notes
Data and results will be made available after publication.
Paths are currently set to myRootPath/NoSeMaze_Experiment.
Please adapt this to your own system before running scripts.
