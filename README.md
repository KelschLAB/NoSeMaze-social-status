**INDIVIDUAL DIFFERENCES DRIVE SOCIAL HIERARCHIES IN MOUSE SOCIETIES (REINWALD, GHANAYEM, ET AL., 2025)**

OVERVIEW

This repository contains the scripts used for the manuscript "Individual differences drive social hierarchies in mouse societies" (Reinwald et al., 2025).
It includes MATLAB and R code for preprocessing, analysis, and visualization.<br>
Note: Data (raw and processed) and results are not included (to be shared after publication).
The repository preserves the folder structure so that scripts can be directly reused when data become available.

REPOSITORY STRUCTURE

NoSeMaze_Experiment_Social_Status/<br>
├── analysis/            # analysis scripts (MATLAB and R)<br>
├── config/              # basic config files (meta information on cohorts)<br>      
├── data/                # (empty) raw and processed data will be placed here<br>
├── results/             # (empty) figures, tables, and output files will be stored here<br>
└── src/                 # source scripts (MATLAB and R)<br>

USAGE

Clone the repository:<br>
git clone https://github.com/USERNAME/NoSeMaze_Experiment.git<br>
Set the root path inside the scripts to match your local setup:<br>
myRootPath = 'C:/your/local/path/NoSeMaze_Experiment_Social_Status';<br>
or in R:<br>
myRootPath <- "C:/your/local/path/NoSeMaze_Experiment_Social_Status"<br>
Place your data in the data/ folder (once available).<br>
Run the scripts from within MATLAB or R.<br>

LICENSE

Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)<br>
This work is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.<br>
You are free to:<br>
  - Share — copy and redistribute the material in any medium or format<br>
Under the following terms:<br>
  - Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made. <br>
    You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.<br>
  - NonCommercial — You may not use the material for commercial purposes.<br>
  - NoDerivatives — If you remix, transform, or build upon the material, you may not distribute the modified material.<br>
No additional restrictions — You may not apply legal terms or technological measures that legally restrict others <br>
from doing anything the license permits.<br>
Full license text available at:<br>
https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode<br>


CITATION

If you use these scripts in your work, please cite:<br>
Title<br>
Individual differences drive social hierarchies in mouse societies <br>
Authors <br>
Jonathan R. Reinwald 1,2,3,$, Sarah Ghanayem 1,2,$, David Wolf 1,2, Julia Lebedeva 1,2, Philipp Lebhardt 4, Oliver Gölz 4, Corentin Nelias 1,2, Wolfgang Kelsch 1,2,%<br>
Affiliations<br>
1 Dept. of Psychiatry and Psychotherapy, University Medical Center Mainz, Johannes-Gutenberg University, Untere Zahlbacher Strasse 8, 55131 Mainz, Germany<br>
2 Dept. of Psychiatry and Psychotherapy, Central Institute of Mental Health, Medical Faculty Mannheim, Heidelberg University, Square J5, 68159 Mannheim, Germany<br>
3 Dept. of Neuroimaging, Central Institute of Mental Health, Medical Faculty Mannheim, Heidelberg University, Square J5, 68159 Mannheim, Germany<br>
4 Dept. of Clinical Health Technologies, Institute for Manufacturing Engineering and Automation, Fraunhofer Society, Theodor-Kutzer-Ufer 1-3, 68167 Mannheim, Germany<br>
% corresponding author: Wolfgang Kelsch, wokelsch@uni-mainz.de<br>
$ these authors shared first authorship <br>
(DOI will be added after publication.)

NOTES

Data and results will be made available after publication.<br>
Paths are currently set to myRootPath/NoSeMaze_Experiment.<br>
Please adapt this to your own system before running scripts.<br>
