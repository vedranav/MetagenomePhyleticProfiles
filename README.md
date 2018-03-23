#Metagenome Phyletic Profiles

This code supports the paper: Vidulin, V., Smuc, T., Dzeroski, S. and Supek, F. Automated gene function prediction using metagenome data. Under review in Microbiome.

Please cite the paper if you find the code useful.

The code is written in Java and is packed as a NetBeans project.

To reproduce the key experiments from the paper, begin with the ExperimentsFromThePaper.java class. At the top of the class choose experiments that you want to run by setting boolean variables in the section "SELECT THE EXAMPLE(S) FROM THE PAPER" to "true" for each experiment you would like to reproduce. Furthermore, set paths in the section "SET PATHS": 1) "dataDir" should contain path to the data packed within this project, which is needed to compute the results (for example, if you put this project into the folder '/home/user/MetagenomePhyleticProfiles', then the path to the data should be '/home/user/MetagenomePhyleticProfiles/src/data'); 2) "outDir" should contain path to the folder where you want to save the results of experiments.

A prerequisite to successfully run this code is to have R installed on your computer. The code generates and runs R scripts using Rscript. In case that you get an error "Cannot run program "Rscript": error=2, No such file or directory", set the path to Rscript in "utils.RUtils" by adding the path to the line "String[] rscript = {"/ADD PATH HERE/Rscript", rscriptFileName};". You can find out the path by running "type -a Rscript" in a terminal (tested on Ubuntu and MacOS).

Some parts of the code depend upon the specific R packages:
- Fig3a and Fig 3c experiments need 'RandomForest' R package
- Fig1_d_e experiments need 'venneuler' R package

Co-evolution networks script produces a .gexf file that is visualized with Gephi (http://gephi.org) software.

