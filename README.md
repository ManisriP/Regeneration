#### Github repository

##### DREM analysis
Input files: 
Gene expression file: gse95135_phdrem_5000.txt;
TF-Gene association file:  mouse_predicted.txt.gz;
Arguments setttings file for DREM: gse95135_drem_arg_ph5000.txt;

Running DREM in batch mode: java -mx1024M -jar drem.jar -b gse95135_drem_arg_ph5000.txt out_ph_5000.txt

Visualizing output model: Run the following command 
java -mx1024M -jar drem.jar -d gse95135_drem_arg_ph5000.txt

DREM GUI opens. In the Saved Model File menu, enter out_ph_5000.txt(i.e the output file from the batch mode)

###### Scripts for WGCNA and input file preparation for DREM analysis are available within the script workflow.R

###### "model1.ode" and "model2.ode" correspond to input XPP ode files for mathematical models.



