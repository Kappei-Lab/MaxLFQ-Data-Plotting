# LFQ Script 

LFQ Script for Kappei Lab

Follow the steps given below to run the analysis (python3 and other packages given in the requirements.txt should be installed)

1. Copy the code file (MaxLFQ_script.py) and MaxLFQ_plot_info.json file to the folder with proteinGroups.txt file.(Make sure the folder name has some details on the experiment)

2. Open terminal at the folder with the proteinGroups.txt file.

(OR) 

Open terminal and use the command cd to go to the folder with the proteinGroups.txt file. 

cd path/to/folder

Note: In case of windows the path is path\to\folder


3. Edit the json file copied into the folder for plot title, labels etc. 

4. Type the below line on the terminal and press enter. 

python MaxLFQ_script.py proteinGroups.txt MaxLFQ_plot_info.json


5. Once the code is running it will print a few lines including column names etc. Next, it will ask you to input the order of treatment columns for which you have to type in the numbers (see "no" column next to "LFQ" column that is displayed) with a space.(Enter y if you selected the right columns or else press n) 

6. Follow the same steps for selecting the control samples. 

7. Next, enter the minimum number of replicates that should have a value (not zero) in either the control or treatment samples for every protein. (Recommended to keep the number of imputed values to <10% of the total data points). Press enter. 

8. Wait for the analysis to run. After it has finished check the Analysis folder to view the plot and excel file generated. 
