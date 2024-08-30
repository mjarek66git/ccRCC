
The sca program was created to generate heatmaps for gene expression. It can be used in both Windows and Linux systems. 
For Linux, it is necessary to install Mono. 

The program can be run in the GUI or console version. In the GUI version, all options are available on the screen and can be selected by clicking the mouse. 
In the console version, the options must be saved to a configuration file. A sample configuration file called setup.dat is included with the scd distribution.
The attached setup.dat configuration file contains a set of options that allow you to generate the heatmap included in the article
"Copper drives remodeling of metabolic state and progression of ccRCC." (Fig. 6E). To run the console version with the setup.dat configuration file execute the following command:

Windows OS:
scaTerminal.exe -f setup.dat

Linux OS:

mono scaTerminal.exe -f setup.dat


GUI version.

Windows OS:

sca.exe

Linux OS:

mono sca.exe

