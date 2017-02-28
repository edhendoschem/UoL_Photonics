This program will read the parameters in config.txt and calculate the first 100 modes of the transversal cross section
of a rectangular waveguide using Marcatili's method. The modes will start at 0,0 (fundamental mode).

It will create two folders:
- mode_data: will contain a .csv file with the mode information
- mode_graphs: Will create a plot of the intensity of the EM field (requires gnuplot installed and added to PATH)

Note on gnuplot installation:
go to:
http://www.gnuplot.info/download.html

double click to install, and keep clicking next until you come to a tab that says 'Select Additional Tasks',
scroll down and select the option 'Add application directory to your PATH environment variable', and then keep
clicking next to complete the installation
