import os
import tqdm
import numpy as np
import re
import openpyxl
import time
import scipy.optimize as optimize
import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy.optimize import curve_fit


#Initial parameters
start = 10.1   #miliseconds (can be changed by tenths of miliseconds 0.1 to 0.9, more decimal places will be ignored)
end = 100      #miliseconds
plots = False  #True, plots all graphs, False doesn't plot any graph but runs faster

blank = OrderedDict()
finger = OrderedDict()
blank_asc = []
finger_asc = []

#Defining single & double exponential functions + initial guess and boundaries
def single_exp(x, a, b, c):
    """
    This function is the single exponential used for curve fitting
    a = Pre exponential factor, initial intensity of the pulse
    b = Lifetime
    c = Curve height adjustment factor
    """
    x = np.float64(x)
    a = np.float64(a)
    b = np.float64(b)
    c = np.float64(c)
    return a*np.exp(-x/b) + c


def double_exp(x, a1, a2, b1, b2, c):
    """
    This function is the single exponential used for curve fitting
    a1 = Pre exponential factor of first lifetime
    a2 = Pre exponential factor of second lifetime
    b1 = Lifetime 1
    b2 = Lifetime 2
    c = curve height adjustment 
    """
    x = np.float64(x)
    a1 = np.float64(a1)
    a2 = np.float64(a2)
    b1 = np.float64(b1)
    b2 = np.float64(b2)
    c = np.float64(c)
    return a1 * np.exp(-x/b1) + a2* np.exp(-x/b2) + c




#Information gatherers
blank_file = re.compile(r"b.+\.asc",re.IGNORECASE) #This extracts the name of blank files
finger_file = re.compile(r"f.+\.asc",re.IGNORECASE) #This extracts the name of finger files
date_regex = re.compile(r"(\d+).(\d+).(\d+)") #This extracts the date and separates it in month, day, year
time_regex = re.compile(r"\d+:\d+:\d+") #This extracts the time
day_loc = re.compile(r".+\\day(\d{1,2})\\.+") #this extracts the day number

start_time = time.time() #Stores the time at the beginning of this stage of the program
#Building the file locations
#The command os.walk basically displays all the folders, subfolders and files of the current directory. We iterate over
#this information to store the file locations and names and store it in two separate lists: finger_asc and blank_asc and
#then writes them in a dictionary named blank or finger depending on the type of file. 
#tqdm.tqdm is basically the module used to display the loading bar
for folder, subfolder, file in tqdm.tqdm(os.walk(os.getcwd()), desc = "Searching for file locations"):
    for elements in file:
        if blank_file.search(elements) != None:
            blank_asc.append(elements)    
            blank[folder] = blank_asc
        elif finger_file.search(elements) != None:
            finger_asc.append(elements)    
            finger[folder] = finger_asc
    finger_asc = []        
    blank_asc = []


#Creating the excel spreadsheet   
wb = openpyxl.Workbook() #Creation of excel workbook
blank_sheet = wb.active #Selecting the active spreadsheet in workbook and assigning it to a variable
blank_sheet.title = 'blanks' #Renaming active spreadsheet, the default name is sheet1
wb.create_sheet(index = 1, title = 'finger') #Creating another spreadsheet with the name 'finger'
finger_sheet = wb.get_sheet_by_name('finger') #Setting a variable to refer to the excel spreadsheet called 'finger
blank_sheet.freeze_panes = 'A2' #Freezing the first row in blank page
finger_sheet.freeze_panes = 'A2' #Freezing the first row in finger page
blank_sheet['A1'] = 'File Name' #This is equivalent to writing in cell A1 in excel sheet 'blanks' the phrase 'File name'
blank_sheet['B1'] = 'Date'
blank_sheet['C1'] = 'Day'
blank_sheet['D1'] = 'Time'
blank_sheet['E1'] = 'a'
blank_sheet['F1'] = 'a_stdev'
blank_sheet['G1'] = 'b (lifetime)'
blank_sheet['H1'] = 'b_stdev'
blank_sheet['I1'] = 'c'
blank_sheet['J1'] = 'c_stdev'
blank_sheet['K1'] = 'a1_dou'
blank_sheet['L1'] = 'a1_dou_stdev'
blank_sheet['M1'] = 'a2_dou'
blank_sheet['N1'] = 'a2_dou_stdev'
blank_sheet['O1'] = 'b1_dou (lifetime1)'
blank_sheet['P1'] = 'b1_dou_stdev'
blank_sheet['Q1'] = 'b2_dou (lifetime2)'
blank_sheet['R1'] = 'b2_dou_stdev'
blank_sheet['S1'] = 'c_dou'
blank_sheet['T1'] = 'c_dou_stdev'
finger_sheet['A1'] = 'File Name'
finger_sheet['B1'] = 'Date'
finger_sheet['C1'] = 'Day'
finger_sheet['D1'] = 'Time'
finger_sheet['E1'] = 'a'
finger_sheet['F1'] = 'a_stdev'
finger_sheet['G1'] = 'b (lifetime)'
finger_sheet['H1'] = 'b_stdev'
finger_sheet['I1'] = 'c'
finger_sheet['J1'] = 'c_stdev'
finger_sheet['K1'] = 'a1_dou'
finger_sheet['L1'] = 'a1_dou_stdev'
finger_sheet['M1'] = 'a2_dou'
finger_sheet['N1'] = 'a2_dou_stdev'
finger_sheet['O1'] = 'b1_dou (lifetime1)'
finger_sheet['P1'] = 'b1_dou_stdev'
finger_sheet['Q1'] = 'b2_dou (lifetime2)'
finger_sheet['R1'] = 'b2_dou_stdev'
finger_sheet['S1'] = 'c_dou'
finger_sheet['T1'] = 'c_dou_stdev'



#Analyzing and writing the information
i=2 #Position setter, this is used to start writing the information in the excel spreadsheet. 
#Starts at 2 because row 1 is used for the description

#Blanks
csvfile = open("blank_results.csv","w") #Creates a file with that name in write mode (hence the "w")
csvfile.write('File Name'+','
              +'Date'+','
              +'Day'+','
              +'Time'+','
              +'a'+','
              +'a_stdev'+','
              +'b'+','
              +'b_stdev'+','
              +'c'+','
              +'c_stdev'+','
              +'a1_dou'+','
              +'a1_dou_stdev'+','
              +'a2_dou'+','
              +'a2_dou_stdev'+','
              +'b1_dou'+','
              +'b1_dou_stdev'+','
              +'b2_dou'+','
              +'b2_dou_stdev'+','
              +'c_dou'+','
              +'c_dou_stdev'+','
              +'\n')
#The previous was all a command, separated in several lines to improve readability, and its basically writing in the text
#file those words, separated by a comma (this is basically a comma separated value file), the \n at the end is used to 
#tell the cursor to move to the next line in a text file.

#We start to iterate over the 'keys' of the dictionary (in this case the file paths) and inside each file path we iterate
#the file names contained, opening each, reading their contents and assigning them to the y variable.
for key,value in tqdm.tqdm(blank.items(), desc = "Fitting and plotting blank measurements", ncols = 110):
    for files in value:
        element = str(key)+"\\"+str(files)
        file = open(element) #Blank Reading
        info = file.readlines()  #File information stored as lines
        file.close()
        x = np.linspace(start,end,len(info[11+int(start*10):10+int(end*10)]), dtype = np.float64)
        y = np.int64(info[11+int(start*10):10+int(end*10)])
        date = date_regex.search(info[3]).group(2)+"/"+date_regex.search(info[3]).group(1)\
        +"/"+date_regex.search(info[3]).group(3)
        
        #Curve fitting single exponential
        #This uses Trust Region Reflective algorithm to find the solution which we will then put as starting point 
        #in another algorithm for increased accuracy. Same procedure for double exponential
        popt_0, pcov_0 = curve_fit(single_exp, x, y, 
                               p0 = (np.mean(y[:10]), 9, np.mean(y[:-10])),
                               method = 'trf',
                               loss = 'linear',
                               bounds = ([np.amin(y)*0.8, 0, np.amin(y)*0.8],[np.amax(y)*1.2, 15, np.amax(y)*1.2]),
                               max_nfev = 1000000)
        
        #This uses Levenberg-Marquardt algorithm using as starting point the previous value found by the Trust Region 
        #Reflective algorithm.
        popt, pcov = curve_fit(single_exp, x, y, 
                               p0 = (popt_0[0], popt_0[1], popt_0[2]), 
                               method = 'lm', maxfev = 1000000)

        #This uses dogleg algorithm with rectangular trust regions (Currently not in use)
        #popt, pcov = curve_fit(single_exp, x, y, 
        #                       p0 = (np.mean(y[10:20]), 9, np.mean(y[:-10])), 
        #                       method = 'trf',
        #                       loss = 'linear',
        #                       max_nfev = 10000000)
        stdev = np.sqrt(np.diag(pcov)) #This calculates the standard deviation
        
        #Curve fitting double exponential
        #This uses Trust Region Reflective algorithm to find the starting point
        popt2_0, pcov2_0 = curve_fit(double_exp, x, y, 
                                 p0 = (np.mean(y[:10]) * 0.5, np.mean(y[:10]) * 0.5, 9, 9, np.mean(y[:-10])),
                                 method = 'trf',
                                 loss = 'linear',
                                 bounds = ([0, 0, 0, 0, np.amin(y)*0.8],[np.amax(y)*5, np.amax(y)*5, 15, 15, np.amax(y)*1.2]),
                                 max_nfev = 1000000)
        
        #This uses Levenberg-Marquardt algorithm to find the solution
        popt2, pcov2 = curve_fit(double_exp, x, y, 
                                 p0 = (popt2_0[0], popt2_0[1], popt2_0[2], popt2_0[3], popt2_0[4]), 
                                 method = 'lm', maxfev = 1000000)
        
        #This uses dogleg algorithm with rectangular trust regions
        #popt2, pcov2 = curve_fit(double_exp, x, y, 
        #                         p0 = (np.mean(y[10:20]), 0.5, 9, 9, np.mean(y[:-10])),
        #                         bounds = ([np.amin(y), 0, 0, 0, np.amin(y)],[np.amax(y), 1, 15, 15, np.amax(y)]),
        #                         method = 'dogbox',
        #                         loss = 'linear',
        #                         max_nfev = 10000000)
        stdev2 = np.sqrt(np.diag(pcov2))
        
        #Excel spreadsheet results
        #This basically writes the data in the appropriate column
        blank_sheet['A{}'.format(i)] = files
        blank_sheet['B{}'.format(i)] = date
        blank_sheet['C{}'.format(i)] = 'Day {}'.format(day_loc.search(key).group(1))
        blank_sheet['D{}'.format(i)] = info[4]
        blank_sheet['E{}'.format(i)] = popt[0]
        #This if clause is used to catch infinity values or divisions by 0 which excel doesnt like and converts them to
        #strings 'inf' or 'nan'. 
        if np.isinf(stdev[0]) or np.isnan(stdev[0]):
            blank_sheet['F{}'.format(i)] = str(stdev[0])
        else:
            blank_sheet['F{}'.format(i)] = (stdev[0])
            
        blank_sheet['G{}'.format(i)] = popt[1]
        if np.isinf(stdev[1]) or np.isnan(stdev[1]):
            blank_sheet['H{}'.format(i)] = str(stdev[1])
        else:
            blank_sheet['H{}'.format(i)] = (stdev[1])
            
        blank_sheet['I{}'.format(i)] = popt[2]
        if np.isinf(stdev[2]) or np.isnan(stdev[2]):
            blank_sheet['J{}'.format(i)] = str(stdev[2])
        else:
            blank_sheet['J{}'.format(i)] = (stdev[2])

        blank_sheet['K{}'.format(i)] = popt2[0]
        if np.isinf(stdev2[0]) or np.isnan(stdev2[0]):
            blank_sheet['L{}'.format(i)] = str(stdev2[0])
        else:
            blank_sheet['L{}'.format(i)] = (stdev2[0])
        
        blank_sheet['M{}'.format(i)] = popt2[1]
        if np.isinf(stdev2[1]) or np.isnan(stdev2[1]):
            blank_sheet['N{}'.format(i)] = str(stdev2[1])
        else:
            blank_sheet['N{}'.format(i)] = (stdev2[1])

        blank_sheet['O{}'.format(i)] = popt2[2]
        if np.isinf(stdev2[2]) or np.isnan(stdev2[2]):
            blank_sheet['P{}'.format(i)] = str(stdev2[2])
        else:
            blank_sheet['P{}'.format(i)] = (stdev2[2])
            
        blank_sheet['Q{}'.format(i)] = popt2[3]
        if np.isinf(stdev2[3]) or np.isnan(stdev2[3]):
            blank_sheet['R{}'.format(i)] = str(stdev2[3])
        else:
            blank_sheet['R{}'.format(i)] = (stdev2[3])

        blank_sheet['S{}'.format(i)] = popt2[4]
        if np.isinf(stdev2[4]) or np.isnan(stdev2[4]):
            blank_sheet['T{}'.format(i)] = str(stdev2[4])
        else:
            blank_sheet['T{}'.format(i)] = (stdev2[4])


        #CSV results
        #This converts the results to strings and writes them to the csv file, adding a comma between each value and the
        #endline character '\n' at the end to tell the cursor to move to the next line for the next iteration.
        csvfile.write(str(files)+','+str(date)+','+str('Day {}'.format(day_loc.search(key).group(1)))+','
                          +str(info[4]).rstrip("\n")+','
                          +str(popt[0])+','
                          +str(stdev[0])+','
                          +str(popt[1])+','
                          +str(stdev[1])+','
                          +str(popt[2])+','
                          +str(stdev[2])+','
                          +str(popt2[0])+','
                          +str(stdev2[0])+','
                          +str(popt2[1])+','
                          +str(stdev2[1])+','
                          +str(popt2[2])+','
                          +str(stdev2[2])+','
                          +str(popt2[3])+','
                          +str(stdev2[3])+','
                          +str(popt2[4])+','
                          +str(stdev2[4])+','
                          +'\n')
        #Plots
        if plots == True:
            #Single Exp plot
            yresult = single_exp(x,popt[0],popt[1],popt[2])
            fig1 = plt.figure(figsize = (10,10), dpi = 100) #Specifies the dots per inch in the plot or 
            #how far can you zoom in in the plot's image
            plt.hold(True)
            plt.subplot(211) #creates another plot in the same window, this will be the residuals plot.
            plt.plot(x,y,"ro", label = "Original Data")
            plt.plot(x, yresult,color = 'k', label = "Fitted data", linewidth = 3.0)
            plt.ylabel("Signal")
            plt.xlabel("time in ms")
            plt.legend()
            plt.subplot(212)
            max_val = np.amax(np.absolute(y-yresult)) #Finds the maximum value of the residual
            normalized_res = (y-yresult)/max_val #Divides all residual values by the maximum values in order to obtain
            #a normalized result between -1 to 1
            plt.plot(x,normalized_res)
            plt.ylabel("Residuals")
            plt.ylim(-1,1)
            plt.xlabel("time in ms")
            modfile = files[:-4]+"_"+"single_exp_blank"+".png" #Saves the file, using the measurement filename plus 
            #'single_exp_blank' and the file extension '.png'
            plt.savefig(str(key)+"\\"+modfile)
        
            #Double Exp plot
            yresult = double_exp(x, popt2[0], popt2[1], popt2[2], popt2[3], popt2[4])
            fig2 = plt.figure(figsize = (10,10), dpi = 100)
            plt.subplot(211)
            plt.plot(x,y,"ro", label = "Original Data")
            plt.plot(x, yresult,color = 'k', label = "Fitted data", linewidth = 3.0)
            plt.ylabel("Signal")
            plt.xlabel("time in ms")
            plt.legend()
            plt.subplot(212)
            max_val = np.amax(np.absolute(y-yresult))
            normalized_res = (y-yresult)/max_val
            plt.plot(x,normalized_res)
            plt.ylabel("Residuals")
            plt.ylim(-1,1)
            plt.xlabel("time in ms")
            modfile = files[:-4]+"_"+"double_exp_blank"+".png"
            plt.savefig(str(key)+"\\"+modfile)
            plt.close("all")
        i+=1
csvfile.close() #Closes the csv file to prevent any 'file in use' type of errors when trying to open it with other
#programs

#Finger measurements
i = 2
csvfile = open("finger_results.csv","w")
csvfile.write('File Name'+','
              +'Date'+','
              +'Day'+','
              +'Time'+','
              +'a'+','
              +'a_stdev'+','
              +'b'+','
              +'b_stdev'+','
              +'c'+','
              +'c_stdev'+','
              +'a1_dou'+','
              +'a1_dou_stdev'+','
              +'a2_dou'+','
              +'a2_dou_stdev'+','
              +'b1_dou'+','
              +'b1_dou_stdev'+','
              +'b2_dou'+','
              +'b2_dou_stdev'+','
              +'c_dou'+','
              +'c_dou_stdev'+','
              +'\n')

for key,value in tqdm.tqdm(finger.items(), desc = "Fitting and plotting finger measurements", ncols = 110):
    for files in value:
        element = str(key)+"\\"+str(files)
        file = open(element) #Finger Reading
        info = file.readlines()  #File information stored as lines
        file.close()
        x = np.linspace(start,end,len(info[11+int(start*10):10+int(end*10)]), dtype = np.float64)
        y = np.int64(info[11+int(start*10):10+int(end*10)])
        date = date_regex.search(info[3]).group(2)+"/"+date_regex.search(info[3]).group(1)\
        +"/"+date_regex.search(info[3]).group(3)
        
        #Curve fitting single exponential
        #This uses Trust Region Reflective algorithm to find the starting point
        popt_0, pcov_0 = curve_fit(single_exp, x, y, 
                               p0 = (np.mean(y[:10]), 9, np.mean(y[:-10])),
                               method = 'trf',
                               loss = 'linear',
                               bounds = ([np.amin(y)*0.8, 0, np.amin(y)*0.8],[np.amax(y)*1.2, 15, np.amax(y)*1.2]),
                               max_nfev = 1000000)
        
        #This uses Levenberg-Marquardt algorithm to find the solution
        popt, pcov = curve_fit(single_exp, x, y, 
                               p0 = (popt_0[0], popt_0[1], popt_0[2]), 
                               method = 'lm', maxfev = 1000000)

        #This uses dogleg algorithm with rectangular trust regions
        #popt, pcov = curve_fit(single_exp, x, y, 
        #                       p0 = (np.mean(y[10:20]), 9, np.mean(y[:-10])), 
        #                       method = 'trf',
        #                       loss = 'linear',
        #                       max_nfev = 10000000)
        stdev = np.sqrt(np.diag(pcov))
        
        #Curve fitting double exponential
        #This uses Trust Region Reflective algorithm to find the starting point
        popt2_0, pcov2_0 = curve_fit(double_exp, x, y, 
                                 p0 = (np.mean(y[:10]) * 0.5, np.mean(y[:10]) * 0.5, 9, 9, np.mean(y[:-10])),
                                 method = 'trf',
                                 loss = 'linear',
                                 bounds = ([0, 0, 0, 0, np.amin(y)*0.8],[np.amax(y)*5, np.amax(y)*5, 15, 15, np.amax(y)*1.2]),
                                 max_nfev = 1000000)
        
        #This uses Levenberg-Marquardt algorithm to find the solution
        popt2, pcov2 = curve_fit(double_exp, x, y, 
                                 p0 = (popt2_0[0], popt2_0[1], popt2_0[2], popt2_0[3], popt2_0[4]), 
                                 method = 'lm', maxfev = 1000000)
        
        #This uses dogleg algorithm with rectangular trust regions
        #popt2, pcov2 = curve_fit(double_exp, x, y, 
        #                         p0 = (np.mean(y[10:20]), 0.5, 9, 9, np.mean(y[:-10])),
        #                         bounds = ([np.amin(y), 0, 0, 0, np.amin(y)],[np.amax(y), 1, 15, 15, np.amax(y)]),
        #                         method = 'dogbox',
        #                         loss = 'linear',
        #                         max_nfev = 10000000)
        stdev2 = np.sqrt(np.diag(pcov2))
        
        #Excel spreadsheet results
        finger_sheet['A{}'.format(i)] = files
        finger_sheet['B{}'.format(i)] = date
        finger_sheet['C{}'.format(i)] = 'Day {}'.format(day_loc.search(key).group(1))
        finger_sheet['D{}'.format(i)] = info[4]
        finger_sheet['E{}'.format(i)] = popt[0]
        if np.isinf(stdev[0]) or np.isnan(stdev[0]):
            finger_sheet['F{}'.format(i)] = str(stdev[0])
        else:
            finger_sheet['F{}'.format(i)] = (stdev[0])
            
        finger_sheet['G{}'.format(i)] = popt[1]
        if np.isinf(stdev[1]) or np.isnan(stdev[1]):
            finger_sheet['H{}'.format(i)] = str(stdev[1])
        else:
            finger_sheet['H{}'.format(i)] = (stdev[1])
            
        finger_sheet['I{}'.format(i)] = popt[2]
        if np.isinf(stdev[2]) or np.isnan(stdev[2]):
            finger_sheet['J{}'.format(i)] = str(stdev[2])
        else:
            finger_sheet['J{}'.format(i)] = (stdev[2])

        finger_sheet['K{}'.format(i)] = popt2[0]
        if np.isinf(stdev2[0]) or np.isnan(stdev2[0]):
            finger_sheet['L{}'.format(i)] = str(stdev2[0])
        else:
            finger_sheet['L{}'.format(i)] = (stdev2[0])
        
        finger_sheet['M{}'.format(i)] = popt2[1]
        if np.isinf(stdev2[1]) or np.isnan(stdev2[1]):
            finger_sheet['N{}'.format(i)] = str(stdev2[1])
        else:
            finger_sheet['N{}'.format(i)] = (stdev2[1])

        finger_sheet['O{}'.format(i)] = popt2[2]
        if np.isinf(stdev2[2]) or np.isnan(stdev2[2]):
            finger_sheet['P{}'.format(i)] = str(stdev2[2])
        else:
            finger_sheet['P{}'.format(i)] = (stdev2[2])
            
        finger_sheet['Q{}'.format(i)] = popt2[3]
        if np.isinf(stdev2[3]) or np.isnan(stdev2[3]):
            finger_sheet['R{}'.format(i)] = str(stdev2[3])
        else:
            finger_sheet['R{}'.format(i)] = (stdev2[3])

        finger_sheet['S{}'.format(i)] = popt2[4]
        if np.isinf(stdev2[4]) or np.isnan(stdev2[4]):
            finger_sheet['T{}'.format(i)] = str(stdev2[4])
        else:
            finger_sheet['T{}'.format(i)] = (stdev2[4])


        #CSV results
        csvfile.write(str(files)+','+str(date)+','+str('Day {}'.format(day_loc.search(key).group(1)))+','
                          +str(info[4]).rstrip("\n")+','
                          +str(popt[0])+','
                          +str(stdev[0])+','
                          +str(popt[1])+','
                          +str(stdev[1])+','
                          +str(popt[2])+','
                          +str(stdev[2])+','
                          +str(popt2[0])+','
                          +str(stdev2[0])+','
                          +str(popt2[1])+','
                          +str(stdev2[1])+','
                          +str(popt2[2])+','
                          +str(stdev2[2])+','
                          +str(popt2[3])+','
                          +str(stdev2[3])+','
                          +str(popt2[4])+','
                          +str(stdev2[4])+','
                          +'\n')
        #Plots
        if plots == True:
            #Single Exp plot
            yresult = single_exp(x,popt[0],popt[1],popt[2])
            fig1 = plt.figure(figsize = (10,10), dpi = 100)
            plt.hold(True)
            plt.subplot(211)
            plt.plot(x,y,"ro", label = "Original Data")
            plt.plot(x, yresult,color = 'k', label = "Fitted data", linewidth = 3.0)
            plt.ylabel("Signal")
            plt.xlabel("time in ms")
            plt.legend()
            plt.subplot(212)
            max_val = np.amax(np.absolute(y-yresult))
            normalized_res = (y-yresult)/max_val
            plt.plot(x,normalized_res)
            plt.ylabel("Residuals")
            plt.ylim(-1,1)
            plt.xlabel("time in ms")
            modfile = files[:-4]+"_"+"single_exp_finger"+".png"
            plt.savefig(str(key)+"\\"+modfile)
        
            #Double Exp plot
            yresult = double_exp(x, popt2[0], popt2[1], popt2[2], popt2[3], popt2[4])
            fig2 = plt.figure(figsize = (10,10), dpi = 100)
            plt.subplot(211)
            plt.plot(x,y,"ro", label = "Original Data")
            plt.plot(x, yresult,color = 'k', label = "Fitted data", linewidth = 3.0)
            plt.ylabel("Signal")
            plt.xlabel("time in ms")
            plt.legend()
            plt.subplot(212)
            max_val = np.amax(np.absolute(y-yresult))
            normalized_res = (y-yresult)/max_val
            plt.plot(x,normalized_res)
            plt.ylabel("Residuals")
            plt.ylim(-1,1)
            plt.xlabel("time in ms")
            modfile = files[:-4]+"_"+"double_exp_finger"+".png"
            plt.savefig(str(key)+"\\"+modfile)
            plt.close("all")
        i+=1
csvfile.close()

wb.save('blanks_finger_results.xlsx') #Saves the excel spreadsheet
end_time = time.time() #Saves the time at the end of the program

print("Done")
print("Elapsed time: {} seconds".format(end_time-start_time)) #Prints the time by substracting time at the end minus time
#at the start