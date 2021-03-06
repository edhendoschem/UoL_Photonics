# Table of contents

## 1. [TCP-IP control](#item1)
## 2. [direct_att_control](#item2)
## 3. [lifetime_processor](#item3)


## TCP-IP control <a name = "item1"></a>
This is a client and server jupyter notebook pair, created to continuosly read data from a power log file as it is being
written and then interpolate the laser attenuator position and send it to the server machine to adjust the attenuator.

The attenuator in the server machine is controlled using the mouse and keyboard via pyautogui instead of directly using
device ports. This is done to prevent problems arising from the control program attempting to acquire resources which
are already in use by this python script

<b>Notes:</b>
- The power log file will be read from the end, if the initial lines need to be skipped, change in line_get():
```python
file_handle.seek(0,2)
```
To
```python
file_handle.seek(0,1)
```

## direct_att_control <a name = "item2"></a>
Similarly to TCP-IP, this is a program designed to continuously read from a datafile, extract the power and attenuator 
values from an internal power meter and adjust the attenuator accordingly to maintain the target value. There is no data
transmitted by TCP-IP and no need to manually measure data as the program will self calibrate.

## lifetime_processor <a name = "item3"></a>
This is a simple program to extract the lifetime data from the .txt files created by the Steady State Time Resolved 
Fluorescence Spectrometer. It uses the Trust Region Reflective algorithm. If the curve fitter fails to properly fit the
curve, you may be in the presence of a double exponential or some other behaviour. Lifetime Processor-double_exp will 
attempt to fit double exponentials.

In order to use it, place the lifetime measurement files in the folder with the .ipynb file and click on "run cell,
select below" twice or alternatively press shift+enter twice.

If you want to save the plots as well, change
```python
save_data("results.csv", "\t", False, False)
```
To
```python
save_data("results.csv", "\t", False, True)
```
 Before pressing shift+enter twice or clicking "run cell, select below" button

