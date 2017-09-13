# Table of contents

## 1. [TCP-IP control](#item1)

## TCP-IP control <a name = "item1"></a>
This is a client and server jupyter notebook pair, created to continuosly read data from a power log file as it is being
written and then interpolate the laser attenuator position and send it to the server machine to adjust the attenuator.

The attenuator in the server machine is controlled using the mouse and keyboard via pyautogui instead of directly using
device ports. This is done to prevent problems arising from the control program attempting to acquire resources which
are already in use by this python script

<b>Notes:</b>
- The power log file will be read from the beginning, if the initial lines need to be skipped, change in line_get():
```python
file_handle.seek(0,1)
```
To
```python
file_handle.seek(0,2)
```