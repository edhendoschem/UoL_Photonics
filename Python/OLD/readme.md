# Table of contents
## 1. [0-Planar waveguide.ipynb](#0-Planar)
## 2. [0-Planar-Rectangular waveguide 0.45](#planar0.45)
## 3. [placeholder](#item3)


##<b>Important note:</b>
Files in the old folder may or may not be correct and may contain several iterations of the same solution 

## 0-Planar waveguide.ipynb <a name = "0-Planar"></a>		
Calculates the parameters of a planar waveguide, splits the domain in 4 equal parts and assigns each to a processing
thread. 

**Note:** This could have been solved orders of magnitude faster using numpy wihtout the need for parallel processing

##0-Planar-Rectangular waveguide 0.45.ipynb <a name = "planar0.45"></a>
Uses classes to solve planar waveguides and rectangular waveguides using Marcatili's method, provides a range of options
to customize plot output and also displays progress bars.