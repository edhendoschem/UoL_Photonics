# Cython: A quick tutorial

## What is Cython?

Cython is a way to compile C extensions using syntax similar to Python, this allows you to have C levels of speed
in your Python code. 

More info on their website:
http://cython.org/

## When to use Cython?

If you have a problem that requires extensive use of loops, or computations using the same type of elements (i.e.
numbers only) and you have tried optimizing with numpy and it still is too slow.

## How do I install Cython?

If you have installed Python Anaconda then you already have Cython installed, otherwise type:
```
pip install cython
```

## Before starting

We will be using jupyter notebook in this tutorial. In order to open jupyter notebook follow this procedure:

In Windows:
1) Go to start -> run and type cmd

2) Go to the directory where you keep your Python scripts using `cd yourdirectoryname`

3) type `jupyter notebook` and press enter

In linux:
1) Open your terminal in the directory you keep your Python scripts

2) type `jupyter notebook` and press enter

### jupyter notebook and cell magic

The reason why we will use jupyter notebook is due to the convenience of the commands you can execute from the 
notebook called 'cell magic'. Here are some examples:

1) Create a file from jupyter notebook:
```
%%file myscript.py
program body
```
The previous command creates a file in the current directory with the .py extension containing whatever is written in
that particular cell

2) Run a Python script from jupyter notebook:
```
%run myscript.py
```
Runs the script using Python as if you had gone to command console or cmd and typed 'python myscript.py' and 
pressed enter

3) Run a system command from jupyter notebook:
In windows:
```
%%cmd
system command (such as dir, cd, etc..)
```
In Linux:
```
!command (such as ls, etc...)
```
Note: You can also use '!' in Windows but it doesnt work as well as '%%cmd'

4) Measure execution time of a statement:
```
%timeit -n 1000 -r 4  somestatement
```
runs 'somestatement' 1000 times and repeats it 4 times. Displays the speed of the fastest run

```
%timeit -n 1000 somestatement
```
runs 'somestatement' 1000 times and repeats it 3 times (the default number). Displays the speed of the fastest run

### Useful jupyter notebook shortcuts
'ctrl + enter' Runs the code in the current cell
'shift + enter' Runs the code in the current cell and creates a new cell

## Cython files: .pyd and .pyx

Cython files have a different extension than Python files. They are **.pyd** to store declarations (we will ignore
these for now) and **.pyx** to store the definitions. If we want to compile our program to C we need to place the
parts to optimize inside a .pyx file, this can be done very easily with the cell magic '%%file ourextname.pyx' and
then just write the body of the file.

## Defining variables in Cython

In order to use Cython we must first statically type our variables, this means we must tell the compiler what sort
of value to expect. Of all the types, the most commonly used are:

'long' represents a 32 bit precision integer
'long long' represents a 64 bit precision integer
'float' represents a 32 bit precision floating point
'double' represents a 64 bit precision floating point

To define a variable in Cython we can write:
```
%%file mymod.pyx
cdef:
	#Notice the indentation, it tells the program when the cdef ends
	double var1 #Defines a 64 bit floating point called var1
	long int1, int2 #Defines two 32 bit integers, int1 and int2
	double array_test[10][20]  #Defines a stack allocated array of 10 x 20 elements, the var name is 'array_test'
	double var2 = 0.0, var3 = 0.0 #Defines and assigns a value to var2 and var3. 
				      #Initializing a variable by assigning a value is strongly recommended

```

continued

## Defining loops in Cython

In Cython loops are defined in exactly the same way as python, the only difference is that any integer used for range
need to be statically defined:

```
cdef:
	int i, n = 10;

for i in range(n): #This is the same as a python 'for i in range(10)' loop
	do stuff
```

## Defining functions in Cython

Functions in Cython are defined in almost the same way as in python, the difference is that you need to specify the
input type (if any) and you use cpdef instead of def. Example:

```
#This function takes an integer 'a' and returns a * a which is also an integer
cpdef some_func(int a):
	return a * a

#This function takes no arguments and returns a double
cpdef some_func2():
	double a = 3.3
	return a 
```

**Note:** You can also use def and cdef when defining functions in Cython which have different properties and 
advantages/disadvantages, but for now let us stick with cpdef

## Compiling the code
Before we can actually use our Cython code we need to compile it. This can be done by creating a Python script
setup.py with the module name and running it in Python. You can copy and paste this code and just change the relevant
names: 
```
%%file setup.py
from distutils.core import setup
from Cython.Build import cythonize

setup (
	name = "mymodulename",
	ext_files = cythonize("mymodulename.pyx")
)
```
The previous code tells which file to compile. Now to actually compile we use:
```
%run setup.py build_ext --inplace
```
If successful you can the simply import your code using:
```
import mymodulename
```

## Full Examples
In these examples we create a .pyx file, compile it and then compare it with an equivalent python loop using %timeit

Example 1: Loop comparison
```
%%file busy_loop.pyx
cpdef busy_loop():
    cdef:
        long i, total = 0, n = 10000
    for i in range(n):
        total += i
    return total
```
Press shift+enter to execute and create another cell
```
%%file setup.py
from distutils.core import setup
from Cython.Build import cythonize

setup (
    name = "busy_loop",
    ext_modules = cythonize("busy_loop.pyx")
)
```
Press shift+enter to execute and create another cell
```
%run setup.py build_ext --inplace
```
Press shift+enter to execute and create another cell
```
def py_busy_loop(): #This is the same as busy_loop(), except it's defined in Python
    n = 100000
    total = 0
    for i in range(n):
        total += i
    return total
```
Press shift+enter to execute and create another cell
```
import busy_loop #Import your compiled c extension before using it!
%timeit -n 1000 busy_loop.busy_loop()
%timeit -n 1000 py_busy_loop()
```
Press shift enter to run this cell and wait
In my computer Cython busy_loop() took 333 nanoseconds to run while py_busy_loop() took 8.45 miliseconds, in other
words the cython implementation was **25000 times** faster




