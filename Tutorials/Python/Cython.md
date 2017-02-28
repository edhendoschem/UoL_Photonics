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

## Cython files: .pyd and .pyx

Cython files have a different extension than Python files. They are **.pyd** to store declarations (we will ignore
these for now) and **.pyx** to store the definitions. If we want to compile our program to C we need to place the
parts to optimize inside a .pyx file, this can be done very easily with the cell magic '%%file ourextname.pyx' and
then just write the body of the file.

## Defining variables in Cython

In order to use Cython we must first statically type our variables, this means we must tell the compiler what sort
of value to expect. Of all the types, the most commonly used are:

'long' represents a 32 bit precision integer
'long long' represents a 64 bit precision int
'float' represents a 32 bit precision integer
'double' represents a 64 bit precision integer

To define a variable in Cython we can write:
```
%%file mymod.pyx
cdef:
	#Notice the indentation, it tells the program when the cdef ends
	double var1 #Defines a 64 bit floating point called var1
	long int1, int2 #Defines two 32 bit integers, in1 and int2
	double array_test[10][20]  #Defines a stack allocated array of 10 x 20 elements, the var name is 'array_test'
	double var2 = 0.0, var3 = 0.0 #Defines and assigns a value to var2 and var3. 
				      #Initializing a variable by assigning a value is strongly recommended

```





