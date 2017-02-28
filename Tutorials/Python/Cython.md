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
1. Go to start -> run and type cmd
2. Go to the directory where you keep your Python scripts using `cd yourdirectoryname`
3. type `jupyter notebook` and press enter

In linux:
1. Open your terminal in the directory you keep your Python scripts
2. type `jupyter notebook` and press enter

### jupyter notebook and cell magic

The reason why we will use jupyter notebook is due to the convenience of the commands you can execute from the 
notebook called 'cell magic'. Here are some examples:

1. Create a file from jupyter notebook:
```
%%file myscript.py
program body
```
The previous command creates a file in the current directory with the .py extension containing whatever is written in
that particular cell

2. Run a Python script from jupyter notebook:
```
%run myscript.py
```
Runs the script using Python as if you had gone to command console or cmd and typed 'python myscript.py' and 
pressed enter

3. Run a system command from jupyter notebook:
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

## Defining a variable in Cython

