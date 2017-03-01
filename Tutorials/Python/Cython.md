# Cython: A quick tutorial

## Table of contents
1. [What is Cython](#cython)
2. [When to use Cython](#when)
3. [How do I install Cython?](#how)
4. [Before starting](#before)
5. [Cython files: .pyd and .pyx](#files)
6. [Defining variables in Cython](#vars)
7. [Defining loops in Cython](#loops)
8. [Defining functions in Cython](#functions)
9. [Full Examples](#examples)
10. [Working with line_profiler to find bottlenecks](#line_profiler)




## What is Cython? <a name="cython"></a>

Cython is a way to compile C extensions using syntax similar to Python, this allows you to have C levels of speed
in your Python code. 

More info on their website:
http://cython.org/

## When to use Cython? <a name="when"></a>
When your code is running slowly the first step is to check which part can be improved, whether you can avoid loops
and change your arrays to numpy arrays. In order to help you determine where the code might be slow you can install
line_profiler with `pip install line_profiler`

If you have a problem that requires extensive use of loops, or computations using the same type of elements (i.e.
numbers only).

## How do I install Cython?<a name="how"></a>

If you have installed Python Anaconda then you already have Cython installed, otherwise type:
```
pip install cython
```

## Before starting<a name="before"></a>

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

## Cython files: .pyd and .pyx <a name="files"></a>

Cython files have a different extension than Python files. They are **.pyd** to store declarations (we will ignore
these for now) and **.pyx** to store the definitions. If we want to compile our program to C we need to place the
parts to optimize inside a .pyx file, this can be done very easily with the cell magic '%%file ourextname.pyx' and
then just write the body of the file.

## Defining variables in Cython <a name="vars"></a>

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
```

***Warning:*** Python by default has unlimited precision, this is not the case in statically typed Cython, if for
example we have defined `long i` The maximum possible integer it can store is -2,147,483,647 to 2,147,483,647. If we
exceed that number it will "wrap around": `2,147,483,647 + 1 = -2,147,483,647` or `-2,147,483,647 - 1 = 2,147,483,647`

If you are unsure whether an integer number will fit inside a long long integer (+/-9,223,372,036,854,775,807) then
it is better to store it in a double floating point 

## Defining loops in Cython <a name="loops"></a>

In Cython loops are defined in exactly the same way as python, the only difference is that any integer used for range
need to be statically defined:

```
cdef:
	int i, n = 10;

for i in range(n): #This is the same as a python 'for i in range(10)' loop
	do stuff
```

## Defining functions in Cython <a name="functions"></a>

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

## Compiling the code <a name="compiling"></a>
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

## Full Examples <a name="examples"></a>
In these examples we create a .pyx file, compile it and then compare it with an equivalent python loop using %timeit

### Example 1: Loop comparison
```
%%file busy_loop.pyx
cpdef busy_loop():
    cdef:
        long i, n = 10000
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
Press ctrl+enter to run this cell and wait. This will run each statement a thousand times on three ocassions and
display the fastest runtime. You should get some output like this:
```
1000 loops, best of 3: 333 ns per loop
1000 loops, best of 3: 8.45 ms per loop
```
In my computer Cython busy_loop() took 333 nanoseconds to run while py_busy_loop() took 8.45 miliseconds, in other
words the Cython implementation was **~25000 times** faster than pure Python.

### Example 2: Fibonacci sequence
In this example the function will take an argument n which is the element number: 
f1 = 1, f2 = 1, f3 = 2, f4 = 3, f5 = 5, f6 = 8, f7 = 13 ... fn = f_(n-1) + f_(n-2)
```
%%file fibonacci.pyx
cpdef fibonacci(long long n):
    cdef:
        long long current = 1, previous = 1, temp = 0, i
    if n > 2:
        for i in range(2, n):
            temp = current
            current += previous
            previous = temp
    else:
        return current
    return current
```
Run and then
```
%%file setup.py
from distutils.core import setup
from Cython.Build import cythonize

setup (
    name = "fibonacci",
    ext_modules = cythonize("fibonacci.pyx")
)
```
Run and then
```
%run setup.py build_ext --inplace
```
In another cell:
```
def py_fibonacci(n):
    current = 1
    previous = 1
    temp = 0
    if n > 2:
        for i in range(2, n):
            temp = current
            current += previous
            previous = temp
    else:
        return current
    return current
```
Finally:
```
import fibonacci as fib
%timeit -n 10000 fib.fibonacci(92)
%timeit -n 10000 py_fibonacci(92)
```
The result:
```
10000 loops, best of 3: 497 ns per loop
10000 loops, best of 3: 7.79 µs per loop
```
In my computer the Cython version ran ***~15 times*** faster than the pure Python version. There is a caveat however;
If we print the values:
```
print(fib.fibonacci(92))
print("_________________________________________________")
print(py_fibonacci(92))
```
The output:
```
7540113804746346429
_________________________________________________
7540113804746346429
```
We notice they are the same as they should, but if we try with n = 93, we exceed the limit of what the long long
integer can store and obtain a wrap around result:
```
print(fib.fibonacci(93))
print("_________________________________________________")
print(py_fibonacci(93))
```
We get:
```
-6246583658587674878
_________________________________________________
12200160415121876738
```

If we are unsure whether an integer value will fit inside an integer variable, it is better to change it to a double.

## Working with line_profiler to find bottlenecks <a name="line_profiler"></a>
pending
