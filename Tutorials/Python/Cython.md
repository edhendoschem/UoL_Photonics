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