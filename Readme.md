# Algorithms for numerical computation.

* Author: LuXu(Oliver Lew)
* email: oliver_lew@outlook.com

just an exercise implementing the algorithms from 《数值计算方法》，林成森，科学出版社。

Requirements:

* make
* gcc(>= c99)
* sed(for makefile dependencies)
* python2

## Project structure

    Myrecipes
    |
    |-- (created) debug: intermediate files
    |
    |-- include: header files
    |
    |-- (created) lib: static library files
    |
    |-- test: test files
    |
    |-- other(start with capital letter): source files

## How test files are compiled

Test folder structure:

    test
    |-- Test.c
    |-- Test.h
    |-- GenerateTest.py
    |-- test*.c

- All the test functions are written in **test\*.c** files. Every test function 
must return `int` type, have no argumemts and be named as 'test*' in order to be
collected automatically.

- **GenerateTest.py** is a python script that scans through **test\*.c** files, 
looks for functions that meet the conditions mentioned above and generate a 
header file **Test.h** which contains declarations, (macro)list of function 
variables, (macro)list of function names and (macro)the number of all the test 
functions.

- **Test.c** is a short C file that includes **Test.h**, uses the macros in it 
to call every test function and print the summary infomation of the test.

To show the flow in a chart:

               GenerateTest.py
    test*.c -----------------> Test.h \
                                      |---> test binary
                               Test.c /
