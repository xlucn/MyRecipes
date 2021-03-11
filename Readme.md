# Algorithms for numerical computation.

Just an exercise implementing the algorithms from 《数值计算方法》，林成森，科学出版社, do not guarantee the programs run perfectly.

## Requirements:

- make
- gcc(>= c99)
- sed(for makefile dependencies)
- python

## Project structure

```
Myrecipes
|-- include: header files
|-- test: test files
    |-- Test.c
    |-- Test.h
    |-- GenerateTest.py
    |-- test*.c
|-- [A-Z]*: source files
|-- (created) debug: intermediate files
|-- (created) lib: static library files
```

## How test files are compiled

- Test files follows the pattern in order to automatically parse and generate a test header file:
  - All the test functions are written in **test&ast;.c** files.
  - Every test function must return `int` type, have no argumemts and be named as 'test&ast;'.

- **GenerateTest.py** is a python script that scans through **test&ast;.c** files, 
  looks for functions that meet the conditions mentioned above and generate a 
  header file **Test.h** which contains declarations, (macro)list of function 
  variables, (macro)list of function names and (macro)the number of all the test 
  functions.

- **Test.c** is a short C file that includes **Test.h**, uses the macros in it 
  to call every test function and print the summary infomation of the test.

To show the flow in a chart:

```
         GenerateTest.py
test*.c -----------------> Test.h \
                                  |---> test binary
                           Test.c /
```
