#! /usr/env/bin python
# Generate the header files including all the declarations of test functions,
# a array of pointers pointing to these functions and a array of function names.

testc = './src/Test.c'
testh = './include/Test.h'
testinfo = './src/Testinfo.c'

def readfuncs(testfile):
    '''
    read source file
    '''
    funcs = []
    source = open(testfile, 'r')
    for line in source.readlines():
        words = line.split()
        if len(words) > 1 and words[0] == 'int':
            if len(words[1]) > 4 and words[1][:4] == 'test':
                funcs.append(words[1].split('(')[0])
    funcs = list(set(funcs))
    funcs.remove('testall')
    source.close()
    return funcs

def gentest(funcs):
    '''
    generate test header file
    '''
    header = open(testh, 'w')
    header.write(
'''/*
** This file is automatically generated by GenerateTest.py.
*/

#ifndef _TESTDECLARATION_H_
#define _TESTDECLARATION_H_

#define PASSED 0
#define FAILED 1

extern int (*tests[])();
extern char *names[];

'''
    )
    for func in funcs:
        header.write('int ' + func + '();\n')
    header.write('int testall();\n')
    header.write('\n#endif\n')
    header.close()

def gentestinfo(funcs):
    '''
    generate testinfo file
    '''
    cfile = open(testinfo, 'w')
    cfile.write(
'''/*
** This file is automatically generated by GenerateTest.py.
*/
#ifndef _TESTINFO_H_
#define _TESTINFO_H_

#include <Test.h>

int (*tests[])() = {
'''
    )
    for func in funcs:
        cfile.write('\t' + func)
        if func != funcs[-1]:
            cfile.write(',')
        cfile.write('\n')
    cfile.write('};\n')
    cfile.write('\nchar *names[] = {\n')
    for func in funcs:
        cfile.write('\t\"' + func + '\"')
        if func != funcs[-1]:
            cfile.write(',')
        cfile.write('\n')
    cfile.write('};\n')
    cfile.write('\n#endif\n')
    cfile.close()

def main():
    funclist = readfuncs(testc)
    funclist2 = readfuncs(testh)
    if cmp(funclist, funclist2):
        gentest(funclist)
        gentestinfo(funclist)
        print funclist, funclist2
        print 'File ' + testh + ' and ' + testinfo + ' regenerated.'
    
if __name__ == '__main__':
    main()
