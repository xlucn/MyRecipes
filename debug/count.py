#!/usr/bin/env python
# count the words of files in the current project and write to a file
import datetime, os, sys

num = os.popen('cat ./include/* ./src/* 2> /dev/null | wc').read().split()

file = open(sys.path[0] + '/count','a')
file.write(datetime.datetime.today().ctime())
file.write('\t%8s lines %8s words %8s bytes\n' % (num[0], num[1], num[2]))
file.close()
