#------------- custom names, something need to be customized -------------------

## header dir, where the header files are
MY_INC_DIR=include

## dir containing the source dirs
MY_SRC_DIR=.

## (create) lib name, and the lib will be named lib$(LIBNAME).a
LIBNAME=NR

## (create) library file dir, library file will be created here
MY_LIB_DIR=lib

## (create) debug files dir, dependency files and object files will be here
MY_DBG_DIR=debug


#--------------- those things are not necessarily needed to change -------------
#
# names of directories
SRC_DIR=$(shell find $(MY_SRC_DIR) -maxdepth 1 -type d -name "[A-Z]*")
INC_DIR=$(MY_INC_DIR)
# new dirs below
LIB_DIR=$(MY_LIB_DIR)
DBG_DIR=$(MY_DBG_DIR)
OBJ_DIR=$(DBG_DIR)/obj
DEP_DIR=$(DBG_DIR)/dep

# names of files
SRC=$(foreach dir,$(SRC_DIR),$(wildcard $(dir)/*.c))
# new files below
OBJ=$(addprefix $(OBJ_DIR)/,$(SRC:.c=.o))
DEP=$(addprefix $(DEP_DIR)/,$(SRC:.c=.d))
LIB=$(LIB_DIR)/lib$(LIBNAME).a

# names of test-related files
TEST_DIR=test
TESTSRC=$(wildcard $(TEST_DIR)/*.c)
GENTEST=$(TEST_DIR)/GenerateTest.py
# new files below
TESTH=$(TEST_DIR)/Test.h
TESTBIN=$(TEST_DIR)/test
TESTOBJ=$(addprefix $(OBJ_DIR)/,$(TESTSRC:.c=.o))
TESTDEP=$(addprefix $(DEP_DIR)/,$(TESTSRC:.c=.d))

# shell config
SHELL=/bin/bash
.SHELLFLAGS = -c -e

# executables and parameters
CC=gcc
PYTHON=python2
CFLAGS=-Wall -g -std=c99
IFLAGS=-I $(INC_DIR)
DFLAGS=-MM
LFLAGS=-lm -lNR -L $(LIB_DIR)

# directories for all new created files
DIRS+=$(LIB_DIR) $(DBG_DIR) $(OBJ_DIR) $(DEP_DIR)
DIRS+=$(foreach dir,$(SRC_DIR) $(TEST_DIR),$(OBJ_DIR)/$(dir))
DIRS+=$(foreach dir,$(SRC_DIR) $(TEST_DIR),$(DEP_DIR)/$(dir))

.PHONY:all help clean remove cleanall rebuild doc test


all:$(DEP) $(TESTDEP) $(TESTBIN)


doc:
	doxygen Doxyfile

# other PHONY targets

help:
	@echo \
"Usage:\n\
(all)	:	build the whole project.\n\
clean	:	remove the object files and the dependancy files.\n\
remove	:	remove the binary files and the library files.\n\
cleanall:	clean and remove.\n\
rebuild	:	cleanall and make.\n\
help    :	show this message"

clean:
	rm -rf $(DBG_DIR)

remove:
	rm -rf $(LIB_DIR) $(TESTBIN)

cleanall:clean remove


rebuild:cleanall
	make

test:
	@echo $(wildcard ^[A-Z]*)

-include $(foreach dir, $(DIRS), $(wildcard $(dir)/*.d))

$(DEP) $(TESTDEP) $(TESTBIN): | $(DIRS)

$(DIRS):
	mkdir -p $(DIRS)

$(TESTH):$(TESTSRC) $(GENTEST)
	$(PYTHON) $(GENTEST)

$(DEP_DIR)/%.d:%.c
	$(CC) $(DFLAGS) $(IFLAGS) $< | sed 's,$(*F).o[ :]*,$(OBJ_DIR)/$*.o $@ : ,g' > $@

$(OBJ_DIR)/%.o:%.c
	$(CC) $(CFLAGS) -c $< -o $@ $(IFLAGS)

$(LIB):$(OBJ)
	ar cr $@ $?

$(TESTBIN):$(LIB) $(TESTOBJ)
	$(CC) $(CFLAGS) -o $(TESTBIN) $(TESTOBJ) $(LFLAGS)
