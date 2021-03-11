#----------------- custom names, which can be customized -----------------------

## (existing) header dir, where the header files are
MY_INC_DIR=include

## (existing) dir containing the source dirs
MY_SRC_DIR=.

## (created) lib name, and the lib will be named lib$(LIBNAME).a
LIBNAME=NR

## (created) library file dir, library file will be created here
MY_LIB_DIR=lib

## (created) debug files dir, dependency files and object files will be here
MY_DBG_DIR=debug


#--------------------------Names of files and folders---------------------------

# Names of existing directories
SRC_DIR=$(shell find $(MY_SRC_DIR) -maxdepth 1 -type d -name "[A-Z]*")
INC_DIR=$(MY_INC_DIR)

# Names of new dirs
LIB_DIR=$(MY_LIB_DIR)
DBG_DIR=$(MY_DBG_DIR)

# Names of source files
SRC=$(foreach dir,$(SRC_DIR),$(wildcard $(dir)/*.c))

# New files
OBJ=$(addprefix $(DBG_DIR)/,$(SRC:.c=.o))
DEP=$(addprefix $(DBG_DIR)/,$(SRC:.c=.d))
LIB=$(LIB_DIR)/lib$(LIBNAME).a

# Directories for all newly created files
DIRS+=$(LIB_DIR) $(DBG_DIR)


#---------------------------------Tests-----------------------------------------

# Names of test-related files
TEST_DIR=test
TESTSRC=$(wildcard $(TEST_DIR)/*.c)
GENTEST=$(TEST_DIR)/GenerateTest.py

# New files
TESTH=$(TEST_DIR)/Test.h
TESTBIN=$(TEST_DIR)/test
TESTOBJ=$(addprefix $(DBG_DIR)/,$(TESTSRC:.c=.o))
TESTDEP=$(addprefix $(DBG_DIR)/,$(TESTSRC:.c=.d))

# Directories for all newly created files
DIRS+=$(foreach dir,$(SRC_DIR) $(TEST_DIR),$(DBG_DIR)/$(dir))

#---------------------------Commands and flags----------------------------------

# Executables
CC=gcc
PYTHON=python
# Shell config
SHELL=/bin/bash
# GCC flags
CFLAGS=-Wall -g -std=c99 -fpic
IFLAGS=-I $(INC_DIR)
DFLAGS=-MM
LFLAGS=-lm -lNR -L $(LIB_DIR)
# Shell flags
.SHELLFLAGS = -c -e


#---------------------------Targets---------------------------------------------

.PHONY:all help clean remove cleanall rebuild doc test


all:$(DEP) $(TESTDEP) $(TESTBIN)


doc:
	doxygen Doxyfile

tags:
	ctags --exclude=.ccls-cache/* --exclude=docs/* --exclude=compile_commands.json -R .

help:
	@echo \
"Usage:\n\
(all)	:	build the whole project.\n\
clean	:	remove the object files and the dependancy files.\n\
remove	:	remove all generated files.\n\
rebuild	:	remove and make.\n\
help    :	show this message"

clean:
	rm -rf $(DBG_DIR)

remove:clean
	rm -rf $(LIB_DIR) $(TESTBIN)

rebuild:remove
	make

test:
	@echo $(SRC_DIR)
	@echo $(SRC)


#------------------------------Dependencies-------------------------------------

-include $(foreach dir, $(DIRS), $(wildcard $(dir)/*.d))

$(DEP) $(TESTDEP) $(TESTBIN): | $(DIRS)

$(DIRS):
	mkdir -p $(DIRS)

$(TESTH):$(TESTSRC) $(GENTEST)
	$(PYTHON) $(GENTEST)

$(DBG_DIR)/%.d:%.c
	$(CC) $(DFLAGS) $(IFLAGS) $< | sed 's,$(*F).o[ :]*,$(DBG_DIR)/$*.o $@ : ,g' > $@

$(DBG_DIR)/%.o:%.c
	$(CC) $(CFLAGS) -c $< -o $@ $(IFLAGS)

$(LIB):$(OBJ)
	ar cr $@ $?

$(TESTBIN):$(LIB) $(TESTOBJ)
	$(CC) $(CFLAGS) -o $(TESTBIN) $(TESTOBJ) $(LFLAGS)
