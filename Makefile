# custom names, something need to be customized
## lib name, and the lib will be named lib$(LIBNAME).a
LIBNAME=NR
## binary file name
BINNAME=MyRecipes
## source dirs, dirs directly contain the source files
MY_SRC_DIR=ODE Basic Integral Interpolation LeastSq LibFunction LinearEquations Solve

#those things are not necessarily needed to change
#
# names of directories
SRC_DIR=$(MY_SRC_DIR)
INC_DIR=include
BIN_DIR=bin
LIB_DIR=lib
DBG_DIR=debug
OBJ_DIR=$(DBG_DIR)/obj
DEP_DIR=$(DBG_DIR)/dep

# names of files
SRC=$(foreach dir,$(SRC_DIR),$(wildcard $(dir)/*.c))
OBJ=$(addprefix $(OBJ_DIR)/,$(SRC:.c=.o))
DEP=$(addprefix $(DEP_DIR)/,$(SRC:.c=.d))
BIN=$(BIN_DIR)/$(BINNAME)
LIB=$(LIB_DIR)/lib$(LIBNAME).a

# names of test-related files
TEST_DIR=test
TESTH=$(TEST_DIR)/Test.h
TESTSRC=$(wildcard $(TEST_DIR)/*.c)
TESTOBJ=$(addprefix $(OBJ_DIR)/,$(TESTSRC:.c=.o))
TESTDEP=$(addprefix $(DEP_DIR)/,$(TESTSRC:.c=.d))
GENTEST=$(TEST_DIR)/GenerateTest.py
TESTIFLAGS=-I $(TEST_DIR)

# compiler and parameters
CC=gcc
CFLAGS=-Wall -g -std=c99
IFLAGS=-I $(INC_DIR)
DFLAGS=-MM
LFLAGS=-lm -lNR -L $(LIB_DIR)
PYTHON=python

# directories for new created files
DIRS+=$(BIN_DIR) $(LIB_DIR) $(DBG_DIR)
DIRS+=$(foreach dir,$(SRC_DIR) $(TEST_DIR),$(OBJ_DIR)/$(dir))
DIRS+=$(foreach dir,$(SRC_DIR) $(TEST_DIR),$(DEP_DIR)/$(dir))

.PHONY:all help clean remove cleanall rebuild test

all:$(DEP) $(TESTDEP) $(BIN)
 
include $(foreach dir, $(DIRS), $(wildcard $(dir)/*.d))

$(DEP) $(TESTDEP) $(BIN): | $(DIRS)

$(DIRS):
	mkdir -p $(DIRS)

$(TESTH):$(TESTSRC)
	$(PYTHON) $(GENTEST)

$(DEP):$(DEP_DIR)/%.d:%.c
	set -e; rm -f $@; $(CC) $(DFLAGS) $(IFLAGS) $< | \
	sed 's,\($(*F)\)\.o[ :]*,$(OBJ_DIR)/\1.o $@ : ,g' > $@

$(TESTDEP):$(DEP_DIR)/%.d:%.c
	set -e; rm -f $@; $(CC) $(DFLAGS) $(IFLAGS) $< | \
	sed 's,\($(*F)\)\.o[ :]*,$(OBJ_DIR)/$(*D)/\1.o $@ : ,g' > $@

$(OBJ_DIR)/%.o:%.c
	$(CC) $(CFLAGS) -c $< -o $@ $(IFLAGS)

$(LIB):$(OBJ)
	ar crv $@ $?

$(BIN):$(LIB) $(TESTOBJ)
	$(CC) $(CFLAGS) -o $(BIN) $(TESTOBJ) $(LFLAGS)
	@echo \
"**************************************************\n\
Thank you for using MyRecipes!\n\
You can see other 'make' usage by \"make help\"\n\
**************************************************"


# PHONY targets
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
	rm -rf $(BIN_DIR) $(LIB_DIR)

cleanall:clean remove


rebuild:cleanall
	make

test:
	@echo $(TESTSRC)
	@echo $(foreach dir, $(DIRS), $(wildcard $(dir)/*.d))
