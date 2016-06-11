# custom names, something need to be customized by user
## lib name, and the lib will be named lib$(LIBNAME).a
LIBNAME=NR
## binary file name
BINNAME=MyRecipes
## source dirs
MY_SRC_DIR=ODE Basic Integral Interpolation LeastSq LibFunction LinearEquations Solve

#those things are not necessarily needed to change
#
# names of directories
SRC_DIR=$(MY_SRC_DIR)
BIN_DIR=bin
INC_DIR=include
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
GENTEST=GenerateTest.py
TESTIFLAGS=-I $(TEST_DIR)

# compiler and parameters
CC=gcc
CFLAGS=-Wall -g -std=c11
IFLAGS=-I $(INC_DIR)
DFLAGS=-MM
LFLAGS=-lm
PYTHON=python


.PHONY:all clean remove help dirs count files test

all:dirs files

DIRS+=$(BIN_DIR)
DIRS+=$(LIB_DIR)
DIRS+=$(foreach dir,$(SRC_DIR) $(TEST_DIR),$(OBJ_DIR)/$(dir))
DIRS+=$(foreach dir,$(SRC_DIR) $(TEST_DIR),$(DEP_DIR)/$(dir))

dirs:
	@mkdir -p $(DIRS)

files:$(DEP) $(TESTDEP) $(BIN)

$(BIN):$(LIB) $(TESTOBJ) $(OBJ)
	$(CC) $(CFLAGS) -o $(BIN) $(TESTOBJ) $(OBJ) $(LFLAGS)
	@echo \
"**************************************************\
\nThank you for using MyRecipes!\n\
You can see other 'make' usage by \"make help\"\n\
**************************************************"

$(LIB):$(OBJ)
	ar crv $@ $?

$(OBJ):$(OBJ_DIR)/%.o:%.c
	$(CC) $(CFLAGS) -c $< -o $@ $(IFLAGS)

$(DEP):$(DEP_DIR)/%.d:%.c
	rm -f $@; $(CC) $(DFLAGS) $(IFLAGS) $< | sed 's,\($*\)\.o[ :]*,$(OBJ_DIR)/\1.o $@ : ,g' > $@

$(TESTOBJ):$(OBJ_DIR)/%.o:%.c $(TESTH)
	$(CC) $(CFLAGS) -c $< -o $@ $(IFLAGS)

$(TESTDEP):$(DEP_DIR)/%.d:%.c
	rm -f $@; $(CC) $(DFLAGS) $(IFLAGS) $< | sed 's,\($*\)\.o[ :]*,$(OBJ_DIR)/\1.o $@ : ,g' > $@

$(TESTH):$(TESTSRC)
	cd $(TEST_DIR);$(PYTHON) $(GENTEST)

-include $(wildcard $(DEP_DIR)/*.d)

help:
	@echo \
"Usage:\n\
(all)	:	build the whole project.\n\
count	:	count the lines, words and bytes of source files.\n\
clean	:	remove the object files and the dependancy files.\n\
remove	:	remove the binary file.\n\
cleanall:	remove all the files created by make.\n\
rebuild	:	clean all the files and rebuild the whole project."

count:
	echo -n `date`"\t" >> ./.count
	cat $(SRC_DIR)/* $(INC_DIR)/* $(TEST_DIR)/* | wc >> ./.count

remove:
	rm -rf $(BIN_DIR) $(LIB_DIR)

clean:
	rm -rf $(OBJ_DIR) $(DEP_DIR)

cleanall:clean remove


rebuild:cleanall all


test:
	@echo $(SRC)
	@echo $(OBJ)
	@echo $(DEP)
