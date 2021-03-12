#--------------------------Names of files and folders---------------------------

# Names of existing directories
SRC_DIR=$(shell find . -maxdepth 1 -type d -name "[A-Z]*")
INC_DIR=include

# Names of new dirs
LIB_DIR=lib
DBG_DIR=debug

# Names of source files
SRC=$(foreach dir,$(SRC_DIR),$(wildcard $(dir)/*.c))

# New files
OBJ=$(addprefix $(DBG_DIR)/,$(SRC:.c=.o))
DEP=$(addprefix $(DBG_DIR)/,$(SRC:.c=.d))
LIB=$(LIB_DIR)/libNR.a

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
SHELL=/bin/sh
# GCC flags
CFLAGS=-Wall -g -std=c99 -fpic
IFLAGS=-I $(INC_DIR)
DFLAGS=-MM
LFLAGS=-lm -lNR -L $(LIB_DIR)
# Shell flags
.SHELLFLAGS = -c -e


#---------------------------Targets---------------------------------------------

.PHONY: all lib help clean remove cleanall rebuild doc test

all: $(LIB) $(DBG) $(TESTDEP) $(TESTBIN)

doc:
	doxygen Doxyfile

tags:
	ctags --exclude=.ccls-cache/* --exclude=docs/* --exclude=compile_commands.json -R .

help:
	@echo "Usage:\n\
	(all)	:	build the whole project.\n\
	clean	:	remove the object files and the dependancy files.\n\
	remove	:	remove all generated files.\n\
	rebuild	:	remove and make.\n\
	help    :	show this message"

clean:
	rm -rf $(DBG_DIR)

remove: clean
	rm -rf $(LIB_DIR) $(TESTBIN)

rebuild: remove
	make

test: $(TESTDEP) $(TESTBIN)
	$(TESTBIN)

_test:
	@echo $(SRC_DIR)
	@echo $(SRC)


#------------------------------Dependencies-------------------------------------

-include $(foreach dir, $(DIRS), $(wildcard $(dir)/*.d))

$(DIRS):
	mkdir -p $(DIRS)

$(TESTH): $(TESTSRC) $(GENTEST)
	$(PYTHON) $(GENTEST)

$(DBG_DIR)/%.d: %.c | $(DIRS)
	$(CC) $(DFLAGS) $(IFLAGS) $< | sed 's,$(*F).o[ :]*,$(DBG_DIR)/$*.o $@ : ,g' > $@

$(DBG_DIR)/%.o: %.c | $(DIRS)
	$(CC) $(CFLAGS) -c $< -o $@ $(IFLAGS)

$(LIB): $(OBJ) | $(DIRS)
	ar cr $@ $?

$(TESTBIN): $(LIB) $(TESTOBJ) | $(DIRS)
	$(CC) $(CFLAGS) -o $(TESTBIN) $(TESTOBJ) $(LFLAGS)
