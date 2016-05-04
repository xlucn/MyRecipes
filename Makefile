# function to remove files in $(1) list
RM=@for i in $(1);do                                                 \
    if [ -e $$i ];                                                   \
    then                                                             \
        rm -f $$i > /dev/null 2>&1; echo "[removing]"$$i;            \
    fi;                                                              \
done

# function to delete the outdated dependance files.
RMDEP=@                                                              \
for i in `ls $(DEP_DIR)`;do                                          \
    if ! [ -e $(SRC_DIR)/`basename $$i .d`.c ];                      \
    then                                                             \
        rm $(DEP_DIR)/$$i; echo "[removing]"$(DEP_DIR)/$$i;          \
    fi;                                                              \
done

# names of directories
SRC_DIR=src
BIN_DIR=bin
INC_DIR=include
LIB_DIR=lib
DBG_DIR=debug
OBJ_DIR=$(DBG_DIR)/obj
DEP_DIR=$(DBG_DIR)/dep
DIRS=$(BIN_DIR) $(LIB_DIR) $(OBJ_DIR) $(DEP_DIR)

# names of files
SRC=$(wildcard $(SRC_DIR)/*.c)
OBJ=$(addprefix $(OBJ_DIR)/,$(notdir $(SRC:.c=.o)))
DEP=$(addprefix $(DEP_DIR)/,$(notdir $(SRC:.c=.d)))
BIN=$(BIN_DIR)/MyRecipes
LIB=$(LIB_DIR)/libNR.a

# names of test-related files
TEST_DIR=./test
TESTH=$(TEST_DIR)/Test.h
TESTC=$(TEST_DIR)/Test.c
TESTSRC=$(wildcard $(TEST_DIR)/*.c)
TESTOBJ=$(addprefix $(OBJ_DIR)/,$(notdir $(TESTSRC:.c=.o)))
TESTDEP=$(addprefix $(DEP_DIR)/,$(notdir $(TESTSRC:.c=.d)))
GENTEST=GenerateTest.py
TESTIFLAGS=-I $(TEST_DIR)

# compiler and parameters
CC=gcc
CFLAGS=-Wall -g -std=c11
IFLAGS=-I $(INC_DIR)
DFLAGS=-MM
LFLAGS=-lm
PYTHON=python


.PHONY:all clean remove help dirs count files

.DEFAULT:
	@echo "see usage by \"make help\""

all:dirs files

dirs:$(DIRS)
	mkdir -p $(DIRS)

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

$(OBJ):$(OBJ_DIR)/%.o:$(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@ $(IFLAGS)

$(DEP):$(DEP_DIR)/%.d:$(SRC_DIR)/%.c
	rm -f $@; $(CC) $(DFLAGS) $(IFLAGS) $< | sed 's,\($*\)\.o[ :]*,$(OBJ_DIR)/\1.o $@ : ,g' > $@

$(TESTOBJ):$(OBJ_DIR)/%.o:$(TEST_DIR)/%.c $(TESTH)
	$(CC) $(CFLAGS) -c $< -o $@ $(IFLAGS)

$(TESTDEP):$(DEP_DIR)/%.d:$(TEST_DIR)/%.c $(TESTH)
	rm -f $@; $(CC) $(DFLAGS) $(IFLAGS) $< | sed 's,\($*\)\.o[ :]*,$(OBJ_DIR)/\1.o $@ : ,g' > $@

$(TESTH):$(TESTC)
	cd $(TEST_DIR);$(PYTHON) $(GENTEST)

-include $(wildcard $(DEP_DIR)/*.d)

help:
	@echo \
"Usage:\n\
(all)	:	build the whole project.\n\
count	:	count the lines, words and bytes of source files.\n\
clean	:	remove the object files and the dependancy files.\n\
cleandep:	clean the useless dependancy files.\n\
remove	:	remove the binary file.\n\
cleanall:	remove all the files created by make.\n\
rebuild	:	clean all the files and rebuild the whole project."

count:
	echo -n `date`"\t" >> ./.count && cat $(SRC_DIR)/* $(INC_DIR)/* $(TEST_DIR)/* | wc >> ./.count

remove:
	$(call RM,$(BIN))

clean:
	$(call RM,$(OBJ) $(DEP) $(TESTDEP) $(TESTOBJ))

cleanall:clean remove

rebuild:cleanall all

cleandep:
	$(call RMDEP)
