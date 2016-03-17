# function to remove files in $(1) list
RM=@for i in $(1);do \
	if [ -e $$i ];\
	then\
		rm -f $$i > /dev/null 2>&1; echo "[removing]"$$i;\
	fi;\
done
# function to delete the outdated dependance files.
RMDEP=@\
for i in `ls $(DEP_DIR)`;do \
	if ! [ -e $(SRC_DIR)/`basename $$i .d`.c ];\
	then\
		rm $(DEP_DIR)/$$i; echo "[removing]"$(DEP_DIR)/$$i;\
	fi;\
done

# names of files and directories
SRC_DIR=./src
OBJ_DIR=./debug/obj
BIN_DIR=./bin
DEP_DIR=./debug/dep
INC_DIR=./include
SRC=$(wildcard $(SRC_DIR)/*.c)
OBJ=$(addprefix $(OBJ_DIR)/,$(notdir $(SRC:.c=.o)))
DEP=$(addprefix $(DEP_DIR)/,$(notdir $(SRC:.c=.d)))
BIN=$(BIN_DIR)/MyRecipes
TESTH=$(INC_DIR)/Test.h
TESTC=$(SRC_DIR)/Test.c
GENTEST=$(SRC_DIR)/GenerateTest.py

# compiler and parameters
CC=gcc
CFLAGS=-Wall -g -std=c11
IFLAGS=-I $(INC_DIR)
DFLAGS=-MM
LFLAGS=-lm
PYTHON=python

.PHONY:all clean remove help dir count

all:dir $(DEP) $(TESTH) $(BIN)

$(BIN):$(OBJ)
	$(CC) $(CFLAGS) -o $(BIN) $(OBJ) $(LFLAGS); chmod a+x $(BIN)

$(OBJ):$(OBJ_DIR)/%.o:$(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@ $(IFLAGS)

$(DEP):$(DEP_DIR)/%.d:$(SRC_DIR)/%.c
	rm -f $@; $(CC) $(DFLAGS) $(IFLAGS) $< | sed 's,\($*\)\.o[ :]*,$(OBJ_DIR)/\1.o $@ : ,g' > $@

$(TESTH):$(TESTC)
	$(PYTHON) $(GENTEST)

dir:
	@mkdir -p $(BIN_DIR) $(DEP_DIR) $(OBJ_DIR)

-include $(DEP)

help:
	@echo \
"Usage:\n\
(all)	:	build the whole project.\n\
count   :	count the lines, words and bytes of source files.\n\
clean	:	remove the object files and the dependancy files.\n\
cleandep:	clean the useless dependancy files.\n\
remove	:	remove the binary file.\n\
cleanall:	remove all the files created by make.\n\
rebuild :	clean all the files and rebuild the whole project.\n\
backup	:	tar the source file into a tar file named src.tar.gz"

count:
	echo -n `date`"\t" >> ./.count && cat $(SRC_DIR)/* $(INC_DIR)/* | wc >> ./.count

remove:
	$(call RM,$(BIN))

clean:
	$(call RM,$(OBJ) $(DEP))

cleanall:clean remove

rebuild: cleanall all

cleandep:
	$(call RMDEP)

test:
	@echo 'a'
