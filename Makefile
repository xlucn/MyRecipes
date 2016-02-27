# function to remove files in $(1) list
RM=@for i in $(1);do \
	if test -e $$i;then rm $$i; echo "[removing]"$$i; fi;\
done
# function to delete the outdated dependance files.
RMDEP=@for i in $(DEP_DIR)/*.d	;do \
	if ! test -e $(SRC_DIR)/`basename $$i .d`.c;then rm $$i;fi;\
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
TEST=$(SRC_DIR)/Test.h
GENTEST=$(SRC_DIR)/GenerateTest.py

# compiler and parameters
CC=/usr/bin/gcc
CFLAGS=-Wall -g -std=c11
IFLAGS=-I $(INC_DIR)
DFLAGS=-MM
LFLAGS=-lm
PYTHON=/usr/bin/python

.PHONY:all clean remove help backup

all:$(TEST) $(DEP) $(BIN)

$(BIN):$(OBJ)
	$(CC) $(CFLAGS) -o $(BIN) $(OBJ) $(LFLAGS)

$(OBJ):$(OBJ_DIR)/%.o:$(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@ $(IFLAGS)

$(DEP):$(DEP_DIR)/%.d:$(SRC_DIR)/%.c
	$(call RM,$@)
	$(CC) $(DFLAGS) $(IFLAGS) $< | sed 's,\($*\)\.o[ :]*,$(OBJ_DIR)/\1.o $@ : ,g' > $@

$(TEST):$(GENTEST) $(SRC_DIR)/Test.c
	$(PYTHON) $<

-include $(DEP)

help:
	@echo -e \
"Usage:\n\
(all)	:	build the whole project.\n\
count   :   \n\
clean	:	remove the object files.\n\
remove	:	remove the binary file.\n\
backup	:	tar the source file into a tar file named src.tar.gz"

remove:
	$(call RM,$(wildcard $(BIN_DIR)/*))

clean:
	$(call RM,$(wildcard $(OBJ_DIR)/*))

count:
	python ./debug/count.py
