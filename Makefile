MOSFHET_DIR = ./src/mosfhet
include $(MOSFHET_DIR)/Makefile.def
INCLUDE_DIRS += ./include/
KEY=BINARY
PARAM=SET_2_3

SRC = sparse_amortized_bootstrap.c
SRC_SELF = $(addprefix ./src/, $(SRC))


all: main

main: $(SRC_MOSFHET) $(SRC_SELF) main.c
	$(CC) -g -o main $^ $(OPT_FLAGS) $(LIBS) -D$(KEY) -D$(PARAM)
