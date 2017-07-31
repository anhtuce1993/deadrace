TOP_DIR:=/home/anhtu/deadrace
SW_DIR:=$(TOP_DIR)/stackwalk/ver0 
SW_INC:=-I$(SW_DIR)
SW_LIB:=-L$(SW_DIR) -lstackwalk
CPPFLAGS:=-I$(TOP_DIR)/inc -g
TRACES_DIR:=$(TOP_DIR)/traces

NPROCS:= 3


VPATH=$(TOP_DIR)/src

C_SRCS:= Loop.c Iter.c Values.c Misc.c PMPI.c Controller.c Memory.c
SUMMARY_SRC:= src/summary.c

OBJS:=$(C_SRCS:.c=.o)

SRC_TEST:=tests/race.c

%.o: %.c
	mpicxx -Wall $(CPPFLAGS) -c $< 
# Loop.o: Loop.c
# 	mpicxx -Wall $(CPPFLAGS) -c $< 
# Iter.o: Iter.c
# 	mpicxx -Wall $(CPPFLAGS) -c $< 
# Values.o: Values.c
# 	mpicxx -Wall $(CPPFLAGS) -c $< 
# Misc.o: Misc.c
# 	mpicxx -Wall $(CPPFLAGS) -c $< 
# PMPI.o: PMPI.c
# 	mpicxx -Wall $(CPPFLAGS) -c $< 

libdeadrace.a: $(OBJS)
	ar -cvq $(TOP_DIR)/$@ $(OBJS)
	#rm -f $(OBJS)

test:  $(SRC_TEST) libdeadrace.a
	mpicxx $(CPPFLAGS) -o test $(SRC_TEST) libdeadrace.a
	#rm -f libdeadrace.a 

run: test
	mpirun -np $(NPROCS) ./test

summary: 
	gcc -o $(TRACES_DIR)/summary $(SUMMARY_SRC)
	cd traces && ./summary $(NPROCS)

clean:
	rm -rf *.o libdeadrace.a test traces/* result

cleantraces:
	rm -rf traces/* 
