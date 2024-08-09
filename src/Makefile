include ../Makefile.inc

LIBexp = libtdma_pdd_cpu.a
INC = -I ../include

OBJDIR = ./obj

OBJS = pdd_host.o								\
       tdma_pdd.o              			\
       pdd_initialize.o                      	\
       communication.o        				\
       tdma_main.o	\
		 
all: $(LIBexp) inc

$(LIBexp): $(OBJS)
	ar qc $@ $(OBJS)
	mv $@ ../lib/generated_by_sourceCode

inc: $(LIBexp)
	cp *.mod ../include/generated_by_sourceCode

%.o: %.f90
	${FC} ${FLAG} -c $< -o $@ $(INC)

${EXE}: ${OBJS}
	${FC} -o ${EXE} ${OBJS} -lstdc++
	cp ${EXE} exam_jcp

clean:
	rm -rf ../lib/generated_by_sourceCode/$(LIBexp) 
	rm -rf ../include/generated_by_sourceCode*.mod
	rm -rf *.o *.mod

