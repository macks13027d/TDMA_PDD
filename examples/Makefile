include ../Makefile.inc

VERSION	= pdd_cpu
EXE		= tdma_${VERSION}

INC	   		= -I ../include/hypre-2.11.1\
	   			  -I ../include/generated_by_sourceCode \
	   			  -I ../include
LIB   	      = -L ../lib/hypre-2.11.1 -lHYPRE \
	   	        -L ../lib/generated_by_sourceCode -ltdma_pdd_cpu \
					  -L /home/mack/project/kisti/MPI/mvapich3/lib

OBJS = global_host.o \
		 global_initialize.o \
		 main.o          
		 




all: $(EXE) 
 
${EXE}: ${OBJS}
	${FC} -o ${EXE} ${OBJS} ${LFLAGS} -lstdc++ ${INC} ${LIB} 
	mv ${EXE} ../run

%.o: %.f90
	${FC} ${FLAG} -c $< -o $@ ${INC}


clean:
	rm *.mod *.o

