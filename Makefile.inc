###############################################################################
# Compile setting with intel compiler and intel MPI.
# Modify flags according to the system and environment
###############################################################################

FC = /home/mack/project/kisti/MPI/mvapich3/bin/mpif90

FLAG  = -O2 -r8 -i4 -shared-intel
#FLAG	= -O3 -r8 -i4 -g -traceback -shared-intel

