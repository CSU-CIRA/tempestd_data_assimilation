VERSION=2.2.6
###########################################
#                                         #
#                                         #
#         Version 2.2.6 Section           #
#               HERA                      #
#                                         #
###########################################
ifeq ($(VERSION),2.2.6)
CRTM_DIR =/scratch2/NCEPDEV/nwprod/NCEPLIBS/compilers/intel/18.0.5.274
FLAGS = -O3 -free -assume byterecl -fp-model source -traceback -xHOST
INC = -I$(CRTM_DIR)/include/crtm_v$(VERSION)
LIB = -L$(CRTM_DIR)/lib -lcrtm_v$(VERSION)
endif
###########################################
#                                         #
#                                         #
#          Additional setup               #
#                                         #
#                                         #
###########################################
FC = ifort
FL = $(FC)

#
# Target definitions
# ------------------
# The main target

all: Get_CRTM

Get_CRTM: Get_CRTM.o
	${FC} -o Get_CRTM ${FLAGS} ${INC} Get_CRTM.o ${LIB}

.SUFFIXES : .f90  .o

.f90.o :
	${FC} ${FLAGS} ${INC} -c $<


# A clean target
clean:
	rm -f *.o *.exe



