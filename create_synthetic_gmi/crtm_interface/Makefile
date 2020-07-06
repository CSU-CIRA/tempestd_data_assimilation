VERSION=2.2.6
#VERSION=2.1.3.ifort
#VERSION=2.2.3
###########################################
#                                         #
#                                         #
#         Version 2.1.3 Section           #
#                                         #
#                                         #
###########################################
#
ifeq ($(VERSION),2.1.3.ifort)
CRTM_DIR = /home/grasso/rtm/crtm_v2.1.3/REL-$(VERSION)
#
FC_FLAGS =\
 -c \
 -O2 \
 -convert little_endian \
 -free \
 -assume byterecl \
 -CB \
 -I$(CRTM_DIR)/include
#
FL_FLAGS = -L$(CRTM_DIR)/lib -lCRTM # use for 2.1.3
endif
###########################################
#                                         #
#                                         #
#         Version 2.2.3 Section           #
#                                         #
#                                         #
###########################################
ifeq ($(VERSION),2.2.3)
CRTM_DIR =/home/grasso/rtm/crtm_v2.2.3/REL-2.2.3/crtm_v2.2.3
#
#*** FC_FLAGS are from rcm02:~/rtm/crtm_v2.2.3/config_setup/ifort.setup.cira
#
FC_FLAGS =\
-c \
-O3 \
-convert little_endian \
-fp-model source \
-free \
-assume byterecl \
-I$(CRTM_DIR)/include
#
FL_FLAGS = -L$(CRTM_DIR)/lib -lcrtm
endif
###########################################
#                                         #
#                                         #
#         Version 2.2.6 Section           #
#               HERA                      #
#                                         #
###########################################
ifeq ($(VERSION),2.2.6)
CRTM_DIR =/scratch2/NCEPDEV/nwprod/NCEPLIBS/compilers/intel/18.0.5.274
#
#*** FC_FLAGS are from rcm02:~/rtm/crtm_v2.2.3/config_setup/ifort.setup.cira
#
#FC_FLAGS = -O3 -convert big_endian -free -assume byterecl -fp-model source -traceback -I$(CRTM_DIR)/include/crtm_v$(VERSION)
#
#FL_FLAGS = -L$(CRTM_DIR)/lib -lcrtm_v$(VERSION)
#FLAGS = -O3 -convert big_endian -free -assume byterecl -fp-model source -traceback -xHOST
FLAGS = -O3 -free -assume byterecl -fp-model source -traceback -xHOST
INC = -I$(CRTM_DIR)/include/crtm_v$(VERSION)
LIB = -L$(CRTM_DIR)/lib -lcrtm_v$(VERSION)
endif
#
###########################################
#                                         #
#                                         #
#          Additional setup               #
#                                         #
#                                         #
###########################################
FC = ifort
FL = $(FC)
# The executable and object files
#SOURCES = Get_CRTM.f90
#OBJECTS = ${SOURCES:.f90=.o}
#EXE_FILE = Get_CRTM
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


#$(EXE_FILE): $(OBJECTS)
#	$(FL) $(OBJECTS) $(FL_FLAGS) -o $(EXE_FILE)
#
#Get_CRTM.o: Get_CRTM.f90
#	$(FC) $(FC_FLAGS) Get_CRTM.f90 -o Get_CRTM.o $(FL_FLAGS)
#
# A clean target
#
#clean:
#	rm -f $(OBJECTS) $(EXE_FILE)
