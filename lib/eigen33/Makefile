# FORTRAN compiler
# ----------------
F77           = gfortran
F77FLAGS      = -fPIC -Wall -O3

# Compiler and linker options
# ---------------------------
LIBFLAGS      =
INCFLAGS      =
LOCAL_LDFLAGS =

# Project files
# -------------
SONAME = lib3x3.so.1
BIN    = $(SONAME).0.0
OBJ    = sqrabs.o slvsec3.o dsyev2.o dsytrd3.o dsyevc3.o dsyevj3.o dsyevq3.o dsyevv3.o \
	 dsyevh3.o dsyevd3.o zhetrd3.o zheevc3.o zheevj3.o zheevq3.o zheevv3.o zheevh3.o zheevd3.o
HEADER =


# Rules
# -----
$(BIN): Makefile $(OBJ)
	$(F77)  -shared -Wl,-soname,$(SONAME) $(F77FLAGS) $(INCFLAGS) $(OBJ) $(LIBFLAGS) \
                $(LDFLAGS) $(LOCAL_LDFLAGS) -o $(BIN)

%.o : %.f Makefile $(HEADER)
	$(F77)  $(F77FLAGS) -c $(INCFLAGS) $<

.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ)



