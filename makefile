#c compiler command
CC = g++-6.2.0

#path for c compiler
CPATH = /usr/local/gcc-6.2.0/bin/

#path for eigen libs
EIGENPATH = /usr/local/include

#other options
COPTIONS = -O2 -g -Wall

#local build directories (shouldn't need to be changed)
IDIR = include
CDIR = src
ODIR = obj

CFLAGS = -I $(IDIR) -I $(EIGENPATH)
LIBS =

_DEPS = lung_model_discrete_branching.hpp;
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

SRC:= $(wildcard $(CDIR)/*.cpp)
_OBJ:= $(notdir $(SRC:%.cpp=%.o))
OBJ:= $(patsubst %, $(ODIR)/%, $(_OBJ))
BIN= PULMsim_linux

$(warning $(RESULT3))

$(BIN): $(OBJ) $(DEPS)
	$(CPATH)$(CC) $^ -o $@ $(CFLAGS) $(COPTIONS) $(LIBS)

$(ODIR)/%.o: $(CDIR)/%.cpp
	$(CPATH)$(CC) -c $^ -o $@ $(CFLAGS) $(COPTIONS) $(LIBS)

