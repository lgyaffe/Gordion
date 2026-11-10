#
# Makefile for multiple executables with common code base.
#
# Include on command line:
#   DBG=1 for debugging compilation with -g
#   NUM64=1 for 64-bit observable numbering
#   MASSSCALE=<exponent> to scale fermion masses with a power of gauge coupling
#
# LGY 11/2026

PRGM	:= gordion

SRC_DIR := Src
OBJ_DIR := Obj
DEP_DIR := Dep
BIN_DIR := Bin
PROGRAM := $(shell echo "$(PRGM)" | tr '[:lower:]' '[:upper:]')
PREVTHY := $(shell test -h $(PRGM) && readlink $(PRGM) | sed 's/.*$(PRGM)-//')

# Set theory target to PREVTHY if defined:
ifndef THY
THY := $(PREVTHY)
endif

TARGET  := $(BIN_DIR)/$(PRGM)-$(THY)
SRCS	:= $(wildcard $(SRC_DIR)/*.cpp)
OBJS	:= $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
DEPS	:= $(SRCS:$(SRC_DIR)/%.cpp=$(DEP_DIR)/%.d)

# Clean object files & replace target symlink if new theory target
ifneq ($(THY),$(PREVTHY))
ifneq ($(wildcard $(PRGM)),)
    CLEAN   = clean
endif
    RELINK  = relink
else
    CLEAN   =
    RELINK  =
endif

CXX      := c++
CXXFLAGS := -std=c++20 -I. -D$(PROGRAM) -DARMA_DONT_USE_WRAPPER
CPPFLAGS  = -MMD -MP -MF $(DEP_DIR)/$*.d
OS	 := $(shell uname)

ifeq ($(OS), Darwin)
LDFLAGS   = -framework Accelerate
else
LDFLAGS   = -lopenblas
endif

ifdef DBG
CXXFLAGS += -O0 -g
else
CXXFLAGS += -O3 -DPARALLEL
LDFLAGS  += -ltbb
endif

ifdef NUM64
CXXFLAGS += -DNUM64
endif

ifdef MASSSCALE
CXXFLAGS += -DMASSSCALE=$(MASSSCALE)
endif

# Add theory target, if set, to preprocessor defines, else gripe:
ifdef THY
CXXFLAGS += -D$(THY)
OOPS =
else
OOPS = oops
endif

all: $(TARGET)

$(TARGET): $(OOPS) $(CLEAN) $(OBJS) $(RELINK)
	$(CXX) $(OBJS) $(LDFLAGS) -o $(TARGET)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(DEP_DIR)/%.d | $(OBJ_DIR) $(DEP_DIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(DEPS):

$(OBJ_DIR) $(DEP_DIR):
	@mkdir -p $@

-include $(DEPS)

clean:
	rm -rf $(OBJ_DIR) $(DEP_DIR) $(TARGET)

relink:
	ln -f -s $(TARGET) $(PRGM)

oops:
	@echo "!!! Unknown theory: must give THY=<thy> command arg, for some <thy> defined in Theory.h !!!"
	@false

.PHONY: clean all relink oops
