#
# Makefile for multiple gordion-* executables with common code base.
# Somewhat based on https://spin.atomicobject.com/makefile-c-projects/
#
# LGY 4/2025

TARGET   = gordion
CXX      = c++

# Run "make DBG=1 ..." for compilation with -g:
ifdef DBG
CXXFLAGS = -O0 -g
LDFLAGS  = -larmadillo
else
CXXFLAGS = -O3 -DPARALLEL
LDFLAGS  = -ltbb -larmadillo
endif
CXXFLAGS += -std=c++20 -I.
# add -DMASSSCALE=<numeric_value> to CXXFLAGS if desired

# Extract previous theory target:
PREVTHY := $(shell test -h $(TARGET) && readlink $(TARGET) | sed 's/$(TARGET)-//')

# Set theory target to PREVTHY if defined:
ifndef THY
THY := $(PREVTHY)
endif

# Clean object files & replace target symlink if new theory target
ifneq ($(THY),$(PREVTHY))
CLEAN   = clean
RELINK  = relink
else
CLEAN   =
RELINK  =
endif

# Add theory target, if set, to preprocessor defines, else gripe:
ifdef THY
CXXFLAGS += -D$(THY)
OOPS =
else
OOPS = oops
endif

# Build list of source files:
SRCS  += Assess.cpp
SRCS  += Blab.cpp
SRCS  += Build.cpp
SRCS  += Canon.cpp
SRCS  += Commute.cpp
SRCS  += Coupling.cpp
SRCS  += Gen.cpp
SRCS  += Global.cpp
SRCS  += Init.cpp
SRCS  += Numerics.cpp
SRCS  += Obs.cpp
SRCS  += Ode.cpp
SRCS  += Op.cpp
SRCS  += Parse.cpp
SRCS  += Poly.cpp
SRCS  += Print.cpp
SRCS  += Rep.cpp
SRCS  += RK.cpp
SRCS  += Save.cpp
SRCS  += Symb.cpp
SRCS  += Symm.cpp
SRCS  += Term.cpp
SRCS  += Test.cpp
SRCS  += Theory.cpp
SRCS  += Gordion.cpp 

OBJS   = $(SRCS:.cpp=.o)
DEPS   = $(SRCS:.cpp=.d)

.PHONY: clean all relink oops

$(TARGET)-$(THY): $(OOPS) $(CLEAN) $(OBJS) $(DEPS) $(RELINK)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) -o $(TARGET)-$(THY)
#	codesign -s - -f --entitlements $(TARGET).entitlements $(TARGET)

$(TARGET).eneitlements:
	/usr/libexec/PlistBuddy -c "Add :com.apple.security.get-task-allow bool true" \
	$(TARGET).entitlements

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.d: %.cpp
	$(CXX) -M $(CXXFLAGS) $< > $@

clean:
	rm -f $(OBJS) $(TARGET)-$(THY) $(DEPS)

relink:
	ln -f -s $(TARGET)-$(THY) $(TARGET)

oops:
	@echo "!!! Unknown theory: must give THY=<thy> command arg, for some <thy> defined in Theory.h !!!"

	@false

all: $(TARGET)-$(THY)

-include $(DEPS)
