ObjSuf        = o
SrcSuf        = cc
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o
HeadSuf       = h

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lMinuit
ROOTGLIBS     = $(shell root-config --glibs) -lMinuit

# Linux with egcs
DEFINES       = -DNO_ORCA_CLASSES
CXX           = g++
CXXFLAGS	= -O -Wall -fPIC $(DEFINES)
LD			= g++
LDFLAGS		= -g
SOFLAGS		= -shared

LDDIR = $(shell echo "-L" $$LD_LIBRARY_PATH |sed -e 's/:/ -L/g')

CXXFLAGS	+= $(ROOTCFLAGS)
LIBS		= $(LDDIR) $(ROOTLIBS) -lRFIO -lEG
GLIBS		= $(ROOTGLIBS)
#------------------------------------------------------------------------------
SOURCES		= TH5.$(SrcSuf)
HEADERS		= ../interface/TH5.$(HeadSuf)
OBJECTS		= $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS		= $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS	= $(SOURCES:.$(SrcSuf)=.$(DllSuf))

all: libTH5.so clean

clean:
	@echo "Cleaning..."
	@rm -f $(OBJECTS) $(DEPENDS) *Dict.* core

.SUFFIXES: .$(SrcSuf) .C .o .so

###

Dict.$(SrcSuf): $(HEADERS) ./LinkDef.h
	@echo "Generating dictionary Dict..."
	@$(ROOTSYS)/bin/rootcint -f Dict.$(SrcSuf) -c $(DEFINES) $(HEADERS) ./LinkDef.h

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<

libTH5.so: $(OBJECTS) Dict.o
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@
	@echo "$@ done"
