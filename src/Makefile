# subdir and objects in current dir
SUBDIRS = 

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,${SOURCES})

ifdef MACRO
  CXXFLAGS += -D${MACRO}
endif

all:subdirs ${OBJECTS}

clean:cleansubdirs
	rm -f ${OBJECTS} ${OBJECTS:.o=.d}

include ${MAKEINCLUDE}
include $(OBJECTS:.o=.d)
