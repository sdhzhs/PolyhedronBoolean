# target, subdir, objects in current dir
TARGET = cgaltest
SUBDIRS = src
SOURCES = $(wildcard ${SUBDIRS}/*.cpp)
OBJECTS := $(patsubst %.cpp,%.o,${SOURCES})
LDFLAGS = -lgmp -lmpfr

all:subdirs ${OBJECTS}
	${CC} ${OBJECTS} ${LDFLAGS} -o ${TARGET}

clean:cleansubdirs
	rm -f ${TARGET} ${OBJECTS}

# path of "make global scripts"
# NOTE, use absolute path. export once, use in all subdirs
export PROJECTPATH=${PWD}
export MAKEINCLUDE=${PROJECTPATH}/make/make.global

# include "make global scripts"
include ${MAKEINCLUDE}
