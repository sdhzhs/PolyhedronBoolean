# target, subdir, objects in current dir
TARGET = CGALNefOp
LIB = libNefOp.so
SUBDIRS = src
SOURCES = $(wildcard ${SUBDIRS}/*.cpp)
OBJECTS := $(patsubst %.cpp,%.o,${SOURCES})
MAINOBJ := ${SUBDIRS}/Main.o
LIBOBJS := $(subst ${MAINOBJ},,${OBJECTS})
LDFLAGS := -lgmp -lmpfr
LIBFLAGS := -L. -lNefOp
LDD := --shared -fPIC -O3

all:subdirs ${LIB} ${OBJECTS}
	${CC} ${MAINOBJ} ${LIBFLAGS} ${LDFLAGS} -o ${TARGET}

${LIB}: ${LIBOBJS}
	${CC} ${LDD} ${LIBOBJS} -o ${LIB}

clean:cleansubdirs
	rm -f ${TARGET} ${LIB}

# path of "make global scripts"
# NOTE, use absolute path. export once, use in all subdirs
export PROJECTPATH=${PWD}
export MAKEINCLUDE=${PROJECTPATH}/make/make.global

# include "make global scripts"
include ${MAKEINCLUDE}
