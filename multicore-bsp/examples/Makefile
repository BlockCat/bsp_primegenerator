include ../include.mk

RELPATH = ../

INCLUDE_DIR = ../include/
LFLAGS := ../lib/libmcbsp${VERSION}.a -lm ${LFLAGS}

CFLAGS := ${CFLAGS} -I${INCLUDE_DIR}
CPPFLAGS := ${CPPFLAGS} -I${INCLUDE_DIR}

all: ../lib/libmcbsp${VERSION}.a multistart parallel_loop hello_world hello_world_cpp hierarchical

../lib/libmcbsp${VERSION}.a:
	(cd ..; make)

multistart: multistart.o
	${CC} ${CFLAGS} -o $@ $^ ${LFLAGS} 

parallel_loop: parallel_loop.o
	${CC} ${CFLAGS} -o $@ $^ ${LFLAGS} 

hello_world: hello_world.o
	${CC} ${CPPFLAGS} -o $@ $^ ${LFLAGS} 

hello_world_cpp: hello_world_cpp.opp
	${CPP} ${CPPFLAGS} -o $@ $^ ${LFLAGS} 

hierarchical: hierarchical.o
	${CC} ${CPPFLAGS} -o $@ $^ ${LFLAGS} 

clean:
	rm -f multistart.o parallel_loop.o hello_world.o hello_world_cpp.opp hierarchical.o multistart parallel_loop hello_world hello_world_cpp hierarchical

