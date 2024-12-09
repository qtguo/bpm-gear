TARGET=bpm
all: ${TARGET}

CFLAGS+=-march=native -std=c++11 -O3 -DBATCH_SIZE=16 


${TARGET}: *.cpp *.h
	g++ ${CFLAGS} bpm.cpp -o ${TARGET} ${LDFLAGS} -ldl -lpthread

clean:
	rm ${TARGET}
