OBJS = src/fast_misc.o

CC = gcc
ifdef notshared
ifeq ($(notshared),yes)
EXEC = main.exe
CFLAGS =
OFLAGS =
endif
else
EXEC = src/c_fast_misc.so
CFLAGS = -fPIC
OFLAGS = -shared
endif

INCL = -I/home/tom/code/gsl/include/ -fopenmp -O3
LIBS = -lgsl -lgslcblas -L/home/tom/code/gsl/lib -lm -fopenmp -O3
.SUFFIXES : .c .o
%.o: %.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

#$(OBJS): $(INCL)

.PHONY : clean

clean:
	rm -f $(OBJS) main.exe src/c_fast_misc.so
	rm -f *~ src/*~ src/*.pyc
