OBJS = fast_miscentering.o

CC = gcc
ifdef notshared
ifeq ($(notshared),yes)
EXEC = main.exe
CFLAGS =
OFLAGS =
endif
else
EXEC = fast_miscentering.so
CFLAGS = -fPIC
OFLAGS = -shared
endif

INCL = -I/${GSLI} -fopenmp -O3
LIBS = -lgsl -lgslcblas -L/${GSLL} -lm -fopenmp -O3
.SUFFIXES : .c .o
%.o: %.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

#$(OBJS): $(INCL)

.PHONY : clean

clean:
	rm -f $(OBJS) main.exe fast_miscentering.so
	rm -f *~ src/*~ *.pyc
