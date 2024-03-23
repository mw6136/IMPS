CC = g++
CFLAGS  = -flto -O3 -Wall

# define the c++ source files
SRCS = mesh.cpp pressure.cpp
LIBS = -lhdf5 -lhdf5_cpp

OBJS = $(SRCS:.cpp=.o)

MAIN = IMPS

.PHONY: depend clean

all:    $(MAIN)
		@echo  Compiled IMPS

$(MAIN): $(OBJS) 
		$(CC) $(CFLAGS) -o $(MAIN) main.cpp $(OBJS) $(LIBS)

.c.o:
		$(CC) $(CFLAGS) -c $<  -o $@

clean:
		$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
		makedepend $(INCLUDES) $^