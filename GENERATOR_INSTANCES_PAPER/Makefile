# COMPILER

CXX= g++ 

# FILES

OBJS =	src/Main.o   
		
DBG= -O3
#DBG= -g

CXXFLAGS =  $(DBG)

.c.o:
	gcc -c $(CXXFLAGS) $< -o $@


all:GEN

GEN: $(OBJS)
		$(CXX) $(CXXFLAGS) -o $@ $(OBJS) 

$(OBJS): Makefile

clean:
	rm -f $(OBJS) rm GEN
