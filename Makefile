include ./arch.gnu
# OPTIMIZATION = -fast
# OPTIMIZATION = -O3
# DEBUG += -g
#mpicxx -I/kuacc/apps/openmpi/4.0.1/include -c cardiacsim.C -o cardiacsim.o


app:		cardiacsim

OBJECTS = cardiacsim.o splot.o cmdLine.o

cardiacsim:	        $(OBJECTS) 
		$(C++LINK) $(LDFLAGS) -o $@ $(OBJECTS)  $(LDLIBS)

clean:	
	$(RM) *.o cardiacsim *~;
	$sdsa(RM) core;
