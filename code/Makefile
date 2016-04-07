include Makefile.edison

SRC=../src
CXXFLAGS+=-I$(SRC)
CFLAGS+=-I$(SRC)

CPPEXAMPLES= matrixfree schur solve power product communicator selection

all: $(CPPEXAMPLES) c_example f90_example mex_example

cpp_examples: $(CPPEXAMPLES)

include $(SRC)/Makefile.SDP

$(CPPEXAMPLES): % : %.cpp $(RRQR) $(headers)
	$(CXX) $(CXXFLAGS) -o $@ $< $(RRQR) $(LIB)

c_example: % : %.c $(RRQR) $(headers) $(c_interface)
	$(CC) $(CFLAGS) -c $<
	$(CXX) -o $@ $@.o $(c_interface) $(RRQR) $(LIB)

f90_example: % : %.F90 $(RRQR) $(headers) $(c_interface) $(f90_interface)
	$(FC) $(FFLAGS) -c $<
	$(FC) $(FFLAGS) -o $@ $@.o $(c_interface) $(f90_interface) $(RRQR) $(LIB) $(LIBCXX)

mex_example: $(RRQR) $(headers) $(mex_interface)
	echo "Done"
clean:
	-rm $(SRC)/*.o *.o *.mex *.mod
	-rm $(RRQR)
	-rm $(CPPEXAMPLES) c_example f90_example

realclean: clean
