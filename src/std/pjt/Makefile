

CXXFLAGS = -g

ALL = Makefile Makefile.pjt  README.pjt 	program.h $(SRC) 	test1.cc test2.cc


SRC = 	getparam.cc main.cc 	program.cc \
	p_description.cc p_examples.cc p_keywords.cc p_usage.cc p_version.cc \


OBJ = $(SRC:.cc=.o)

help:
	@echo OBJ = $(OBJ)

lib: $(OBJ)
	ar ruv libp.a $(OBJ)

test1: lib test1.cc
	g++ -g -o test1 test1.cc libp.a

test2: lib test2.cc
	g++ -g -o test2 test2.cc libp.a


clean:
	rm -f *.o test1 test2 libp.a

tar:
	tar cvf pjt.tar $(ALL)
