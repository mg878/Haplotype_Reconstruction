CC	      = g++
CC_FLAGS	= -g3 -O3 -Wall -I/usr/local/include/ -std=c++11
LD_FLAGS	= -L/usr/local/lib -lm  -lgsl
OBJECTS    = test.o 
TEST_2   = test_2.o
TEST_3   = test_3.o
HIV_TEST   = hiv_test.o
BOTTLENECK   = bottleneck.o
REAL_BOTTLENECK   = real_bottleneck.o
KNOWN_HAPS   = known_haps.o
GOOD_BAD   = good_bad.o

test:  $(OBJECTS)
	$(CC) $(CC_FLAGS) $(OBJECTS) -o run_code $(LD_FLAGS)
test.o: test.cpp
	$(CC) $(CC_FLAGS) -c test.cpp

test_2:  $(TEST_2)
	$(CC) $(CC_FLAGS) $(TEST_2) -o run_test_2 $(LD_FLAGS)
test_2.o: test_2.cpp
	$(CC) $(CC_FLAGS) -c test_2.cpp

test_3:  $(TEST_3)
	$(CC) $(CC_FLAGS) $(TEST_3) -o run_test_3 $(LD_FLAGS)
test_3.o: test_3.cpp
	$(CC) $(CC_FLAGS) -c test_3.cpp
		
hiv_test:  $(HIV_TEST)
	$(CC) $(CC_FLAGS) $(HIV_TEST) -o run_hiv $(LD_FLAGS)
hiv_test.o: hiv_test.cpp
	$(CC) $(CC_FLAGS) -c hiv_test.cpp
	
bottleneck:	$(BOTTLENECK)
	$(CC) $(CC_FLAGS) $(BOTTLENECK) -o run_bottleneck $(LD_FLAGS)
bottleneck.o: bottleneck.cpp
	$(CC) $(CC_FLAGS) -c bottleneck.cpp

real_bottleneck:	$(REAL_BOTTLENECK)
	$(CC) $(CC_FLAGS) $(REAL_BOTTLENECK) -o run_real_bottleneck $(LD_FLAGS)
real_bottleneck.o: real_bottleneck.cpp
	$(CC) $(CC_FLAGS) -c real_bottleneck.cpp
	
known_haps:	$(KNOWN_HAPS)
	$(CC) $(CC_FLAGS) $(KNOWN_HAPS) -o run_known_haps $(LD_FLAGS)
known_haps.o: known_haps.cpp
	$(CC) $(CC_FLAGS) -c known_haps.cpp
	
good_bad:	$(GOOD_BAD)
	$(CC) $(CC_FLAGS) $(GOOD_BAD) -o run_good_bad $(LD_FLAGS)
good_bad.o: good_bad.cpp
	$(CC) $(CC_FLAGS) -c good_bad.cpp


