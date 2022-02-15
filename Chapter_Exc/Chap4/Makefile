default: Exc_4_3

libs = -lblas -lpython3.8 -lboost_program_options

Include_Paths = -I/usr/include/python3.8  -I/usr/lib/python3/dist-packages/numpy/core/include/numpy/ 

Exc_4_3: Exc_4_3.o Exc_4_3_funcs.o
	$(CXX) -o $@ $^ $(libs) $(Include_Paths)

%.o: %.cpp Exc_4_3_funcs.hpp
	$(CXX) -o $@ -c $< $(Include_Paths)

.PHONY: clean
clean:
	rm -f Exc_4_3*.o
