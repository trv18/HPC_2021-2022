default: myprog

libs = -lblas

myprog: compiler-main.o compiler-norm.o
	$(CXX) -o $@ $^ $(libs)

%.o: %.cpp compiler-norm.h
	$(CXX) -o $@ -c $<

.PHONY: clean
clean: 
	rm -f *.o myprog

#compiler-main.o: compiler-main.cpp compiler-norm.h
# 	g++ -o $@ -c $<
