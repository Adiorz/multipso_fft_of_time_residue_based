main: additional.cpp  data_preprosessor.cpp  fft.cpp  main.cpp  pso.cpp  scheduler.cpp  worker.cpp barrier.cpp
	g++ additional.cpp  data_preprosessor.cpp  fft.cpp  main.cpp  pso.cpp  scheduler.cpp  worker.cpp barrier.cpp -std=c++14 -pedantic -lpthread -o main -O3

debug: additional.cpp  data_preprosessor.cpp  fft.cpp  main.cpp  pso.cpp  scheduler.cpp  worker.cpp barrier.cpp
	g++ additional.cpp  data_preprosessor.cpp  fft.cpp  main.cpp  pso.cpp  scheduler.cpp  worker.cpp barrier.cpp -std=c++14 -pedantic -lpthread -o main

profile: additional.cpp  data_preprosessor.cpp  fft.cpp  main.cpp  pso.cpp  scheduler.cpp  worker.cpp barrier.cpp
	g++ additional.cpp  data_preprosessor.cpp  fft.cpp  main.cpp  pso.cpp  scheduler.cpp  worker.cpp barrier.cpp -std=c++14 -pedantic -lpthread -o main -pg

test: additional.cpp  data_preprosessor.cpp  fft.cpp  main.cpp  pso.cpp  scheduler.cpp  worker.cpp barrier.cpp
	g++ -std=c++14 -pedantic additional.cpp  data_preprosessor.cpp  fft.cpp  main.cpp  pso.cpp  scheduler.cpp  worker.cpp barrier.cpp -lpthread -o test -O3

after_copernicon: additional.cpp  data_preprosessor.cpp  fft.cpp  main.cpp  pso.cpp  scheduler.cpp  worker.cpp barrier.cpp
	g++ -std=c++14 -pedantic additional.cpp  data_preprosessor.cpp  fft.cpp  main.cpp  pso.cpp  scheduler.cpp  worker.cpp barrier.cpp -lpthread -o after_copernicon -O3

clean:
	rm -f main
