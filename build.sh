g++ -fopenmp -std=c++11 main.o simulation.o -o build/runnable -I/home/akomporday/git-project/include/ -L/home/akomporday/git-project/lib/ -lkokkos -Larray -ldl -lpthread  -O3
