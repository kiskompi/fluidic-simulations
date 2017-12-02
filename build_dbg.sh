clear
g++ -fopenmp -std=c++14 main.cpp Simulation/simulation.cpp -g -O0 -o build/runnable_dbg -I/home/akomporday/git-project/include/ -L /home/akomporday/git-project/lib/ -lkokkos -Larray -ldl -lpthread  
