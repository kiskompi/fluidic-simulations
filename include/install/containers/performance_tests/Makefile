KOKKOS_SETTINGS=KOKKOS_SRC_PATH=/home/akomporday/kokkos GTEST_PATH=/home/akomporday/kokkos/tpls/gtest KOKKOS_PATH=/home/akomporday/kokkos

all:
	$(MAKE) -f /home/akomporday/kokkos/containers/performance_tests/Makefile KOKKOS_SRC_PATH=/home/akomporday/kokkos GTEST_PATH=/home/akomporday/kokkos/tpls/gtest KOKKOS_PATH=/home/akomporday/kokkos

test: all
	$(MAKE) -f /home/akomporday/kokkos/containers/performance_tests/Makefile KOKKOS_SRC_PATH=/home/akomporday/kokkos GTEST_PATH=/home/akomporday/kokkos/tpls/gtest KOKKOS_PATH=/home/akomporday/kokkos test

clean:
	$(MAKE) -f /home/akomporday/kokkos/containers/performance_tests/Makefile KOKKOS_SRC_PATH=/home/akomporday/kokkos GTEST_PATH=/home/akomporday/kokkos/tpls/gtest KOKKOS_PATH=/home/akomporday/kokkos clean
