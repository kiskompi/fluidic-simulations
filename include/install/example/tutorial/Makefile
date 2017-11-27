KOKKOS_SETTINGS=KOKKOS_SRC_PATH=/home/akomporday/kokkos GTEST_PATH=/home/akomporday/kokkos/tpls/gtest KOKKOS_PATH=/home/akomporday/kokkos/install/install

build:
	$(MAKE) -f /home/akomporday/kokkos/example/tutorial/Makefile KOKKOS_SETTINGS='KOKKOS_SRC_PATH=/home/akomporday/kokkos GTEST_PATH=/home/akomporday/kokkos/tpls/gtest KOKKOS_PATH=/home/akomporday/kokkos/install/install' KOKKOS_PATH=/home/akomporday/kokkos build

test: build
	$(MAKE) -f /home/akomporday/kokkos/example/tutorial/Makefile KOKKOS_SETTINGS='KOKKOS_SRC_PATH=/home/akomporday/kokkos GTEST_PATH=/home/akomporday/kokkos/tpls/gtest KOKKOS_PATH=/home/akomporday/kokkos/install/install' KOKKOS_PATH=/home/akomporday/kokkos test

clean:
	$(MAKE) -f /home/akomporday/kokkos/example/tutorial/Makefile KOKKOS_SETTINGS='KOKKOS_SRC_PATH=/home/akomporday/kokkos GTEST_PATH=/home/akomporday/kokkos/tpls/gtest KOKKOS_PATH=/home/akomporday/kokkos/install/install' KOKKOS_PATH=/home/akomporday/kokkos clean
