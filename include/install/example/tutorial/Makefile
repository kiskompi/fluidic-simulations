KOKKOS_SETTINGS=KOKKOS_SRC_PATH=/home/akomporday/git-project/include KOKKOS_DEVICES=OpenMP GTEST_PATH=/home/akomporday/git-project/include/tpls/gtest KOKKOS_PATH=/home/akomporday/git-project/include/install/install

build:
	$(MAKE) -f /home/akomporday/git-project/include/example/tutorial/Makefile KOKKOS_SETTINGS='KOKKOS_SRC_PATH=/home/akomporday/git-project/include KOKKOS_DEVICES=OpenMP GTEST_PATH=/home/akomporday/git-project/include/tpls/gtest KOKKOS_PATH=/home/akomporday/git-project/include/install/install' KOKKOS_PATH=/home/akomporday/git-project/include build

test: build
	$(MAKE) -f /home/akomporday/git-project/include/example/tutorial/Makefile KOKKOS_SETTINGS='KOKKOS_SRC_PATH=/home/akomporday/git-project/include KOKKOS_DEVICES=OpenMP GTEST_PATH=/home/akomporday/git-project/include/tpls/gtest KOKKOS_PATH=/home/akomporday/git-project/include/install/install' KOKKOS_PATH=/home/akomporday/git-project/include test

clean:
	$(MAKE) -f /home/akomporday/git-project/include/example/tutorial/Makefile KOKKOS_SETTINGS='KOKKOS_SRC_PATH=/home/akomporday/git-project/include KOKKOS_DEVICES=OpenMP GTEST_PATH=/home/akomporday/git-project/include/tpls/gtest KOKKOS_PATH=/home/akomporday/git-project/include/install/install' KOKKOS_PATH=/home/akomporday/git-project/include clean
