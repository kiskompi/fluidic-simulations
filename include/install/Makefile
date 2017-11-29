KOKKOS_SETTINGS=KOKKOS_SRC_PATH=/home/akomporday/git-project/include KOKKOS_DEVICES=OpenMP GTEST_PATH=/home/akomporday/git-project/include/tpls/gtest KOKKOS_PATH=/home/akomporday/git-project/include

kokkoslib:
	cd core; \
	$(MAKE) -f /home/akomporday/git-project/include/core/src/Makefile KOKKOS_SRC_PATH=/home/akomporday/git-project/include KOKKOS_DEVICES=OpenMP GTEST_PATH=/home/akomporday/git-project/include/tpls/gtest KOKKOS_PATH=/home/akomporday/git-project/include PREFIX=/home/akomporday/git-project/ build-lib

install: kokkoslib
	cd core; \
	$(MAKE) -f /home/akomporday/git-project/include/core/src/Makefile KOKKOS_SRC_PATH=/home/akomporday/git-project/include KOKKOS_DEVICES=OpenMP GTEST_PATH=/home/akomporday/git-project/include/tpls/gtest KOKKOS_PATH=/home/akomporday/git-project/include PREFIX=/home/akomporday/git-project/ install

kokkoslib-test:
	cd core; \
	$(MAKE) -f /home/akomporday/git-project/include/core/src/Makefile KOKKOS_SRC_PATH=/home/akomporday/git-project/include KOKKOS_DEVICES=OpenMP GTEST_PATH=/home/akomporday/git-project/include/tpls/gtest KOKKOS_PATH=/home/akomporday/git-project/include PREFIX=/home/akomporday/git-project/include/install/install build-lib

install-test: kokkoslib-test
	cd core; \
	$(MAKE) -f /home/akomporday/git-project/include/core/src/Makefile KOKKOS_SRC_PATH=/home/akomporday/git-project/include KOKKOS_DEVICES=OpenMP GTEST_PATH=/home/akomporday/git-project/include/tpls/gtest KOKKOS_PATH=/home/akomporday/git-project/include PREFIX=/home/akomporday/git-project/include/install/install install

build-test: install-test
	$(MAKE) -C core/unit_test
	$(MAKE) -C core/perf_test
	$(MAKE) -C containers/unit_tests
	$(MAKE) -C containers/performance_tests
	$(MAKE) -C algorithms/unit_tests
	$(MAKE) -C example/fixture
	$(MAKE) -C example/feint
	$(MAKE) -C example/fenl
	$(MAKE) -C example/tutorial build

test: build-test
	$(MAKE) -C core/unit_test test
	$(MAKE) -C core/perf_test test
	$(MAKE) -C containers/unit_tests test
	$(MAKE) -C containers/performance_tests test
	$(MAKE) -C algorithms/unit_tests test
	$(MAKE) -C example/fixture test
	$(MAKE) -C example/feint test
	$(MAKE) -C example/fenl test
	$(MAKE) -C example/tutorial test

unit-tests-only:
	$(MAKE) -C core/unit_test test
	$(MAKE) -C containers/unit_tests test
	$(MAKE) -C algorithms/unit_tests test

clean:
	$(MAKE) -C core/unit_test clean
	$(MAKE) -C core/perf_test clean
	$(MAKE) -C containers/unit_tests clean
	$(MAKE) -C containers/performance_tests clean
	$(MAKE) -C algorithms/unit_tests clean
	$(MAKE) -C example/fixture clean
	$(MAKE) -C example/feint clean
	$(MAKE) -C example/fenl clean
	$(MAKE) -C example/tutorial clean
	cd core; \
	$(MAKE) -f /home/akomporday/git-project/include/core/src/Makefile KOKKOS_SRC_PATH=/home/akomporday/git-project/include KOKKOS_DEVICES=OpenMP GTEST_PATH=/home/akomporday/git-project/include/tpls/gtest KOKKOS_PATH=/home/akomporday/git-project/include clean
