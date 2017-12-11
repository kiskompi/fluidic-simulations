/* ---------------------------------------------
Makefile constructed configuration:
Mon Dec 11 11:02:42 CET 2017
----------------------------------------------*/
#if !defined(KOKKOS_MACROS_HPP) || defined(KOKKOS_CORE_CONFIG_H)
#error "Do not include KokkosCore_config.h directly; include Kokkos_Macros.hpp instead."
#else
#define KOKKOS_CORE_CONFIG_H
#endif
/* Execution Spaces */
#define KOKKOS_HAVE_OPENMP 1
#ifndef __CUDA_ARCH__
#define KOKKOS_USE_ISA_X86_64
#endif
/* General Settings */
#define KOKKOS_HAVE_CXX11 1
#define KOKKOS_ENABLE_PROFILING
/* Optimization Settings */
/* Cuda Settings */
#define KOKKOS_ARCH_AVX 1
