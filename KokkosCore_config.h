/* ---------------------------------------------
Makefile constructed configuration:
Mon Dec 11 01:08:33 CET 2017
----------------------------------------------*/
#if !defined(KOKKOS_MACROS_HPP) || defined(KOKKOS_CORE_CONFIG_H)
#error "Do not include KokkosCore_config.h directly; include Kokkos_Macros.hpp instead."
#else
#define KOKKOS_CORE_CONFIG_H
#endif
/* Execution Spaces */
#define KOKKOS_HAVE_CUDA 1
#define KOKKOS_HAVE_SERIAL 1
#ifndef __CUDA_ARCH__
#define KOKKOS_USE_ISA_X86_64
#endif
/* General Settings */
#define KOKKOS_HAVE_CXX11 1
#define KOKKOS_ENABLE_PROFILING
/* Optimization Settings */
/* Cuda Settings */
#define KOKKOS_CUDA_USE_LAMBDA 1
#define KOKKOS_ARCH_AVX 1
#define KOKKOS_ARCH_KEPLER 1
#define KOKKOS_ARCH_KEPLER35 1
