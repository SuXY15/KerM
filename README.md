## KerM

Kernel Mixing Model for Turbulent Combustion and Transport PDF method.

This repo is for 0-D simple validation, for PaSR validation please refer to [PaSR](https://github.com/SuXY15/PaSR)

### 1. Usage

+ C++ Version:
  + [src_cpp/MixingModels.hpp](src_cpp/MixingModels.hpp): implementation of mixing models
  + [src_cpp/main.cpp](src_cpp/main.cpp): run simulation
  + [performance.py](performance.py): show results

  ```shell
  # build executable file
  make
  # run simulation, taking around 10 seconds for N=1000
  ./mix
  # show results
  python performance.py
  ```

  Comparison results of 4k particles for EMST and 100k particles for other models (`KerM` sigma_k=0.25)

  ![](figs/comparison_cpp_PoF_1996_Fig9b_uniform_4000&50000.png)
  
  Performance of mixing models (`EMST-D` do not account aging properties and no `IEM` assisted, for original `EMST` implementation, please refer to the Fortran version)
  
  <img src="figs/performance_cpp_PoF_1996_Fig9b.png" style="width:50%;" />

+ Fortran version:
  + [src_fortran/mixing_test.f90](src_fortran/mixing_test.f90): main program for testing.
  + [src_fortran/mixing_model_kerm.f90](src_fortran/mixing_model_kerm.f90): source code of `KerM` with differential diffusion supported.
  + [src_fortran/mixing_model_iem.f90](src_fortran/mixing_model_iem.f90): source code of `IEM`, used for comparison.
  + [src_fortran/mixing_model_mcurl.f90](src_fortran/mixing_model_mcurl.f90): source code of `MC`, used for comparison.
  + [src_fortran/emst.f](src_fortran/emst.f) and [src_fortran/emst_subs.f](src_fortran/emst_subs.f): source code of `EMST` adopted from [^1], used for comparison.
  ```shell
  # build executable file (single precision float)
  make fortran
  # run simulation, faster than cpp version
  ./mix
  # show results
  python performance.py fortran
  ```
  + Comparison results of 10k particles for EMST and 50k particles for other models (`KerM` sigma_k=0.25)
  ![](figs/comparison_fortran_PoF_1996_Fig9b_uniform_10000&50000.png)
  + Performance of mixing models
  <img src="figs/performance_fortran_PoF_1996_Fig9b.png" style="width:50%;" />
  
+ Python Version:
  + [src_python/models.py](src_python/models.py): implementation of mixing models
  + [src_python/mixing_py.py](src_python/mixing_py.py): run simulation

  ```shell
  # run simulation, EMST taking around 10 minute for N=1000
  python src_python/mixing_py.py
  ```
  Comparison results of 1k particles for EMST and 10k particles for other models (uniform weighted samples, `KerM` sigma_k=0.25)
  
  ![](figs/comparison_python_PoF_1996_Fig9b_uniform_1000&10000.png)
  

### 2. Implementations

Please refer to [TheoryGuide.pdf](TheoryGuide.pdf)

[^1]: Original EMST implementation https://tcg.mae.cornell.edu/emst/
