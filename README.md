## KerM

Kernel Mixing Model for Turbulent Combustion and Transport PDF method.

This repo is for 0-D simple validation, for PaSR validation please refer to [PaSR](https://github.com/SuXY15/PaSR)

### 1. Usage

+ Python Version:
  + [models.py](models.py): implementation of mixing models
  + [mixing_py.py](mixing_py.py): run simulation

  ```shell
  # run simulation, taking around 10 minute for N=1000
  python mixing_py.py
  ```

  Comparison results of 1k particles for EMST and 10k particles for other models (uniform weighted samples, KerM sigma_k=0.3)
  
  ![](figs/py_PoF_1996_Fig9b_comparison_uniform_1000&10000.png)

+ C++ Version:
  + [src/MixingModels.hpp](src/MixingModels.hpp): implementation of mixing models
  + [src/main.cpp](src/main.cpp): run simulation
  + [performance_c.py](performance_c.py): show results

  ```shell
  # build executable file
  make 
  # run simulation, taking around 10 secdons for N=1000
  ./mix
  # show results
  python performance_c.py
  ```

  Comparison results of 4k particles for EMST and 100k particles for other models (KerM sigma_k=0.25)

  ![](figs/c++_PoF_1996_Fig9b_comparison_4000&100000.png)
  
  Performance of mixing models
  
  <img src="figs/c++_PoF_1996_Fig9b_performance.png" alt="c++_PoF_1996_Fig9b_performance" style="width:50%;" />
  
  

### 2. Implementations

Please refer to [UserGuide.pdf](UserGuide.pdf)


