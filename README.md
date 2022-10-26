# ViennaDD
***
This is a library to solve the Drift-Diffusion model. 

## Pre-Requisits
* `conan` (for dependency management)
* `eigen-3.4.0`
* `CMAKE` 
* `C++17 Compiler`
* `python3`
* `python3-pip`

## Supported Platforms
Currently, only Linux is supported and tested for. If you find problems for Windows/Mac please let us know how you solved them!
***

## Building the python-module 

* `mkdir build` -create a build directory
* `cd build` -change directory to build folder
* `cmake .. -D BUILD_PYTHON=ON` -create makefiles
* `make PSSolverPythonWheel` -build python library
* `cd wheel` - change to wheel directory
* `pip install -e .` - install python library

***

File Structure
-----------------

* `/wrapping` - contains the python library package
* `/example` - contains example code
* `/docs`- contains information for doxygen to build the documentation
* `/lib` - contains the C++ sourcecode
* `/include` - contains the C++ header-files
License
--------
ViennaPS is provided under a MIT license that can be found in the [MIT LICENSE](LICENSE) file.


Used Libraries
-----------------

[conan](https://github.com/conan-io/conan) ([MIT](https://raw.githubusercontent.com/conan-io/conan/develop/LICENSE.md))
[pybind11](https://github.com/pybind/pybind11/) ([BSD](https://raw.githubusercontent.com/pybind/pybind11/master/LICENSE))
[eigen](https://gitlab.com/libeigen/eigen) ([MPL2](https://gitlab.com/libeigen/eigen/-/raw/master/COPYING.MPL2))

## Documentation
To build the documentation:
```
cd build
cmake -D BUILD_DOCS=ON ..
make docs
```
--------------------------
Acknowledgment
----------------
This project was done with the help of the [Christian Doppler Forschungsgesellschaft](https://www.cdg.ac.at/)





