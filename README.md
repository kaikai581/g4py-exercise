# g4py-exercise
My exercise in Geant4 made easy with python

## Geant4Py Installation
### Site-package g4py
Before you undertake installation, make sure whether or not your reference scripts import a module named `g4py`. As of writing, this "site-package" is no longer available since `Geant4.10.06`, and the developers have no plan to make it available any time soon.
#### Install g4py alongside Geant4Py
If you decide `g4py` is a necessity of your project, your best bet is to install `Geant4.10.05`. Follow the steps below to compile `Geant4.10.05` first.
```
$ wget https://geant4-data.web.cern.ch/releases/geant4.10.05.p01.tar.gz
$ tar zxf geant4.10.05.p01.tar.gz
$ cd geant4.10.05.p01
$ mkdir build install
$ cd build
$ # You might want to enable OpenGL since Geant4Py by default depends on libG4OpenGL.so
$ cmake -DGEANT4_INSTALL_DATA=ON -DCMAKE_INSTALL_PREFIX=../install -DGEANT4_USE_GDML=ON -DGEANT4_USE_OPENGL_X11=ON ../
$ make -j`nproc`
$ make install
$ source ../install/bin/geant4.sh
```
After Geant4 is installed, stay in `build` directory, and follow the steps below.
```
$ mkdir g4py
$ cmake -DCMAKE_INSTALL_PREFIX=../../install -DBoost_LIBRARY_DIRS:FILEPATH=<your boost_python path> ../../environments/g4py
$ make -j`nproc`
$ make install
```
Lastly, in order for your python environment to find `g4py`, include the following paths to your `PYTHONPATH`.
```
$ export PYTHONPATH=../../install/lib64:../../install/lib64/examples:../../install/lib64/tests:$PYTHONPATH
```
### Install only Geant4Py without superuser privileges
If the modules you will use are all included in `Geant4Py`, then you can install newer Geant4 releases easily.

If you are like me who work with a remote server, it's unlikely you have root access. In this case, [Anaconda](https://www.anaconda.com/products/individual) comes to your rescue. Note that the third-party repository `conda-forge` does have Geant4 binaries. However, as far as I know, those binaries don't have `Geant4Py` available.

Install Anaconda locally, and install the packages necessary to compile Geant4.
```
$ conda create -n g4py
$ conda activate g4py
$ conda config --env --add channels conda-forge
$ conda install --yes root boost xersec-c mesa-libgl-devel-cos6-x86_64 xorg-libxmu
```
Starting from `Geant4.10.06`, `Geant4Py` can be easily built with the main Geant4 by invoking the cmake option `GEANT4_USE_PYTHON` when building the main application. Here are reference steps.
```
$ wget https://geant4-data.web.cern.ch/releases/geant4.10.07.p01.tar.gz
$ tar zxf geant4.10.07.p01.tar.gz
$ cd geant4.10.07.p01
$ mkdir build install
$ cd build
$ cmake -DGEANT4_INSTALL_DATA=ON -DCMAKE_INSTALL_PREFIX=../install -DGEANT4_USE_GDML=ON -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_USE_PYTHON=ON ../
$ make -j`nproc`
$ make install
$ source ../install/bin/geant4.sh
```
After successfully compiled, include the module path to `PYTHONPATH`.
```
export PYTHONPATH=../install/lib64/python3.7/site-packages:$PYTHONPATH
```
Change `python3.7` to whatever version you have in your `install` path.

## An almost minimal Geant4Py simulation
Can be found in `p-on-be/p-on-be.py`.
Usage of the script is as simple as
```
$ ./p-on-be.py -n 1000
```
Use `./p-on-be.py -h` for all available options.

`Geant4Py` is amazing in the sense that a minimal C++ counterpart of this 200-line script could easily involve 4 header files, 5 source files, and many times more lines of code.

### Known issues
This script eventually will encounter segmentation fault that I cannot quite clear yet, which renders this script not very useful in practical purposes.
Here are some notes for the debugging processes I have gone through.

In order to save the core dump file, this command has to be issued.
```
$ ulimit -c unlimited
```
After core is dumped, use the command to examine it.
```
$ gdb python <core_filename>
```
Relevant output is excerpted here.
```
Core was generated by `python ./p_on_be.py -n 2000'.
Program terminated with signal 11, Segmentation fault.
#0  0x00007ff32efab57d in G4CrossSectionDataStore::ComputeCrossSection(G4DynamicParticle const*, G4Material const*) ()
   from /cshare/vol2/users/shihkai/software/geant4/geant4.10.07-install/lib64/python3.7/site-packages/Geant4/../../../libG4processes.so
```
Does the cross section calculation encounter some corner cases? Is my compilation faulty due to dependency setup?

All these might require arduous investigation.
