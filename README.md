# Parallel Adaptive Swarm Search (PASS)

Version 0.10.0
==============

PASS is written in C++14 and uses [Armadillo](http://arma.sourceforge.net/) (developed by Conrad Sanderson et al., NICTA, Australia) for highly efficient linear algebra calculations.

#Prerequirements

- GCC 4.8.1 or later
- CMake 3.2 or later
- Armadillo C++ 8.200.0 or later

# Installation

## Linux Installation
The following commands are based on Ubuntu, using apt-get. If you are using another Linux distributions, the names of each package might also differ.

If you are using any recent Linux distribution, you should already be equipped with a C++11 feature complete compiler, like GCC 4.8.1+

1. Install Cmake

```
sudo apt-get install cmake
```

2. Install Armadillo C++ with OpenBLAS support (visit [Armadillo’s download page](http://arma.sourceforge.net/download.html) to find more information on how to use other implementations of BLAS and LAPACK):

```
sudo apt-get install cmake
```
3. Download and install PASS

```
git clone --depth 1 --branch master https://github.com/SRAhub/PASS.git
cd PASS
cmake .
make
sudo make install
```

## MacOS Installation
The following commands are based on [Homebrew] (https://brew.sh), a package manager for OS X.

1. Install Cmake

```
brew install cmake
```

2. Install GCC

```
brew install gcc
```

3. Install Armadillo

```
brew install armadillo
```

4. Download and install PASS

```
git clone --depth 1 --branch master https://github.com/SRAhub/PASS.git
cd PASS
cmake .
make
sudo make install
```

5. Change default gcc compiler: open .bash_profiles and paste the following code (if you are using GCC 7.*)
```
# Modify system PATH for new GNU compiler install
export PATH="/usr/local/bin:$PATH"

alias gcc='gcc-7'
alias cc='gcc-7'
alias g++='g++-7'
alias c++='c++-7'
```

Implemented Algorithms
-------
- Standard Particle Swarm Optimisation 2011
- Random Searcb
- Parallel Swarm Search (own)

Implemented Problems
-------
- Benchmark Problems
  - Ackley Function
  - Rosenbrock Function
  - Sphere Function
  - Sum of Different Powers

- Astro Problems
  - [GTOC1 C++ implementation by the ESA](http://www.esa.int/gsp/ACT/inf/projects/gtop/gtoc1.html)

- Parallel Kinematic Machine
  - Implements robot simulation models described in [Influence of kinematic redundancy on the singularity-free workspace of parallel kinematic machines.](https://link.springer.com/article/10.1007/s11465-012-0321-8)

License
-------
Distributed under [MIT license](http://opensource.org/licenses/MIT).

Credits
-------
Big thanks goes to [Sebastian Niemann](https://github.com/SebastianNiemann) helping during the developement. 
