# Parallel Adaptive Swarm Search (PASS)

PASS is written in C++14 and uses [Armadillo](http://arma.sourceforge.net/) (developed by Conrad Sanderson et al., NICTA, Australia) for highly efficient linear algebra calculations.

## Prerequirements

- GCC 4.8.1 or later
- CMake 3.2 or later
- Armadillo C++ 8.200.0 or later

## Installation

### Linux Installation

The following commands are based on Ubuntu, using apt-get. If you are using another Linux distributions, the names of each package might also differ.

If you are using any recent Linux distribution, you should already be equipped with a C++11 feature complete compiler, like GCC 4.8.1+

- Install Cmake

```bash
sudo apt-get install cmake
```

- Install Armadillo C++ with OpenBLAS support (visit [Armadillo’s download page](http://arma.sourceforge.net/download.html) to find more information on how to use other implementations of BLAS and LAPACK):

```bash
sudo apt-get install cmake
```

- Download and install PASS

```bash
git clone --depth 1 --branch master https://github.com/SRAhub/PASS.git
cd PASS
cmake .
make
sudo make install
```

### MacOS Installation

The following commands are based on [Homebrew](https://brew.sh), a package manager for OS X.

- Install Cmake

```bash
brew install cmake
```

- Install GCC

```bash
brew install gcc
```

- Install Armadillo

```bash
brew install armadillo
```

- Download and install PASS

```bash
git clone --depth 1 --branch master https://github.com/SRAhub/PASS.git
cd PASS
cmake .
make
sudo make install
```

- If you don't want to use the MacOS CLANG compiler, you have to change it to the gcc compiler: open *.bash_profiles* and paste the following code (if you are using GCC 7.*)

```bash
# Modify system PATH for new GNU compiler install
export PATH="/usr/local/bin:$PATH"

alias gcc='gcc-7'
alias cc='gcc-7'
alias g++='g++-7'
alias c++='c++-7'
```

## Implemented Algorithms

- Hooke Jeeves Algorithm
- Random Search
- Parallel Swarm Search (own)
- Standard Particle Swarm Optimisation 2011

## Implemented Problems

- Optimisation Benchmark Problems
  - Ackley Function
  - De Jong's Function
  - Rastrigin Function
  - Rosenbrock Function
  - Styblinski-Tang Function
  - Sum of Different Powers Function

- Space Mission Problems
  - [GTOC1](http://www.esa.int/gsp/ACT/inf/projects/gtop/gtoc1.html)
  - [Messenger (Full Version)](http://www.esa.int/gsp/ACT/inf/projects/gtop/messenger_full.html)

## License

Distributed under [MIT license](http://opensource.org/licenses/MIT).

## Credits

Big thanks goes to [Sebastian Niemann](https://github.com/SebastianNiemann) helping during the development.
