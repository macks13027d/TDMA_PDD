# TDMA_PDD_CPU: Multi-dimensional Diagonal Dominant Tridiagonal Matrix Algorithm Library
![](https://img.shields.io/badge/Fortran-Fortran_90-blue.svg)
[![](https://img.shields.io/badge/docs-passing-green.svg)](https://xccels.github.io/TDMA_PDD_CPU)
![](https://img.shields.io/badge/license-MIT_License-yellow.svg)
This library provides efficient parallel computation for multidimensional TDMA.

The parallel multi-dimensional TDMA for data of size (nx, ny, nz) using npx*npy*npz processors with this library proceeds through the following steps:

1.	Allocation: The tri-diagonal matrix and the RHS (Right-hand side) vector are allocated to each processor in sizes of (nx/npx, ny/npy, nz/npz) using MPI_Send and MPI_Recv.
2.	Transpose: Transpose the data so that the index of the dimension to be solved (e.g., y-direction) comes first.
(nx/npx, ny/npy, nz/npz) -> (ny/npy, nx/npx, nz/npz)
3.	Computation: Perform TDMA using the PDD algorithm.

(Modules corresponding to steps (1), (2), and (3))

The PDD algorithm, proposed by Sun et al. (1989), is a method for solving TDMA of diagonally dominant tridiagonal matrices. This library applies this method to CPU multi-dimensional TDMA, reducing data communication between processors, thereby improving computation speed and achieving high scalability.
 
# Authors
- Hojun Moon (mhj2013@postech.ac.kr), Mechanical Engineering, POSTECH
- Seungchan Kim (macks1029@postech.ac.kr), Mechanical Engineering, POSTECH
- Jihoo Kim (hugh577@postech.ac.kr), Mechanical Engineering, POSTECH
- Donghyun You (dhyou@postech.ac.kr), Mechanical Engineering, POSTECH

# Usage
## Downloading TDMA_PDD_CPU
The repository can be cloned as follows:

```
git clone https://github.com/MPMC-Lab/TDMA_PDD_CPU.git
```
Alternatively, the source files can be downloaded through github menu 'Download ZIP'.

## Compile
### [Prerequisites](./doc/2_installation.md)
Prerequisites to compile TDMA_PDD_CPU are as follows:
* MPI (MVAPICH 2-2.3.6)
* Fortran compiler (Intel Fortran Compiler 18.0.3)

### Compile and build
* Build TDMA_PDD_CPU
    ```
	make lib
	```
* Build an example problem after build TDMA_PDD_CPU

    ```
	make example
	```
* Build all

    ```
	make all
	```

## Running the example
After building the example file, an executable binary, `tdma_pdd_cpu`, is built in the `run` folder. The `input.dat` file in the `run` folder is a pre-defined input file, and the `tdma_pdd_cpu` can be executed as follows:
    ```
	mpirun -np 8 ./tdma_pdd_cpu >log &
    ```
Output is the error of PDD method and error would be machine error if you choose appropriated number of grid and processors. Note that in this example, diagonal dominance was forced by setting a=0.1*rand(), b=1, and c=0.1*rand(), where a, b, and c are entries in the tridiagonal matrix.

# Folder structure
* `src` : source files of TDMA_PDD_CPU
* `example` : validation and performance checking of TDMA_PDD_CPU
* `include` : contains header files created after building and required codes
* `lib` : a static library of TDMA_PDD_CPU is created after building
and header files of required codes are contained.
* `doc` : documentation
* `run` : an executable binary file for the example problem is created after building

# Cite
Please use the following bibtex, when you refer to this project.


@article{moon2020application,
  title={Application of the parallel diagonal dominant algorithm for the incompressible Navier-Stokes equations},
  author={Moon, Hojun and Hong, Seungpyo and You, Donghyun},
  journal={Journal of Computational Physics},
  volume={423},
  pages={109795},
  year={2020},
  publisher={Elsevier}
}

