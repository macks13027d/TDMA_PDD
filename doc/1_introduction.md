# TDMA_PDD: Multi-dimensional Diagonal Dominant Tridiagonal Matrix Algorithm Library
![](https://img.shields.io/badge/Fortran-Fortran_90-blue.svg)
[![](https://img.shields.io/badge/docs-passing-green.svg)](https://xccels.github.io/TDMA_PDD)
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

# Cite
Please use the following bibtex, when you refer to this project.

```bibtex
@article{moon2020application,
  title={Application of the parallel diagonal dominant algorithm for the incompressible Navier-Stokes equations},
  author={Moon, Hojun and Hong, Seungpyo and You, Donghyun},
  journal={Journal of Computational Physics},
  volume={423},
  pages={109795},
  year={2020},
  publisher={Elsevier}
}

<div class="section_buttons">

|                        Read Next |
|---------------------------------:|
| [Installation](2_installation.md) |

</div>
