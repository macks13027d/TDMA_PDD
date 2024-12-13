Installation                       
============

## Downloads
The repository can be cloned as follows:

```
git clone https://github.com/MPMC-Lab/TDMA_PDD.git
```
Alternatively, the source files can be downloaded through github menu 'Download ZIP'.




## Compile
### Prerequisites
Prerequisites to compile TDMA_PDD are as follows:
* MPI (IntelMPI)
* Fortran compiler (Intel Fortran Compiler)


### Compile and build
#  ![](https://img.shields.io/badge/Tested-IntelMPI_21.2_and_Intel_Fortran_Compiler_21.2-silver.svg?logo=cachet)
* Build TDMA_PDD
    ```
   make lib
   ```
* Build an example problem after build TDMA_PDD

    ```
   make example
   ```
* Build all

    ```
   make all
    ```

### Running the example
After building the example file, an executable binary, `tdma_pdd`, is built in the `run` folder. The `input.dat` file in the `run` folder is a pre-defined input file, and the `tdma_pdd` can be executed as follows:
    ```
   mpirun -np 8 ./tdma_pdd_cpu >log &
    ```
Output is the error of PDD method and error would be machine error if you choose appropriated number of grid and processors. Note that in this example, diagonal dominance was forced by setting a=0.1*rand(), b=1, and c=0.1*rand(), where a, b, and c are entries in the tridiagonal matrix.


## Folder structure
* `src` : source files of TDMA_PDD
* `example` : validation and performance checking of TDMA_PDD
* `include` : contains header files created after building and required codes
* `lib` : a static library of TDMA_PDD is created after building
and header files of required codes are contained.
* `doc` : documentation
* `run` : an executable binary file for the example problem is created after building

<div class="section_buttons">

| Previous          |                              Next |
|:------------------|----------------------------------:|
| [Introduction](1_introduction.md) | [Performance](3_performance.md) |
</div>
