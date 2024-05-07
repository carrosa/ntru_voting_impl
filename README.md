# More Efficient Post-Quantum Electronic Voting from NTRU 

Code accompannying the paper "More Efficient Post-Quantum Electronic Voting from NTRU" (Submission #29) for SAC 2024. 

Depedencies are the [NFLlib](https://github.com/quarkslab/NFLlib) and [FLINT](https://flintlib.org/doc/) 2.9 libraries.
NFLLib is already included in this repository, but instructions for installing its dependencies can be found in the link above.
FLINT is usually included in package managers and can be easily installed in most systems out there.

### Building dependencies

To build NFLLib, run the following inside a cloned version of this repository:

```
$ mkdir deps
$ cd deps
$ cmake ../NFLlib -DCMAKE_BUILD_TYPE=Release -DNFL_OPTIMIZED=ON
$ make
$ make test
```
### Building and running the code

For building the actual code, run `make` inside the source directory. This will build the binaries for `ntru_bdlop`, `ntru`, `ntru_pismall` and `ntru_shuffle` to test and benchmark different modules of the code.

The binaries respectively implement the BDLOP commitment scheme, the distributed NTRU cryptosystem, the zero-knowledge proof and the shuffle itself. Tests and benchmarks are included for each of them, such that they can be used independently. NFLlib is quite memory-hungry due to being a template library, so we recommend to adjust the stack size with `ulimit -s unlimited` to avoid crashing in the largest benchmarks.

__WARNING__: This is an academic proof of concept, and in particular has not received code review. This implementation is NOT ready for any type of production use.
