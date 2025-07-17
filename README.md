# Fast-Amortized-Bootstrapping

## Paper

[eprint 2025/686](https://eprint.iacr.org/2025/686)

```
@misc{guimaraes_fast_2025,
 author = {Guimarães, Antonio and Pereira, Hilder V. L.},
 keywords = {Amortized bootstrapping, Bootstrapping, Fully homomorphic encryption, Lattice-based Cryptography},
 note = {Publication info: Preprint - To appear at CCS 2025},
 title = {Fast amortized bootstrapping with small keys and polynomial noise overhead},
 url = {https://eprint.iacr.org/2025/686},
 urldate = {2025-04-19},
 year = {2025}
}


```


## Build

For processors with `AVX-512` and `VAES`:

``` make [parameters]```

For processors with only `AVX2/FMA`:

``` make FFT_LIB=spqlios A_PRNG=none ENABLE_VAES=false [parameters] ```

> [!WARNING]
> The implementation will be slower without [AVX-512](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions). Results in the paper are for an [`r7i.metal-24xl`](https://instances.vantage.sh/aws/ec2/r7i.metal-24xl) instance on AWS.


For other compiling options, see [MOSFHET](https://github.com/antoniocgj/MOSFHET). 

## Available parameters

For each of the parameters below, SET_X_Y_Z refers to the bootstrapping of X-bit messages if the function is arbitrary or Y-bit messages if the function is negacyclic (See Remark 7.1.3 in the paper) for Z messages.

> [!WARNING]
> Don't forget to add `FFT_LIB=spqlios A_PRNG=none ENABLE_VAES=false` to the `make` command if you don't have AVX512.

### Binary keys
``` make -B PARAM=SET_2_3_2048 && ./main``` 

``` make -B PARAM=SET_2_3_4096 && ./main``` 

``` make -B PARAM=SET_2_3_8192 && ./main``` 

``` make -B PARAM=SET_4_5_2048 && ./main``` 

``` make -B PARAM=SET_4_5_4096 && ./main``` 

``` make -B PARAM=SET_4_5_8192 && ./main``` 

``` make -B PARAM=SET_6_7_4096 && ./main``` 

``` make -B PARAM=SET_6_7_8192 && ./main``` 

``` make -B PARAM=SET_8_9_4096 && ./main``` 

``` make -B PARAM=SET_8_9_8192 && ./main``` 

Result with very high Failure Probability:

``` make -B PARAM=SET_8_9_HIGH_FR && ./main``` 

### Ternary keys

``` make -B KEY=TERNARY PARAM=SET_2_3_2048 && ./main``` 

``` make -B KEY=TERNARY PARAM=SET_2_3_4096 && ./main``` 

``` make -B KEY=TERNARY PARAM=SET_2_3_8192 && ./main``` 

``` make -B KEY=TERNARY PARAM=SET_4_5_2048 && ./main``` 

``` make -B KEY=TERNARY PARAM=SET_4_5_4096 && ./main``` 

``` make -B KEY=TERNARY PARAM=SET_4_5_8192 && ./main``` 

``` make -B KEY=TERNARY PARAM=SET_6_7_4096 && ./main``` 

``` make -B KEY=TERNARY PARAM=SET_6_7_8192 && ./main``` 

``` make -B KEY=TERNARY PARAM=SET_8_9_4096 && ./main``` 

``` make -B KEY=TERNARY PARAM=SET_8_9_8192 && ./main``` 


### Arbitrary keys

``` make -B KEY=ARBITRARY PARAM=SET_A2 && ./main``` 

``` make -B KEY=ARBITRARY PARAM=SET_A3 && ./main``` 

``` make -B KEY=ARBITRARY PARAM=SET_A4 && ./main``` 

``` make -B KEY=ARBITRARY PARAM=SET_A5 && ./main``` 

## Measuring noise

To measure noise, uncomment line 101 of `src/sparse_amortized_bootstrapping.c`.

Performance measurements are unreliable, and correctness checks may fail when measuring noise.

## License

[Apache License Version 2.0](LICENSE)

This repository contains code from:

- [MOSFHET](https://github.com/antoniocgj/MOSFHET): [Apache License Version 2.0](https://github.com/antoniocgj/MOSFHET/blob/main/LICENSE) - Copyright Antonio Guimarães et al. - See their [detailed copyright information](https://github.com/antoniocgj/MOSFHET/tree/main?tab=readme-ov-file#license).
