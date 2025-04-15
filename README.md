# Fast-Amortized-Bootstrapping


## Build

For processors with `AVX-512` and `VAES`:

``` make [parameters]```

For processors with only `AVX2/FMA`:

``` make FFT_LIB=spqlios A_PRNG=none ENABLE_VAES=false [parameters] ```

> [!WARNING]
> The implementation will be slower without [AVX-512](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions). Results in the paper are for an [`r7i.metal-24xl`](https://instances.vantage.sh/aws/ec2/r7i.metal-24xl) instance on AWS.


For other compiling options, see [MOSFHET](https://github.com/antoniocgj/MOSFHET). 

## Available parameters

For each of the parameters below, SET_X_Y refers to the bootstrapping of X-bit messages if function is abitrary or Y-bit messages if the function is negacyclic (See Remark 7.1.3 in the paper).

> [!WARNING]
> Don't forget to add `FFT_LIB=spqlios A_PRNG=none ENABLE_VAES=false` to the `make` command if you don't have AVX512.

### Binary keys
``` make -B PARAM=SET_2_3 && ./main``` 

``` make -B PARAM=SET_4_5 && ./main``` 

``` make -B PARAM=SET_6_7 && ./main``` 

``` make -B PARAM=SET_8_9 && ./main``` 

``` make -B PARAM=SET_8_9_b && ./main``` 

Result with very high Failure Probability:

``` make -B PARAM=SET_8_9_HIGH_FR && ./main``` 

### Ternary keys

``` make -B KEY=TERNARY PARAM=SET_2_3 && ./main``` 

``` make -B KEY=TERNARY PARAM=SET_4_5 && ./main``` 

``` make -B KEY=TERNARY PARAM=SET_6_7 && ./main``` 

``` make -B KEY=TERNARY PARAM=SET_8_9 && ./main``` 

Result with very high Failure Probability:

``` make -B KEY=TERNARY PARAM=SET_8_9_HIGH_FR && ./main``` 

### Arbitrary keys

``` make -B KEY=ARBITRARY PARAM=SET_A2 && ./main``` 

``` make -B KEY=ARBITRARY PARAM=SET_A4 && ./main``` 

``` make -B KEY=ARBITRARY PARAM=SET_A6 && ./main``` 

### Larger parameters

``` make -B KEY=TERNARY PARAM=SET_L0 && ./main``` 

``` make -B KEY=TERNARY PARAM=SET_L1 && ./main```

``` make -B KEY=TERNARY PARAM=SET_L2 && ./main```

``` make -B KEY=TERNARY PARAM=SET_L3 && ./main```

``` make -B KEY=ARBITRARY PARAM=SET_LA1 && ./main```

``` make -B KEY=ARBITRARY PARAM=SET_LA2 && ./main```

``` make -B KEY=ARBITRARY PARAM=SET_LA3 && ./main```

``` make -B KEY=ARBITRARY PARAM=SET_LA4 && ./main```

## Measuring noise

To measure noise, uncomment line 278 of `src/sparse_amortized_bootstrapping.c`.

Performance measurement are unreliable and correctness checks may fail when measuring noise.
