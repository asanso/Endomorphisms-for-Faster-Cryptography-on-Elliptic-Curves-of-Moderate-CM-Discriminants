# Benchmark Results Overview

This document presents benchmark results for scalar multiplication on three elliptic curves: **GC256C (Russian curve)**, **Lollipop-489-201**, and **Lollipop-574-261**. The benchmarks compare the performance of the optimized **GLV (Gallant–Lambert–Vanstone)** method against the standard **double-and-add** scalar multiplication. For more details, refer to the paper and especially to its Section 3.


## Usage

```shell
cd python
sage -python bench.py 
================================================================================
                                 GC256C (Russian curve)                                    
================================================================================
scalar multiplication by 2^l':    0.898ms
endomorphism ϕ:                   0.330ms (172% faster)

standard scalar multiplication:   2.710ms
GLV based on ϕ:                   1.929ms (41% faster)
================================================================================
                                 Lollipop-489-201                                 
================================================================================
scalar multiplication by 2^l':    0.703ms
endomorphism ϕ:                   0.334ms (111% faster)

standard scalar multiplication:   2.029ms
GLV based on ϕ:                   1.551ms (31% faster)
================================================================================
                                 Lollipop-574-261                                 
================================================================================
scalar multiplication by 2^l':    0.960ms
endomorphism ϕ:                   0.373ms (157% faster)

standard scalar multiplication:   2.903ms
GLV based on ϕ:                   2.064ms (41% faster)
```

**N.B.** The benchmarks are averaged over 1000 instances of the routine. Besides, the endomorphism $\phi$ is not evaluated via (the homogeneous version of) Horner's scheme as proposed in the paper, but more elementarily. Therefore, (GLV based on) $\phi$ is even faster when implemented more properly.

## How to Interpret Results

The benchmark results provide valuable insights into the efficiency of scalar multiplication optimizations. Below are guidelines to help interpret the results:

1. **Execution Time**
  * The results are displayed in milliseconds (ms), indicating the time taken for each operation.
  * Lower values signify better performance.
2. **Percentage Improvement**
 * The percentage improvements demonstrate how much more efficient the optimized solutions (the endomorphism $\phi$ and GLV based on $\phi$) are compared to the double(-and-add) scalar multiplication.
3. **Comparative Context**
 * In the first two lines of each curve’s results, the time of evaluating the endomorphism $\phi$ is compared against that taken to double a point $\ell^\prime = \lceil \ell/2 \rceil$ times, where the value $\ell$ is indicated in the third column of Table 1 from the paper.
 * The second entry compares the standard scalar multiplication against the GLV optimization using $\phi$.
 
