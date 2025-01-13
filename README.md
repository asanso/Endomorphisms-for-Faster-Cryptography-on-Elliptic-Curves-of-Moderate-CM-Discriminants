# Benchmark Results Overview

This document presents benchmark results for scalar multiplication optimizations on three elliptic curves: **Lollipop489201**, **Lollipop574261**, and **GC256C**. The benchmarks compare the performance of standard scalar multiplication against optimized methods, such as **endomorphism** and **GLV (Gallant–Lambert–Vanstone)**.

For more details, refer to **Section 3** of the paper.


## Usage

```shell
cd python
sage -python bench.py 
================================================================================
                                 Lollipop489201                                 
================================================================================
scalar multiplication:                  0.703ms
endomorphism:           0.334ms (111% faster)

scalar multiplication:                  2.029ms
GLV:            1.551ms (31% faster)
================================================================================
                                 Lollipop574261                                 
================================================================================
scalar multiplication:                  0.960ms
endomorphism:           0.373ms (157% faster)

scalar multiplication:                  2.903ms
GLV:            2.064ms (41% faster)
================================================================================
                                     GC256C                                     
================================================================================
scalar multiplication:                  0.898ms
endomorphism:           0.330ms (172% faster)

scalar multiplication:                  2.710ms
GLV:            1.929ms (41% faster)
```
## How to Interpret Results

The benchmark results provide valuable insights into the efficiency of scalar multiplication optimizations. Below are guidelines to help interpret the results:

1. **Execution Time**
  * The results are displayed in milliseconds (ms), indicating the time taken for each operation.
  * Lower values signify better performance.
2. **Percentage Improvement**
 * The percentage improvements (e.g., "111% faster") demonstrate how much more efficient the optimized methods (endomorphism and GLV) are compared to standard scalar multiplication.
3. **Comparative Context**
 * In the first line of each curve’s results, the endomorphism $\phi$ optimization is compared against the time taken to double a point l/2 times, as indicated in the third column of Table 1 for the corresponding curve.
 * The second entry compares standard Double-and-Add scalar multiplication against scalar multiplication combined with the GLV optimization based on $\phi$.
 
