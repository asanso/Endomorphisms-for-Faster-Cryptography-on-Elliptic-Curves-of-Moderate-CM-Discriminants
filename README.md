# Benchmark Results Overview

This document presents benchmark results for scalar multiplication optimizations on three elliptic curves: **Lollipop489201**, **Lollipop574261**, and **GC256C**. The benchmarks compare the performance of standard scalar multiplication against optimized methods, such as **endomorphism** and **GLV (Gallant–Lambert–Vanstone)**.

For more details, refer to **Section 3** of the paper.


### Usage

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
