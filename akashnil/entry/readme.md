## Chia network VDF competition track 1 submission

Author: Akashnil Dutta

The install.sh script is run by the server to install any dependencies, and/or compile the code. Here we are installing the GMP library as a dependency.

The run.sh file is what executed to run the VDF. It takes two arguments:
* A discriminant in hex
* The number of iterations, in decimal

```
sh ./install.sh
sh ./run.sh -0xdc2a335cd2b355c99d3d8d92850122b3d8fe20d0f5360e7aaaecb448960d57bcddfee12a229bbd8d370feda5a17466fc725158ebb78a2a7d37d0a226d89b54434db9c3be9a9bb6ba2c2cd079221d873a17933ceb81a37b0665b9b7e247e8df66bdd45eb15ada12326db01e26c861adf0233666c01dec92bbb547df7369aed3b1fbdff867cfc670511cc270964fbd98e5c55fbe0947ac2b9803acbfd935f3abb8d9be6f938aa4b4cc6203f53c928a979a2f18a1ff501b2587a93e95a428a107545e451f0ac6c7f520a7e99bf77336b1659a2cb3dd1b60e0c6fcfffc05f74cfa763a1d0af7de9994b6e35a9682c4543ae991b3a39839230ef84dae63e88d90f457 10000
```

The script should output the result of the VDF.

## Description

This is an improvement from the sample submission provided by Chia network. The overall squaring algorithm is the same one from the Binary Quadratic Forms paper [1].

When the sample entry is profiled, it is clear that 90% of the time is taken in the reduction steps, about 10% in the GCD computation. The rest of the time is insignificant. For 2048 bit discriminants, the reduction loop runs about about 200 times. In each iteration, large integer division / remainder operations are performed which take time. What I have done is approximate the large integers using small int64s. Then perform the reduction steps on these small approximants.

During reduction on int64s, we can update the matrix by which (x, y) is multiplied in the quadratic form (a x2 + b xy + c y2). The matrix is (u, v, w, x) in the code and (alpha, beta, gamma, delta) in [2]. From these values we can find the matrix by which (a, b, c) is multiplied (See formulas in [2]). The matrix is M = (aa, ab, ac, ba, bb, bc, ca, cb, cc) in the code. Ater that we apply M * (a, b, c) in each stage of the reduction. The size of (a, b, c) is reduced significantly after each reduction. In second stage we approximate the new (a, b, c) again following the same process. This way we can reduce the total number of stages to 20-30 from 200. The whole algorithm improves the reduction by almost 5x.

We were unable to improve extended GCD calculation of GMP which currently takes up 35% time in the final code. Reduction now takes up 60% time of the overall time budget.

[1] https://github.com/Chia-Network/vdf-competition/blob/master/classgroups.pdf

[2] http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.590.2666&rep=rep1&type=pdf
	(Section 2 : Representation and Equivalence)