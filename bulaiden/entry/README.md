## Description

This implementation is written in C/C++ and it uses the GMP and FLINT library for arithmetic and GCD.

It uses NUDUPL to square a form, and Chia's reference implementation to reduce a form.
A NUDUPL implementation is described in the paper Computational Aspects of NUCOMP
https://www.researchgate.net/publication/221451638_Computational_aspects_of_NUCOMP

The initial element is (2, 1, c), where c is calculated using the discriminant passed in. A form is represented by (a, b, c) and the discriminant. Here, reduction is performed after every composition/multiplication.