# Ellpitic-Curve-Cryptography
Simple C++ code for ecc. 

The code uses MPIR (forked from the GMP bignum library). It is not intended to be extremely fast, but rather for educational and learning purposes.
Basic operations on E(F_q) are implemented, with q a prime >= 5, and the curve being described in Weierstrass reduced form.
Cardinality of a curve can be determined with a O(q^(1/4)) algorithm with a fine use of baby step giant step. With a few GB of memory 
you work on curves of order approx. 2^95.
Implemented polynomials on finite field with the goal of implementing schoof algorithm, which I had not enough time to do.
You can find a pdf explaining the theory behind all this and a powerpoint presentation which is more centered around the actual implementation (written in french).
