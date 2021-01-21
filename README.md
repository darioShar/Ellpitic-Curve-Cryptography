# Ellpitic-Curve-Cryptography
Simple C++ code for ecc. 

The code uses MPIR (forked from the GMP bignum library) for number manipulation, and Vtk (Visual ToolKit) for displaying images. It is not intended to be extremely fast, but rather made for educational and learning purposes.
Basic operations on E(F_q) are implemented, with q a prime >= 5, and the curve is being described in Weierstrass reduced form.
Cardinality of a curve can be determined with a O(q^(1/4)) algorithm thanks to a fine use of baby step giant step. With a few GB of memory 
you can determine curves order with q approx 2^95.

Implemented polynomials on finite field with the goal of  writing schoof's algorithm, which I had not enough time to finish.
You can find a pdf explaining the theory behind all this and a powerpoint presentation which is more centered around the actual implementation (written in french).
