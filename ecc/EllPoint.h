#pragma once
#include <stdarg.h>
#include <stdio.h>
#include <time.h>
#include <mpir.h>
#include "utils.h"
#include "ectypes.h"


#define MPZ_PRINT(x) mpz_out_str(stdout, 10, x); printf("\n");

/* Elliptic Curve Point*/
class EllPoint {
public:
	/*Constructs an EllPoint*/
	EllPoint(const ECParam& param);
	EllPoint(const EllPoint& p);
	~EllPoint();

	void operator=(const EllPoint& p);
	bool operator==(const EllPoint& p) const;
	bool operator+=(const EllPoint& p);
	void operator*=(mpir_ui n);
	void operator*=(const mpz_t& n);
	void inverse();

	const ECCoord& getCoord() const;
	const ECParam& getECParam() const;

	/*Returns true and sets given coordinates if they correspond to 
	a point on the curve. Returns false otherwise*/
	bool setCoord(const ECCoord& coord);

	/*Sets new coordinates by randomly choosing x until x^3 + ax + b
	is a quadratic residue mod p and setting y accordingly (Shanks-Tonelli)*/
	void setRandomCoord();

	/*Returns true if the point is at infinite (neutral element)*/
	bool isInf() const;

	/*Bruteforces by calculating nP for all n in Fp. Can be expected to be very slow*/
	void order(mpz_t& dest) const;

	/*If curve order prime factorization is known, we can compute point's order much faster
	thanks to Lagrange's theorem */
	void order_co(mpz_t& dest, prime_factorization curveOrderFactorization) const;

	/*Prints the point coordinates plus its name *name* if given*/
	void print(const char* name = "") const;

	/*Internal function to determine whether a point is on the curve or not*/
	bool isOnCurve(const ECCoord& coord) const;

	/*Sets point to neutral element (point at infinity)*/
	void setInf();

private:
	ECParam			m_param;
	mpz_t			m_lc; // leading coefficient, used when adding etc.
	ECCoord			m_coord;
};