#pragma once
#include <stdarg.h>
#include <stdio.h>
#include <type_traits>
#include <time.h>
#include <utility>
#include <mpir.h>


#define DEFAULT_A 0
#define DEFAULT_B 0
#define DEFAULT_P 5
#define DEFAULT_ORDER 0

#define PT_INFTY 1
#define NOT_INFTY 0

void extern mpz_print(const mpz_t print);

/* Parameters of an elliptic curve in its reduced Weierstrass form
y^2 = x^3 + ax + b*/
struct ECParam {
	ECParam(mpir_ui _a = DEFAULT_A, mpir_ui _b = DEFAULT_B, mpir_ui _p = DEFAULT_P, mpir_ui _order = DEFAULT_ORDER);
	ECParam(const char* _a, const char* _b, const char* _p, const char* _order);
	ECParam(mpz_t _a, mpz_t _b, mpz_t _p, mpz_t _order = NULL);

	~ECParam();
	ECParam(const ECParam& p);
	void operator=(const ECParam& p);
	bool operator==(const ECParam& p) const;
	mpz_t a;
	mpz_t b;
	mpz_t p; // Z / pZ
	mpz_t order;
};

/* Coordinate of a point on an ec*/
struct ECCoord {
	ECCoord(int isInf = PT_INFTY, const char* x = "0", const char* y = "0");
	ECCoord(int _isInf, mpz_t _x, mpz_t _y);
	~ECCoord();
	ECCoord(const ECCoord& p);
	void operator=(const ECCoord& p);
	bool operator==(const ECCoord& p) const;
	int isInf;
	mpz_t x;
	mpz_t y;
};


/* Helper struct for hash function */
struct ECPair {
	const ECCoord coord;
	const mpir_ui val;
	ECPair(const ECCoord& c, mpir_ui v) : coord(c), val(v) {}
	bool operator==(const ECPair& ecpair) const {
		return coord == ecpair.coord;
	}
};

/* (ECCoord, mpz_t) Hash Function */
struct BSGSHasher {
	mpir_ui operator()(const ECPair& elem) const
	{
		using std::hash;

		return ((hash<mpir_ui>()(mpz_get_ui(elem.coord.x))
				^ (hash<mpir_ui>()(mpz_get_ui(elem.coord.y)) << 1)) >> 1);
	}
};

