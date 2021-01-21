#pragma once
#include <stdarg.h>
#include <stdio.h>
#include <time.h>
#include <mpir.h>
#include <vector>
#include <algorithm>
#include <mutex>

//#include <boost/math/tools/polynomial.hpp>
//#include "ZpZ.h"

//using namespace boost::math;
//using namespace boost::math::tools; // for polynomial
//using boost::lexical_cast;

//polynomial<ZpZ> gcd(const polynomial<ZpZ>& a, const polynomial<ZpZ>& b);

static class RandomClass{
public :
	RandomClass() { gmp_randinit_default(m_rand_state);}
	~RandomClass() { gmp_randclear(m_rand_state); }

	bool is_likely_prime(const mpz_t& prime)			{ return mpz_likely_prime_p(prime, m_rand_state, 0); }
	void next_prime_candidate(mpz_t& dest, mpz_t src)	{ mpz_next_prime_candidate(dest, src, m_rand_state); }
	void randomb(mpz_t& dest, mpir_ui bitcnt)			{ mpz_urandomb(dest, m_rand_state, bitcnt); }
	void randomm(mpz_t& dest, const mpz_t& mod)				{ mpz_urandomm(dest, m_rand_state, mod); }
private :
	gmp_randstate_t m_rand_state;
} Rand;


void get_primes(std::vector<mpir_ui>& dest, unsigned int n);

bool is_product_geq(const std::vector<mpir_ui>& src, mpz_t N);

struct PrimeFactComp {
	mpz_t p;
	mpir_ui exp;
};

typedef std::vector<PrimeFactComp> prime_factorization;

/* Sets dest to the number described by primeFact under its prime factorization*/
void mpz_set_from_prime_fact(mpz_t dest, prime_factorization primeFact);


/* Solve the modular equatioon x^2 = n (mod p) using the Shanks-Tonelli
* algorihm. x will be placed in q and 1 returned if the algorithm is
* successful. Otherwise 0 is returned (currently in case n is not a quadratic
* residue mod p). A check is done if p = 3 (mod 4), in which case the root is
* calculated as n ^ ((p+1) / 4) (mod p).
*
* Note that currently mpz_legendre is called to make sure that n really is a
* quadratic residue. The check can be skipped, at the price of going into an
* eternal loop if called with a non-residue.
*/
int mpz_sqrtm(mpz_t q, const mpz_t n, const mpz_t p);

/* Sets dest to rhs of the weierstrass equation for an ecc : x^3 + ax +b */
void set_weierstrass(mpz_t &dest, const mpz_t &x, const mpz_t &a, const mpz_t &b, const mpz_t &p);

/* Sets dest to the discriminant of the curve y^2 = x^3 + ax + b, and j to its j_invariant :
discr = 4a^3 + 27b^2 and j = -1728 * 64 * A^3 / discr*/
void set_weierstrass_discriminant_j_invariant(mpz_t& dest, mpz_t& j,const mpz_t &a, const mpz_t &b, const mpz_t& p);

/* Returns a^(-1) % b, assuming gcd(a, b) = 1 */ 
mpir_ui mod_inv(mpir_ui a, mpir_ui b);

/* Returns x such that x equiv a[i] mod n[i]*/
void chinese_remainder(mpz_t& dest, const std::vector<mpir_ui>& a, const std::vector<mpir_ui>& n);