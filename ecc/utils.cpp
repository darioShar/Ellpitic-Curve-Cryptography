#include "utils.h"
/*
polynomial<ZpZ> gcd(const polynomial<ZpZ>& a, const polynomial<ZpZ>& b) {
	if (a.degree() == 0 && b.degree() == 0) return polynomial<ZpZ>{ZpZ(1, 5)};
	polynomial<ZpZ> x, y;
	if (a.degree() > b.degree()) {
		x = b;
		y = a;
	}
	else {
		x = a;
		y = b;
	}
	while (x.degree() != 0) {
		auto tmp = x;
		x = y % x;
		y = tmp;
	}
	return y / y[y.degree()];
}
*/

void get_primes(std::vector<mpir_ui>& dest, unsigned int n) {
	std::vector<bool> vec_primes(n, true);
	vec_primes[0] = false;
	vec_primes[1] = false;

	for (auto i = 4; i < n; i += 2) { // even numbers are handled first for efficiency
		vec_primes[i] = false;
	}

	for (auto i = 3; i < n; i += 2) {
		if (vec_primes[i]) {
			//the usual is for(int j = i * 2; j < MAXIMUM; j += i) { but I am using a bit more optimized loop
			for (auto j = i * i; j < n; j += 2 * i) {
				vec_primes[j] = false;
			}
		}
	}
	dest.clear();
	dest.reserve(n / log(n));
	for (auto i = 2; i < n; i++) {
		if (vec_primes[i])
			dest.push_back(i);
	}
}

/*Checks if >=*/
bool is_product_geq(const std::vector<mpir_ui>& src, mpz_t N) {
	mpz_t tmp; mpz_init_set_ui(tmp, 1);
	for (auto& x : src) {
		mpz_mul_ui(tmp, tmp, x);
		if (mpz_cmp(tmp, N) >= 0) {
			mpz_clear(tmp);
			return true;
		}
	}
	mpz_clear(tmp);
	return false;
}

void mpz_set_from_prime_fact(mpz_t dest, prime_factorization primeFact) {
	mpz_set_ui(dest, 0);
	mpz_t op;
	mpz_init(op);
	for (auto& x : primeFact) {
		mpz_pow_ui(op, x.p, x.exp);
		mpz_mul(dest, dest, op);
	}
	mpz_clear(op);
}

/* Sets x to rhs of the weierstrass equation for an ecc : x^3 + ax +b */
void set_weierstrass(mpz_t &dest, const mpz_t &x, const mpz_t &a, const mpz_t &b, const mpz_t &p) {
	mpz_t op1, op2;
	mpz_inits(op1, op2, NULL);
	mpz_mul(op1, x, a);
	mpz_powm_ui(op2, x, 2, p);
	mpz_addmul(op1, op2, x);
	mpz_add(op1, op1, b);
	mpz_mod(op1, op1, p); // op1 = x^3 + ax + b
	mpz_set(dest, op1);
	mpz_clears(op1, op2, NULL);
}


int mpz_sqrtm(mpz_t q, const mpz_t n, const mpz_t p) {
	mpz_t w, n_inv, y;
	mpir_ui i, s;
	//TMP_DECL;
	//TMP_MARK;

	if (mpz_divisible_p(n, p)) {         /* Is n a multiple of p?            */
		mpz_set_ui(q, 0);               /* Yes, then the square root is 0.  */
		return 1;                       /* (special case, not caught        */
	}                                   /* otherwise)                       */
	if (mpz_legendre(n, p) != 1)         /* Not a quadratic residue?         */
		return 0;                       /* No, so return error              */
	if (mpz_tstbit(p, 1) == 1) {         /* p = 3 (mod 4) ?                  */
		mpz_set(q, p);
		mpz_add_ui(q, q, 1);
		mpz_fdiv_q_2exp(q, q, 2);
		mpz_powm(q, n, q, p);           /* q = n ^ ((p+1) / 4) (mod p)      */
		return 1;
	}
	/*
	MPZ_TMP_INIT(y, 2 * SIZ(p));
	MPZ_TMP_INIT(w, 2 * SIZ(p));
	MPZ_TMP_INIT(n_inv, 2 * SIZ(p));
	*/
	mpz_init2(y, mpz_sizeinbase(p, 2));
	mpz_init2(w, mpz_sizeinbase(p, 2));
	mpz_init2(n_inv, mpz_sizeinbase(p, 2));

	mpz_set(q, p);
	mpz_sub_ui(q, q, 1);                /* q = p-1                          */
	s = 0;                              /* Factor out 2^s from q            */
	while (mpz_tstbit(q, s) == 0) s++;
	mpz_fdiv_q_2exp(q, q, s);           /* q = q / 2^s                      */
	mpz_set_ui(w, 2);                   /* Search for a non-residue mod p   */
	while (mpz_legendre(w, p) != -1)     /* by picking the first w such that */ 
		mpz_add_ui(w, w, 1);            /* (w/p) is -1                      */
	mpz_powm(w, w, q, p);               /* w = w^q (mod p)                  */
	mpz_add_ui(q, q, 1);
	mpz_fdiv_q_2exp(q, q, 1);           /* q = (q+1) / 2                    */
	mpz_powm(q, n, q, p);               /* q = n^q (mod p)                  */
	mpz_invert(n_inv, n, p);
	for (;;) {
		mpz_powm_ui(y, q, 2, p);        /* y = q^2 (mod p)                  */
		mpz_mul(y, y, n_inv);
		mpz_mod(y, y, p);               /* y = y * n^-1 (mod p)             */
		i = 0;
		while (mpz_cmp_ui(y, 1) != 0) {
			i++;
			mpz_powm_ui(y, y, 2, p);    /*  y = y ^ 2 (mod p)               */
		}
		if (i == 0) {                    /* q^2 * n^-1 = 1 (mod p), return   */
			mpz_clear(w);
			mpz_clear(n_inv);
			mpz_clear(y);
			return 1;
		}
		if (s - i == 1) {                  /* In case the exponent to w is 1,  */
			mpz_mul(q, q, w);           /* Don't bother exponentiating      */
		}
		else {
			mpz_powm_ui(y, w, (mpir_ui)1 << (s - i - 1), p);
			mpz_mul(q, q, y);
		}
		mpz_mod(q, q, p);               /* r = r * w^(2^(s-i-1)) (mod p)    */
	}

	mpz_clear(w);
	mpz_clear(n_inv);
	mpz_clear(y);
}

/* Sets dest to the discriminant of the curve y^2 = x^3 + ax + b, and j to its j_invariant :
discr = -16*(4a^3 + 27b^2) and j = 1728 * 4 * A^3 / discr*/
void set_weierstrass_discriminant_j_invariant(mpz_t& discr, mpz_t& j, const mpz_t &a, const mpz_t &b, const mpz_t& p) {
	mpz_t op1, op2;
	mpz_init(op2);
	mpz_init_set(op1, b);
	mpz_mul_ui(op1, b, 27);
	mpz_mul(op1, op1, b); // op1 = 27b^2
	mpz_powm_ui(op2, a, 3, p);
	mpz_mul_ui(op2, op2, 4); // op2 = 4 * A^3
	mpz_add(discr, op1, op2);
	mpz_mul_ui(discr, discr, 16);
	mpz_neg(discr, discr);	
	mpz_mod(discr, discr, p); // discr ok
	if (mpz_sgn(discr) == 0) return; // discr is null
	mpz_mul_ui(op2, op2, 1728); // op2 = 1728 * 4 * A^3
	mpz_invert(op1, discr, p);
	mpz_mul(j, op1, op2); 
	mpz_mod(j, j, p);// j ok
	mpz_clears(op1, op2, NULL);
}

/* Returns a^(-1) % b, assuming gcd(a, b) = 1 */
mpir_ui mod_inv(mpir_ui a, mpir_ui b) {
	mpir_ui b0 = b, t, q;
	mpir_ui x0 = 0, x1 = 1;
	if (b == 1) return 1;
	while (a > 1) {
		q = a / b;
		t = b, b = a % b, a = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}

/* Returns x such that x equiv a[i] mod n[i]*/
void chinese_remainder(mpz_t &dest, const std::vector<mpir_ui>& a, const std::vector<mpir_ui>& n) {
	if (a.size() != n.size()) {
		printf("Dimension error in chinese remainder calculation");
		return;
	}
	mpz_t prod; mpz_init_set_ui(prod, 1);
	mpz_t sum; mpz_init_set_ui(sum, 0);
	mpz_t p; mpz_init(p);
	mpz_t tmp_mod; mpz_init(tmp_mod);

	for (size_t i = 0; i < n.size(); i++) 
		mpz_mul_ui(prod, prod, n[i]);

	for (size_t i = 0; i < n.size(); i++) {
		mpz_divexact_ui(p, prod, n[i]);
		// MAYBE PROBLEM HERE WITH mpz_mod_ui : 
		mpz_mod_ui(tmp_mod, p, n[i]);
		mpz_addmul_ui(sum, p, a[i] * mod_inv(mpz_get_ui(tmp_mod), n[i]));
	}

	mpz_mod(dest, sum, prod);
	mpz_clears(prod, sum, p, tmp_mod, NULL);
}

