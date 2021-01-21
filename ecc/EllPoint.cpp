#include "EllPoint.h"


EllPoint::EllPoint(const ECParam& param)
	:	m_param(param)
{
	mpz_init2(m_lc, mpz_sizeinbase(m_param.p, 2)); // Exepcting same size as p
}

EllPoint::EllPoint(const EllPoint& p)
	:	m_param(p.m_param)
	,	m_coord(p.m_coord)
{
	mpz_init2(m_lc, mpz_sizeinbase(m_param.p, 2)); // Exepcting same size as p
}

EllPoint::~EllPoint()
{
	mpz_clear(m_lc);
}

void EllPoint::operator=(const EllPoint& p) {
	m_param = p.m_param;
	m_coord = p.m_coord;
}

bool EllPoint::operator==(const EllPoint& p) const {
	return (m_coord == p.m_coord);
}

void EllPoint::inverse() {
	mpz_neg(m_coord.y, m_coord.y);
}


bool EllPoint::operator+=(const EllPoint& p) {
	// Adding 0
	if (p.m_coord.isInf)
		return true;
	if (m_coord.isInf) {
		*this = p;
		return true;
	}
	//Else
	if (mpz_cmp(m_coord.x, p.m_coord.x)) {	//mpz_cmp returns sign of op1 - op2
		mpz_t op1, op2; // temp variables
		mpz_sub(m_lc, p.m_coord.y, m_coord.y);		// m = y2 - y1
		mpz_init_set(op1, m_lc); // op1 = m
		mpz_sub(m_lc, p.m_coord.x, m_coord.x);		// m = x2 - x1
		if(mpz_invert(m_lc, m_lc, m_param.p) == 0)	// m = 1 / m [p]
			return false; // inversion failed
		mpz_mul(m_lc, m_lc, op1);					// m = m * op1
		// m_lc is now set to (y2 - y1) / (x2 - x1) ;
		mpz_powm_ui(op1, m_lc, 2, m_param.p);		// op1 = m^2 [p]
		mpz_sub(op1, op1, p.m_coord.x);				// op1 = op1 - x2
		mpz_sub(op1, op1, m_coord.x);				// op1 = op1 - x1
		// x3 is now set in op1 = m^2 - x1 - x2;
		mpz_init_set(op2, m_coord.x);				// op2 = x1
		mpz_sub(op2, op2, op1);						// op2 = op2 - op1 = x1 - x3
		mpz_mul(op2, op2, m_lc);					// op2 = op2 * m
		mpz_sub(m_coord.y, op2, m_coord.y);			// y3 = op2 - y1

		mpz_mod(m_coord.x, op1, m_param.p);			// x3 = op1 % p
		mpz_mod(m_coord.y, m_coord.y, m_param.p);	// y3 = y3 % p
		//Done
		mpz_clear(op1);
		mpz_clear(op2);
	}
	else {
		if (mpz_cmp(m_coord.y, p.m_coord.y)) {
			m_coord.isInf = PT_INFTY;
		}
		else {
			// POURQOI CA NE MARCHE BORDEL DE BORDEL DE PAS ?????!!!!!
			if (mpz_sgn(m_coord.y)) { // != 0
				mpz_t op1, op2; // temp variables
				mpz_init_set(op1, m_coord.x);			// op1 = x1
				mpz_init_set(op2, m_coord.x);			// op2 = x1
				mpz_powm_ui(op1, op1, 2, m_param.p);	// op1 = op1^2
				mpz_mul_ui(op1, op1, 3);				// op1 *= 3
				mpz_add(op1, op1, m_param.a);			// op1 += a
				mpz_mod(op1, op1, m_param.p);			// op1 = op1%p
				// op1 = 3 * x1 * x1 + a

				mpz_add(op2, m_coord.y, m_coord.y);		// op2 = 2* y1
				if(mpz_invert(op2, op2, m_param.p) == 0)// op2 = 1/op2 % p
					return false; // inversion failed
				mpz_mul(m_lc, op1, op2);				// m = op1 * op2
				mpz_mod(m_lc, m_lc, m_param.p);			// m = m % p
				//m_lc now set
				mpz_powm_ui(op1, m_lc, 2, m_param.p);	// op1 = m^2 % p
				mpz_sub(op1, op1, m_coord.x); 
				mpz_sub(op1, op1, m_coord.x);			// op1 = op1 - 2*x1
				// x3 now in op1
				mpz_sub(op2, m_coord.x, op1);			// op2 = x1 - x3
				mpz_mul(op2, op2, m_lc);				// op2 = m * op2
				mpz_sub(op2, op2, m_coord.y);			// y3 = op2 - y1

				mpz_mod(m_coord.x, op1, m_param.p);		// x3 = x3 % p
				mpz_mod(m_coord.y, op2, m_param.p);		// y3 = y3 % p
				//Done

				mpz_clears(op1, op2, NULL);
			}
			else {
				m_coord.isInf = PT_INFTY;
			}
		}
	}
	return true;
}

void EllPoint::operator*=(mpir_ui n) {
	if (this->m_coord.isInf) return;
	if (n == 0) {
		this->m_coord.isInf = PT_INFTY;
		return;
	}
	EllPoint tmp(*this);
	// set *this to neutral element
	this->setCoord(ECCoord());
	while (n != 0) {
		if (n & 1) *this += tmp;
		tmp += tmp;
		n >>= 1;
	}
}

void EllPoint::operator*=(const mpz_t& n) {
	if (this->m_coord.isInf) return;
	if (mpz_sgn(n) == 0) {
		this->m_coord.isInf = PT_INFTY;
		return;
	}
	EllPoint tmp(*this);
	mpz_t t;
	mpz_init_set(t, n);
	this->setCoord(ECCoord());
	while (mpz_sgn(t) != 0) { // While t != 0
		if(mpz_odd_p(t)) *this += tmp;
		tmp += tmp;
		mpz_tdiv_q_2exp(t, t, 1);
	}
	mpz_clear(t);
}

const ECCoord& EllPoint::getCoord() const {
	return m_coord;
}

const ECParam& EllPoint::getECParam() const {
	return m_param;
}


bool EllPoint::setCoord(const ECCoord& coord) {
	if (isOnCurve(coord)) {
		this->m_coord = coord;
		return true;
	}
	return false;
}

void EllPoint::setRandomCoord() {
	mp_bitcnt_t n = mpz_sizeinbase(m_param.p, 2);
	mpz_t op1;
	mpz_init(op1);
	while (true) {
		Rand.randomb(m_coord.x, n);
		mpz_mod(m_coord.x, m_coord.x, m_param.p); // Even if m_coord.x has same bitcount than m_param.p, it could come out larger
		set_weierstrass(op1, m_coord.x, m_param.a, m_param.b, m_param.p); // op1 = x^3 + ax + b
		if (mpz_sqrtm(m_coord.y, op1, m_param.p))
			break;
	}
	m_coord.isInf = NOT_INFTY;
	mpz_clear(op1);
}

bool EllPoint::isInf() const {
	return m_coord.isInf;
}

void EllPoint::setInf() {
	m_coord.isInf = PT_INFTY;
}

void EllPoint::order(mpz_t& dest) const {
	// Naively brute-forcing is so unefficient that we won't need smthng larger than mpir_ui
	mpir_ui n = 0;
	EllPoint Q(*this);
	while (!Q.isInf()) {
		Q += *this; n++;
	}
	mpz_set_ui(dest, n);
}

void EllPoint::order_co(mpz_t& dest, prime_factorization curveOrderFactorization) const{
	mpz_t m, op1;
	mpz_set_from_prime_fact(m, curveOrderFactorization);
	mpz_init(op1);
	for (auto& comp : curveOrderFactorization) {
		mpz_pow_ui(op1, comp.p, comp.exp);
		mpz_divexact(m, m, op1);
		EllPoint Q(*this);
		Q *= m;
		while (!Q.isInf()) {
			Q *= comp.p;
			mpz_mul(m, m, comp.p);
		}
	}
	mpz_clear(op1);
	mpz_set(dest, m);
	mpz_clear(m);
}

void EllPoint::print(const char* name) const {
	if (m_coord.isInf) {
		printf("[%s] Point is at infinty (neutral element)\n", name);
		return;
	}
	printf("[%s] Point coordinates are \n\tx : ", name);
	mpz_out_str(stdout, 10, m_coord.x);
	printf("\n\ty : ");
	mpz_out_str(stdout, 10, m_coord.y);
	printf("\nwith \tq : ");
	mpz_out_str(stdout, 10, m_param.p);
	printf("\n");
}


bool EllPoint::isOnCurve(const ECCoord& coord) const {
	if (coord.isInf) return true;
	mpz_t op1, op2, op3;
	mpz_init_set(op1, coord.y);
	mpz_powm_ui(op1, op1, 2, m_param.p);
	//op1 = y^2
	mpz_init_set(op2, coord.x);
	mpz_init_set(op3, coord.x);
	mpz_mul(op2, coord.x, m_param.a);
	mpz_mul(op3, coord.x, coord.x);
	mpz_addmul(op2, op3,coord.x);
	mpz_add(op2, op2, m_param.b);
	mpz_mod(op2, op2, m_param.p);
	//op2 = x^3 + ax + b

	int res = mpz_cmp(op1, op2);
	mpz_clear(op1);
	mpz_clear(op2);
	mpz_clear(op3);
	return res == 0;
}
