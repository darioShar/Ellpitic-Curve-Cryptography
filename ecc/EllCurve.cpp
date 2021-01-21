#include "EllCurve.h"

EllCurve::EllCurve() : m_gen(ECParam())
{
	mpz_init(m_genOrder);
}

EllCurve::EllCurve(const ECParam& param) : m_gen(ECParam())
{
	mpz_init(m_genOrder);
	setECParam(param);
}

EllCurve::EllCurve(ECParam&& param) : m_gen(ECParam())
{
	mpz_init(m_genOrder);
	setECParam(param);
}

EllCurve::EllCurve(const EllCurve& curve) 
	: m_gen(curve.getECParam()) // Generator is initialiazed to 0 (inty)
	, m_param(curve.getECParam())
{
	mpz_init(m_genOrder);
}


EllCurve::~EllCurve() {
	mpz_clear(m_genOrder);
}

/*Returns i*P where P is the generator in use*/
EllPoint EllCurve::operator[](const mpz_t& i) const{
	EllPoint tmp(m_gen);
	tmp *= i;
	return tmp;
}

/*Returns i*P where P is the generator in use*/
EllPoint EllCurve::operator[](mpir_ui i) const {
	EllPoint tmp(m_gen);
	tmp *= i;
	return tmp;
}

/*Returns the parameters of the curve in use*/
const ECParam& EllCurve::getECParam() const {
	return m_param;
}

/*Sets new curve parameters, finds new generator and then returns true
if parameters are usable (i.e no problems such as curve singularity,
null discriminant or weak choice of p (?)). Returns false otherwise.*/
bool EllCurve::setECParam(ECParam&& param, bool verbose) {
	return setECParam(param, verbose);
}

/*Sets new curve parameters, sets geenrator to 0 and then returns true
if parameters are usable (i.e no problems such as curve singularity,
null discriminant or weak choice of p (?)). Returns false otherwise.*/
bool EllCurve::setECParam(const ECParam& param, bool verbose) {
	if (verbose) {
		if (mpz_divisible_ui_p(param.p, 2) || mpz_divisible_ui_p(param.p, 3)) {
			printf("Error : elliptic curves defined on field of characteristic 2 or 3 aren't supported.\n");
			return false;
		}
		mpz_t discr, j; mpz_inits(discr, j, NULL); // Also treating them as temp var
		if (!Rand.is_likely_prime(param.p)) {
			printf("Warning : field order is not prime.\n");
			return false;
		}
		else {
			int shifts = 0;
			mpz_set(discr, param.p);
			while (mpz_sgn(discr)) {
				mpz_fdiv_q_2exp(discr, discr, 1);
				shifts++;
			}
			printf("Curve security level : %d bits\n", shifts / 2);
		}
		set_weierstrass_discriminant_j_invariant(discr, j, param.a, param.b, param.p);
		if (mpz_sgn(discr) == 0) {
			printf("Error : new curve parameters yield null discriminant : curve is singular and cannot be used.\n");
			mpz_clears(discr, j, NULL);
			return false;
		}
		if (mpz_sgn(j) == 0) {
			printf("Warning : new curve parameters yield null j_invariant : curve is super-singular, potentially weak\n");
		}
		else if (mpz_cmp_ui(j, 1728) == 0) {
			printf("Warning : new curve parameters yield j_invariant = 1728 : curve is super-singular, potentially weak\n");
		}
		m_param = param;
		if (mpz_sgn(param.order) == 0) {
			printf("Order not specified.\n");
		}
		// if order is given, assume it is right
		printf("New curve set with :\n\tDiscriminant = "); mpz_print(discr);
		printf("\tj_invariant = "); mpz_print(j);
		m_gen = EllPoint(param); // Resetting gen point to neutral element on new curve
		mpz_set_ui(m_genOrder, 1);
		mpz_clears(discr, j, NULL);
		print("new curve", true);
		return true;
	}

	//else
	if (mpz_divisible_ui_p(param.p, 2) || mpz_divisible_ui_p(param.p, 3))
		return false;
	mpz_t discr, j; mpz_inits(discr, j, NULL); // Also treating them as temp var
	if (!Rand.is_likely_prime(param.p))
		return false;

	set_weierstrass_discriminant_j_invariant(discr, j, param.a, param.b, param.p);
	if (mpz_sgn(discr) == 0) {
		mpz_clears(discr, j, NULL);
		return false;
	}
	m_param = param;

	m_gen = EllPoint(param); // Resetting gen point to neutral element on new curve
	mpz_set_ui(m_genOrder, 1);
	mpz_clears(discr, j, NULL);
	return true;
}

/* Finds new random parameter A, B and P for the curve. If an argument 
	P is passed, the most probable prime greater than P will be used and 
	new A and B will be found accordingly. Gen point is reset to neutral point*/
void EllCurve::setRandomParam(mpir_ui p, bool findOrder, bool verbose) {
	mpz_t tmp; mpz_init_set_ui(tmp, p);
	setRandomParam(tmp, findOrder, verbose);
	mpz_clear(tmp);
}

/*	Finds new random parameter A, B and P for the curve. If an argument 
	P is passed, the most probable prime greater than P will be used and 
	new A and B will be found accordingly. Gen point is reset to neutral point. */
void EllCurve::setRandomParam(const mpz_t p, bool findOrder, bool verbose) {
	if (p != NULL) {
		mpz_set(m_param.p, p);
		while (Rand.is_likely_prime(m_param.p) == 0)
			Rand.next_prime_candidate(m_param.p, m_param.p);
	}
	else {
		do {
			Rand.randomb(m_param.p, rand() % MAX_RANDOM_BITCOUNT);
		} while (!Rand.is_likely_prime(m_param.p));
	}
	ECParam param(m_param);
	//mpz_set_ui(param.order, 0); // so curve order is found
	while (true) {
		if (verbose) {
			printf("Constructing curve in E(Fp) where p = "); 
			mpz_print(param.p);
		}
		for (int i = 0; i < 3; i++) {
			Rand.randomm(param.a, param.p);
			Rand.randomm(param.b, param.p);
			if (setECParam(param)) {
				if (verbose)
					printf("Finding curve order...\n");
				if (!findCurveOrderHasseBSGS())
					continue;
				return;
			}
		}
		do { Rand.next_prime_candidate(param.p, param.p); }
			while (Rand.is_likely_prime(param.p) == 0);
	}

}

/*Returns generator in use*/
const EllPoint& EllCurve::getGen() const {
	return m_gen;
}

bool EllCurve::setGen(const ECCoord& coord, const mpz_t& order) {
	mpz_set(m_genOrder, m_param.order);
	if (m_gen.setCoord(coord)) {
		if (order) {
			mpz_set(m_genOrder, order);
			return true;
		}
		else {
			find_exact_order(m_genOrder, getGen());
			return true;
		}
	}
	return false;
}

/*Finds a new generator by efficiently randomly finding a new point on
the elliptic curve in use.*/
void EllCurve::findNewGen() {
	m_gen.setRandomCoord();
}

void EllCurve::findGenOrder() {
	mpz_set(m_genOrder, m_param.order);
	find_exact_order(m_genOrder, m_gen);
}

/*Returns curve order.*/
const mpz_t& EllCurve::getCurveOrder() const {
	return m_param.order;
}

/*Returns order of the subgroup generated by generator. For now, since order is prime, same as curve order.*/
const mpz_t& EllCurve::getGenOrder() const {
	return m_genOrder;
}

/*Finds a point on the elliptic curve in use, by random trials.
Because this can be very slow for large Fp, user can specify a time threshold
in ms after which the funciton will fail, returning false. 0 means no threshold,
which may lead to freezing the program.
On success returns true with found point stored in dest*/
bool EllCurve::getRandomPointRandomTrials(EllPoint& dest, mpir_ui time_threshold_ms) {
	mpz_t op1, op2;
	mpz_inits(op1, op2, NULL);
	clock_t begin = clock();
	while (true) {
		for (int i = 0; i < 1024; i++) {
			Rand.randomm(op1, m_param.p);
			Rand.randomm(op2, m_param.p);
			if (dest.setCoord(ECCoord(NOT_INFTY, op1, op2))) {
				mpz_clears(op1, op2, NULL);
				return true;
			}
		}
		if (1000.0 * double(clock() - begin) / (double)CLOCKS_PER_SEC > time_threshold_ms) break;
	}
	mpz_clears(op1, op2, NULL);
	return false;
}

/*Returns kP with random k which will then be stored in random_k parameter*/
EllPoint EllCurve::getRandomPoint(mpz_t& random_k) {
	Rand.randomb(random_k, mpz_sizeinbase(m_param.p, 2));
	EllPoint tmp(m_gen);
	tmp *= random_k;
	return tmp;
}

/*Finds curve order using simple formula |E(Fq)| = 1 + sum for x in Fq (1  + legendre(x^3 + ax + b, q) )
Result will be stored in internal ECParam.order*/
void EllCurve::findCurveOrderNaive() {
	// using simple formula |E(Fq)| = 1 + sum for x in Fq (1  + legendre(x^3 + ax + b, q) )
	// or |E(Fq)| = 1 + q + sum for x in Fq legendre(x^3 + ax + b, q) but we'll use the former 
	// because we can only easily add unsigned integers with mpir.
	mpz_set_ui(m_param.order, 1);
	mpz_t x, weier_x;
	mpz_init(weier_x);
	mpz_init_set_ui(x, 0);
	while (mpz_cmp(x, m_param.p) < 0) {
		set_weierstrass(weier_x, x, m_param.a, m_param.b, m_param.p);
		mpz_add_ui(m_param.order, m_param.order, 1 + mpz_legendre(weier_x, m_param.p));
		mpz_add_ui(x, x, 1);
	}
	mpz_clears(weier_x, x, NULL);

}

/*Finds curve order naively using Hasse Theorem
Result will be stored in internal ECParam.order*/
void EllCurve::findCurveOrderHasseNaive() {
	std::vector<mpir_ui> orders;
	mpz_t ppcm; mpz_init_set_ui(ppcm, 1);

	mpz_t inf; mpz_init(inf);
	mpz_sqrt(inf, m_param.p);
	mpz_mul_ui(inf, inf, 2); // inf = 2sqrt(p)
	mpz_sub(inf, m_param.p, inf);
	mpz_add_ui(inf, inf, 1);
	// inf = p + 1 - 2sqrt(p)

	mpz_t width; mpz_init_set_ui(width, 0);
	mpz_sqrt(width, m_param.p);
	mpz_mul_ui(width, width, 4);
	// width = 4sqrt(p)
	mpir_ui w = mpz_get_ui(width) + 1;

	mpz_t k, m; mpz_inits(k, m, NULL);
	 // Finding k such that k*gen = 0, with p + 1 - 2sqrt(p) <= k <= p + 1 + 2sqrt(p)
	do {
		do {
			findNewGen();
		} while (getGen().isInf());
		EllPoint P(getGen());
		EllPoint Q(getGen());
		Q *= inf;

		mpz_set(k, inf);
		while(!Q.isInf()) {
			Q += P;
			mpz_add_ui(k, k, 1);
		}

		// prime factorization of k now

		// if k is prime we're done
		if (Rand.is_likely_prime(k)) {
			if (mpz_cmp(k, inf) < 0 // otherwise we have found the curve order, which is then prime (k > p + 1 - 2sqrt(p))
				&& std::find(orders.begin(), orders.end(), mpz_get_ui(k)) == orders.end()) { // hoping k is sufficiently small
				orders.push_back(mpz_get_ui(k));
			}
			mpz_set(m_genOrder, k);
			mpz_lcm(ppcm, ppcm, k);
			continue;
		}

		// k = p1^a1 * .. * pr^ar. Finding k' such that ord(gen) = k' = p1^a1' * ... * pr^ar'
		find_exact_order(k, getGen());

		if (std::find(orders.begin(), orders.end(), mpz_get_ui(k)) == orders.end()) { // so we are sure it is not an element already seen
			orders.push_back(mpz_get_ui(k));
			mpz_lcm(ppcm, ppcm, k);
		}
		mpz_set(m_genOrder, k);

	} while (mpz_cmp(ppcm, width) <= 0);

	// Finding unique N in right range : N = n*ppcm, ppcm > 4sqrt(p), p + 1 - 2sqrt(p) <= n*ppcm <= p + 1 + 2sqrt(p)
	// so inf <= ppcm*n <= inf + 4sqrt(p) and inf / ppcm <= n <= inf / ppcm + e where e < 1 : n = ceil(inf / ppcm)
	mpz_t n; mpz_init(n);
	mpz_cdiv_q(n, inf, ppcm); // n = ceil(inf / ppcm)
	mpz_mul(m_param.order, n, ppcm); // ok : N = n * ppcm

	mpz_clears(n, ppcm, inf, width, k, m, NULL);
}

/*Finds curve order using baby-steps giant-steps. Result will be stored in internal ECParam.order*/
bool EllCurve::findCurveOrderHasseBSGS() {

	std::vector<mpir_ui> orders;
	mpz_t ppcm; mpz_init_set_ui(ppcm, 1);

	mpz_t sup; mpz_init(sup);
	mpz_sqrt(sup, m_param.p);
	mpz_mul_ui(sup, sup, 2); 
	mpz_add(sup, sup, m_param.p);
	mpz_add_ui(sup, sup, 1); // sup =  p + 1 + 2sqrt(p)

	mpz_t width; mpz_init_set_ui(width, 0);
	mpz_sqrt(width, m_param.p);
	mpz_mul_ui(width, width, 4); // width =  4sqrt(p)

	mpz_t m, k, l, tmp;
	mpz_inits(m, k, l, tmp, NULL);
	mpz_sqrt(m, width); // m = sqrt(4sqrt(p))

	do {
		do {
			findNewGen();
		} while (getGen().isInf());

		// If space needed is too big for memory allocation (let's say 4*sizeof(ECPair) , if hash function isn't too good)
		mpz_mul_ui(tmp, m, 4*sizeof(ECPair));
		if (mpz_cmp_ui(tmp, BSGS_MEMORY) > 0) {
			printf("Cannot run bsgs : too much memory needed : "); mpz_print(tmp);
			mpz_clears(m, k, l, tmp, ppcm, sup, width, NULL);
			return false;
		}

		// Finding k such that k*gen = 0, with p + 1 - 2sqrt(p) <= k <= p + 1 + 2sqrt(p)
		// k =(p + 1 -2sqrt(p)) + am +b where m = sqrt(4sqrt(p)) and -(p + 1)gen = (am + b)gen

		// Pre-computing j*gen fo j in 0, m
		std::unordered_set<ECPair, BSGSHasher> data;
		pre_compute_bsgs(data, getGen(), mpz_get_ui(m)+1);

		// Now finding k
		EllPoint Q(getGen());
		mpz_divexact_ui(tmp, width, 2); // tmp = 2sqrt(p)
		mpz_sub(tmp, m_param.p, tmp);
		mpz_add_ui(tmp, tmp, 1); // tmp = p + 1 - 2sqrt(p)
		Q *= tmp;
		Q.inverse();
		mpz_set_ui(k, 0);
		if (!Q.isInf())
			find_k_bsgs(k, data, getGen(), Q, mpz_get_ui(m));
		
		mpz_add(k, k, tmp); // ok
		data.clear();

		EllPoint P = getGen();
		P *= k;
		if (!P.isInf()) {
			printf("problem, with k = "); mpz_print(k);
			P.print();
			mpz_clears(m, k, l, tmp, ppcm, sup, width, NULL);
			return false;
		}

		// prime fact of k now
		std::vector<mpir_ui> primes;
		EllPoint T(getGen());

		// If k is big enough
		if (mpz_cmp(k, width) > 0) {
			// gen order is still undefined
			mpz_lcm(ppcm, ppcm, k);
			continue;
		}
		// if k is prime we're done
		if (Rand.is_likely_prime(k)) {
			mpz_sub(tmp, sup, width);
			if (mpz_cmp(k, tmp) < 0 // otherwise we have found the curve order, which is prime (k > p + 1 - 2sqrt(p))
				&& std::find(orders.begin(), orders.end(), mpz_get_ui(k)) == orders.end()) { // hoping k is sufficiently small
				orders.push_back(mpz_get_ui(k));
			}
			mpz_set(m_genOrder, k);
			mpz_lcm(ppcm, ppcm, k);
			continue;
		}

		// k = p1^a1 * .. * pr^ar. Finding k' such that ord(gen) = k' = p1^a1' * ... * pr^ar'
		// this function can take a lot of memory (prime table) // only about sqrt(sqrt(p)) since k < width. ok
		find_exact_order(k, getGen());
		
		if (std::find(orders.begin(), orders.end(), mpz_get_ui(k)) == orders.end()){ // so we are sure it is not an element already seen
			orders.push_back(mpz_get_ui(k));
			mpz_lcm(ppcm, ppcm, k);
		}
		mpz_set(m_genOrder, k);

	} while (mpz_cmp(ppcm, width) <= 0);

	// Finding unique N in right range : N = n*ppcm, ppcm > 4sqrt(p), p + 1 - 2sqrt(p) <= n*ppcm <= p + 1 + 2sqrt(p)
	// so sup - 4sqrt(p) <= ppcm*n <= sup and sup / ppcm - e<= n <= sup / ppcm where e < 1 : n = floor(sup / ppcm)
	mpz_t n; mpz_init(n);
	mpz_fdiv_q(n, sup, ppcm); // n = floor(sup / ppcm)
	mpz_mul(m_param.order, n, ppcm); // ok : N = n * ppcm

	mpz_clears(n, m, k, l, tmp, ppcm, sup, width, NULL);
	return true;
}

/*Finds k satisfying K = k*P where P is the generator in use.
Cracking the discrete logarithm problem naively.
Returns true if found before time threshold, false otherwise.
On failure, all number up to updated k will have been tested.*/
bool EllCurve::crackDiscreteLogNaive(mpz_t& k, const EllPoint& K, mpir_ui time_threshold_ms) const {
	mpz_set_ui(k, 1);
	clock_t begin = clock();
	EllPoint tmp(this->getGen());
	while (true) {
		for (int i = 0; i < 1000; i++) {
			tmp += this->getGen();
			if (tmp == K) {
				mpz_add_ui(k, k, i + 1);
				return true;
			}
		}
		mpz_add_ui(k, k, 1000);
		if (1000.0 * double(clock() - begin) / (double)CLOCKS_PER_SEC > time_threshold_ms) break;
	}
	return false;
}

/*Finds k satisfying K = k*P where P is the generator in use.
Cracking the discrete logarithm problem with baby step giant step method.
max_memory_usage is the maximum amount of memory the function 
will be able to use (default being BSGS_MEMEORY bytes). sqrt(order) must be < (unsigned long long)~0
Returns true if found before time threshold and with less memory than max to allocate,
false otherwise*/
bool EllCurve::crackDiscreteLogBSGS(mpz_t& k, const EllPoint& K, mpir_ui time_threshold_ms, mpir_ui max_memory_usage) const {
	// If order is 0
	if (mpz_sgn(this->m_param.order) == 0) {
		printf("Curve order is null\n");
		return false;
	}
	// assumes order is correct
	mpz_t m, tmp;
	mpz_inits(m, tmp, NULL);
	mpz_sqrtrem(m, tmp, this->m_param.order);
	if (mpz_sgn(tmp) != 0) mpz_add_ui(m, m, 1); // now m is ceil(sqrt(order))
	

	// If space needed is superior to max_memory_usage
	mpz_mul_ui(tmp, m, 4*sizeof(ECCoord));
	if (mpz_cmp_ui(tmp, max_memory_usage) > 0) {
		printf("Cannot run bsgs : too much memory needed : "); mpz_print(tmp);
		mpz_clears(m, tmp, NULL);
		return false;
	}

	mpir_ui s = mpz_get_ui(m);
	mpz_clears(m, tmp, NULL);

	// Pre-computing jP fo j in 0, m-1
	std::unordered_set<ECPair, BSGSHasher> data;
	pre_compute_bsgs(data, getGen(), s);

	// Now finding k
	find_k_bsgs(k, data, getGen(), K, s);

	return true;
}

void EllCurve::print(const char* name, bool printGen) const {
	printf("[%s] Curve is defined by ", name);
	printf("\n\tE(Fq) : y^2 = x^3 + "); mpz_out_str(stdout, 10, m_param.a);
	printf(" * x + ");					mpz_out_str(stdout, 10, m_param.b);
	printf("\n\twhere q = ");			mpz_out_str(stdout, 10, m_param.p);
	printf("\n\torder is = ");			mpz_out_str(stdout, 10, m_param.order);
	if (printGen) {
		printf("\n");
		m_gen.print("Generator");
	}
	printf("\n");
}

void EllCurve::find_exact_order(mpz_t& k, const EllPoint& G) {
	std::vector<mpir_ui> primes;
	EllPoint T(G);

	mpz_t m; mpz_init(m);
	mpz_sqrt(m, k);
	mpir_ui s = mpz_get_ui(m);
	get_primes(primes, s);

	// k = p1^a1 * .. * pr^ar. Finding k' such that ord(gen) = k' = p1^a1' * ... * pr^ar'
	for (auto& i : primes) {
		if (!mpz_divisible_ui_p(k, i)) continue;
		while (mpz_divisible_ui_p(k, i)) {
			mpz_divexact_ui(k, k, i);
			T = G;
			T *= k;
			if (!T.isInf()) {
				mpz_mul_ui(k, k, i);
				break;
			}
		}
	}
	mpz_clear(m);
}

void EllCurve::pre_compute_bsgs(std::unordered_set<ECPair, BSGSHasher>& data, const EllPoint& G, mpir_ui m) {
	EllPoint P(G);
	data.reserve(m + 1); // so there is no rehash during computation
	EllPoint InvGen(P); InvGen.inverse();
	P *= m;
	while (m > 0) {
		data.emplace(ECPair(P.getCoord(), m));
		P += InvGen;
		m--;
	}
	data.emplace(ECPair(P.getCoord(), m)); // m = 0
}

// find k != such that k * G = K using baby step giant step. m is such that it is possible to write k = (am + b), 0<= a,b <= m.
// data must hold j*gen for 0<= j <= m
void EllCurve::find_k_bsgs(mpz_t& k, const std::unordered_set<ECPair, BSGSHasher>& data, const EllPoint& gen, const EllPoint& K, mpir_ui m) {
	// assumes m can be stored in mpir_ui but m^2 might not
	mpz_set_ui(k, 1);
	mpir_ui l=0; // l <= m can be stored in mpir_ui 
	EllPoint _mPoint(gen);
	_mPoint *= m;
	_mPoint.inverse();
	EllPoint R(K);
	R.isInf() ? l = 1 : l = 0; // beginning at l = 1 so k is not set to 0 if R.isInf()

	while (true) {
		std::unordered_set<ECPair, BSGSHasher>::const_iterator it = data.find(ECPair(R.getCoord(), 0)); // Since second value doesn't count in hash
		if (it != data.end()) {
			mpir_ui j_value = (*it).val;
			mpz_mul_ui(k, k, m);
			mpz_mul_ui(k, k, l);
			mpz_add_ui(k, k, j_value); // k = ml + j
			break;
		}
		R += _mPoint;
		l++;
	}
}


/* For BSGS 
// Since m can be very big, we have to split the task (not enough memory)
// we'll make chunks of few hundred megabytes (of size max_memory_usage bytes)

mpir_ui chunk_size;
mpir_ui remains;
mpir_ui num_chunks = 1;
mpz_mul_ui(tmp, m, sizeof(ECCoord)); // do not forget sizeof
while (mpz_cmp_ui(tmp, max_memory_usage) > 0) {
mpz_tdiv_q_2exp(tmp, tmp, 1);
num_chunks *= 2;
}
chunk_size = mpz_get_ui(tmp);
//Do not forget the possibly "forgotten" values
mpz_mul_ui(tmp, tmp, num_chunks*sizeof(ECCoord));
mpz_sub(tmp, m, tmp);
remains = mpz_get_ui(tmp);
*/


//BSGS in while(true) loop, because of possibly expansive time check
/*
	for (int i = 0; i < 1000; i++, point += _mPoint) {
			auto iter = data.find(ECPair(point.getCoord(), 0)); // Since second value doesn't count in hash
			if (iter != data.end()) {
				mpir_ui j_value = (*iter).val;
				mpz_add_ui(l, l, i);
				mpz_mul(l, l, m);
				mpz_add_ui(l, l, j_value);
				mpz_set(k, l);
				mpz_clears(m, l, j, tmp, NULL);
				data.clear();
				return true;
			}
		}
		mpz_add_ui(l, l, 1000);
*/



/*Finds curve order using basic Schoof algorithm. Result will be stored in internal ECParam.order*/
void EllCurve::findCurveOrderSchoof() {
	/*
	mpz_t N; mpz_init(N);
	mpz_sqrt(N, m_param.p); mpz_add_ui(N, N, 1);
	mpz_mul_ui(N, N, 4); // N now is approx 4 * sqrt(p)
	// First construct set of primes. First 100 should be sufficient, but well.
	std::vector<mpir_ui> primes;
	mpir_ui nbPrimes = 0;
	do {
	nbPrimes += 100;
	get_primes(primes, nbPrimes);
	} while (!is_product_geq(primes, N));

	// Preparing array of remainders
	std::vector<mpir_ui> remainders(primes.size());

	polynomial<ZpZ> X;
	polynomial<ZpZ> unit;
	// For l = primes[0] = 2
	X = polynomial<ZpZ>{ { ZpZ(0, 2), ZpZ(1, 2) } };
	unit = polynomial<ZpZ>{ { ZpZ(1, 2) } };
	// FUUUUUUUUUUUUUUUUUUUCK FUUUUUUUUUUUUUUUUUUUCK FUUUUUUUUUUUUUUUUUUUCK FUUUUUUUUUUUUUUUUUUUCK
	polynomial<ZpZ> XP;
	polynomial<ZpZ> tmp = X;
	mpz_t op1; mpz_init_set(op1, m_param.p);
	while (mpz_sgn(op1) > 0) {
	if (mpz_odd_p(op1)) XP *= tmp;
	tmp *= tmp;
	mpz_tdiv_q_2exp(op1, op1, 1);
	}
	// XP = X^p
	// gcd(Xp - X, X^3 + AX + B) = 1 <=> remainders[0] equiv 0 [2]
	gcd(XP - X, X*X*X + X*ZpZ(m_param.a, 2) + unit*ZpZ(m_param.b, 2)).degree() == 0 ? remainders[0] = 1 : remainders[0] = 0;
	*/
	/**************************************************************************************************************************************************************/


	//mpz_clears(N, op1, NULL);
}

/*Finds curve order using improved Schoof algorithm (Schoof-Elkies-Atkin aka SEA). Result will be stored in internal ECParam.order*/
void EllCurve::findCurveOrderSEA() {

}