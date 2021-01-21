#pragma once
#include <thread>
#include <algorithm>
#include <mutex>
#include "Constants.h"
#include "plot.h"
#include "ECDH.h"
#include "ECDSA.h"
#include "Poly.h"

#define THREADS 4

std::mutex order_mutex;

void test_bitcoin() {
	EllCurve curve;
	curve_init_bitcoin(curve);

	EllPoint Q(curve.getGen());
	curve.print("TEST", true);
	Q += Q;
	Q.print("Doubling");
	Q += curve.getGen();
	Q.print("Addition");
	Q = curve[curve.getECParam().order];
	if (!Q.isInf()) printf("Error in Big Multiplication !!");
	else Q.print("Big Multiplication");
	Q = curve[1000];
	Q.print("Small Multiplication");
	curve.findNewGen();
	curve.print("New Gen", true);
	Q = curve[1 << 10];
	Q.print("Key to crack : (1 << 10) * Gen");
	mpz_t n; mpz_init(n);
	if (curve.crackDiscreteLogNaive(n, Q, 5000)) {
		printf("\nCracked discrete log problem naively : ");
		mpz_out_str(stdout, 10, n);
		printf("\n");
	}
	else {
		printf("\nFailed to crack discrete log problem naively under 5000 ms. Tested up to : ");
		mpz_out_str(stdout, 10, n);
		printf("\n");
	}
	if (curve.crackDiscreteLogBSGS(n, Q, 5000)) {
		printf("Cracked discrete log problem with bsgs : ");
		mpz_out_str(stdout, 10, n);
		printf("\n");
	}
	else
		printf("\nFailed to crack discrete log problem with bsgs under 5000 ms\n");

	curve.print("TEST, New Generator", true);
	if (curve.getRandomPointRandomTrials(Q, 100)) Q.print("Random Point");
	else printf("Failed to find random point under 100 ms.\n");
	printf("---------------------------------------------------\n");

	mpz_clears(n, NULL);
}

void test_default() {
	EllCurve curve;
	curve_init_default(curve);

	EllPoint Q(curve.getGen());
	curve.print("TEST", true);
	Q += Q;
	Q.print("Doubling");
	Q += curve.getGen();
	Q.print("Addition");
	Q = curve[curve.getECParam().order];
	if (!Q.isInf()) printf("Error in Big Multiplication !!\n"); // Oui effectivement on ne connait pas l'ordre de la courbe
	else Q.print("Big Multiplication");
	Q = curve[1000];
	Q.print("Small Multiplication");
	curve.findNewGen();
	curve.print("New Gen", true);
	Q = curve[1 << 10];
	Q.print("Key to crack : (1 << 10) * Gen");
	mpz_t n; mpz_init(n);
	if (curve.crackDiscreteLogNaive(n, Q, 5000)) {
		printf("\nCracked discrete log problem naively : ");
		mpz_out_str(stdout, 10, n);
		printf("\n");
	}
	else {
		printf("\nFailed to crack discrete log problem naively under 5000 ms. Tested up to : ");
		mpz_out_str(stdout, 10, n);
		printf("\n");
	}
	if (curve.crackDiscreteLogBSGS(n, Q, 5000)) {
		printf("Cracked discrete log problem with bsgs : ");
		mpz_out_str(stdout, 10, n);
		printf("\n");
	}
	else
		printf("\nFailed to crack discrete log problem with bsgs under 5000 ms\n");

	curve.print("TEST, New Generator", true);
	if (curve.getRandomPointRandomTrials(Q, 100)) Q.print("Random Point");
	else printf("Failed to find random point under 100 ms.\n");
	printf("---------------------------------------------------\n");

	mpz_clears(n, NULL);
}

void test_bignum() {
	EllCurve curve;
	curve_init_bignum(curve);

	EllPoint Q(curve.getGen());
	curve.print("TEST", true);
	Q += Q;
	Q.print("Doubling");
	Q += curve.getGen();
	Q.print("Addition");
	mpz_t p; mpz_init_set_str(p, BIG_P, 0);
	Q = curve[p];
	Q.print("Big Multiplication");
	Q = curve[1000];
	Q.print("Small Multiplication");
	curve.findNewGen();
	curve.print("New Gen", true);

	curve.print("TEST, New Generator", true);
	if (curve.getRandomPointRandomTrials(Q, 100)) Q.print("Random Point");
	else printf("Failed to find random point under 100 ms.\n");
	printf("---------------------------------------------------\n");
}

void long_test() {
	mpz_t n;
	mpz_init(n);
	EllCurve curve;
	printf("Testing with defaut parameters.\n");
	curve_init_default(curve);
	printf("Finding curve order 5 times\n");
	for (int i = 0; i < 5; i++)
		curve.findCurveOrderNaive();
	printf("Found order : "); mpz_print(curve.getECParam().order);
	printf("Should be approx. "); mpz_print(curve.getECParam().p);
	curve.findNewGen();
	curve.crackDiscreteLogNaive(n, curve[5000], 500);
	curve.crackDiscreteLogBSGS(n, curve[5000], 500);
	printf("Cracked discrete log problem : ");
	mpz_print(n);

	printf("Testing with big nums.\n");
	curve_init_bignum(curve);
	printf("Finding generator 1000 times\n");
	for (int i = 0; i < 1000; i++)
		curve.findNewGen();
	EllPoint Q(curve.getGen());
	printf("1000 Doubling\n");
	for (int i = 0; i < 1000; i++)
		Q += Q;
	printf("1000 multiple point operation\n");
	for (int i = 0; i < 1000; i++) {
		Q += curve.getGen();
		Q += Q;
		Q *= 5198109484198019801;
	}
}

void helper_func(std::vector<float>& x_values, std::vector<float>& y_values, mpir_ui begin, mpir_ui end, mpir_ui incr, EllCurve curve) {
	/*
	curve.setRandomParam(begin, true);
	mpz_t tmp; mpz_init(tmp);
	std::vector<float> x_tmp; x_tmp.reserve((end / log(end)) - (begin / log(begin)));
	std::vector<float> y_tmp; y_tmp.reserve((end / log(end)) - (begin / log(begin)));
	std::size_t counter = 1;
	do {
		if (counter % 100 == 0) printf("Did %d iterations on this thread [%d : %d]\n", counter, begin, end);
		x_tmp.push_back(static_cast<float>(mpz_get_ui(curve.getECParam().p)));
		float order = static_cast<float>(mpz_get_ui(curve.getCurveOrder()));
		while (order > 1.5f*(float)(mpz_get_ui(curve.getECParam().p))) { // ????? WHY this error ? threads ?
			curve.findCurveOrderHasseNaive();
			order = static_cast<float>(mpz_get_ui(curve.getCurveOrder()));
		}
		y_tmp.push_back(order);
		mpz_add_ui(tmp, curve.getECParam().p, 1);
		curve.setRandomParam(tmp, true);
		counter++;
	} while (mpz_get_ui(curve.getECParam().p) < end);
	order_mutex.lock();
	x_values.insert(x_values.end(), x_tmp.begin(), x_tmp.end());
	y_values.insert(y_values.end(), y_tmp.begin(), y_tmp.end());
	order_mutex.unlock();
	printf("Thread calculation done : \
		x thread : %d, y thread : %d ; x total : %d, y total : %d\n",
		x_tmp.size(), y_tmp.size(), x_values.size(), y_values.size());
	mpz_clear(tmp);*/
}

/*Displays curves order for random curve defined over Fq with begin <= q <= end and q prime.*/
void display_random_curve_order(mpir_ui begin, mpir_ui end, mpir_ui incr) {
	EllCurve curve;
	std::vector<std::vector<float>> x_values(3);
	for (auto& x : x_values)
		x.reserve((end / log(end)) - (begin / log(begin)));
	std::vector<std::vector<float>> y_values(3);
	for (auto& y : y_values)
		y.reserve((end / log(end)) - (begin / log(begin)));
	for (mpir_ui i = begin; i <= end; i += incr) {
		curve.setRandomParam(i, true);
		int p = mpz_get_ui(curve.getECParam().p);
		x_values[0].push_back(p);
		if (!mpz_get_ui(curve.getCurveOrder())) {
			curve.print("problem", true);
			curve.findCurveOrderHasseNaive();
		}
		y_values[0].push_back(mpz_get_ui(curve.getCurveOrder()));
		x_values[1].push_back(p);
		x_values[2].push_back(p);
		y_values[1].emplace_back(p + 1 + 2 * std::sqrtf(p)); // Hasse Theorem upper bound
		y_values[2].emplace_back(p + 1 - 2 * std::sqrtf(p)); // Hasse Theorem lower bound
	}
	plot(x_values, y_values, std::vector<Vtk::Color> {Vtk::Color(1.0, 0.0, 0.0), Vtk::Color(1.0, 0.5, 0.0), Vtk::Color(1.0, 0.5, 0.0)}, ChartType::CTLINE, "Curve Order");
}

void test_ecdh(void init_curve(EllCurve&) = curve_init_secp112r1) {
	EllCurve curve;
	init_curve(curve);
	ECDH(curve);
}

void test_ecdsa(void init_curve(EllCurve&) = curve_init_secp112r1) {
	EllCurve curve;
	init_curve(curve);
	// Message sender private and public key
	mpz_t priv_key; mpz_init(priv_key);
	EllPoint publicKey = curve.getRandomPoint(priv_key);

	// Message to send
	// ...
	// ...

	// Hash of the message
	mpz_t hash; mpz_init_set_str(hash, "0xF1561A561F6A16151D1F20019BFFA90110A", 0);

	// Signing message
	mpz_t r, s;
	mpz_inits(r, s, NULL);
	SignMessage(r, s, curve, priv_key, hash);

	// Checking Signature
	if (CheckSignature(r, s, curve, publicKey, hash))
		printf("ECDSA test : success !\n");
	else
		printf("ECDSA test : failure !\n");

	mpz_clears(priv_key, r, s, hash, NULL);
}

void compute_formula_1(std::vector<Poly> &P, const Poly& Q, mpir_ui i) {
	Poly R;
	//5
	P[i] += P[(i / 2) + 2];
	P[i] *= P[i / 2];
	P[i] *= P[i / 2];
	P[i] *= P[i / 2];
	P[i] *= Q;
	R += P[(i / 2) - 1];
	R *= P[(i / 2) + 1];
	R *= P[(i / 2) + 1];
	R *= P[(i / 2) + 1];
	P[i] -= R;
	//P[i] = (Q * Q) * (P[(i / 2) + 2] * P[i / 2] * P[i / 2] * P[i / 2])
	//- (P[(i / 2) - 1] * P[(i / 2) + 1] * P[(i / 2) + 1] * P[(i / 2) + 1]);
	//m = 2n, P[2m + 1] = Y^4 * P[2n+2]P[2n]^3 - P[m-1]P[m+1]^3
}

void compute_formula_2(std::vector<Poly> &P, const Poly& Q, mpir_ui i) {
	Poly R;
	P[i + 1] += P[((i + 1) / 2) + 2];
	P[i + 1] *= P[((i + 1) / 2) - 1];
	P[i + 1] *= P[((i + 1) / 2) - 1];
	R += P[((i + 1) / 2) - 2];
	R *= P[((i + 1) / 2) + 1];
	R *= P[((i + 1) / 2) + 1];
	P[i + 1] -= R;
	P[i + 1] *= P[(i + 1) / 2];
	//P[i + 1] = P[(i + 1) / 2] * (P[((i + 1) / 2) + 2] * P[((i + 1) / 2) - 1] * P[((i + 1) / 2) - 1]
	//- P[((i + 1) / 2) - 2] * P[((i + 1) / 2) + 1] * P[((i + 1) / 2) + 1]);
	for (int j = 0; j < P[i + 1].size(); j++)
		mpz_divexact_ui(P[i + 1][j], P[i + 1][j], 2);
}

void compute_division_polynomials(int n, const mpz_t& a, const mpz_t& b, const std::string& name) {
	n = max(5, n);
	n += 4 - (n % 4);
	std::vector<Poly> P(n+1);

	P[0].resize(1);
	mpz_init_set_ui(P[0][0], 0);

	P[1].resize(1);
	mpz_init_set_ui(P[1][0], 1);

	P[2].resize(1);
	mpz_init_set_ui(P[2][0], 2);

	P[3].resize(5);
	mpz_init_set_ui(P[3][0], 0);
	mpz_submul(P[3][0], a, a);
	mpz_init_set_ui(P[3][1], 0);
	mpz_addmul_ui(P[3][1], b, 12);
	mpz_init_set_ui(P[3][2], 0);
	mpz_addmul_ui(P[3][2], a, 6);
	mpz_init_set_ui(P[3][3], 0);
	mpz_init_set_ui(P[3][4], 3);
	//P[3] =  -a^2 + 12bx + 6ax^2 + 3x^4

	mpz_t tmp; mpz_init_set_ui(tmp, 0);
	P[4].resize(7);
	mpz_init_set_ui(P[4][0], 0);
	mpz_submul(P[4][0], a, a);
	mpz_submul_ui(P[4][0], a, 4);
	mpz_submul(P[4][0], a, a);
	mpz_submul(tmp, b, b);
	mpz_submul_ui(P[4][0], tmp, 32);
	mpz_init_set_ui(P[4][1], 0);
	mpz_mul(tmp, a, b);
	mpz_submul_ui(P[4][1], tmp, 16);
	mpz_init_set_ui(P[4][2], 0);
	mpz_mul(tmp, a, a);
	mpz_submul_ui(P[4][1], tmp, 20);
	mpz_init_set_ui(P[4][3], 0);
	mpz_addmul_ui(P[4][3], b, 80);
	mpz_init_set_ui(P[4][4], 0);
	mpz_addmul_ui(P[4][4], a, 20);
	mpz_init_set_ui(P[4][5], 0);
	mpz_init_set_ui(P[4][6], 4);
	// P[4] = -4a^3-32b^2 -16abx -20a^2x^2 -80bx^3+20ax^4 +4x^6
	mpz_clear(tmp);

	Poly Q(4); // Q = Y^2 = X^3 + aX + b
	mpz_set(Q[0], b);
	mpz_set(Q[1], a);
	mpz_set_ui(Q[3], 1);
	Q *= Q; // Q = Q^2
	// Now computing, 4 different formulas depending on i mod 4
	std::thread t[4];
	for (int i = 5; i < n;) {
		t[0] = std::thread(compute_formula_1, std::ref(P), std::ref(Q), i);
		t[1] = std::thread(compute_formula_2, std::ref(P), std::ref(Q), i);
		i += 2;
		t[2] = std::thread(compute_formula_1, std::ref(P), std::ref(Q), i);
		t[3] = std::thread(compute_formula_2, std::ref(P), std::ref(Q), i);
		
		t[0].join();
		t[1].join();
		t[2].join();
		t[3].join();
		/*compute_formula_1(P, Q, i);
		compute_formula_2(P, Q, i);
		i += 2;
		compute_formula_1(P, Q, i);
		compute_formula_2(P, Q, i);*/

		i += 2;
	}

	ofstream f(name, ios::out);
	std::vector<char> vec(1 << 20);
	f.rdbuf()->pubsetbuf(&vec.front(), vec.size());
	for (int i = 0; i <= n; i++) {
		f << "\nP_" << i << " = "; P[i].writeOut(f);
	}
	f.close();
}

void test_1() {
	testVtk();
	display_random_curve_order(1000, 2500, 100);
	test_ecdh(curve_init_secp521r1);
	test_ecdsa(curve_init_secp521r1);
}

void test_2() {
	/*EllCurve curve(ECParam("11567", "34619", "100003", "99823"));
	curve.setGen(ECCoord(0, "30027", "93977"));*/
	EllCurve curve(ECParam("79282401", "92328317", "134217757", "134226857"));
	curve.setGen(ECCoord(0, "29141079", "26194183"), curve.getCurveOrder());
	curve.print("", true);
	mpz_t k; mpz_init(k);
	curve.crackDiscreteLogBSGS(k, curve[60000000]);
	mpz_print(k);
}

void test_3(int n) {
	// secp112r1
	ECParam p_s(secp112_A, secp112_B, secp112_P, secp112_n);
	//bitcoin
	ECParam p_b(BITCOIN_A, BITCOIN_B, "7", BITCOIN_N);

	//compute_division_polynomials(n, p_b.a, p_b.b, "div_pol_secp121r1.txt");
	compute_division_polynomials(n, p_s.a, p_s.b, "div_pol_bitcoin.txt");
}

void test_4(int n) {
	Poly Q(4); // Q = X^3 + X + 2
	mpz_set_ui(Q[3], 1);
	mpz_set_ui(Q[1], 1);
	mpz_set_ui(Q[0], 2);
	/*mpz_set_ui(Q[3], 254242040540);
	mpz_set_ui(Q[1], 540540510545);
	mpz_set_ui(Q[0], 525354040055);*/
	Poly P(Q);
	for (int i = 0; i < n; i++) {
		P *= P;
	}
	//P.print();
}

void test_5(int n) {
	mpz_t p; mpz_init_set_str(p , "0xDB7C2ABF62E35E668076BEAD208B", 0);
	Poly_p Q(p, 4); // Q = X^3 + X + 2
	mpz_set_ui(Q[3], 1);
	mpz_set_ui(Q[1], 1);
	mpz_set_ui(Q[0], 2);
	/*mpz_set_ui(Q[3], 254242040540);
	mpz_set_ui(Q[1], 540540510545);
	mpz_set_ui(Q[0], 525354040055);*/
	Poly_p R(Q);
	for (int i = 0; i < n; i++) {
		//R.mul(R);
		R *= R;
	}
	//R.print();
}

void test_hasse_bsgs_order() {
	EllCurve curve(ECParam("968113241544", "19232197347131", "23640121541677", "23640123950704"));
	curve.findCurveOrderHasseBSGS();
	mpz_print(curve.getCurveOrder());
}

void test_hasse_naive_order() {
	EllCurve curve(ECParam("968113241544", "19232197347131", "23640121541677", "23640123950704"));
	curve.setGen(ECCoord(0, "29141079", "26194183"), curve.getCurveOrder());
	curve.findCurveOrderHasseNaive();
	mpz_print(curve.getCurveOrder());
}

/*mpz_t p; mpz_init_set_str(p, "0xDB7C2ABF62E35E668076BEAD208B", 0);
Poly_p Q(p, 2); // Q = X + 1
mpz_set_ui(Q[1], 1);
mpz_set_ui(Q[0], 1);
Poly_p P(Q);

for (int j = 0; j < 9; j++)
P *= P;
Q = P;
Poly_p tmp(P);

std::vector<float> runTimes;
for (int i = 0; i < 7; i++) {
P *= P;
tmp = P;
QueryPerformanceCounter(&counterBegin);
tmp.mul (Q);
QueryPerformanceCounter(&counterEnd);
runTimes.push_back((double)(counterEnd.QuadPart - counterBegin.QuadPart) / (double)frequency.QuadPart);
}
std::for_each(runTimes.begin(), runTimes.end(), [](double time) {printf("%lf\n", time); });*/

/*EllCurve curve;
mpz_t p, t; mpz_inits(p,t,NULL);
mpz_set_ui(t, 2);
mpz_mul_2exp(p, t, 92);
do {
curve.setRandomParam(p, true);
mpz_add_ui(p, curve.getCurveOrder(), 1);
} while (!Rand.is_likely_prime(curve.getCurveOrder()));
curve.print("", true);*/


void compute_formula_1_p(std::vector<Poly_p>& P, const Poly_p& Q, const mpz_t& p, mpir_ui i) {
	Poly_p R(p);
	//5
	P[i] += P[(i / 2) + 2];
	P[i] *= P[i / 2];
	P[i] *= P[i / 2];
	P[i] *= P[i / 2];
	P[i] *= Q;
	R += P[(i / 2) - 1];
	R *= P[(i / 2) + 1];
	R *= P[(i / 2) + 1];
	R *= P[(i / 2) + 1];
	P[i] -= R;
	//P[i] = (Q * Q) * (P[(i / 2) + 2] * P[i / 2] * P[i / 2] * P[i / 2])
	//- (P[(i / 2) - 1] * P[(i / 2) + 1] * P[(i / 2) + 1] * P[(i / 2) + 1]);
	//m = 2n, P[2m + 1] = Y^4 * P[2n+2]P[2n]^3 - P[m-1]P[m+1]^3
}

void compute_formula_2_p(std::vector<Poly_p>& P, const Poly_p& Q, const mpz_t& p,  mpir_ui i) {
	Poly_p R(p);
	P[i + 1] += P[((i + 1) / 2) + 2];
	P[i + 1] *= P[((i + 1) / 2) - 1];
	P[i + 1] *= P[((i + 1) / 2) - 1];
	R += P[((i + 1) / 2) - 2];
	R *= P[((i + 1) / 2) + 1];
	R *= P[((i + 1) / 2) + 1];
	P[i + 1] -= R;
	P[i + 1] *= P[(i + 1) / 2];

	//P[i + 1] = P[(i + 1) / 2] * (P[((i + 1) / 2) + 2] * P[((i + 1) / 2) - 1] * P[((i + 1) / 2) - 1]
	//- P[((i + 1) / 2) - 2] * P[((i + 1) / 2) + 1] * P[((i + 1) / 2) + 1]);
	for (int j = 0; j < P[i + 1].size(); j++)
		mpz_divexact_ui(P[i + 1][j], P[i + 1][j], 2);
}

void compute_division_polynomials_p(int n, const ECParam& p, const std::string& name) {
	n = max(5, n);
	n += 4 - (n % 4);
	std::vector<Poly_p> P;
	Poly_p R(p.p);
	for (int i = 0; i < n+1; i++) P.push_back(Poly_p(R));

	P[0].resize(1);
	mpz_init_set_ui(P[0][0], 0);

	P[1].resize(1);
	mpz_init_set_ui(P[1][0], 1);

	P[2].resize(1);
	mpz_init_set_ui(P[2][0], 2);

	P[3].resize(5);
	mpz_init_set_ui(P[3][0], 0);
	mpz_submul(P[3][0], p.a, p.a);
	mpz_init_set_ui(P[3][1], 0);
	mpz_addmul_ui(P[3][1], p.b, 12);
	mpz_init_set_ui(P[3][2], 0);
	mpz_addmul_ui(P[3][2], p.a, 6);
	mpz_init_set_ui(P[3][3], 0);
	mpz_init_set_ui(P[3][4], 3);
	//P[3] =  -a^2 + 12bx + 6ax^2 + 3x^4

	mpz_t tmp; mpz_init_set_ui(tmp, 0);
	P[4].resize(7);
	mpz_init_set_ui(P[4][0], 0);
	mpz_submul(P[4][0], p.a, p.a);
	mpz_submul_ui(P[4][0], p.a, 4);
	mpz_submul(P[4][0], p.a, p.a);
	mpz_submul(tmp, p.b, p.b);
	mpz_submul_ui(P[4][0], tmp, 32);
	mpz_init_set_ui(P[4][1], 0);
	mpz_mul(tmp, p.a, p.b);
	mpz_submul_ui(P[4][1], tmp, 16);
	mpz_init_set_ui(P[4][2], 0);
	mpz_mul(tmp, p.a, p.a);
	mpz_submul_ui(P[4][1], tmp, 20);
	mpz_init_set_ui(P[4][3], 0);
	mpz_addmul_ui(P[4][3], p.b, 80);
	mpz_init_set_ui(P[4][4], 0);
	mpz_addmul_ui(P[4][4], p.a, 20);
	mpz_init_set_ui(P[4][5], 0);
	mpz_init_set_ui(P[4][6], 4);
	// P[4] = -4a^3-32b^2 -16abx -20a^2x^2 -80bx^3+20ax^4 +4x^6
	mpz_clear(tmp);

	Poly_p Q(p.p, 4); // Q = Y^2 = X^3 + aX + b
	mpz_set(Q[0], p.b);
	mpz_set(Q[1], p.a);
	mpz_set_ui(Q[3], 1);
	Q *= Q; // Q = Q^2
			// Now computing, 4 different formulas depending on i mod 4
	std::thread t[4];

	for (int i = 5; i < n;) {
		compute_formula_1_p(P, Q, p.p,  i);
		compute_formula_2_p(P, Q, p.p, i);
		i += 2;
		compute_formula_1_p(P, Q, p.p, i);
		compute_formula_2_p(P, Q, p.p, i);
		i += 2;
	}
	
	ofstream f(name, ios::out);
	std::vector<char> vec(1 << 20);
	f.rdbuf()->pubsetbuf(&vec.front(), vec.size());
	for (int i = 0; i <= n; i++) {
		f << "\nP_" << i << " = "; P[i].writeOut(f);
	}
	f.close();
}



/*t[0] = std::thread(compute_formula_1_p, std::ref(P), std::ref(Q), p.p, i);
t[1] = std::thread(compute_formula_2_p, std::ref(P), std::ref(Q), p.p, i);
i += 2;
t[2] = std::thread(compute_formula_1_p, std::ref(P), std::ref(Q), p.p, i);
t[3] = std::thread(compute_formula_2_p, std::ref(P), std::ref(Q), p.p, i);

t[0].join();
t[1].join();
t[2].join();
t[3].join();*/