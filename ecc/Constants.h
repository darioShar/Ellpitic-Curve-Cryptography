#pragma once
#include "EllCurve.h"

#define BIG_A "3626830188090040170470538223632857507340287510495561807368105032792956752466460499091754944292670519"
#define BIG_B "6504080042008251533437753757878856267866982998831103001251197043718345565337995842176672637937899949"
#define BIG_P "35265075282301222022179652071920843985804547105655617470846478149269603109205365857799131187032150 \
238941117870812985999909923480249871256893147030937370668858201140634615843266527731318894294518024331"

#define DEFAULT_A "11"
#define DEFAULT_B "17"
#define DEFAULT_P "1000003"

// Bitcoin chain parameters
#define BITCOIN_A "0"
#define BITCOIN_B "7"
// P to be defined later
// Gen Point
#define BITCOIN_Gx "0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798"
#define BITCOIN_Gy "0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8"
// Number of points in the field
#define BITCOIN_N "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141"

// secp112r1
#define secp112_A  "0xDB7C2ABF62E35E668076BEAD2088"
#define secp112_B  "0x659EF8BA043916EEDE8911702B22"
#define secp112_Gx "0x09487239995A5EE76B55F9C2F098"
#define secp112_Gy "0xA89CE5AF8724C0A23E0E0FF77500"
#define secp112_P  "0xDB7C2ABF62E35E668076BEAD208B"
#define secp112_n  "0xDB7C2ABF62E35E7628DFAC6561C5"

// secp521r1
#define secp521_A  "0x01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC"
#define secp521_B  "0x0051953EB9618E1C9A1F929A21A0B68540EEA2DA725B99B315F3B8B489918EF109E156193951EC7E937B1652C0BD3BB1BF073573DF883D2C34F1EF451FD46B503F00"
#define secp521_Gx "0x00C6858E06B70404E9CD9E3ECB662395B4429C648139053FB521F828AF606B4D3DBAA14B5E77EFE75928FE1DC127A2FFA8DE3348B3C1856A429BF97E7E31C2E5BD66"
#define secp521_Gy "0x011839296A789A3BC0045C8A5FB42C7D1BD998F54449579B446817AFBD17273E662C97EE72995EF42640C550B9013FAD0761353C7086A272C24088BE94769FD16650"
#define secp521_P  "0x01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
#define secp521_n  "0x01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA51868783BF2F966B7FCC0148F709A5D03BB5C9B8899C47AEBB6FB71E91386409"


/* Initializing with bitcoin chain parameters */
static void curve_init_bitcoin(EllCurve& curve) {
	mpz_t p; mpz_init(p);
	mpz_t op1;
	mpz_set_ui(p, 1);
	mpz_mul_2exp(p, p, 256);
	mpz_init_set_ui(op1, 1);
	mpz_mul_2exp(op1, op1, 4); mpz_sub(p, p, op1);
	mpz_mul_2exp(op1, op1, 2); mpz_sub(p, p, op1);
	mpz_mul_2exp(op1, op1, 1); mpz_sub(p, p, op1);
	mpz_mul_2exp(op1, op1, 1); mpz_sub(p, p, op1);
	mpz_mul_2exp(op1, op1, 1); mpz_sub(p, p, op1);
	mpz_mul_2exp(op1, op1, 23); mpz_sub(p, p, op1);
	mpz_sub_ui(p, p, 1);
	// p == 2^256 - 2^32 - 2^9 - 2^8 - 2^7 - 2^6 - 2^4 - 1

	char* s = mpz_get_str(NULL, 10, p);
	curve.setECParam(ECParam(BITCOIN_A, BITCOIN_B, s, BITCOIN_N));

	if (!curve.setGen(ECCoord(NOT_INFTY, BITCOIN_Gx, BITCOIN_Gy), curve.getCurveOrder())) {
		printf("Generator given is not on curve.\n");
		curve.findNewGen();
		return;
	}
	mpz_clears(p, op1, NULL);
}

/* Initializing with some parameters and a generator point*/
static void curve_init_default(EllCurve& curve) {
	ECParam p("968113241544", "19232197347131", "23640121541677", "23640123950704");	
	curve.setECParam(p);
	curve.findNewGen();
}

/* Initializing with some Big parameters */
static void curve_init_bignum(EllCurve& curve) {
	curve.setECParam(ECParam(BIG_A, BIG_B, BIG_P, "0"));
	curve.findNewGen();
}

static void curve_init_secp112r1(EllCurve& curve) {
	curve.setECParam(ECParam(secp112_A, secp112_B, secp112_P, secp112_n));
	if (!curve.setGen(ECCoord(NOT_INFTY, secp112_Gx, secp112_Gy), curve.getCurveOrder())) {
		printf("Generator given is not on curve.\n");
		curve.findNewGen(); // oupsi
		return;
	}
}

static void curve_init_secp521r1(EllCurve& curve) {
	curve.setECParam(ECParam(secp521_A, secp521_B, secp521_P, secp521_n));
	if (!curve.setGen(ECCoord(NOT_INFTY, secp521_Gx, secp521_Gy), curve.getCurveOrder())) {
		printf("Generator given is not on curve.\n");
		curve.findNewGen(); // oupsi
		return;
	}
}


// Some other curve i've found
#define ds1a "118916092566490565748229184"
#define ds1b "66274645112519462583913069"
#define ds1p "2475880078570760549798248507"
#define ds1n "2475880078570690266226895996"

#define ds2a "2972362603913098775400772217"
#define ds2b "5191751986946129553737567891"
#define ds2p "9903520314283042199192993897"
#define ds2n "9903520314282870095649001203"