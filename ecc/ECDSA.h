#pragma once
#include "EllCurve.h"
#include "Constants.h"

/*
The scenario is the following : Alice wants to sign a message 
with her private key (dA), and Bob wants to validate the signature
using Alice's public key (HA = dA*G). Nobody but Alice should be able to
produce valid signatures. Everyone should be able to check 
signatures. Some pre-requisite : 
- Order of the subgroup must be prime (<G> is a subgroup of prime order ok)
- The hash of the message to be signed should be the same bit length as n
- The random k used to generate signature must be changed for each signature (c.f Sony hack)*/

/*Signs messages given : 
- Publicly knwon elliptic curve 
- priv_key the private key of the sender (that is dA defined earlier)
- Hash hash of the message to be sent
Returns (r, s), r the x-coord of P = kG with random k 
and s the variable where are "melted" r, the message hash, priv_key and k*/
void SignMessage(mpz_t& r, mpz_t& s, 
	EllCurve& curve, const mpz_t& priv_key, const mpz_t& hash) {
	printf("Signing message of hash "); mpz_print(hash);
	curve.print("ECDSA Curve", true);
	printf("Private key is "); mpz_print(priv_key);
	printf("Corresponding public key is ");
	mpz_t k; mpz_init(k);
	EllPoint P(curve.getGen());
	do {
		do { P = curve.getRandomPoint(k); } while (P.isInf());//mpz_sgn(P.getCoord().x) == 0);
		mpz_mod(r, P.getCoord().x, curve.getGenOrder()); // r ok
		mpz_set(s, hash);
		mpz_addmul(s, r, priv_key);
		mpz_invert(k, k, curve.getGenOrder());
		mpz_mul(s, s, k);
		mpz_mod(s, s, curve.getGenOrder()); // s ok
		mpz_print(curve.getGenOrder());
	} while (mpz_sgn(s) == 0);
	mpz_clear(k);

	printf("Signature (r, s) is : \n\tr = "); mpz_print(r);
	printf("\ts = "); mpz_print(s);
}

bool CheckSignature(const mpz_t& r, const mpz_t& s, 
	const EllCurve& curve, const EllPoint& publicKey, mpz_t hash) {
	if (!(curve.getECParam() == publicKey.getECParam())) {
		printf("(ECDSA)Point given isn't defined over the same curve\n"); return false;
	}
	if (!publicKey.isOnCurve(publicKey.getCoord())) {
		printf("(ECDSA)Point given is not on curve\n"); return false;
	}
	if (mpz_sgn(r) == 0 || mpz_sgn(s) == 0) {
		printf("(ECDSA)Signature given is invalid (at least one parameter is null\n"); return false;
	}
	mpz_t u, v; mpz_inits(u, v, NULL);
	mpz_invert(u, s, curve.getGenOrder());
	mpz_mul(v, u, r);
	mpz_mod(v, v, curve.getGenOrder()); // v ok
	mpz_mul(u, u, hash);
	mpz_mod(u, u, curve.getGenOrder()); // u ok

	EllPoint P(curve.getGen()), tmp(publicKey);
	tmp *= v; P *= u; P += tmp;

	mpz_sub(u, r, P.getCoord().x); // using u as tmp variable
	mpz_mod(u, u, curve.getGenOrder());
	if (mpz_sgn(u) == 0) {
		printf("(ECDSA)Signature is valid.\n"); mpz_clears(u, v, NULL); return true;
	}
	printf("(ECDSA)Signature is invalid.\n"); mpz_clears(u, v, NULL); return false;
}



