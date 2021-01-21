#pragma once
#include "EllCurve.h"
#include "Constants.h"


/*Elliptic Curve Diffie-Hellman lets two users securely generate a secret key
over an insecure channel. The private key obtained can be used to encrypt/decrypt
a file by both parties with AES for instance.*/
void ECDH(EllCurve& curve) {
	// Publicly known parameters
	printf("Public parameters used for ECDH : \n");
	curve.print("ECDH", true);

	//Alice
	mpz_t k_alice;
	mpz_init(k_alice);
	EllPoint P_Alice = curve.getRandomPoint(k_alice);
	printf("Alice private key : "); mpz_print(k_alice);
	
	//Bob
	mpz_t k_bob;
	mpz_init(k_bob);
	EllPoint P_Bob = curve.getRandomPoint(k_bob);
	printf("Bob private key : "); mpz_print(k_bob);

	// Sent over insecure channel
	printf("Alice sends to Bob : \n");
	P_Alice.print("Alice");
	printf("Bob sends to Alice : \n");
	P_Bob.print("Bob");

	// Private key for Alice, received P_Bob
	EllPoint P_Alice_private(P_Bob);
	P_Alice_private *= k_alice;

	// Private key for Bob, received P_Alice
	EllPoint P_Bob_private(P_Alice);
	P_Bob_private *= k_bob;

	// Shared private key is now the same for both Alice and Bob
	if (!(P_Alice_private == P_Bob_private)) {
		// Should never happen
		printf("Error : private keys don't match !\n");
		return;
	}

	P_Alice.print("Shared private key");

	// Now private key can be used for securing communication between
	// and bob, for instance using the x coordinate stripped of its first
	// 16 bytes to encrypt/ decrypt file using AES/DES etc.

	mpz_clears(k_alice, k_bob, NULL);
}