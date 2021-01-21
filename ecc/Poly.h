#pragma once
#include "EllPoint.h"
#include <vector>
#include <fstream>
#include <iostream>

#define NEW_ALLOC_SIZE(x) (1.5*x) // soft
//#define NEW_ALLOC_SIZE(x) (1.5*x*x + 7*x + 1) // HARD

#define BASE_n 10

struct Poly {
	Poly(mpir_ui size = 1);
	Poly(const Poly& p);
	Poly(Poly&& p);

	~Poly();

	void resize(mpir_ui n);
	mpir_ui size() const;

	void operator+=(const Poly& p);
	void operator-=(const Poly& p);
	void operator*=(const Poly& p);
	void operator=(Poly&& p);
	void operator=(const Poly& p);
	mpz_t& operator[](mpir_ui i);
	const mpz_t& operator[](mpir_ui i) const;

	void print() const;
	void writeOut(std::ofstream& out);

private :
	// nb_coef must be <= size.
	void replace(mpir_ui size, mpir_ui nb_coef = 0, mpz_t* coef = nullptr);

	mpir_ui		m_size;
	mpz_t*		m_coef;
	mpir_ui		m_alloc_size;
};

struct Poly_p {
	Poly_p(const mpz_t& p, mpir_ui size = 1);
	Poly_p(const Poly_p& p);
	~Poly_p();

	void resize(mpir_ui n);
	mpir_ui size() const;

	void operator+=(const Poly_p& p);
	void operator-=(const Poly_p& p);
	void operator*=(const Poly_p& p);
	void operator<<=(mpir_ui n);
	void operator=(Poly_p&& p);
	void operator=(const Poly_p& p);
	mpz_t& operator[](mpir_ui i);
	const mpz_t& operator[](mpir_ui i) const;

	void print() const;
	void writeOut(std::ofstream& out);

	// returns false fi there has been no cut
	bool cut(Poly_p& A, Poly_p& B, mpir_ui deg) const;
private:
	// nb_coef must be <= size.
	void replace(mpir_ui size, mpir_ui nb_coef = 0, mpz_t* coef = nullptr);

	mpir_ui		m_size;
	mpz_t*		m_coef;
	mpir_ui		m_alloc_size;
	mpz_t		m_p;
};


