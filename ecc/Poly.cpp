#include "Poly.h"

// div pol : huge deg jump, hence bigger alloc size
Poly::Poly(mpir_ui size) : m_size(size), m_alloc_size(NEW_ALLOC_SIZE(size)), m_coef(nullptr)
{
	m_coef = new mpz_t[m_alloc_size];
	for (int i = 0; i < m_alloc_size; i++)
		mpz_init_set_ui(m_coef[i], 0);
}

Poly::Poly(const Poly& p) : m_size(p.m_size), m_alloc_size(p.m_alloc_size), m_coef(nullptr)
{
	m_coef = new mpz_t[m_alloc_size];
	for (int i = 0; i < m_alloc_size; i++)
		mpz_init_set(m_coef[i], p[i]);
}

Poly::Poly(Poly&& p) : m_size(p.m_size), m_alloc_size(p.m_alloc_size), m_coef(p.m_coef) {

}


Poly::~Poly() {
	for (int i = 0; i < m_alloc_size; i++)
		mpz_clear(m_coef[i]);
	delete[] m_coef;
}

void Poly::resize(mpir_ui n) {
	replace(n, 0, nullptr);
}

mpir_ui Poly::size() const {
	return m_size;
}

void Poly::operator+=(const Poly& p){
	mpir_ui m = std::min(m_size, p.size());
	mpir_ui M = std::max(m_size, p.size());
	replace(M, m_size, m_coef);

	for (mpir_ui i = 0; i < m; i++)
		mpz_add(m_coef[i], m_coef[i], p[i]);
	if (p.size() == M) {
		for (mpir_ui i = m; i < M; i++)
			mpz_set(m_coef[i], p[i]);
	}
}

void Poly::operator-=(const Poly& p) {
	mpir_ui m = std::min(m_size, p.size());
	mpir_ui M = std::max(m_size, p.size());
	replace(M, m_size, m_coef);

	for (mpir_ui i = 0; i < m; i++)
		mpz_sub(m_coef[i], m_coef[i], p[i]);
	if (p.size() == M) {
		for (mpir_ui i = m; i < M; i++)
			mpz_neg(m_coef[i], p[i]);
	}
}

void Poly::operator*=(const Poly& p){
	mpir_ui m = m_size; mpir_ui M = p.size();
	mpir_ui d = m + M - 1 ;
	mpz_t* coef = new mpz_t[d];
	for (mpir_ui i = 0; i < d; i++)
		mpz_init_set_ui(coef[i], 0);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < M; j++)
			mpz_addmul(coef[i + j], m_coef[i], p[j]);
	
	replace(d, d, coef);
	for (mpir_ui i = 0; i < d; i++)
		mpz_clear(coef[i]);
	delete[] coef;
}

void Poly::operator=(Poly&& p) {
	replace(p.size(), p.size(), p.m_coef);
	/*mpir_ui m = m_size; mpir_ui M = p.size();
	if (M <= m) {
		for (int i = 0; i < m; i++)
			mpz_set(m_coef[i], p[i]);
	}
	else {
		resize(M);
		for (int i = 0; i < M; i++)
			mpz_set(m_coef[i], p[i]);
	}*/
}

void Poly::operator=(const Poly& p) {
	replace(p.size(), p.size(), p.m_coef);
	/*mpir_ui m = m_size; mpir_ui M = p.size();
	if (M <= m) {
		for (int i = 0; i < M; i++)
			mpz_set(m_coef[i], p[i]);
		for (int i = M; i < m; i++)
			mpz_set_ui(m_coef[i], 0);
	}
	else {
		resize(M);
		for (int i = 0; i < M; i++)
			mpz_set(m_coef[i], p[i]);
	}*/
}

mpz_t& Poly::operator[](mpir_ui i) { return m_coef[i]; }
const mpz_t& Poly::operator[](mpir_ui i) const { return m_coef[i]; }

void Poly::print() const{
	char str[16384];
	mpz_get_str(str, BASE_n, m_coef[0]);
	printf("%s", str);
	for (int i = 1; i < m_size; i++) {
		mpz_get_str(str, BASE_n, m_coef[i]);
		if (mpz_sgn(m_coef[i]) >= 0)
			printf(" + %s X^%d", str, i);
		else
			printf(" %sX^%d", str, i);
	}
	printf("\n");
}

void Poly::writeOut(std::ofstream& out) {
	char str[1 << 18];
	mpz_get_str(str, 16, m_coef[0]);
	out << str;
	for (int i = 1; i < m_size; i++) {
		mpz_get_str(str, 16, m_coef[i]);
		if (mpz_sgn(m_coef[i]) >= 0)
			out << " + " << str << "X^" << i;
		else
			out << " " << str << "X^" << i;
	}
	out << "\n";
}


void Poly::replace(mpir_ui size, mpir_ui nb_coeff, mpz_t* coef) {
	if (size > m_alloc_size) {
		mpir_ui m_alloc_size_new = NEW_ALLOC_SIZE(size);
		mpz_t* m_coef_new = new mpz_t[m_alloc_size_new];
		//printf("Reallocation with : size = %d, m_alloc_size = %d, m_alloc_size_new = %d, coef = %p\n", size, m_alloc_size, m_alloc_size_new, coef);
		if (coef) {
			for (mpir_ui i = 0; i < nb_coeff; i++)
				mpz_init_set(m_coef_new[i], coef[i]);
			for (mpir_ui i = nb_coeff; i < size; i++)
				mpz_init_set_ui(m_coef_new[i], 0);
		}
		else {
			for (mpir_ui i = 0; i < size; i++)
				mpz_init_set_ui(m_coef_new[i], 0);
		}
		for (mpir_ui i = size; i < m_alloc_size_new; i++)
			mpz_init_set_ui(m_coef_new[i], 0);

		for (mpir_ui i = 0; i < m_alloc_size; i++)
			mpz_clear(m_coef[i]);
		delete[] m_coef;
		m_coef = m_coef_new;
		m_size = size;
		m_alloc_size = m_alloc_size_new;
	}
	else {
		if (coef && (m_coef != coef)) {
			for (mpir_ui i = 0; i < size; i++)
				mpz_set(m_coef[i], coef[i]);
			for (mpir_ui i = size; i < m_size; i++)
				mpz_set_ui(m_coef[i], 0);
		}
		else if (!coef) { // reset
			for (mpir_ui i = 0; i < m_size; i++)
				mpz_set_ui(m_coef[i], 0);
		}
		else {
			for (mpir_ui i = size; i < m_size; i++)
				mpz_set_ui(m_coef[i], 0);
		}
		m_size = size;
	}
}




Poly_p::Poly_p(const mpz_t& p, mpir_ui size) : m_size(size), m_alloc_size(NEW_ALLOC_SIZE(size)), m_coef(nullptr)
{
	mpz_init_set(m_p, p);
	m_coef = new mpz_t[m_alloc_size];
	for (int i = 0; i < m_alloc_size; i++)
		mpz_init_set_ui(m_coef[i], 0);
}

Poly_p::Poly_p(const Poly_p& p) : m_size(p.m_size), m_alloc_size(p.m_alloc_size), m_coef(nullptr)
{
	mpz_init_set(m_p, p.m_p);
	m_coef = new mpz_t[m_alloc_size];
	for (int i = 0; i < m_alloc_size; i++)
		mpz_init_set(m_coef[i], p[i]);
}

/*Poly_p::Poly_p(Poly_p&& p) : m_size(p.m_size), m_alloc_size(p.m_alloc_size), m_coef(p.m_coef) {
	mpz_init_set(m_p, p.m_p);
	mpz_clear(p.m_p);
}*/


Poly_p::~Poly_p() {
	for (int i = 0; i < m_alloc_size; i++)
		mpz_clear(m_coef[i]);
	delete[] m_coef;
	mpz_clear(m_p);
}

void Poly_p::resize(mpir_ui n) {
	replace(n, 0, nullptr);
}

mpir_ui Poly_p::size() const {
	return m_size;
}

void Poly_p::operator+=(const Poly_p& p) {
	mpir_ui m = std::min(m_size, p.size());
	mpir_ui M = std::max(m_size, p.size());
	replace(M, m_size, m_coef);

	for (mpir_ui i = 0; i < m; i++) {
		mpz_add(m_coef[i], m_coef[i], p[i]);
		mpz_mod(m_coef[i], m_coef[i], m_p);
	}
	if (p.size() == M) {
		for (mpir_ui i = m; i < M; i++)
			mpz_set(m_coef[i], p[i]);
	}
}

void Poly_p::operator-=(const Poly_p& p) {
	mpir_ui m = std::min(m_size, p.size());
	mpir_ui M = std::max(m_size, p.size());
	replace(M, m_size, m_coef);

	for (mpir_ui i = 0; i < m; i++) {
		mpz_sub(m_coef[i], m_coef[i], p[i]);
		mpz_mod(m_coef[i], m_coef[i], m_p);
	}
	if (p.size() == M) {
		for (mpir_ui i = m; i < M; i++) {
			mpz_neg(m_coef[i], p[i]);
			mpz_mod(m_coef[i], m_coef[i], m_p);
		}
	}
}

void Poly_p::operator<<=(mpir_ui n) {
	mpz_t* coef_new = new mpz_t[m_alloc_size + n];
	for (mpir_ui i = 0; i < n; i++)
		mpz_init_set_ui(coef_new[i],0);
	for (mpir_ui i = 0; i < m_alloc_size; i++)
		mpz_init_set(coef_new[i + n],m_coef[i]);
	for (mpir_ui i = 0; i < m_alloc_size; i++)
		mpz_clear(m_coef[i]);
	delete[] m_coef;
	m_coef = coef_new;
	m_size += n;
	m_alloc_size += n;
}


void Poly_p::operator*=(const Poly_p& p) {
	if (m_size >= 500 || p.size() >= 500) {
		Poly_p A(m_p), B(m_p), C(m_p), D(m_p);
		mpir_ui deg = m_size / 2;
		cut(A, B, deg); //*this = A; // HEY HEY
		if (!p.cut(C, D, deg)) { // C = 0
			//printf("\n");
			//print(); p.print();
			//printf("m_size = %ld, deg = %ld\n", m_size, deg);
			//A.print(); B.print(); C.print(); D.print();
			A *= D; B *= D;
			*this = A;
			*this <<= deg;
			*this += B;
			return;
		};
		//printf("\n");
		//print(); p.print();
		//printf("m_size = %ld, deg = %ld\n", m_size, deg);
		//A.print(); B.print(); C.print(); D.print();
		Poly_p tmp_1(A);
		Poly_p tmp_2(C);
		tmp_1 -= B;
		tmp_2 -= D;
		A *= C;
		B *= D;
		tmp_2 *= tmp_1;
		tmp_1 = A;
		tmp_1 += B;
		tmp_1 -= tmp_2;
		// *this = AC, B = BD, tmp_2 = (A - B)(C - D), tmp_1 = AC + BD - (A-B)(C-D) = AD + BC
		A <<= m_size;
		tmp_1 <<= deg;
		A += tmp_1;
		A += B;
		*this = A;
		return;
	}
	
	mpir_ui m = m_size; mpir_ui M = p.size();
	mpir_ui d = m + M - 1;
	mpz_t* coef = new mpz_t[d];
	for (mpir_ui i = 0; i < d; i++)
		mpz_init_set_ui(coef[i], 0);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < M; j++) {
			mpz_addmul(coef[i + j], m_coef[i], p[j]);
			mpz_mod(coef[i + j], coef[i + j], m_p);
		}

	replace(d, d, coef);
	for (mpir_ui i = 0; i < d; i++)
		mpz_clear(coef[i]);
	delete[] coef;
}

void Poly_p::operator=(Poly_p&& p)					{replace(p.size(), p.size(), p.m_coef);}
void Poly_p::operator=(const Poly_p& p)				{replace(p.size(), p.size(), p.m_coef);}
mpz_t& Poly_p::operator[](mpir_ui i)				{ return m_coef[i]; }
const mpz_t& Poly_p::operator[](mpir_ui i) const	{ return m_coef[i]; }

void Poly_p::print() const {
	char str[16384];
	mpz_get_str(str, BASE_n, m_coef[0]);
	printf("%s", str);
	for (int i = 1; i < m_size; i++) {
		mpz_get_str(str, BASE_n, m_coef[i]);
		if (mpz_sgn(m_coef[i]) >= 0)
			printf(" + %s X^%d", str, i);
		else
			printf(" %sX^%d", str, i);
	}
	printf("\n");
}

void Poly_p::writeOut(std::ofstream& out) {
	char str[1 << 18];
	mpz_get_str(str, 10, m_coef[0]);
	out << str;
	for (int i = 1; i < m_size; i++) {
		mpz_get_str(str, 10, m_coef[i]);
		if (mpz_sgn(m_coef[i]) >= 0)
			out << " + " << str << "X^" << i;
		else
			out << " " << str << "X^" << i;
	}
	out << "\n";
}

void Poly_p::replace(mpir_ui size, mpir_ui nb_coeff, mpz_t* coef) {
	if (size > m_alloc_size) {
		mpir_ui m_alloc_size_new = NEW_ALLOC_SIZE(size);
		mpz_t* m_coef_new = new mpz_t[m_alloc_size_new];
		//printf("Reallocation with : size = %d, m_alloc_size = %d, m_alloc_size_new = %d, coef = %p\n", size, m_alloc_size, m_alloc_size_new, coef);
		if (coef) {
			for (mpir_ui i = 0; i < nb_coeff; i++)
				mpz_init_set(m_coef_new[i], coef[i]);
			for (mpir_ui i = nb_coeff; i < size; i++)
				mpz_init_set_ui(m_coef_new[i], 0);
		}
		else {
			for (mpir_ui i = 0; i < size; i++)
				mpz_init_set_ui(m_coef_new[i], 0);
		}
		for (mpir_ui i = size; i < m_alloc_size_new; i++)
			mpz_init_set_ui(m_coef_new[i], 0);

		for (mpir_ui i = 0; i < m_alloc_size; i++)
			mpz_clear(m_coef[i]);
		delete[] m_coef;
		m_coef = m_coef_new;
		m_size = size;
		m_alloc_size = m_alloc_size_new;
	}
	else {
		if (coef && (m_coef != coef)) {
			for (mpir_ui i = 0; i < size; i++)
				mpz_set(m_coef[i], coef[i]);
			for (mpir_ui i = size; i < m_size; i++)
				mpz_set_ui(m_coef[i], 0);
		}
		else if (!coef) { // reset
			for (mpir_ui i = 0; i < m_size; i++)
				mpz_set_ui(m_coef[i], 0);
		}
		else {
			for (mpir_ui i = size; i < m_size; i++)
				mpz_set_ui(m_coef[i], 0);
		}
		m_size = size;
	}
}

bool Poly_p::cut(Poly_p& A, Poly_p& B, mpir_ui deg) const {
	signed long long size_a = m_size - deg;
	if (size_a <= 0) {
		A.replace(1);
		B = *this;
		return false;
	}
	A.replace(size_a, size_a, &m_coef[deg]);
	B.replace(deg, deg, m_coef);
	return true;
}
