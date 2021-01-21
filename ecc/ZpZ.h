#pragma once
//#include <mpirxx.h>
//#include <cstdio>

//#define FAIL printf("Trying to do operation on element of distinct field. Exiting..."); exit(EXIT_FAILURE)

// We'll consider that p = 0 yields a ZpZ object comparable to all others ZpZ objects.
/*
class ZpZ {
public :
	template <class T = int, class S = int>
	ZpZ(T value = 0, S p = 5) : m_value(value), m_p(p)
	{
		m_value = m_value % m_p;
	}

	ZpZ(ZpZ&& x) : m_value(x.m_value), m_p(x.m_p)
	{}
	
	ZpZ&& operator+(const ZpZ& x) const {
		//if (m_p != x.m_p) FAIL;
		return ZpZ((x.m_value + m_value) % m_p, m_p);
	}

	ZpZ&& operator-(const ZpZ& x) const {
		//if (m_p != x.m_p) FAIL;
		return ZpZ((x.m_value - m_value) % m_p, m_p);
	}

	ZpZ&& operator=(const ZpZ& x) {
		m_value = x.m_value;
		m_p = x.m_p;
	}

	ZpZ&& operator*(const ZpZ& x) {
		//if (m_p != x.m_p) FAIL;
		return ZpZ((x.m_value * m_value) % m_p, m_p);
	}

	ZpZ&& operator/(const ZpZ& x) {
		//if (m_p != x.m_p) FAIL;
		mpz_class tmp;
		mpz_invert(tmp.get_mpz_t(), x.m_value.get_mpz_t(), m_p.get_mpz_t());
		return ZpZ((m_value * tmp) % m_p, m_p);
	}

private :
	mpz_class m_value;
	mpz_class m_p;
};
*/