#include "ectypes.h"

void mpz_print(const mpz_t print) {
	mpz_out_str(stdout, 10, print);
	printf("\n");
}

ECParam::ECParam(mpir_ui _a, mpir_ui _b , mpir_ui _p, mpir_ui _order) {
	mpz_init_set_ui(a, _a);
	mpz_init_set_ui(b, _b);
	mpz_init_set_ui(p, _p);
	mpz_init_set_ui(order, _order);
}

ECParam::ECParam(const char* _a, const char* _b, const char* _p, const char* _order) {
	mpz_init_set_str(a, _a, 0);
	mpz_init_set_str(b, _b, 0);
	mpz_init_set_str(p, _p, 0);
	mpz_init_set_str(order, _order, 0);
}

ECParam::ECParam(mpz_t _a, mpz_t _b, mpz_t _p, mpz_t _order) {
	mpz_init_set(a, _a);
	mpz_init_set(b, _b);
	mpz_init_set(p, _p);
	mpz_init_set(order, _order);
}

ECParam::~ECParam() {
	mpz_clears(a, b, p, order, NULL);
}

ECParam::ECParam(const ECParam& param) {
	mpz_init_set(a, param.a);
	mpz_init_set(b, param.b);
	mpz_init_set(p, param.p);
	mpz_init_set(order, param.order);
}

void ECParam::operator=(const ECParam& param) {
	mpz_set(a, param.a);
	mpz_set(b, param.b);
	mpz_set(p, param.p);
	mpz_set(order, param.order);
}

bool ECParam::operator==(const ECParam& param) const {
	return	(mpz_cmp(a, param.a) == 0)
		&&	(mpz_cmp(b, param.b) == 0)
		&&	(mpz_cmp(p, param.p) == 0)
		&&	(mpz_cmp(order, param.order) == 0);
}

ECCoord::ECCoord(int isInf, const char* x, const char* y) : isInf(isInf) {
	mpz_init_set_str(this->x, x, 0);
	mpz_init_set_str(this->y, y, 0);

}

ECCoord::ECCoord(int _isInf, mpz_t _x, mpz_t _y) {
	isInf = _isInf;
	mpz_init_set(x, _x);
	mpz_init_set(y, _y);
}

ECCoord::~ECCoord() {
	mpz_clear(x);
	mpz_clear(y);
}

ECCoord::ECCoord(const ECCoord& p) : isInf(p.isInf) {
	mpz_init_set(x, p.x);
	mpz_init_set(y, p.y);
}

void ECCoord::operator=(const ECCoord& p) {
	isInf = p.isInf;
	mpz_set(x, p.x);
	mpz_set(y, p.y);
}

bool ECCoord::operator==(const ECCoord& p) const {
	// We have to be carefull with isInf because it is an integer and not a boolean
	if (isInf && p.isInf) return true;
	if (isInf ^ p.isInf) return false;
	return (mpz_cmp(x, p.x) == 0) && (mpz_cmp(y, p.y) == 0);
}