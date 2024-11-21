#include "mimalloc-new-delete.h"
#include "dicer.h"
#include <iostream>
#include <mutex>
#include <ctime>
#include "omp.h"
#include <complex>

using namespace spaceless;

static void PolynomialTest() {
	{
		polynomial<DYNAMIC> a, b;
		a += polynomial(monomial<DYNAMIC>({ 2 }), 1);
		a += polynomial(monomial<DYNAMIC>({ 2, 1 }), 2);
		a += polynomial(monomial<DYNAMIC>({ 1, 1 }), -2);
		a += polynomial(monomial<DYNAMIC>({ 1, 2 }), 2);
		a += polynomial(monomial<DYNAMIC>({ 0, 2 }), -3);

		b += polynomial(monomial<DYNAMIC>({ 1 }), 1);
		b += polynomial(monomial<DYNAMIC>({ 0, 1 }), 1);

		auto c = a / b, d = a % b;
		std::cout << a << '\n';
		std::cout << b << '\n';
		std::cout << c << '\n';
		std::cout << d << '\n';
		auto e = c * b;
		std::cout << e << '\n';
		auto f = e + d;
		std::cout << f << '\n';
		auto g = f - a;
		g.trim();
		std::cout << g << '\n';
	}
	{
		polynomial<DYNAMIC> a, b;
		a += polynomial(monomial<DYNAMIC>({ 3 }), 1);
		a += polynomial(monomial<DYNAMIC>({ 2 }), -2);
		a -= 4;

		b += polynomial(monomial<DYNAMIC>({ 1 }), 1);
		b -= 3;

		auto c = a / b, d = a % b;
		std::cout << a << '\n';
		std::cout << b << '\n';
		std::cout << c << '\n';
		std::cout << d << '\n';
		auto e = c * b;
		std::cout << e << '\n';
		auto f = e + d;
		std::cout << f << '\n';
		auto g = f - a;
		std::cout << g << '\n';
	}
}

int main(int argc, char **argv) {
	int begin = clock(), end = 0;

	PolynomialTest();

	end = clock();
	std::cout << "cost time : " << (end - begin) / 1000.0 << "s" << std::endl;
	return 0;
}
