/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace cel�ch ��sel z pohledu l�m�n� RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    utils.c                                                    **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

#include <cstdio>
#include <gmpxx.h>
#include <iostream>

#include "utils.h"

/* TODO: Dokoncit!!! */
progData_t *ParseArgs(int argCnt, char **args) {
	return NULL;
}

/* TODO: Dokoncit!!! */
unsigned char *CheckNumber(char *str) {
	return NULL;
}

/***
*/
uint16_t DecNumToBinNum(uint8_t *decNum, uint16_t digitCnt, uint64_t *binNum) {
	// Bude fungovat tak, ze se vezme pole dekadickych cislic a od nejvyssi
	// cislice po nejmensi se bude delit 2. Pokud bude zbytek po deleni 1, pak
	// se musi k nasledujici cislici pricis 10. Zbytek zaroven vyjadruje hodnotu
	// bitu, ktery se musi ulozit
	return 0;
}

/** This function represents polynom that is used by Pollard Rho Method
 *  Most common polynom is g(x) = x^2 + 1
 *
 *  Function is private
 *
 * @param x Input
 * @param N Z_N
 * @return  Result
 */
void Pollard_polynom(mpz_t x, mpz_t N) {

	// x^2
	mpz_mul(x, x, x);

	// +1
	mpz_add_ui(x, x, 1);

	// mod N
	mpz_mod(x,x,N);
}

/** Factorization method: Pollar Rho Method
 * This method is used for smaller numbers.
 *
 * @param N Number to factor
 * @return  Result
 */
void Pollard_rho_method_GMP(mpz_t res, mpz_t N) {
	mpz_t b_i;
	mpz_t b_2i;

	using namespace std;

#ifdef UTILS_DEBUG_PRINT_POLLARD_RHO
	cout << "Going to print Pollard Rho method debug msgs" << endl;
#endif

	// Initialize the number
	mpz_init(b_i);
	mpz_init(b_2i);
	mpz_set_ui(b_i, 2);
	mpz_set_ui(b_2i, 2);

	do {
		Pollard_polynom(b_i, N);
		Pollard_polynom(b_2i, N);
		Pollard_polynom(b_2i, N);

#ifdef UTILS_DEBUG_PRINT_POLLARD_RHO
		cout << " b_i = ";
		mpz_out_str(stdout, 10, b_i);
		cout << endl;

		cout << " b_2i = ";
		mpz_out_str(stdout, 10, b_2i);
		cout << endl;
#endif

		mpz_sub(res, b_i, b_2i);

#ifdef UTILS_DEBUG_PRINT_POLLARD_RHO
		cout << " b_i - b_2i = ";
		mpz_out_str(stdout, 10, res);
		cout << endl;
#endif

		mpz_abs(res,res);

#ifdef UTILS_DEBUG_PRINT_POLLARD_RHO
		cout << " abs = ";
		mpz_out_str(stdout, 10, res);
		cout << endl;
#endif

		mpz_gcd(res, res, N);

#ifdef UTILS_DEBUG_PRINT_POLLARD_RHO
		cout << " gcd(res, N) = ";
		mpz_out_str(stdout, 10, res);
		cout << endl;
#endif
	} while(mpz_cmp_ui(res,1) <= 0);
}
