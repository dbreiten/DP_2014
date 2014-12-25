/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace cel�ch ��sel z pohledu l�m�n� RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    main.c                                                     **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

/*********************************************************************************
** TODO List
**
** SIQS
** Single Large Prime Variation
** Lanczos matrix
** Pro test B-smooth detekovat kandidaty na zaklade velkych prvocisel
** z faktorizacni baze. Tim se vyeliminuje mnoho relaci a na tech, co zbydou,
** pouzit bernsteina (primes 139)
**
*********************************************************************************/

// Standard libraries
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cassert>

#include <string>

// Implemented libraries
#include "siqs.h"
#include "utils.h"

// Big Numbers library
#include <gmpxx.h>

using namespace std;

int main(int argc, char **argv) {
	progData_t    *progData;
	bool primalityCheck = true; // TODO: udelat jako argument
	bool verbose        = true; // TODO: udelat jako argument

	uint8_t primalityResult;

	string inputStr;


	// mpz_t is the type defined for GMP integers.
	// It is a pointer to the internals of the GMP integer data structure
	mpz_t n;
	mpz_t f1;     // Only for tests
	mpz_t f2;     // Only for tests
	mpz_t fake_n; // Only for tests
	mpz_t res;
	int flag;

	cout << "Enter your number: " << endl;
	cin >> inputStr;

	// Check how long the inputed number is
	if(inputStr.length() > MAX_DIGITS) {
		cerr << "Input is too long! (" << inputStr.length() << ") digits" << endl;
		return INPUT_IS_TOO_LONG;
	}

	// Initialize the number
	mpz_init(n);
	mpz_init(f1);
	mpz_init(f2);
	mpz_init(fake_n);
	mpz_init(res);
	mpz_set_ui(n, 0);
	mpz_set_ui(f1, 961748941);
	mpz_set_ui(f2, 982451653);
	mpz_set_ui(fake_n, 0);
	mpz_set_ui(res, 0);

	// Parse the input string as a base 10 number
	flag = mpz_set_str(n, inputStr.c_str(), 10);
	assert(flag == 0); /* If flag is not 0 then the operation failed */

	// Print n
	cout << "N = ";
	mpz_out_str(stdout, 10, n);
	cout << endl;

	// Primality check
	if(primalityCheck) {
		primalityResult = mpz_probab_prime_p(n,25);
		if(primalityResult > 0) {
			cerr << "N is prime => Nothing to factor" << endl;
			return INPUT_IS_PRIME;
		}
	}

	// If the input is too small, factor with Pollard rho method else SIQS
	if(inputStr.length() <= 30) {
		Pollard_rho_method_GMP(res, n);
	} else {
		;
	}

	// Print result
	cout << "Result = ";
	mpz_out_str(stdout, 10, res);
	cout << endl;

	/* Clean up the mpz_t handles or else we will leak memory */
	mpz_clear(n);
	mpz_clear(f1);
	mpz_clear(f2);
	mpz_clear(fake_n);
	mpz_clear(res);

	return OK;
}
