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
** C++ !!!!!!!!!!!!!!!!
** SIQS
** Single Large Prime Variation
** Lanczos matrix
** Test prvociselnosti
** Pro mala pouzit jinou faktorizaci
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
	uint16_t      digitCnt;
	bool primalityCheck = true; // TODO: udelat jako argument

	uint8_t primalityResult;
	//result = FermatPrimalityTest(15485761);

	//uint64_t result = Euklid_algorithm(50, 15);

	// Ziskani nastaveni zadane uzivatelem
	/* TODO: Dokoncit!!! */

	//progData = ParseArgs(argc, argv);

	//Assuming input is correct
	//digitCnt = strlen(progData->inputNumStr);

	/******************************************************************************
	** GMP                                                                       **
	******************************************************************************/

	//char inputStr[1025];
	string inputStr;


	// mpz_t is the type defined for GMP integers.
	// It is a pointer to the internals of the GMP integer data structure
	mpz_t n;
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
	mpz_set_ui(n, 0);

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
		;
	} else {
		;
	}

	/* 3. Add one to the number */

	mpz_add_ui(n, n, 1); /* n = n + 1 */

	/* 4. Print the result */

	cout << " n +1 = ";
	mpz_out_str(stdout, 10, n);
	cout << endl;


	/* 5. Square n+1 */

	mpz_mul(n, n, n); /* n = n * n */


	cout << " (n +1)^2 = ";
	mpz_out_str(stdout, 10, n);
	cout << "\n";


	/* 6. Clean up the mpz_t handles or else we will leak memory */
	mpz_clear(n);

	/******************************************************************************
	** GMP                                                                       **
	******************************************************************************/

	return 0;
}
