/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace celých èísel z pohledu lámání RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    siqs.cpp                                                   **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "siqs.h"
#include "utils.h"

/** Function performs convertion from string to num_t structure
***
*** @string Input to convert
*** @return Converted string
*/
num_t *ConvertStringToNum(char *string) {
	int32_t      i;
	int32_t strLen = strlen(string);
	int32_t   last = strLen - 1;
	num_t   *numSt;

	// Empty string check
	if (!strLen) {
		fprintf(stderr, "Nothing to factor\n");
		return NULL;
	}

	// Init struct for factored number
	numSt = (num_t *)malloc(sizeof(num_t));
	if (!numSt) {
		fprintf(stderr, "ConvertStringToNum() - Malloc error\n");
		return NULL;
	}
	numSt->digitCnt = 0;

	// Init digit array, assumed that input is correct
	numSt->N = (uint8_t *)malloc(sizeof(uint8_t)* strLen);

	// Convert string to number array
	// First digit is at position 0
	// Last digit is at position last
	for (i = strLen - 1; i >= 0; --i) {
		if (isdigit(string[i])) {
			numSt->N[last - i] = NUM(string[i]);
			i++;
		}
		else {
			fprintf(stderr, "Bad Input\n");
			return NULL;
		}

	}

	return numSt;
}

/** Implementation of Fermat Primality test
***
*** n is input
***
*** Repeat k times:
***	    a = k
***	    If a^(n-1) is not congruent with 1 (mod n), then return composite
***
*** If composite is never returned: return probably prime
**/
bool FermatPrimalityTest(uint64_t input) {
	uint16_t i;
	uint32_t j;
	uint32_t a_orig = 2;
	//uint64_t a_orig = 275;
	uint64_t a = a_orig;
	uint64_t n = input;
	uint64_t exp = n - 1;     //exponent
	uint64_t k = 1024;        //repeat Cnt

	// Input num is smaller than repeat Cnt
	if (n < k)
		k = n - a; // We have to set number of repeat right so it goes [2..n-1]

	for (i = 0; i < k; ++i) {
		for (j = 2; j <= exp; ++j) {
			a = (a * a_orig) % n;
			//fprintf(stderr, "J: %u EXP: %u A: %u\n", j, exp, a);
		}
		if (a != 1)
			return false; // Input is composite number

		fprintf(stdout, "I: %u A: %u\n", i, a_orig);

		a = ++a_orig;

	}
	return true; // Input is probably prime
}

/** Implementation of great common divisor
*/
uint64_t Euklid_algorithm(uint64_t u, uint64_t w)
{
	uint64_t r;

	while (w != 0) {
		r = u % w;
		u = w;
		w = r;
	}

	return u;
}