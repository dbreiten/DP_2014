/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace celych cisel z pohledu lamani RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    mpqs.h                                                     **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

#ifndef _SIQS_H_
#define _SIQS_H_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <cstdbool>

#include <string>
#include <vector>
#include <map>

// Big Numbers library
#if defined (_WIN32) || defined (_WIN64)
#include <mpir.h>
#else
#include <gmpxx.h>
#endif

// OpenMP
#include <omp.h>

#include "c401.h"

using namespace std;

/******************************************************************************
** MACROS                                                                    **
******************************************************************************/

#define MAX_DIGITS 256
//#define NO_PARALLEL
//#define LOG_ENABLE
//#define LOG_VALUES_ENABLE
//#define LOG_VECTOR_ENABLE
//#define LOG_GAUSS_ENABLE

/******************************************************************************
** STRUCTS                                                                   **
******************************************************************************/

typedef struct {
	mpz_t a;
	mpz_t b;
	mpz_t c;
	int64_t x;
	mpz_t a2;    // If relation is made from partials
	mpz_t b2;    // If relation is made from partials
	mpz_t c2;    // If relation is made from partials
	int64_t x2; // If relation is made from partials
	uint32_t *relation;
	bool composed;
} relInfoUInt_t;

typedef struct {
	mpz_t N;                         // Zadane cislo
	uint32_t digits;
	mpz_t origN;
	uint32_t M;                      // Pocet iteraci provedenych nad jednim polynomem
	uint32_t F;                      // Pocet prvocisel, ktere se budou zkoumat a pripadne i ukladat do FB
	uint32_t k;                      // Pocet delitelu, ktere je treba najit pro a
	uint32_t h;                      // Promenna pro NEXTKSB
	uint32_t m;                      // Promenna pro NEXTKSB
	uint32_t multiplier;             // kN mod 8 = 1
	vector<uint64_t> primes;         // Vector prvocisel, pro ktere plati, ze Legendreuv symbol = 1 pro zadane N
	mpz_t *primesMPZ;
	mpz_t *quadraticResidue;
	vector<uint64_t> nexksbPrimes;  // Primes with odd indicies
	vector<uint64_t> treePrimes;     // Primes with even indicies
	map<uint64_t, bool> largePrimes; // Mapa obsahujici vsechna velka prvocisla, ktera mohou byt pouzita pro nalezeni castecne relace
	mpz_t idealA;                    // Idealni hodnota a
	double logIdealA;                // Logaritmus idealni hodnoty a
	tBSTNodePtr log2primes;          // Binarni strom uchovavajici logaritmy vsech dvojic z FB
	tBSTNodePtr log3primes;          // Binarni strom uchovavajici logaritmy vsech trojic z FB
	uint32_t leftBound;              // Spodni hranice prvocisel pro hledani delitele a
	uint32_t rightBound;             // Horni hranice prvocisel pro hledani delitele a
	uint32_t divNum;                 // Optimalni pocet delitelu
	vector<uint32_t> divIndexes;     // Uchovava indexy poslednich delitelu a pro pripadne generovani noveho a
	double *logs;
	double *threshold;
	double logLargestPrime;          // Logaritmus nejvetsiho prvocisla, se kterym se bude pocitat pri hledani castecnych relaci
	uint32_t largePrimeMultiplier;   // largePrimeMultiplier * F => Large FB
	uint32_t relationsNeeded;        // Pocet relaci, kterych je treba nasbirat
	uint32_t relationsNeededPart;
	uint32_t expVectorLength;        // Delka vektoru exponentu
	uint32_t expVectorLengthAlloc;   // Alokovana delka vektoru exponeny
	vector<relInfoUInt_t> relations; // Vektor relaci
	multimap<uint64_t, relInfoUInt_t> partialRelationsUInt; // Mapa castecnych relaci

	vector<uint64_t> outValues;        // Delitele ziskani z NEXTKSB
	vector<uint64_t> outIndexes;       // Indexy pro ziskani dalsich delitelu z NEXTKSB
} data_t;

typedef struct {
	uint32_t bits;       /* size of integer this config info applies to */
	uint32_t fb_size;    /* number of factor base primes */
	uint32_t large_mult; /* the large prime multiplier */
	uint32_t sieve_size; /* the size of the sieve (actual sieve is 2x this) */
	uint32_t div_num;    /* number of divisors */
} sieve_param_t;

/******************************************************************************
** FUNCTIONS                                                                 **
******************************************************************************/

bool SIQS(string primesFileName, data_t &data, uint16_t digits, mpz_t res);

#endif
