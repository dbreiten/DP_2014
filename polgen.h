/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace celych cisel z pohledu lamani RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    polgen.h                                                   **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

#ifndef POLGEN_H_
#define POLGEN_H_

#include <cstdint>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <malloc.h>

// Big Numbers library
#if defined (_WIN32) || defined (_WIN64)
#include <mpir.h>
#else
#include <gmpxx.h>
#endif

// OpenMP
#include <omp.h>

#include "siqs.h"
#include "utils.h"

using namespace std;

/******************************************************************************
** STRUCTS                                                                   **
******************************************************************************/
typedef struct {
	mpz_t a;                           // a koeficient
	mpz_t b;
	uint32_t bsize;                    // Pocet b koeficientu
	mpz_t *B;                          // B hodnoty
	uint32_t Bsize;                    // Pocet B hodnot
	//mpz_t **Bainv2;                    // Predvypoctene hodnoty pro rychle ziskani korenu pro dalsi polynom
	//uint32_t Bainv2Size;               // Velikost pole uchovavajici predvypoctene hodnoty pro nasledujici koreny
	mpz_t c;                           // c koeficient polynomu
	uint32_t acc_pol;                  // Index aktualniho polynomu
	uint32_t s;                        // Pocet delitelu a
	map<uint64_t, uint64_t> divisors;  // Delitele a
	mpz_t *a_inv_prime;                
	//mpz_t *amodp;
	double logDivisors;
} pol_t;

/******************************************************************************
** FUNCTIONS                                                                 **
******************************************************************************/

void InitPolynom(pol_t &Q_ab, data_t &data);

void NextPolynom(pol_t &Q_ab, data_t &data, int32_t &e, uint32_t &v);

void GeneratePolynom(data_t &data, pol_t &Q_a);


#endif /* POLGEN_H_ */
