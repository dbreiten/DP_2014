/*
/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace celych cisel z pohledu lamani RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    sieve.h                                                    **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

#ifndef SIEVE_H_
#define SIEVE_H_

#include <cstdint>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <cmath>

#include "polgen.h"
#include "utils.h"

/******************************************************************************
** STRUCTS                                                                   **
******************************************************************************/
typedef struct {
	mpz_t root1;
	mpz_t root2;
	int32_t root1_int;
	int32_t root2_int;
	uint64_t prime;
} root_t;

typedef struct {
	uint64_t *divisors;
	uint8_t max_divs;
	uint8_t divs_cnt;
} divData_t;

/******************************************************************************
** FUNCTIONS                                                                 **
******************************************************************************/

using namespace std;

bool Sieve(data_t &data, pol_t &Q_ab);


#endif /* SIEVE_H_ */
