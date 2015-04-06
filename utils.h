/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace celych cisel z pohledu lamani RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    utils.h                                                    **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

#ifndef _UTILS_H_
#define _UTILS_H_

#include <cstdio>
#include <iostream>
#include <cstdint>
#include <vector>

// Big Numbers library
#if defined (_WIN32) || defined (_WIN64)
#include <mpir.h>
#else
#include <gmpxx.h>
#endif

#include "siqs.h"
#include "c401.h"

/******************************************************************************
** MACROS                                                                    **
******************************************************************************/

#define NUM(chr) (chr - 48)

// Error codes
#define OK        0
#define ERR       1
#define INPUT_IS_TOO_LONG 1
#define INPUT_IS_PRIME 2

// Debug
#define TIME_MEASUREMENT
//#define DEBUG_PRINT
#define UTILS_DEBUG_PRINT
//#define UTILS_DEBUG_PRINT_POLLARD_RHO

/******************************************************************************
** ENUMS                                                                     **
******************************************************************************/


/******************************************************************************
** STRUCTS                                                                   **
******************************************************************************/
typedef struct {
	uint32_t parsedArgs;
	char     *inputNumStr;
} progData_t;

/******************************************************************************
** FUNCTIONS                                                                 **
******************************************************************************/

using namespace std;

progData_t *ParseArgs(int argCnt, char **args);

unsigned char *CheckNumber(char *str);

uint16_t DecNumToBinNum(uint8_t *decNum, uint16_t digitCnt, uint64_t *binNum);

uint64_t Euklid_algorithm(uint64_t u, uint64_t w);

void Euklid_algorithm_BIG(mpz_t u, mpz_t w, mpz_t res);

void Pollard_rho_method_GMP(mpz_t res, mpz_t N);

uint32_t NEXTKSB(vector<uint64_t> &values, vector<uint64_t> &outIndexes, vector<uint64_t> &outValues, uint32_t n, uint32_t k, uint32_t &h, uint32_t &m);

void Build2PrimesLogTree(data_t &data);

void Build3PrimesLogTree(data_t &data);

uint32_t mp_modsqrt_1(uint32_t a, uint32_t p);

int quadratic_residue(mpz_t x, mpz_t q, mpz_t n);

#endif
