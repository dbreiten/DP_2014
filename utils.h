/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace cel�ch ��sel z pohledu l�m�n� RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    utils.h                                                    **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdint.h>
#include <gmpxx.h>

/******************************************************************************
** MACROS                                                                    **
******************************************************************************/

#define NUM(chr) (chr - 48)

// Error codes
#define OK        0
#define INPUT_IS_TOO_LONG 1
#define INPUT_IS_PRIME 2

// Debug
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

typedef struct {
	uint16_t bitCnt;
	uint8_t  *bits;
} binNum_t;

/******************************************************************************
** FUNCTIONS                                                                 **
******************************************************************************/

progData_t *ParseArgs(int argCnt, char **args);

unsigned char *CheckNumber(char *str);

uint16_t DecNumToBinNum(uint8_t *decNum, uint16_t digitCnt, uint64_t *binNum);

void Pollard_rho_method_GMP(mpz_t res, mpz_t N);

#endif
