/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace cel�ch ��sel z pohledu l�m�n� RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    mpqs.h                                                     **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

#ifndef _MPQS_H_
#define _MPQS_H_

#include <stdint.h>
#include <stdbool.h>

/******************************************************************************
** MACROS                                                                    **
******************************************************************************/

#define MAX_DIGITS 256
#define MAX_BITS   768

/******************************************************************************
** STRUCTS                                                                   **
******************************************************************************/
typedef struct {
	uint16_t digitCnt; // Pocet cislic
	uint16_t bitCnt;   // Pocet bitu
	uint8_t  *N;       // Faktorovane cislo v dekadicke podobe
	uint64_t *binN;    // Faktorovane cislo v ninarni podobe
} num_t;

/******************************************************************************
** FUNCTIONS                                                                 **
******************************************************************************/

bool FermatPrimalityTest(uint64_t n);

uint64_t Euklid_algorithm(uint64_t u, uint64_t w);

//TODO: pro num_t
//bool FermatPrimalityTest(uint64_t n);

//TODO: Trial division
//TODO: Trial divisionLongNum
//TODO: Pollard rho method
//TODO: MPQS {
//          Overit zda cislo neni 0 ci jedna, pak neni co faktorizovat
//          Pouzit zkusme deleni a bud vyhledat vsechny faktory anebo aspon faktorovane cislo snizit
//          Overit, ze cislo neni prvocislo pomoci FermatPrimalityTest()
//          Pokud se nepodarilo jeste plno cislo rozlozit a cislo je dostatecne male, pouzit Pollard rho metodu
//          Jinak nasadit MPQS
//      }
//TODO: Mit vse pripravene tak, at muzu porovnavat QS (seriove) vs MPQS (paralelne)
#endif
