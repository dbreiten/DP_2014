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
