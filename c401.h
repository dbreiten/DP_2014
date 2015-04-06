
/* Hlavièkový soubor pro c401.c - rekurzívní implementace operaci nad BVS,
obsahuje jednak nutné knihovny a externí promìnné, ale rovnì¾ 
definici datových typù, se kterými se pracuje v jdenotlivých
funkcích. Nemodifikovat! */

/* ********************** SOUBOR S HLAVIÈKOU ********************** */
/* ********************** ------------------ ********************** */

/*  vytvoèil: Martin Tuèek
    pøedmìt: Algoritmy (IAL) - FIT (Fakulta Informacnich Technologií)
    hlavicka pro soubor: c401.c
    datum: záøí 2005
    upravil: Bohuslav Køena, listopad 2009                           
    upravil: Masarik Karel,  listopad 2010                           */
/* ***************************************************************** */

#ifndef C401_H_
#define C401_H_

#include<cstdio>
#include<iostream>
#include<cstdlib>
#include<cstdarg>
#include<cstdint>
#include<cmath>

#define TRUE 1
#define FALSE 0

extern int solved;                        /* indikace, zda byla funkce øe¹ena */

/* uzel stromu */
                                                                                                            
typedef struct tBSTNode {
	double Key;			                                                      /* klíè */
	uint64_t *BSTNodeCont;                                     /* u¾iteèný obsah uzlu */
	struct tBSTNode * LPtr;                                    /* levý podstrom */
	struct tBSTNode * RPtr;                                   /* pravý podstrom */
} *tBSTNodePtr;	

/* prototypy funkcí */

void BSTInit   (tBSTNodePtr *);
uint64_t *BSTSearch(tBSTNodePtr *, double, uint64_t *);
void BSTInsert(tBSTNodePtr *, double, uint64_t *);
void BSTDelete (tBSTNodePtr *, double);
void BSTDispose(tBSTNodePtr *);

/* konec c401.h */
#endif
