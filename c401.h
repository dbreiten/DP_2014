
/* Hlavi�kov� soubor pro c401.c - rekurz�vn� implementace operaci nad BVS,
obsahuje jednak nutn� knihovny a extern� prom�nn�, ale rovn� 
definici datov�ch typ�, se kter�mi se pracuje v jdenotliv�ch
funkc�ch. Nemodifikovat! */

/* ********************** SOUBOR S HLAVI�KOU ********************** */
/* ********************** ------------------ ********************** */

/*  vytvo�il: Martin Tu�ek
    p�edm�t: Algoritmy (IAL) - FIT (Fakulta Informacnich Technologi�)
    hlavicka pro soubor: c401.c
    datum: z��� 2005
    upravil: Bohuslav K�ena, listopad 2009                           
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

extern int solved;                        /* indikace, zda byla funkce �e�ena */

/* uzel stromu */
                                                                                                            
typedef struct tBSTNode {
	double Key;			                                                      /* kl�� */
	uint64_t *BSTNodeCont;                                     /* u�ite�n� obsah uzlu */
	struct tBSTNode * LPtr;                                    /* lev� podstrom */
	struct tBSTNode * RPtr;                                   /* prav� podstrom */
} *tBSTNodePtr;	

/* prototypy funkc� */

void BSTInit   (tBSTNodePtr *);
uint64_t *BSTSearch(tBSTNodePtr *, double, uint64_t *);
void BSTInsert(tBSTNodePtr *, double, uint64_t *);
void BSTDelete (tBSTNodePtr *, double);
void BSTDispose(tBSTNodePtr *);

/* konec c401.h */
#endif
