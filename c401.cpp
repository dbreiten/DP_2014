
/* c401.c: **********************************************************}
{* Téma: Rekurzivní implementace operací nad BVS
**                                         Vytvoøil: Petr Pøikryl, listopad 1994
**                                         Úpravy: Andrea Nìmcová, prosinec 1995
**                                                      Petr Pøikryl, duben 1996
**                                                   Petr Pøikryl, listopad 1997
**                                  Pøevod do jazyka C: Martin Tuèek, øíjen 2005
**                                         Úpravy: Bohuslav Køena, listopad 2009
**                                         Úpravy: Masarik Karel,  listopad 2010
**
** Implementujte rekurzivním zpùsobem operace nad binárním vyhledávacím
** stromem (BVS; v angliètinì BST - Binary Search Tree).
**
** Klíèem uzlu stromu je jeden znak (obecnì jím mù¾e být cokoliv, podle
** èeho se vyhledává). U¾iteèným (vyhledávaným) obsahem je zde integer.
** Uzly s men¹ím klíèem le¾í vlevo, uzly s vìt¹ím klíèem le¾í ve stromu
** vpravo. Vyu¾ijte dynamického pøidìlování pamìti.
** Rekurzivním zpùsobem implementujte následující funkce:
**
**   BSTInit ...... inicializace vyhledávacího stromu
**   BSTSearch .... vyhledávání hodnoty uzlu zadaného klíèem
**   BSTInsert .... vkládání nové hodnoty
**   BSTDelete .... zru¹ení uzlu se zadaným klíèem
**   BSTDispose ... zru¹ení celého stromu
**
** ADT BVS je reprezentován koøenovým ukazatelem stromu (typ tBSTNodePtr).
** Uzel stromu (struktura typu tBSTNode) obsahuje klíè (typu char), podle
** kterého se ve stromu vyhledává, vlastní obsah uzlu (pro jednoduchost
** typu int) a ukazatel na levý a pravý podstrom (LPtr a RPtr). 
** Pøesnou definici typù naleznete v souboru c401.h.
**
** Pozor! Je tøeba správnì rozli¹ovat, kdy pou¾ít dereferenèní operátor *
** (typicky pøi modifikaci) a kdy budeme pracovat pouze se samotným ukazatelem 
** (napø. pøi vyhledávání). V tomto pøíkladu vám napoví prototypy funkcí.
** Pokud pracujeme s ukazatelem na ukazatel, pou¾ijeme dereferenci.
**/

#include "c401.h"

using namespace std;

int solved;

void BSTInit (tBSTNodePtr *RootPtr) {
/*   -------
** Funkce provede poèáteèní inicializaci stromu pøed jeho prvním pou¾itím.
**
** Ovìøit, zda byl ji¾ strom pøedaný pøes RootPtr inicializován, nelze,
** proto¾e pøed první inicializací má ukazatel nedefinovanou (tedy libovolnou)
** hodnotu. Programátor vyu¾ívající ADT BVS tedy musí zajistit, aby inicializace
** byla volána pouze jednou, a to pøed vlastní prací s BVS. Provedení
** inicializace nad neprázdným stromem by toti¾ mohlo vést ke ztrátì pøístupu
** k dynamicky alokované pamìti (tzv. "memory leak").
**	
** V¹imnìte si, ¾e se v hlavièce objevuje typ ukazatel na ukazatel.	
** Proto je tøeba pøi pøiøazení pøes RootPtr pou¾ít dereferenèní operátor *.
** Ten bude pou¾it i ve funkcích BSTDelete, BSTInsert a BSTDispose.
**/
	
	/* Inicializace stromu */
	*RootPtr = NULL;	
}	

uint64_t *BSTSearch(tBSTNodePtr *RootPtr, double K, uint64_t *Content)	{
/*  ---------
** Funkce vyhledá uzel v BVS s klíèem K.
**
** Pokud je takový nalezen, vrací funkce hodnotu TRUE a v promìnné Content se
** vrací obsah pøíslu¹ného uzlu.´Pokud pøíslu¹ný uzel není nalezen, vrací funkce
** hodnotu FALSE a obsah promìnné Content není definován (nic do ní proto
** nepøiøazujte).
**
** Pøi vyhledávání v binárním stromu bychom typicky pou¾ili cyklus ukonèený
** testem dosa¾ení listu nebo nalezení uzlu s klíèem K. V tomto pøípadì ale
** problém øe¹te rekurzivním volání této funkce, pøièem¾ nedeklarujte ¾ádnou
** pomocnou funkci.
**/

	double rootKeyDiff;
	double lPtrKeyDiff;
	double rPtrKeyDiff;

	rootKeyDiff = fabs(K - (*RootPtr)->Key);

	if ((*RootPtr)->LPtr != NULL)
		lPtrKeyDiff = fabs(K - (*RootPtr)->LPtr->Key);
	else
		lPtrKeyDiff = 1000.0;

	if ((*RootPtr)->RPtr != NULL)
		rPtrKeyDiff = fabs(K - (*RootPtr)->RPtr->Key);
	else
		rPtrKeyDiff = 1000.0;

	/* Nejdrive musime overit, zda neni RootPtr roven NULL, kdyz ano, nic se 
	** nedeje. Jinak porovnavame hledany klic s klicem aktualniho uzlu. Pokud
	** se rovna, nasli jsme hledany klic a do promenne Content ulozime hodnotu
	** uzlu. Pokud je klic mensi, pokracuje hledani do leveho potomka, pokud je
	** vetsi, pak pokracuje hledani v pravem potomkovi.
	*/
	if( RootPtr != NULL )
	{
		if (rootKeyDiff < lPtrKeyDiff && rootKeyDiff < rPtrKeyDiff) {
			//Content = (*RootPtr)->BSTNodeCont;
			return (*RootPtr)->BSTNodeCont;
		}
		else if (lPtrKeyDiff < rPtrKeyDiff)
			return BSTSearch(&((*RootPtr)->LPtr), K, Content);
		else
			return BSTSearch(&((*RootPtr)->RPtr), K, Content);
	}
	else 
		return NULL;
} 


void BSTInsert (tBSTNodePtr* RootPtr, double K, uint64_t *Content)	{	
/*   ---------
** Vlo¾í do stromu RootPtr hodnotu Content s klíèem K.
**
** Pokud ji¾ uzel se zadaným klíèem ve stromu existuje, bude obsah uzlu
** s klíèem K nahrazen novou hodnotou. Pokud bude do stromu vlo¾en nový
** uzel, bude vlo¾en v¾dy jako list stromu.
**
** Funkci implementujte rekurzivnì. Nedeklarujte ¾ádnou pomocnou funkci.
**
** Rekurzivní implementace je ménì efektivní, proto¾e se pøi ka¾dém
** rekurzivním zanoøení ukládá na zásobník obsah uzlu (zde integer).
** Nerekurzivní varianta by v tomto pøípadì byla efektivnìj¹í jak z hlediska
** rychlosti, tak z hlediska pamì»ových nárokù. Zde jde ale o ¹kolní
** pøíklad, na kterém si chceme ukázat eleganci rekurzivního zápisu.
**/
	
	/* Overeni, zda uzel jeste neni alokovany. Provede se prislusna alokace.
	** Jinak se porovnava zadany klic s klicem aktualniho uzlu. Pokud se 
	** rovnaji, provede se pouze jeho aktualizace. Je-li klic mensi,
	** klic aktualniho, pak se pokracuje levym potomkem. Naopak je-li vetsi,
	** pokracuje se pravym potomkem.
	*/
	if( *RootPtr == NULL )
	{
		*RootPtr = (struct tBSTNode*)malloc (sizeof (struct tBSTNode));
		(*RootPtr)->Key = K;
		//(*RootPtr)->BSTNodeCont = Content;
		(*RootPtr)->BSTNodeCont = Content;
		(*RootPtr)->LPtr = NULL;
		(*RootPtr)->RPtr = NULL;
	}
	else
	{
		if( K == (*RootPtr)->Key )
			(*RootPtr)->BSTNodeCont = Content;
		else if( K < (*RootPtr)->Key )
			BSTInsert(&(*RootPtr)->LPtr, K, Content);
		else
			BSTInsert(&(*RootPtr)->RPtr, K, Content);
	}
}

void ReplaceByRightmost (tBSTNodePtr PtrReplaced, tBSTNodePtr *RootPtr) {
/*   ------------------
** Pomocná funkce pro vyhledání, pøesun a uvolnìní nejpravìj¹ího uzlu.
**
** Ukazatel PtrReplaced ukazuje na uzel, do kterého bude pøesunuta hodnota
** nejpravìj¹ího uzlu v podstromu, který je urèen ukazatelem RootPtr.
** Pøedpokládá se, ¾e hodnota ukazatele RootPtr nebude NULL (zajistìte to
** testováním pøed volání této funkce). Tuto funkci implementujte rekurzivnì. 
**
** Tato pomocná funkce bude pou¾ita dále. Ne¾ ji zaènete implementovat,
** pøeètìte si komentáø k funkci BSTDelete(). 
**/
	
	tBSTNodePtr Righmost = NULL;

	if( (*RootPtr) != NULL )
	{
		/* Presun a uvolneni nejpravejsiho uzlu */
		if( (*RootPtr)->RPtr == NULL )
		{
			PtrReplaced->Key = (*RootPtr)->Key;
			PtrReplaced->BSTNodeCont = (*RootPtr)->BSTNodeCont;
			Righmost = (*RootPtr);
			(*RootPtr) = (*RootPtr)->LPtr;
			free(Righmost->BSTNodeCont);
			free(Righmost);
		}
		else
			/* Vyhledani nejpravejsiho uzlu */
			ReplaceByRightmost(PtrReplaced, &(*RootPtr)->RPtr);
	}	
}

void BSTDelete (tBSTNodePtr *RootPtr, double K) 		{
/*   ---------
** Zru¹í uzel stromu, který obsahuje klíè K.
**
** Pokud uzel se zadaným klíèem neexistuje, nedìlá funkce nic. 
** Pokud má ru¹ený uzel jen jeden podstrom, pak jej zdìdí otec ru¹eného uzlu.
** Pokud má ru¹ený uzel oba podstromy, pak je ru¹ený uzel nahrazen nejpravìj¹ím
** uzlem levého podstromu. Pozor! Nejpravìj¹í uzel nemusí být listem.
**
** Tuto funkci implementujte rekurzivnì s vyu¾itím døíve deklarované
** pomocné funkce ReplaceByRightmost.
**/
	
	tBSTNodePtr PomPtr = NULL;

	if( (*RootPtr) != NULL )
	{ 
		if( K == (*RootPtr)->Key )
		{ /* Pokud jsme nalezli hledany uzel, ulozime si ukazatel do pomocne
		  ** promenne.
		  */
			PomPtr = (*RootPtr);
			if(PomPtr->RPtr == NULL)
			{ /* Ruseny uzel ma jeden podstrom a to levy. Dedi jej otec
			  ** ruseneho prvku.
			  */
				(*RootPtr) = PomPtr->LPtr;
				free(PomPtr);
			}
			else if( PomPtr->LPtr == NULL )
			{ /* Ruzeny uzel ma jeden podstrom a to pravy. Dedi jej otec
			  ** ruseneho prvku.
			  */
				(*RootPtr) = PomPtr->RPtr;
				free(PomPtr);
			}
			else
				/* ruseny uzel ma oba podstromy. Ruseny uzel je nahrazen
				** nejpravejsim leveho podstromu.
				*/
				ReplaceByRightmost(*RootPtr, &(*RootPtr)->LPtr);
		}
		else if( K < (*RootPtr)->Key )
			BSTDelete(&(*RootPtr)->LPtr, K);
		else
			BSTDelete(&(*RootPtr)->RPtr, K);
	}	
} 

void BSTDispose (tBSTNodePtr *RootPtr) {	
/*   ----------
** Zru¹í celý binární vyhledávací strom a korektnì uvolní pamì».
**
** Po zru¹ení se bude BVS nacházet ve stejném stavu, jako se nacházel po
** inicializaci. Tuto funkci implementujte rekurzivnì bez deklarování pomocné
** funkce.
**/
	
	if( (*RootPtr) != NULL ) 
	{ /* Pokud koren neni prazdny. */
		BSTDispose(&(*RootPtr)->LPtr);
		BSTDispose(&(*RootPtr)->RPtr);
		free((*RootPtr)->BSTNodeCont);
		free(*RootPtr);
		(*RootPtr) = NULL;
	}	
}

/* konec c401.c */

