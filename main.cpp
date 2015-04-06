/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace celych cisel z pohledu lamani RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    main.cpp                                                   **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

// Standard libraries
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <ctime>

#include <string>

// Big Numbers library
#if defined (_WIN32) || defined (_WIN64)
#include <mpir.h>
#else
#include <gmpxx.h>
#endif

// Implemented libraries
#include "siqs.h"
#include "utils.h"
#include "todo.h"

using namespace std;

int main(int argc, char **argv) {
	bool primalityCheck = true; // TODO: udelat jako argument
	bool verbose        = true; // TODO: udelat jako argument
	bool resErr         = false;

	uint8_t primalityResult;

	double elapsed_secs = 0.0;

#ifdef TIME_MEASUREMENT
	clock_t begin = clock();
#endif

	// 5915587277 * 3267000013 - FACTORIZED
	//string inputStr = "19326223710861634601";

	// 5754853343 * 4093082899 - FACTORIZED
	//string inputStr = "23555091804486281357";

	// 5754853343 * 3628273133 - FACTORIZED
	//string inputStr = "20880179768762133619";

	// 3628273133 * 48112959837082048697 - UNFACTORIZABLE BECAUSE OF SETTING
	//string inputStr = "174566959525992854403922757701";

	// 3628273133 * 29497513910652490397 - FACTORIZED
	//string inputStr = "107025037212314193406975603801";

	// 40 dec
	// 40206835204840513073  * 29497513910652490397 - FACTORIZED
	//string inputStr = "1186001680758095307567274308126685459981";

	// 48112959837082048697 * 54673257461630679457 - FACTORIZED
	//string inputStr = "2630492240413883318777134293253671517529";

	// 50 dec
	// 48112959837082048697 * 521419622856657689423872613771 - FACTORIZED
	//string inputStr = "25087041372768840419953622916560595543980894806387";

	// 48112959837082048697 * 521419622856657689423872613771 - FACTORIZED
	//string inputStr = "13511532193249947664972616118268332877291635311699";

	// 60 dec
	// 590872612825179551336102196593 * 416064700201658306196320137931  - FACTORIZED
	//string inputStr = "245841236512478852752909734912575581815967630033049838269083";

	// 70 dec
	// 521419622856657689423872613771 * 4384165182867240584805930970951575013697
	//string inputStr = "2285989756191926277520207738684272017787417494118059342357922515821387";

	// 80 dec
	// 4384165182867240584805930970951575013697 * 5991810554633396517767024967580894321153
	//string inputStr = "26269087215960187077040851650618265220601648677509488874447059303016675491832641";

	// 90 dec
	// 58645563317564309847334478714939069495243200674793 * 4384165182867240584805930970951575013697
	//string inputStr = "257111836826501668906341282406241950886425250422968530563869938518309060597672268617639721";

	// 100 dec
	// 58645563317564309847334478714939069495243200674793 * 48705091355238882778842909230056712140813460157899 
	//string inputStr = "2856337518961416002334059434204644692528968865631657222194681271598884399851871623337343257129139907";
	string inputStr;

	if (argc > 1)
		inputStr = argv[1];
	else {
		cerr << "No number to factor" << endl;
		return ERR;
	}

	string primesFileName = "primesDelta.txt";


	// mpz_t is the type defined for GMP integers.
	// It is a pointer to the internals of the GMP integer data structure
	mpz_t n;
	mpz_t res;
	int flag;

	// Check how long the inputed number is
	if(inputStr.length() > MAX_DIGITS) {
		cerr << "Input is too long! (" << inputStr.length() << ") digits" << endl;
		return INPUT_IS_TOO_LONG;
	}

	// Initialize the number
	mpz_init(n);
	mpz_init(res);

	// Parse the input string as a base 10 number
	flag = mpz_set_str(n, inputStr.c_str(), 10);
	assert(flag == 0); /* If flag is not 0 then the operation failed */

	// Print n
	cout << "N = ";
	mpz_out_str(stdout, 10, n);
	cout << endl;

	// Primality check
	if(primalityCheck) {
		primalityResult = mpz_probab_prime_p(n,25);
		if(primalityResult > 0) {
			cerr << "N is prime => Nothing to factor" << endl;
			return INPUT_IS_PRIME;
		}
	}

	// If the input is too small, factor with Pollard rho method else SIQS
	data_t data;
	if(inputStr.length() <= 30) {
		// TODO: Can't factor powers of primes!
		Pollard_rho_method_GMP(res, n);
	} else {
		mpz_init_set(data.N, n);
		
		resErr = SIQS(primesFileName, data, inputStr.length(), res);
	}

	// Print result
	if(!resErr) {
		cout << "Result = ";
		mpz_out_str(stdout, 10, res);
		cout << endl;
	}

	// Save the result to the file
	ofstream resFile;
	string fileName = "result_" + inputStr + ".txt";
	resFile.open(fileName);
	resFile << mpz_get_str(NULL, 10, res);
	resFile << endl;
	resFile << "Multiplier: " << data.multiplier << endl;

#ifdef TIME_MEASUREMENT
	clock_t end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "Elapsed time: " << elapsed_secs << "s" << endl;
	resFile << "Elapsed time: " << elapsed_secs << "s" << endl;
#endif

	resFile.close();

	/* Clean up the mpz_t handles or else we will leak memory */
	mpz_clear(n);
	mpz_clear(res);

	if(resErr)
		return ERR;

	return OK;
}
