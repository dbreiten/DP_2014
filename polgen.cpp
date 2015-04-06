/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace celych cisel z pohledu lamani RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    polgen.cpp                                                 **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

#include "polgen.h"

using namespace std;

void InitPolynom(pol_t &Q_ab, data_t &data)
{
	Q_ab.s = data.divNum;
	//mpz_t nmodp;

	mpz_init(Q_ab.a);
	//mpz_init(nmodp);

	// init space for all B's and pre-computed root offsets
	Q_ab.B = (mpz_t *)malloc(sizeof(mpz_t)* Q_ab.s);
	Q_ab.a_inv_prime = (mpz_t *)malloc(sizeof(mpz_t)*data.primes.size());
	//Q_ab.amodp = (mpz_t *)malloc(sizeof(mpz_t)*data.primes.size());
	//Q_ab.Bainv2 = (mpz_t **)malloc(sizeof(mpz_t *)* Q_ab.s);

	for (uint32_t i = 0; i < Q_ab.s; ++i) {
		mpz_init(Q_ab.B[i]);
		//Q_ab.Bainv2[i] = (mpz_t *)malloc(sizeof(mpz_t)* data.primes.size());
		//for (uint32_t j = 0; j < data.primes.size(); ++j) {
		//	mpz_init(Q_ab.Bainv2[i][j]);
		//}
	}

	for (uint32_t i = 0; i < data.primes.size(); ++i) {
		mpz_init(Q_ab.a_inv_prime[i]);
		//mpz_init(Q_ab.amodp[i]);
	}

	// num of b's
	Q_ab.bsize = (uint32_t)pow(2, Q_ab.s - 1);
	
	mpz_init(Q_ab.b);

	mpz_init(Q_ab.c);

	//mpz_clear(nmodp);

	return;
}

void NextPolynom(pol_t &Q_ab, data_t &data, int32_t &e, uint32_t &v)
{
	mpz_t tmp;
	mpz_t B_v;
	mpz_t b;

	uint32_t exp = 0;
	uint32_t dbl_i = 0;

	int32_t x = 0;

	mpz_init(tmp);
	mpz_init(B_v);
	mpz_init(b);

	Q_ab.acc_pol++;

#ifdef DEBUG_PRINT
	cout << "NextPolynom()" << endl;
#endif

	// We just use next b
	if (Q_ab.acc_pol < Q_ab.bsize) {
		//old b
		mpz_set(b, Q_ab.b);

		// exp
		dbl_i = 2 * Q_ab.acc_pol;
		//cout << "dbl_i: " << dbl_i << endl;

		exp = 1;
		dbl_i >>= 1;

		while ((dbl_i & 1) != 1) {
			exp += 1;
			dbl_i >>= 1;
		}
		//cout << "v: " << exp << endl;
		v = exp - 1;
		//cout << "Real v: " << v << endl;

		// B_v
		// exp - 1 because in continis paper it's 1..s
		// there it has to be 0..s-1
		mpz_set(B_v, Q_ab.B[exp-1]);
		//cout << "B_v: ";
		//mpz_out_str(stdout, 10, B_v);
		//cout << endl;

		// 2^exp
		exp = (uint32_t)pow(2, exp);
		//cout << "2^v: " << exp << endl;

		// ceil(i/exp)
		x = (int32_t)ceil(Q_ab.acc_pol / (double)exp);
		//cout << "x = ceil(i/exp): " << x << endl;

		// (-1)^x
		x = (int32_t)pow(-1, x);
		//cout << "(-1)^x: " << x << endl;
		e = x;
		//cout << "e: " << e << endl;

		// 2(-1)^x
		x *= 2;
		//cout << "2(-1)^x: " << x << endl;

		// 2(-1)^xB
		mpz_mul_si(tmp, B_v, x);
		//cout << "2(-1)^xB_v: ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		// b_new = b + 2(-1)^xB
		mpz_add(b, b, tmp);
		//cout << "b_new = b + 2(-1)^xB: ";
		//mpz_out_str(stdout, 10, b);
		//cout << endl;

		mpz_set(Q_ab.b, b);
		
		///////////////////////
		// c
		mpz_mul(tmp, Q_ab.b, Q_ab.b);
		mpz_sub(tmp, tmp, data.N);
		mpz_tdiv_q(Q_ab.c, tmp, Q_ab.a);

		//cout << "a = ";
		//mpz_out_str(stdout, 10, Q_ab.a);
		//cout << endl;

		//cout << "b = ";
		//mpz_out_str(stdout, 10, Q_ab.b[Q_ab.acc_pol]);
		//cout << endl;

		//cout << "c = ";
		//mpz_out_str(stdout, 10, Q_ab.c);
		//cout << endl;

#ifdef LOG_ENABLE
		ofstream polynomFile;

		string fileName = "pol_" + to_string(omp_get_thread_num()) + ".txt";

		polynomFile.open(fileName, ios::out | ios::app);

		polynomFile << "a: " << mpz_get_str(NULL, 10, Q_ab.a);
		polynomFile << "\tdivisors:\t";
		for (map<uint64_t, uint64_t>::iterator it = Q_ab.divisors.begin(); it != Q_ab.divisors.end(); ++it)
			polynomFile << it->first << "\t";
		polynomFile << "\tb: " << mpz_get_str(NULL, 10, Q_ab.b);
		polynomFile << "\tc: " << mpz_get_str(NULL, 10, Q_ab.c);
		polynomFile << endl;

		polynomFile.close();
#endif
	}
	// All b's has been used => we need to generate new A and b's for it
	else {
		GeneratePolynom(data, Q_ab);
	}

	// Cleaning
	mpz_clear(tmp);
	mpz_clear(B_v);
	mpz_clear(b);

	return;
}

void GenerateA(pol_t &Q_ab, data_t &data)
{
	double logNextksb = 0.0;
	double logDiff;
	uint32_t logNum;

	uint64_t *content = NULL;

	mpz_t prime;

	mpz_init(prime);

#ifdef DEBUG_PRINT
	cout << "GenerateA()" << endl;
#endif

	// if acc_pol != 0 => clear divisors for next usage
	Q_ab.divisors.clear();

	mpz_set_ui(Q_ab.a, 1);

	if (Q_ab.s > 3) {
		// Zde by se musel nejak zmenit interval, aby outValues melo smysl, pokud tohle bude fungovat, tak outValues ani nejsou treba
#ifdef NO_PARALLEL
		NEXTKSB(data.nexksbPrimes, data.outIndexes, data.outValues, data.nexksbPrimes.size(), data.k, data.h, data.m);

		for (uint32_t i = 0; i < data.outIndexes.size(); ++i) {
			Q_ab.divisors.insert(pair<uint64_t, uint64_t>(data.nexksbPrimes.at(data.outIndexes.at(i) - 1), data.nexksbPrimes.at(data.outIndexes.at(i) - 1)));
			mpz_mul_ui(Q_ab.a, Q_ab.a, data.nexksbPrimes.at(data.outIndexes.at(i) - 1));
			logNextksb += log(data.nexksbPrimes.at(data.outIndexes.at(i) - 1));
		}
#else
#pragma omp critical
		{
			NEXTKSB(data.nexksbPrimes, data.outIndexes, data.outValues, data.nexksbPrimes.size(), data.k, data.h, data.m);

			for (uint32_t i = 0; i < data.outIndexes.size(); ++i) {
				Q_ab.divisors.insert(pair<uint64_t, uint64_t>(data.nexksbPrimes.at(data.outIndexes.at(i) - 1), data.nexksbPrimes.at(data.outIndexes.at(i) - 1)));
				mpz_mul_ui(Q_ab.a, Q_ab.a, data.nexksbPrimes.at(data.outIndexes.at(i) - 1));
				logNextksb += log(data.nexksbPrimes.at(data.outIndexes.at(i) - 1));
			}
		}
#endif
	}

	logDiff = data.logIdealA - logNextksb;

	if (Q_ab.s == 2) {
		content = BSTSearch(&data.log2primes, logDiff, NULL);
		logNum = 2;
	}
	else {
		content = BSTSearch(&data.log3primes, logDiff, NULL);
		logNum = 3;
	}

	Q_ab.logDivisors = 0;
	for (uint32_t i = 0; i < logNum; ++i) {
		Q_ab.divisors.insert(pair<uint64_t, uint64_t>(content[i], content[i]));
		Q_ab.logDivisors += data.logs[content[i]];
		mpz_mul_ui(Q_ab.a, Q_ab.a, content[i]);
	}

	//cout << "a = ";
	//mpz_out_str(stdout, 10, Q_ab.a);
	//cout << endl;

#ifdef LOG_ENABLE
	ofstream logFile;

	string fileName = "log_" + to_string(omp_get_thread_num()) + ".txt";

	logFile.open(fileName, ios::out | ios::app);

	logFile << data.logIdealA - log(mpz_get_d(Q_ab.a)) << endl;

	logFile.close();
#endif

	for (uint32_t i = 0; i < data.primes.size(); ++i) {
		//mpz_init();
		mpz_set_ui(prime, data.primes.at(i));
		mpz_invert(Q_ab.a_inv_prime[i], Q_ab.a, prime);
		//mpz_mod(Q_ab.amodp[i], Q_ab.a, prime);
	}

	return;
}

void GenerateFakeA(data_t &data, pol_t &Q_ab) {

	string inputStr;
	uint64_t primeInt = 0;

	Q_ab.divisors.clear();

	mpz_set_ui(Q_ab.a, 1);

	for (uint32_t i = 0; i < data.divNum; ++i) {
		cout << "Zadej " << i + 1 << ". prvocislo z " << data.divNum << endl;
		cin >> inputStr;
		cout << "Zadano: " << atoi(inputStr.c_str()) << endl;
		Q_ab.divisors.insert(pair<uint64_t, uint64_t>(atoi(inputStr.c_str()), atoi(inputStr.c_str())));
		mpz_mul_ui(Q_ab.a, Q_ab.a, atoi(inputStr.c_str()));
	}

	cout << "a = ";
	mpz_out_str(stdout, 10, Q_ab.a);
	cout << endl;

	return;
}

void GeneratePolynom(data_t &data, pol_t &Q_ab)
{	
	uint64_t nextUsablePrime;

	mpz_t prime;
	mpz_t res;
	mpz_t a_inv;
	mpz_t a_k;
	mpz_t a_k_inv;
	mpz_t b;
	mpz_t B;
	mpz_t Bainv2;
	mpz_t tmp;
	mpz_t t1;
	mpz_t t2;
	mpz_t gamma1;
	mpz_t gamma2;
	mpz_t nmodp;

	uint32_t exp = 0;
	uint32_t dbl_i = 0;

	int32_t x = 0;

#ifdef DEBUG_PRINT
	cout << "Generate polynom()" << endl;
#endif

	mpz_init(prime);
	mpz_init(res);
	mpz_init(a_inv);
	mpz_init(a_k);
	mpz_init(a_k_inv);
	mpz_init(b);
	mpz_init(B);
	mpz_init(Bainv2);
	mpz_init(tmp);
	mpz_init(t1);
	mpz_init(t2);
	mpz_init(gamma1);
	mpz_init(gamma2);
	mpz_init(nmodp);

	// Generate A
	GenerateA(Q_ab, data);
	//GenerateFakeA(data, Q_ab);

	// Get B_l
	uint32_t i = 0;
	for (map<uint64_t, uint64_t>::iterator it = Q_ab.divisors.begin(); it != Q_ab.divisors.end(); ++it, ++i) {
		mpz_set_ui(prime, it->first);
		//cout << "prime = ";
		//mpz_out_str(stdout, 10, prime);
		//cout << endl;

		quadratic_residue(t1, data.N, prime);
		//cout << "t = ";
		//mpz_out_str(stdout, 10, t1);
		//cout << endl;

		mpz_sub(t2, prime, t1);
		//cout << "t2 = ";
		//mpz_out_str(stdout, 10, t2);
		//cout << endl;

		mpz_tdiv_q(a_k, Q_ab.a, prime);
		//cout << "a_k = ";
		//mpz_out_str(stdout, 10, a_k);
		//cout << endl;

		mpz_invert(a_k_inv, a_k, prime);
		//cout << "a_k_inv = ";
		//mpz_out_str(stdout, 10, a_k_inv);
		//cout << endl;

		mpz_mul(gamma1, t1, a_k_inv);
		mpz_mod(gamma1, gamma1, prime);
		//cout << "gamma1 = ";
		//mpz_out_str(stdout, 10, gamma1);
		//cout << endl;

		mpz_mul(gamma2, t2, a_k_inv);
		mpz_mod(gamma2, gamma2, prime);
		//cout << "gamma2 = ";
		//mpz_out_str(stdout, 10, gamma2);
		//cout << endl;

		if (mpz_cmp(gamma1, gamma2) < 0)
			mpz_mul(res, a_k, gamma1);
		else
			mpz_mul(res, a_k, gamma2);

		//cout << "res = ";
		//mpz_out_str(stdout, 10, res);
		//cout << endl;

		mpz_set(Q_ab.B[i], res);
		//cout << "==========" << endl;

		// Binv2
		/*
		for (uint32_t j = 0; j < data.primes.size(); ++j) {
			if (Q_ab.divisors.find(data.primes.at(j)) != Q_ab.divisors.end()) {
				//cout << "Divisor found: " << data.primes.at(j);
				mpz_init_set_ui(Q_ab.Bainv2[i][j], data.M);
				continue;
			}

			mpz_set_ui(prime, data.primes.at(j));
			//cout << "prime = ";
			//mpz_out_str(stdout, 10, prime);
			//cout << endl;

			mpz_set(Bainv2, Q_ab.B[i]);
			//cout << "B = ";
			//mpz_out_str(stdout, 10, Bainv2);
			//cout << endl;

			mpz_mul_ui(Bainv2, Bainv2, 2);
			//cout << "2B = ";
			//mpz_out_str(stdout, 10, Bainv2);
			//cout << endl;

			mpz_set(tmp, Q_ab.a);
			//cout << "a = ";
			//mpz_out_str(stdout, 10, tmp);
			//cout << endl;

			mpz_invert(a_inv, tmp, prime);
			//cout << "a_inv mod p = ";
			//mpz_out_str(stdout, 10, a_inv);
			//cout << endl;

			mpz_mul(Bainv2, Bainv2, a_inv);
			//cout << "Bainv2 = ";
			//mpz_out_str(stdout, 10, Bainv2);
			//cout << endl;

			mpz_mod(Bainv2, Bainv2, prime);
			//cout << "Bainv2 mod prime = ";
			//mpz_out_str(stdout, 10, Bainv2);
			//cout << endl;

			//mpz_init_set(Q_ab.Bainv2[i][j], Bainv2);
			mpz_set(Q_ab.Bainv2[i][j], Bainv2);
			//cout << "==========" << endl;
		}
		*/
	}

	Q_ab.Bsize = Q_ab.s;
	//Q_ab.Bainv2Size = Q_ab.s;

	// Get all b's
	// b_1 = B_1 ... B_s
	for (uint32_t i = 0; i < Q_ab.s; ++i)
		mpz_add(b, b, Q_ab.B[i]);

	//cout << "b = ";
	//mpz_out_str(stdout, 10, b);
	//cout << endl;
	
	// Save b_1
	mpz_set(Q_ab.b, b);

	Q_ab.acc_pol = 0;

	// compute c
	mpz_mul(tmp, Q_ab.b, Q_ab.b);
	mpz_sub(tmp, tmp, data.N);
	mpz_tdiv_q(Q_ab.c, tmp, Q_ab.a);
	//cout << "c: ";
	//mpz_out_str(stdout, 10, Q_ab.c);
	//cout << endl;

#ifdef LOG_ENABLE
	ofstream polynomFile;

	string fileName = "pol_" + to_string(omp_get_thread_num()) + ".txt";

	polynomFile.open(fileName, ios::out | ios::app);

	polynomFile << "a: " << mpz_get_str(NULL, 10, Q_ab.a);
	polynomFile << "\tdivisors:\t";
	for (map<uint64_t, uint64_t>::iterator it = Q_ab.divisors.begin(); it != Q_ab.divisors.end(); ++it)
		polynomFile << it->first << "\t";
	polynomFile << "\tb: " << mpz_get_str(NULL, 10, Q_ab.b);
	polynomFile << "\tc: " << mpz_get_str(NULL, 10, Q_ab.c);
	polynomFile << "\tk: " << data.k;
	polynomFile << "\th: " << data.h;
	polynomFile << "\tm: " << data.m;
	polynomFile << "\tn: " << data.nexksbPrimes.size();
	polynomFile << endl;

	polynomFile.close();
#endif

	// cleaning
	mpz_clear(prime);
	mpz_clear(res);
	mpz_clear(a_inv);
	mpz_clear(a_k);
	mpz_clear(a_k_inv);
	mpz_clear(b);
	mpz_clear(B);
	mpz_clear(Bainv2);
	mpz_clear(tmp);
	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(gamma1);
	mpz_clear(gamma2);
	mpz_clear(nmodp);

	return;
}
