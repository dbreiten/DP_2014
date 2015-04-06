/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace celych cisel z pohledu lamani RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    sieve.cpp                                                  **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

#include "sieve.h"

using namespace std;

/*
bool ComputeNextRoots(data_t &data, root_t *roots, pol_t &Q_ab, int32_t e, uint32_t v) {
	mpz_t tmp;

	mpz_init(tmp);

	//cout << "a = ";
	//mpz_out_str(stdout, 10, Q_ab.a);
	//cout << endl;

	//cout << "b = ";
	//mpz_out_str(stdout, 10, Q_ab.b);
	//cout << endl;

	//cout << "c = ";
	//mpz_out_str(stdout, 10, Q_ab.c);
	//cout << endl;

	for (uint32_t i = 0; i < data.primes.size(); ++i) {
		//cout << "Prime: " << data.primes.at(i) << endl;
		//cout << "Prime root: " << roots[i].prime << endl;

		if (Q_ab.divisors.find(data.primes.at(i)) != Q_ab.divisors.end()) {
			roots[i].prime = data.primes.at(i);
			roots[i].root1_int = data.M;
			roots[i].root2_int = data.M;
			continue;
		}

		//cout << "v: " << v << endl;
		//cout << "e: " << e << endl;
		//e = (-1)*e;
		//cout << "e: " << e << endl;

		//cout << "B_v = ";
		//mpz_out_str(stdout, 10, Q_ab.B[v]);
		//cout << endl;

		//cout << "Bainv2 = ";
		//mpz_out_str(stdout, 10, Q_ab.Bainv2[v][i]);
		//cout << endl;

		mpz_mul_si(tmp, Q_ab.Bainv2[v][i], e);
		//cout << "e*Bainv2 = ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		mpz_mul_si(tmp, tmp, e);
		//cout << "(-1)*e*Bainv2 = ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		//cout << "root1 = ";
		//mpz_out_str(stdout, 10, roots[i].root1);
		//cout << endl;

		mpz_add(roots[i].root1, roots[i].root1, tmp);
		//cout << "root1 + e*Bainv2 = ";
		//mpz_out_str(stdout, 10, roots[i].root1);
		//cout << endl;

		mpz_mod_ui(roots[i].root1, roots[i].root1, roots[i].prime);
		//cout << "root1 mod prime = ";
		//mpz_out_str(stdout, 10, roots[i].root1);
		//cout << endl;

		roots[i].root1_int = mpz_get_ui(roots[i].root1);

		//cout << "root2 = ";
		//mpz_out_str(stdout, 10, roots[i].root2);
		//cout << endl;

		mpz_add(roots[i].root2, roots[i].root2, tmp);
		//cout << "root2 + e*Bainv2 = ";
		//mpz_out_str(stdout, 10, roots[i].root2);
		//cout << endl;

		mpz_mod_ui(roots[i].root2, roots[i].root2, roots[i].prime);
		//cout << "root2 mod prime = ";
		//mpz_out_str(stdout, 10, roots[i].root2);
		//cout << endl;

		roots[i].root2_int = mpz_get_ui(roots[i].root2);
	}

	//cleaning
	mpz_clear(tmp);

	return true;
}
*/

bool ComputeRoots(data_t &data, root_t *roots, pol_t &Q_ab)
{
	//mpz_t a_inv;
	mpz_t prime;

	//mpz_t a;
	mpz_t b;
	mpz_t tmp;

	bool skip = false;

	uint32_t i = 0;

#ifdef DEBUG_PRINT
	cout << "ComputeRoots()" << endl;
#endif

	//mpz_init(a_inv);
	mpz_init(prime);

	//mpz_init(a);
	mpz_init(b);
	mpz_init(tmp);

	//cout << "a = ";
	//mpz_out_str(stdout, 10, Q_ab.a);
	//cout << endl;

	//cout << "b = ";
	//mpz_out_str(stdout, 10, Q_ab.b[Q_ab.acc_pol]);
	//cout << endl;

	if (data.primes.at(0) == 2) {
		//cout << "Prime: " << 2 << endl;
		roots[i].prime = data.primes.at(0);
		mpz_mul(tmp, Q_ab.b, Q_ab.b);
		mpz_sub(tmp, tmp, data.N);
		mpz_mod_ui(tmp, tmp, 2);

		if (mpz_cmp_ui(tmp, 0) == 0) {
			mpz_set_si(roots[i].root1, 0);
			mpz_set_si(roots[i].root2, -2);
		}
		else {
			mpz_set_si(roots[i].root1, 1);
			mpz_set_si(roots[i].root2, -1);
		}
		roots[i].root1_int = mpz_get_si(roots[i].root1);
		roots[i].root2_int = mpz_get_si(roots[i].root2);

		//cout << "root1 mod prime = ";
		//mpz_out_str(stdout, 10, roots[i].root1);
		//cout << endl;

		//cout << "root2 mod prime = ";
		//mpz_out_str(stdout, 10, roots[i].root2);
		//cout << endl;

		++i;
	}

	for (; i < data.primes.size(); ++i) {
		mpz_set_ui(prime, data.primes.at(i));
		//cout << "Prime = ";
		//mpz_out_str(stdout, 10, prime);
		//cout << endl;

		// It's nesesary to test if prime is divisor of a.
		// If yes there would be divison by zero and thats not allowed
		// so we have to skip this prime.
		// It must be noted that skipping this prime means no smooth
		// numbers with this prime for this polynom.
		if (Q_ab.divisors.find(data.primes.at(i)) != Q_ab.divisors.end()) {
			roots[i].prime = data.primes.at(i);
			roots[i].root1_int = data.M;
			roots[i].root2_int = data.M;
			continue;
		}

		//mpz_mod(a, Q_ab.a, prime);
		//cout << "a = ";
		//mpz_out_str(stdout, 10, a);
		//cout << endl;

		mpz_mod(b, Q_ab.b, prime);
		//cout << "b = ";
		//mpz_out_str(stdout, 10, b);
		//cout << endl;

		//mpz_mod(tmp, data.N, prime);
		//cout << "N = ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		// -t
		//mpz_neg(tmp, t);
		mpz_neg(tmp, data.quadraticResidue[i]);
		//cout << "-t = ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		// a^-1
		//mpz_invert(a_inv, a, prime);
		//cout << "a^-1 = ";
		//mpz_out_str(stdout, 10, a_inv);
		//cout << endl;

		// t - b
		mpz_sub(roots[i].root1, data.quadraticResidue[i], b);
		//cout << "t - b = ";
		//mpz_out_str(stdout, 10, root1);
		//cout << endl;

		// -t - b
		mpz_sub(tmp, tmp, b);
		//cout << "-t - b = ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		// root1
		//mpz_mul(roots[i].root1, a_inv, roots[i].root1);
		mpz_mul(roots[i].root1, Q_ab.a_inv_prime[i], roots[i].root1);
		//cout << "root1 = ";
		//mpz_out_str(stdout, 10, root1);
		//cout << endl;

		mpz_mod(roots[i].root1, roots[i].root1, prime);
		//cout << "root1 mod prime = ";
		//mpz_out_str(stdout, 10, root1);
		//cout << endl;

		// root2
		//mpz_mul(tmp, a_inv, tmp);
		mpz_mul(tmp, Q_ab.a_inv_prime[i], tmp);
		//cout << "root2 = ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		mpz_mod(roots[i].root2, tmp, prime);
		//cout << "root2 mod prime = ";
		//mpz_out_str(stdout, 10, root2);
		//cout << endl;

		roots[i].prime = data.primes.at(i);
		roots[i].root1_int = mpz_get_si(roots[i].root1);
		roots[i].root2_int = mpz_get_si(roots[i].root2);
	}

	// Cleaning
	//mpz_clear(a_inv);
	mpz_clear(prime);
	//mpz_clear(a);
	mpz_clear(b);
	mpz_clear(tmp);

	return true;
}

void SieveValues(data_t &data, root_t *roots, vector<uint32_t> &xCandidates, vector<uint32_t> &xNegCandidates, pol_t &Q_ab, double *xValues, double *xNegValues)
{
	bool skip = false;
	int32_t M = data.M;
	uint32_t rootsSize = data.primes.size();
	double eps;

#ifdef DEBUG_PRINT
	cout << omp_get_thread_num() << " SieveValues()" << endl;
#endif

	if (data.divNum <= 3)
		eps = 2.0;
	else {
		// We are not sieving with the coeficient a divisors, so we need to add some error.
		// Error is equal to the log of the biggest divisor
		map<uint64_t, uint64_t>::iterator it = Q_ab.divisors.end();
		it--;
		eps = data.logs[it->first];
	}

	for (uint32_t i = 0; i < data.M; ++i) {
		xValues[i] = 0.0;
		xNegValues[i] = 0.0;
	}

	int64_t root1;
	int64_t root2;
	int64_t negRoot1;
	int64_t negRoot2;

	uint32_t i = 0;
	if (roots[i].prime == 2)
	{
		double primeLog = data.logs[roots[i].prime];
		root1 = roots[i].root1_int;
		negRoot1 = roots[i].root1_int - roots[i].prime;
		negRoot1 = abs(negRoot1);
		for (; root1 < M && negRoot1 < M; root1 += roots[i].prime, negRoot1 += roots[i].prime) {
			xValues[root1] += primeLog;
			xNegValues[negRoot1] += primeLog;
		}

		if (root1 < M)
			xValues[root1] += primeLog;
		if (negRoot1 < M)
			xNegValues[negRoot1] += primeLog;

		//negRoot1 = roots[i].root1_int - roots[i].prime;
		//for (root1 = abs(negRoot1); root1 < M; root1 += roots[i].prime)
		//	xNegValues[root1] += data.logs[roots[i].prime];

		++i;
	}

	for (; i < rootsSize; ++i) {
		double primeLog = data.logs[roots[i].prime];
		root1 = roots[i].root1_int;
		root2 = roots[i].root2_int;
		negRoot1 = roots[i].root1_int - roots[i].prime;
		negRoot1 = abs(negRoot1);
		negRoot2 = roots[i].root2_int - roots[i].prime;
		negRoot2 = abs(negRoot2);
		for (; root1 < M && root2 < M && negRoot1 < M && negRoot2 < M; root1 += roots[i].prime, root2 += roots[i].prime, negRoot1 += roots[i].prime, negRoot2 += roots[i].prime) {
			xValues[root1] += primeLog;
			xValues[root2] += primeLog;
			xNegValues[negRoot1] += primeLog;
			xNegValues[negRoot2] += primeLog;
		}
		if (root1 < M)
			xValues[root1] += primeLog;
		if (root2 < M)
			xValues[root2] += primeLog;
		if (negRoot1 < M)
			xNegValues[negRoot1] += primeLog;
		if (negRoot2 < M)
			xNegValues[negRoot2] += primeLog;

		//for (; root < M; root += roots[i].prime)
		//	xValues[root] += primeLog;

		//for (root = abs(negRoot1); root < M; root += roots[i].prime)
		//	xNegValues[root] += primeLog;

		//for (root = abs(negRoot2); root < M; root += roots[i].prime)
		//	xNegValues[root] += primeLog;
	}

#ifdef LOG_VALUES_ENABLE
	ofstream xDivisorsFile;

	string fileName = "xDivisors_" + to_string(omp_get_thread_num()) + "_" + mpz_get_str(NULL, 10, Q_ab.a) + ".txt";

	xDivisorsFile.open(fileName, ios::out | ios::app);
	//xDivisorsFile.open("xDivisors.txt");
	xDivisorsFile << "a: " << mpz_get_str(NULL, 10, Q_ab.a);
	xDivisorsFile << "\tb: " << mpz_get_str(NULL, 10, Q_ab.b[Q_ab.acc_pol]);
	xDivisorsFile << "\tc: " << mpz_get_str(NULL, 10, Q_ab.c);
	xDivisorsFile << endl;
	for (map<int64_t, vector<uint64_t>>::iterator it = xDivisors.begin(); it != xDivisors.end(); ++it) {
		xDivisorsFile << "x: " << it->first << "\t log: " << xValues[it->first] << "\t largePrime: " << xValues[it->first] + data.logLargestPrime;
		if (it->first != 0)
			xDivisorsFile << "\t Threshold: " << (log(2 * (it->first)*sqrt(mpz_get_d(data.N))) - eps);
		xDivisorsFile << endl;
		xDivisorsFile << "\t";
		for (vector<uint64_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			xDivisorsFile << *it2 << " ";
		xDivisorsFile << endl;
	}
	xDivisorsFile.close();

	ofstream xNegDivisorsFile;

	string fileNegName = "xNegDivisors_" + to_string(omp_get_thread_num()) + "_" + mpz_get_str(NULL, 10, Q_ab.a) + ".txt";

	xNegDivisorsFile.open(fileNegName, ios::out | ios::app);
	//xDivisorsFile.open("xDivisors.txt");
	xNegDivisorsFile << "a: " << mpz_get_str(NULL, 10, Q_ab.a);
	xNegDivisorsFile << "\tb: " << mpz_get_str(NULL, 10, Q_ab.b[Q_ab.acc_pol]);
	xNegDivisorsFile << "\tc: " << mpz_get_str(NULL, 10, Q_ab.c);
	xNegDivisorsFile << endl;
	for (map<int64_t, vector<uint64_t>>::iterator it = xNegDivisors.begin(); it != xNegDivisors.end(); ++it) {
		xNegDivisorsFile << "x: " << it->first << "\t log: " << xNegValues[it->first] << "\t largePrime: " << xNegValues[it->first] + data.logLargestPrime;
		if (it->first != 0)
			xNegDivisorsFile << "\t Threshold: " << (log(2 * (it->first)*sqrt(mpz_get_d(data.N))) - eps);
		xNegDivisorsFile << endl;
		xNegDivisorsFile << "\t";
		for (vector<uint64_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			xNegDivisorsFile << *it2 << " ";
		xNegDivisorsFile << endl;
	}
	xNegDivisorsFile.close();
#endif

	// Save candidates
	
	for (uint32_t x = 1; x < data.M; ++x) {
		xValues[x] += data.logLargestPrime;
		xNegValues[x] += data.logLargestPrime;

		if (xValues[x] >= (data.threshold[x] - eps))
			xCandidates.push_back(x);

		if (xNegValues[x] >= (data.threshold[x] - eps))
			xNegCandidates.push_back(x);
	}

	return;
}

void TrialDivision(data_t &data, vector<uint32_t> &xCandidates, vector<uint32_t> &xNegCandidates, pol_t &Q_ab) {
	mpz_t tmp;
	mpz_t tmp2;
	mpz_t a_inv;

	mpz_t quotient;
	mpz_t remainder;

	bool negative = false;

	uint64_t largePrime;

	vector<uint64_t> divisors;
	map<uint64_t, uint32_t> relationPrimes;
	
	vector<uint32_t> polynomDivisorsIndexes;

#ifdef DEBUG_PRINT
	cout << "TrialDivison()" << endl;
#endif

	for (map<uint64_t, uint64_t>::iterator it = Q_ab.divisors.begin(); it != Q_ab.divisors.end(); ++it) {
		for (uint32_t i = 0; i < data.primes.size(); ++i) {
			if (data.primes.at(i) == it->first) {
				polynomDivisorsIndexes.push_back(i);
				break;
			}
		}
	}

#ifdef LOG_ENABLE
	ofstream polDivFile;

	string fileName = "polDiv_" + to_string(omp_get_thread_num()) + ".txt";

	polDivFile.open(fileName, ios::out | ios::app);

	for (uint32_t i = 0; i < data.primes.size(); ++i) {
		if (i < 10)
			polDivFile << "    ";
		else if (i < 100)
			polDivFile << "   ";
		else if (i < 1000)
			polDivFile << "  ";
		else if (i < 10000)
			polDivFile << " ";
		polDivFile << i << " ";
	}
	polDivFile << endl;

	for (uint32_t i = 0; i < data.primes.size(); ++i) {
		if (data.primes.at(i) < 10)
			polDivFile << "    ";
		else if (data.primes.at(i) < 100)
			polDivFile << "   ";
		else if (data.primes.at(i) < 1000)
			polDivFile << "  ";
		else if (data.primes.at(i) < 10000)
			polDivFile << " ";
		polDivFile << data.primes.at(i) << " ";
	}
	polDivFile << endl;

	for (uint32_t j = 0; j < polynomDivisorsIndexes.size(); ++j) {
		polDivFile << "Divisor: " << data.primes.at(polynomDivisorsIndexes.at(j)) << endl;
	}
	polDivFile << endl;
	polDivFile.close();
#endif

	mpz_init(tmp);
	mpz_init(tmp2);
	mpz_init(a_inv);

	mpz_init(quotient);
	mpz_init(remainder);

	for (uint32_t i = 0; i < data.primes.size(); ++i)
		relationPrimes.insert(pair<uint64_t, uint32_t>(data.primes.at(i), 0));

	for (uint32_t j = 0; j < polynomDivisorsIndexes.size(); ++j)
		relationPrimes.at(data.primes.at(polynomDivisorsIndexes.at(j))) = 1;

	//cout << "a = ";
	//mpz_out_str(stdout, 10, Q_ab.a);
	//cout << endl;

	//cout << "b = ";
	//mpz_out_str(stdout, 10, Q_ab.b[Q_ab.acc_pol]);
	//cout << endl;

	//cout << "c = ";
	//mpz_out_str(stdout, 10, Q_ab.c);
	//cout << endl;

	for (vector<uint32_t>::iterator it = xCandidates.begin(); it != xCandidates.end(); ++it) {
		//cout << "x: " << it->first << endl;

		mpz_mul_si(tmp, Q_ab.a, *it);
		//mpz_mul_si(tmp, Q_ab.a, -65490);
		//cout << "ax = ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		mpz_mul_ui(tmp2, Q_ab.b, 2);
		//mpz_mul_ui(tmp2, tmp2, 2);
		//cout << "2b = ";
		//mpz_out_str(stdout, 10, tmp2);
		//cout << endl;

		mpz_add(tmp, tmp, tmp2);
		//cout << "ax + 2b = ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		mpz_mul_si(tmp, tmp, *it);
		//mpz_mul_si(tmp, tmp, -65490);
		//cout << "(ax + 2b)x = ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		mpz_add(tmp, tmp, Q_ab.c);
		//cout << "(ax + 2b)x + c = ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		for (uint32_t i = 0; i < data.primes.size(); i++) {
			if (Q_ab.divisors.find(data.primes.at(i)) != Q_ab.divisors.end())
				continue;

			while (1) {
				mpz_tdiv_qr(quotient, remainder, tmp, data.primesMPZ[i]);
				if (mpz_cmp_ui(remainder, 0) == 0) {
					mpz_set(tmp, quotient);
					divisors.push_back(data.primes.at(i));
					relationPrimes.at(data.primes.at(i)) = (relationPrimes.at(data.primes.at(i))) ^ 1;
				}
				else {
					break;
				}
			}
		}

		if (mpz_cmp_ui(tmp, 0) == -1) {
			negative = true;
			mpz_neg(tmp, tmp);
		}

		if (mpz_cmp_ui(tmp, 1) == 0) {
			relInfoUInt_t relInfo;

			uint32_t *relation = (uint32_t *)malloc(sizeof(uint32_t)*data.expVectorLengthAlloc);
			for (uint32_t j = 0; j < data.expVectorLengthAlloc; ++j) {
				relation[j] = 0;
			}

			uint32_t i = 0;
			uint32_t blockIndex = 0;
			uint32_t index = 0;
			for (map<uint64_t, uint32_t>::iterator it2 = relationPrimes.begin(); it2 != relationPrimes.end(); ++it2, ++i) {
				relation[blockIndex] <<= 1;
				if (it2->second > 0)
					relation[blockIndex] |= 1;
				else
					relation[blockIndex] |= 0;

				++index;
				if (index == 32) {
					index = 0;
					++blockIndex;
				}
			}

			// Including negative
			relation[blockIndex] <<= 1;
			if (negative)
				relation[blockIndex] |= 1;
			else
				relation[blockIndex] |= 0;

			++index;
			if (index == 32) {
				index = 0;
				++blockIndex;
			}

			uint32_t diff = 32 - index;
			relation[blockIndex] <<= diff;

#ifdef LOG_ENABLE
			ofstream relFile;

			string file2Name = "rel_" + to_string(omp_get_thread_num()) + ".txt";

			relFile.open(file2Name, ios::out | ios::app);

			for (uint32_t j = 0; j < data.primes.size(); ++j) {
				if (data.primes.at(j) < 10)
					relFile << "    ";
				else if (data.primes.at(j) < 100)
					relFile << "   ";
				else if (data.primes.at(j) < 1000)
					relFile << "  ";
				else if (data.primes.at(j) < 10000)
					relFile << " ";
				relFile << data.primes.at(j) << " ";
			}
			relFile << endl;

			for (uint32_t j = 0; j < relationPrimes.size(); ++j) {
				relFile << "    ";
				if (relation[j])
					relFile << "1 ";
				else
					relFile << "0 ";
			}
			relFile << endl;
			relFile << "X: " << it->first << endl;
			relFile.close();
#endif
			
			// Save info about relation
			relInfo.x = *it;
			relInfo.composed = false;
			relInfo.relation = relation;
			mpz_init_set(relInfo.a, Q_ab.a);
			mpz_init_set(relInfo.b, Q_ab.b);
			mpz_init_set(relInfo.c, Q_ab.c);

#ifdef LOG_ENABLE_SMOOTH
			ofstream smoothFile;

			string fileName = "smooth_" + to_string(omp_get_thread_num()) + ".txt";

			smoothFile.open(fileName, ios::out | ios::app);

			smoothFile << data.relationsNum << ":\t";
			smoothFile << "a: " << mpz_get_str(NULL, 10, relInfo.a);
			smoothFile << "\tb: " << mpz_get_str(NULL, 10, relInfo.b);
			smoothFile << "\tc: " << mpz_get_str(NULL, 10, relInfo.c);
			smoothFile << "\tx: " << it->first;
			smoothFile << "\t Divisors: ";
			for (uint32_t i = 0; i < data.primes.size(); ++i) {
				if (relInfo.relation[i])
					smoothFile << data.primes.at(i) << " ";
			}
			if (relInfo.relation[data.primes.size()])
				smoothFile << "-1";
			smoothFile << endl;

			smoothFile.close();
#endif

#ifdef NO_PARALLEL
			data.relations.push_back(relInfo);
#else
#pragma omp critical
			{
				data.relations.push_back(relInfo);
			}
#endif
		}
		else if (mpz_cmp_ui(tmp, UINT64_MAX) < 0) {
			largePrime = mpz_get_ui(tmp);
			if (data.largePrimes.find(largePrime) != data.largePrimes.end()) {
				relInfoUInt_t parRelInfo;

				uint32_t *parRelation = (uint32_t *)malloc(sizeof(uint32_t)*data.expVectorLengthAlloc);
				for (uint32_t j = 0; j < data.expVectorLengthAlloc; ++j) {
					parRelation[j] = 0;
				}

				uint32_t i = 0;
				uint32_t blockIndex = 0;
				uint32_t index = 0;
				for (map<uint64_t, uint32_t>::iterator it2 = relationPrimes.begin(); it2 != relationPrimes.end(); ++it2, ++i) {
					parRelation[blockIndex] <<= 1;
					if(it2->second > 0)
						parRelation[blockIndex] |= 1;
					else
						parRelation[blockIndex] |= 0;

					++index;
					if (index == 32) {
						index = 0;
						++blockIndex;
					}
				}

				// Including negative
				parRelation[blockIndex] <<= 1;
				if (negative)
					parRelation[blockIndex] |= 1;
				else
					parRelation[blockIndex] |= 0;

				++index;
				if (index == 32) {
					index = 0;
					++blockIndex;
				}

				uint32_t diff = 32 - index;
				parRelation[blockIndex] <<= diff;

				parRelInfo.x = *it;
				parRelInfo.composed = false;
				parRelInfo.relation = parRelation;
				mpz_init_set(parRelInfo.a, Q_ab.a);
				mpz_init_set(parRelInfo.b, Q_ab.b);
				mpz_init_set(parRelInfo.c, Q_ab.c);

#ifdef NO_PARALLEL
				data.partialRelationsUInt.insert(pair<uint64_t, relInfoUInt_t>(largePrime, parRelInfo));
#else
#pragma omp critical
				{
					data.partialRelationsUInt.insert(pair<uint64_t, relInfoUInt_t>(largePrime, parRelInfo));
				}
#endif

#ifdef LOG_ENABLE_PAR
				ofstream parRelFile;

				string file2Name = "parRel_" + to_string(omp_get_thread_num()) + ".txt";

				parRelFile.open(file2Name, ios::out | ios::app);

				for (uint32_t j = 0; j < data.primes.size(); ++j) {
					if (data.primes.at(j) < 10)
						parRelFile << "    ";
					else if (data.primes.at(j) < 100)
						parRelFile << "   ";
					else if (data.primes.at(j) < 1000)
						parRelFile << "  ";
					else if (data.primes.at(j) < 10000)
						parRelFile << " ";
					parRelFile << data.primes.at(j) << " ";
				}
				parRelFile << "   ";
				parRelFile << "-1";
				parRelFile << endl;

				for (uint32_t j = 0; j < relationPrimes.size() + 1; ++j) {
					parRelFile << "    ";
					if (parRelation[j])
						parRelFile << "1 ";
					else
						parRelFile << "0 ";
				}
				parRelFile << endl;

				blockIndex = 0;
				index = 0;
				uint32_t tmp = 0;
				uint32_t block = parRelation[blockIndex];
				for (uint32_t k = 0; k < 32; ++k) {
					tmp <<= 1;
					tmp |= (block & 1);
					block >>= 1;
				}
				for (uint32_t j = 0; j < data.primes.size() + 1; ++j) {
					parRelFile << "    ";
					parRelFile << (tmp & 1) << " ";

					tmp >>= 1;
					++index;
					if (index == 32) {
						++blockIndex;
						index = 0;
						block = parRelation[blockIndex];
						tmp = 0;
						for (uint32_t k = 0; k < 32; ++k) {
							tmp <<= 1;
							tmp |= (block & 1);
							block >>= 1;
						}
					}
				}
				parRelFile << endl;

				parRelFile << "X: " << it->first << endl;
				parRelFile.close();
#endif
			}
		}

		for (vector<uint64_t>::iterator it2 = divisors.begin(); it2 != divisors.end(); ++it2)
			relationPrimes.at(*it2) = 0;

		divisors.clear();
		negative = false;
	}

	for (vector<uint32_t>::iterator it = xNegCandidates.begin(); it != xNegCandidates.end(); ++it) {
		//cout << "x: " << it->first << endl;

		mpz_mul_si(tmp, Q_ab.a, (*it)*(-1));
		//mpz_mul_si(tmp, Q_ab.a, -65490);
		//cout << "ax = ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		mpz_mul_ui(tmp2, Q_ab.b, 2);
		//mpz_mul_ui(tmp2, tmp2, 2);
		//cout << "2b = ";
		//mpz_out_str(stdout, 10, tmp2);
		//cout << endl;

		mpz_add(tmp, tmp, tmp2);
		//cout << "ax + 2b = ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		mpz_mul_si(tmp, tmp, (*it)*(-1));
		//mpz_mul_si(tmp, tmp, -65490);
		//cout << "(ax + 2b)x = ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		mpz_add(tmp, tmp, Q_ab.c);
		//cout << "(ax + 2b)x + c = ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;

		for (uint32_t i = 0; i < data.primes.size(); i++) {
			if (Q_ab.divisors.find(data.primes.at(i)) != Q_ab.divisors.end())
				continue;

			while (1) {
				mpz_tdiv_qr(quotient, remainder, tmp, data.primesMPZ[i]);
				if (mpz_cmp_ui(remainder, 0) == 0) {
					mpz_set(tmp, quotient);
					divisors.push_back(data.primes.at(i));
					relationPrimes.at(data.primes.at(i)) = (relationPrimes.at(data.primes.at(i))) ^ 1;
				}
				else {
					break;
				}
			}
		}

		if (mpz_cmp_ui(tmp, 0) == -1) {
			negative = true;
			mpz_neg(tmp, tmp);
		}

		if (mpz_cmp_ui(tmp, 1) == 0) {
			relInfoUInt_t relInfo;

			uint32_t *relation = (uint32_t *)malloc(sizeof(uint32_t)*data.expVectorLengthAlloc);
			for (uint32_t j = 0; j < data.expVectorLengthAlloc; ++j) {
				relation[j] = 0;
			}

			uint32_t i = 0;
			uint32_t blockIndex = 0;
			uint32_t index = 0;
			for (map<uint64_t, uint32_t>::iterator it2 = relationPrimes.begin(); it2 != relationPrimes.end(); ++it2, ++i) {
				relation[blockIndex] <<= 1;
				if (it2->second > 0)
					relation[blockIndex] |= 1;
				else
					relation[blockIndex] |= 0;

				++index;
				if (index == 32) {
					index = 0;
					++blockIndex;
				}
			}

			// Including negative
			relation[blockIndex] <<= 1;
			if (negative)
				relation[blockIndex] |= 1;
			else
				relation[blockIndex] |= 0;

			++index;
			if (index == 32) {
				index = 0;
				++blockIndex;
			}

			uint32_t diff = 32 - index;
			relation[blockIndex] <<= diff;
			
			// Save info about relation
			relInfo.x = (*it)*(-1);
			relInfo.composed = false;
			relInfo.relation = relation;
			mpz_init_set(relInfo.a, Q_ab.a);
			mpz_init_set(relInfo.b, Q_ab.b);
			mpz_init_set(relInfo.c, Q_ab.c);

#ifdef LOG_ENABLE_SMOOTH
			ofstream smoothFile;

			string fileName = "smooth_" + to_string(omp_get_thread_num()) + ".txt";

			smoothFile.open(fileName, ios::out | ios::app);

			smoothFile << data.relationsNum << ":\t";
			smoothFile << "a: " << mpz_get_str(NULL, 10, relInfo.a);
			smoothFile << "\tb: " << mpz_get_str(NULL, 10, relInfo.b);
			smoothFile << "\tc: " << mpz_get_str(NULL, 10, relInfo.c);
			smoothFile << "\tx: " << it->first*(-1);
			smoothFile << "\t Divisors: ";
			for (uint32_t i = 0; i < data.primes.size(); ++i) {
				if (relInfo.relation[i])
					smoothFile << data.primes.at(i) << " ";
			}
			if (relInfo.relation[data.primes.size()])
				smoothFile << "-1";
			smoothFile << endl;

			smoothFile.close();
#endif

#ifdef NO_PARALLEL
			data.relations.push_back(relInfo);
#else
#pragma omp critical
			{
				data.relations.push_back(relInfo);
			}
#endif
		}
		else if (mpz_cmp_ui(tmp, UINT64_MAX) < 0) {
			largePrime = mpz_get_ui(tmp);
			if (data.largePrimes.find(largePrime) != data.largePrimes.end()) {
				relInfoUInt_t parRelInfo;

				uint32_t *parRelation = (uint32_t *)malloc(sizeof(uint32_t)*data.expVectorLengthAlloc);
				for (uint32_t j = 0; j < data.expVectorLengthAlloc; ++j) {
					parRelation[j] = 0;
				}

				uint32_t i = 0;
				uint32_t blockIndex = 0;
				uint32_t index = 0;
				for (map<uint64_t, uint32_t>::iterator it2 = relationPrimes.begin(); it2 != relationPrimes.end(); ++it2, ++i) {
					parRelation[blockIndex] <<= 1;
					if (it2->second > 0)
						parRelation[blockIndex] |= 1;
					else
						parRelation[blockIndex] |= 0;

					++index;
					if (index == 32) {
						index = 0;
						++blockIndex;
					}
				}

				// Including negative
				parRelation[blockIndex] <<= 1;
				if (negative)
					parRelation[blockIndex] |= 1;
				else
					parRelation[blockIndex] |= 0;

				++index;
				if (index == 32) {
					index = 0;
					++blockIndex;
				}

				uint32_t diff = 32 - index;
				parRelation[blockIndex] <<= diff;

				parRelInfo.x = (*it)*(-1);
				parRelInfo.composed = false;
				parRelInfo.relation = parRelation;
				mpz_init_set(parRelInfo.a, Q_ab.a);
				mpz_init_set(parRelInfo.b, Q_ab.b);
				mpz_init_set(parRelInfo.c, Q_ab.c);

#ifdef NO_PARALLEL
				data.partialRelationsUInt.insert(pair<uint64_t, relInfoUInt_t>(largePrime, parRelInfo));
#else
#pragma omp critical
				{
					data.partialRelationsUInt.insert(pair<uint64_t, relInfoUInt_t>(largePrime, parRelInfo));
				}
#endif
#ifdef LOG_ENABLE_PAR
				ofstream parRelFile;

				string file2Name = "parRel_" + to_string(omp_get_thread_num()) + ".txt";

				parRelFile.open(file2Name, ios::out | ios::app);

				for (uint32_t j = 0; j < data.primes.size(); ++j) {
					if (data.primes.at(j) < 10)
						parRelFile << "    ";
					else if (data.primes.at(j) < 100)
						parRelFile << "   ";
					else if (data.primes.at(j) < 1000)
						parRelFile << "  ";
					else if (data.primes.at(j) < 10000)
						parRelFile << " ";
					parRelFile << data.primes.at(j) << " ";
				}
				parRelFile << "   ";
				parRelFile << "-1";
				parRelFile << endl;

				for (uint32_t j = 0; j < relationPrimes.size() + 1; ++j) {
					parRelFile << "    ";
					if (parRelation[j])
						parRelFile << "1 ";
					else
						parRelFile << "0 ";
				}
				parRelFile << endl;

				blockIndex = 0;
				index = 0;
				uint32_t tmp = 0;
				uint32_t block = parRelation[blockIndex];
				for (uint32_t k = 0; k < 32; ++k) {
					tmp <<= 1;
					tmp |= (block & 1);
					block >>= 1;
				}
				for (uint32_t j = 0; j < data.primes.size() + 1; ++j) {
					parRelFile << "    ";
					parRelFile << (tmp & 1) << " ";

					tmp >>= 1;
					++index;
					if (index == 32) {
						++blockIndex;
						index = 0;
						block = parRelation[blockIndex];
						tmp = 0;
						for (uint32_t k = 0; k < 32; ++k) {
							tmp <<= 1;
							tmp |= (block & 1);
							block >>= 1;
						}
					}
				}
				parRelFile << endl;

				parRelFile << "X: " << (it->first)*(-1) << endl;
				parRelFile.close();
#endif
			}
		}

		for (vector<uint64_t>::iterator it2 = divisors.begin(); it2 != divisors.end(); ++it2)
			relationPrimes.at(*it2) = 0;

		divisors.clear();
		negative = false;
	}

#ifdef NO_PARALLEL
	cout << string(80, '\b');
	cout << "Relations gathered: " << data.relations.size() << " (" << data.partialRelationsUInt.size() << ") / " << data.relationsNeededPart;
#else
	cout << "Relations gathered: " << data.relations.size() << " (" << data.partialRelationsUInt.size() << ") / " << data.relationsNeededPart << endl;
#endif

	relationPrimes.clear();
	polynomDivisorsIndexes.clear();

	// cleaning
	mpz_clear(tmp);
	mpz_clear(tmp2);
	mpz_clear(a_inv);

	mpz_clear(quotient);
	mpz_clear(remainder);

	return;
}

bool Sieve(data_t &data, pol_t &Q_ab)
{
	root_t *roots;

	bool init = true;

	int32_t e;
	uint32_t v;

	vector<uint32_t> xCandidates;
	vector<uint32_t> xNegCandidates;

	// Init roots
	roots = (root_t *)malloc(sizeof(root_t)*data.primes.size());
	for (uint32_t i = 0; i < data.primes.size(); ++i) {
		mpz_init(roots[i].root1);
		mpz_init(roots[i].root2);
	}

	double *xValues = (double *)malloc(sizeof(double)*data.M);
	if (xValues == NULL) {
		cout << omp_get_thread_num() << " SieveValues() - not enough memory for doubles" << endl;
		cout << "M:" << data.M << endl;
	}

	double *xNegValues = (double *)malloc(sizeof(double)*data.M);
	if (xNegValues == NULL) {
		cout << omp_get_thread_num() << " SieveValues() - not enough memory for doubles" << endl;
		cout << "M:" << data.M << endl;
	}

	// We try to get only half of needed relatioes, remaining we get from partial ones
	while (data.relations.size() < data.relationsNeededPart) {

		if (!init)
			NextPolynom(Q_ab, data, e, v);
		else
			init = false;

		// Now everytime new root computation
		//if (Q_ab.acc_pol == 0)
			ComputeRoots(data, roots, Q_ab);
		//else
		//	ComputeNextRoots(data, roots, Q_ab, e, v);

		SieveValues(data, roots, xCandidates, xNegCandidates, Q_ab, xValues, xNegValues);

		TrialDivision(data, xCandidates, xNegCandidates, Q_ab);

		xCandidates.clear();
		xNegCandidates.clear();
	}

	for (uint32_t i = 0; i < data.primes.size(); ++i) {
		mpz_clear(roots[i].root1);
		mpz_clear(roots[i].root2);
	}

	free(roots);
	free(xValues);
	free(xNegValues);
#ifdef NO_PARALLEL
	cout << endl;
#endif

	return true;
}

