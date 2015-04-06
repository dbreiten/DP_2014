/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace celych cisel z pohledu lamani RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    siqs.cpp                                                   **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

#include "siqs.h"
#include "utils.h"
#include "polgen.h"
#include "sieve.h"
#include "matrix.h"

using namespace std;

// Table inspired by msieve
sieve_param_t prebuilt_params[] = {
	{ 64, 100, 40, 1 * 65536, 0 },
	{ 101, 1000, 40, 2 * 65536, 0 },
	{ 128, 450, 40, 1 * 65536, 2 },
	{ 132, 1053, 40, 1 * 65536, 5 },
	{ 166, 3000, 40, 1 * 65536, 6 },
	//	{ 183, 2000, 40, 1 * 65536, 0 },
	{ 200, 6000, 50, 1 * 65536, 8 },
	//	{ 200, 3000, 50, 1 * 65536, 8 },
	{ 212, 5400, 50, 3 * 65536, 0 },
	{ 235, 10000, 100, 3 * 65536, 9 },
	{ 249, 27000, 100, 3 * 65536, 0 },
	{ 266, 100000, 100, 3 * 65536, 10 },
	//	{ 266, 50000, 100, 3 * 65536, 0 },

	{ 283, 55000, 80, 3 * 65536, 0 },
	//	{ 298, 60000, 80, 9 * 65536, 0 },
	{ 299, 130000, 80, 9 * 65536, 12 },
	{ 315, 80000, 150, 9 * 65536, 0 },
	//	{ 332, 100000, 150, 9 * 65536, 0 },
	{ 332, 200000, 150, 9 * 65536, 13 },
	{ 348, 140000, 150, 9 * 65536, 0 },
	{ 363, 210000, 150, 13 * 65536, 0 },
	{ 379, 300000, 150, 17 * 65536, 0 },
	{ 395, 400000, 150, 21 * 65536, 0 },
	{ 415, 500000, 150, 25 * 65536, 0 }, /* beyond this point you're crazy */
	{ 440, 700000, 150, 33 * 65536, 0 },
	{ 465, 900000, 150, 50 * 65536, 0 },
	{ 490, 1100000, 150, 75 * 65536, 0 },
	{ 512, 1300000, 150, 100 * 65536, 0 },
};

/* Compute logs for primes from FB */
void ComputeLogs(data_t &data)
{
	double logres;

	for (vector<uint64_t>::iterator it = data.primes.begin(); it != data.primes.end(); ++it) {
		logres = log(*it);
		//data.logs.insert(pair<uint64_t, double>((*it), logres));
		data.logs[*it] = logres;
	}

	double threshold = 2 * sqrt(mpz_get_d(data.N));
	for (uint32_t i = 0; i < data.M; ++i) {
		data.threshold[i] = log(i * threshold);
	}

	return;
}

/** Function loads deltas between primes and saves them to the vector
 *
 * @param primeFileName File with deltas
 * @param primeDeltas   Vector to store deltas
 */
bool GetPrimeDeltas(string primeFileName, data_t &data)
{
	string line;
	int32_t legendre;
	uint32_t delta = 0;
	uint32_t cnt = 0; 
	uint32_t cnt2 = 0;
	mpz_t prime;
	bool siqs = true;

	mpz_init(prime);

	ifstream primesFile(primeFileName);

	if(primesFile.is_open()) {
		// Save primes to factor base
		while(getline(primesFile,line)) {
			delta = stoi(line);
			// Make a prime from the delta and the last prime
			mpz_add_ui(prime, prime, delta);

			// We tried enough of primes
			if (cnt >= data.F) break;

			legendre = mpz_legendre(data.N, prime);
			// Legendre symbol = 1 => prime might divide the polynom
			if (legendre == 1 || mpz_cmp_ui(prime, data.multiplier) == 0) {
				data.primes.push_back(mpz_get_ui(prime));
				mpz_init_set(data.primesMPZ[cnt2], prime);
				mpz_init(data.quadraticResidue[cnt2]);
				quadratic_residue(data.quadraticResidue[cnt2], data.N, prime);

				cnt2++;
			}
			else if (legendre == 0) {
				cout << "Divisor found: \n" << mpz_get_ui(prime) << endl;
				mpz_set(data.N, prime);
				siqs = false;
				return siqs;
			}

			cnt++;
		}

#ifdef DEBUG_PRINT
		cout << "There were " << cnt << " primes in the file, but only " << data.primes.size() << " will be used" << endl;
#endif

		// Now save large primes for Single Large Prime Variation
		uint32_t largePrimeBound = data.largePrimeMultiplier * data.F;
		cnt = cnt2 = 0;
		legendre = mpz_legendre(data.N, prime);
		if (legendre == 1) {
			data.largePrimes.insert(pair<uint64_t, bool>(mpz_get_ui(prime), true));

			cnt2++;
		}
		else if (legendre == 0) {
			cout << "Divisor found: \n" << mpz_get_ui(prime) << endl;
			mpz_set(data.N, prime);
			siqs = false;
			return siqs;
		}

		cnt++;

		while (getline(primesFile, line)) {
			delta = stoi(line);
			// Make a prime from the delta and the last prime
			mpz_add_ui(prime, prime, delta);

			data.largePrimes.insert(pair<uint64_t, bool>(mpz_get_ui(prime), true));

			legendre = mpz_legendre(data.N, prime);

			if (legendre == 0) {
				cout << "Divisor found: \n" << mpz_get_ui(prime) << endl;
				mpz_set(data.N, prime);
				siqs = false;
				return siqs;
			}

			cnt++;
			if (cnt >= largePrimeBound) break;
		}

		data.logLargestPrime = log(mpz_get_ui(prime));

#ifdef DEBUG_PRINT
		cout << "There were " << cnt << " large primes in the file, but only " << data.largePrimes.size() << " will be used" << endl;
#endif

		primesFile.close();
	} else
		cerr << "Unable to open " << primeFileName << endl;

	// cleaning
	mpz_clear(prime);

	return siqs;
}

void InitRelations(data_t &data)
{
	data.relationsNeeded = data.primes.size() * 105 / 100;
	data.relationsNeededPart = data.relationsNeeded / 2;

	data.expVectorLengthAlloc = ceil((data.primes.size() + 3) / (double)32);
	cout << "Num of uint32_t\'s: " << data.expVectorLengthAlloc << endl;

	return;
}

void ComputeResult(vector<vector<uint32_t>> &S, data_t &data, mpz_t res, uint32_t attempt)
{
	mpz_t t;
	mpz_t t2;
	mpz_t u;
	mpz_t u2;
	mpz_t X;
	mpz_t Y;
	mpz_t tmp;
	mpz_t tmp2;
	mpz_t tmpData;

	mpz_init(t);
	mpz_init(t2);
	mpz_init(u);
	mpz_init(u2);
	mpz_init(tmp);
	mpz_init(tmp2);
	mpz_init_set_ui(X, 1);
	mpz_init_set_ui(Y, 1);
	mpz_init_set(tmpData, data.N);
	uint32_t resCnt = 0;
	for (vector<vector<uint32_t>>::iterator it = S.begin(); it != S.end(); ++it) {
		mpz_set_ui(tmp, 0);
		mpz_set_ui(tmp2, 0);
		mpz_set_ui(X, 1);
		mpz_set_ui(Y, 1);
		mpz_set(tmpData, data.N);

#ifdef LOG_ENABLE
		ofstream resultFile;

		string fileName = "result_"+ to_string(attempt) + "_" + to_string(resCnt) + ".txt";

		resultFile.open(fileName, ios::out | ios::app);
#endif

		for (vector<uint32_t>::iterator it2 = (*it).begin(); it2 != (*it).end(); ++it2) {
#ifdef LOG_ENABLE
			resultFile << "Index: " << *it2 << "\t";
			resultFile.flush();
#endif
			if (data.relations.at((*it2)).composed == false) {
				// a
				mpz_set(tmp, data.relations.at((*it2)).a);
				//cout << "a: ";
				//mpz_out_str(stdout, 10, data.relations.at((*it2)).a);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "a: " << mpz_get_str(NULL, 10, tmp);
				resultFile.flush();
#endif

				// b
				mpz_set(tmp, data.relations.at((*it2)).b);
				//cout << "b: ";
				//mpz_out_str(stdout, 10, data.relations.at((*it2)).b);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "\tb: " << mpz_get_str(NULL, 10, tmp);

				// x
				//cout << "x: " << data.relations.at((*it2)).x << endl;
				resultFile << "\tx: " << data.relations.at((*it2)).x;

				// N
				//cout << "N: ";
				//mpz_out_str(stdout, 10, data.N);
				//cout << endl;
				resultFile << "\tN: " << mpz_get_str(NULL, 10, data.N);
				resultFile << "\tRelationIndex: " << (*it2) << "z " << data.relations.size();
				resultFile.flush();
#endif

				// ax
				mpz_mul_si(tmp, data.relations.at((*it2)).a, data.relations.at((*it2)).x);
				//cout << "ax: ";
				//mpz_out_str(stdout, 10, tmp);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "\tax: " << mpz_get_str(NULL, 10, tmp) << endl;
				resultFile.flush();
#endif

				// t = ax + b
				mpz_add(tmp, tmp, data.relations.at((*it2)).b);
				mpz_set(t, tmp);
				//cout << "t = ax + b: ";
				//mpz_out_str(stdout, 10, tmp);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "t = ax + b: " << mpz_get_str(NULL, 10, t) << endl;
				resultFile.flush();
#endif

				mpz_mul(X, X, t);
				//cout << "X: ";
				//mpz_out_str(stdout, 10, X);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "X: " << mpz_get_str(NULL, 10, X) << endl;
				resultFile.flush();
#endif

				// (ax + b)^2
				mpz_mul(tmp, tmp, tmp);
				//cout << "(ax + b)^2: ";
				//mpz_out_str(stdout, 10, tmp);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "(ax + b)^2: " << mpz_get_str(NULL, 10, tmp) << endl;
				resultFile.flush();
#endif

				// u = (ax + b)^2 - N
				mpz_sub(tmp, tmp, data.N);
				mpz_set(u, tmp);
				//cout << "u = (ax + b)^2 - N: ";
				//mpz_out_str(stdout, 10, tmp);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "u = (ax + b)^2 - N: " << mpz_get_str(NULL, 10, u) << endl;
				resultFile.flush();
#endif

				mpz_mul(Y, Y, u);
				//cout << "Y: ";
				//mpz_out_str(stdout, 10, Y);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "Y: " << mpz_get_str(NULL, 10, Y) << endl << endl;
				resultFile.flush();
#endif
			}
			else {
#ifdef LOG_ENABLE
				resultFile << "Composed -";
				resultFile.flush();
#endif
				// a
				mpz_set(tmp, data.relations.at((*it2)).a);
				//cout << "a: ";
				//mpz_out_str(stdout, 10, data.relations.at((*it2)).a);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "\ta: " << mpz_get_str(NULL, 10, tmp);
				resultFile.flush();
#endif

				// a2
				mpz_set(tmp2, data.relations.at((*it2)).a2);
#ifdef LOG_ENABLE
				resultFile << "\ta2: " << mpz_get_str(NULL, 10, tmp2);
				resultFile.flush();
#endif

				// b
				mpz_set(tmp, data.relations.at((*it2)).b);
				//cout << "b: ";
				//mpz_out_str(stdout, 10, data.relations.at((*it2)).b);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "\tb: " << mpz_get_str(NULL, 10, tmp);
				resultFile.flush();
#endif

				// b2
				mpz_set(tmp2, data.relations.at((*it2)).b2);
#ifdef LOG_ENABLE
				resultFile << "\tb2: " << mpz_get_str(NULL, 10, tmp2);

				// x
				//cout << "x: " << data.relations.at((*it2)).x << endl;
				resultFile << "\tx: " << data.relations.at((*it2)).x;

				// x2
				resultFile << "\tx2: " << data.relations.at((*it2)).x2;

				// N
				//cout << "N: ";
				//mpz_out_str(stdout, 10, data.N);
				//cout << endl;
				resultFile << "\tN: " << mpz_get_str(NULL, 10, data.N);
				resultFile << "\tRelationIndex: " << (*it2) << "z " << data.relations.size();
				resultFile.flush();
#endif

				// ax
				mpz_mul_si(tmp, data.relations.at((*it2)).a, data.relations.at((*it2)).x);
				//cout << "ax: ";
				//mpz_out_str(stdout, 10, tmp);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "\tax: " << mpz_get_str(NULL, 10, tmp);
				resultFile.flush();
#endif

				// a2x2
				mpz_mul_si(tmp2, data.relations.at((*it2)).a2, data.relations.at((*it2)).x2);
#ifdef LOG_ENABLE
				resultFile << "\ta2x2: " << mpz_get_str(NULL, 10, tmp2) << endl;
				resultFile.flush();
#endif

				// t = ax + b
				mpz_add(tmp, tmp, data.relations.at((*it2)).b);
				mpz_set(t, tmp);
				//cout << "t = ax + b: ";
				//mpz_out_str(stdout, 10, tmp);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "t = ax + b: " << mpz_get_str(NULL, 10, t) << endl;
				resultFile.flush();
#endif

				// t2 = a2x2 + b2
				mpz_add(tmp2, tmp2, data.relations.at((*it2)).b2);
				mpz_set(t2, tmp2);
#ifdef LOG_ENABLE
				resultFile << "t2 = a2x2 + b2: " << mpz_get_str(NULL, 10, t2) << endl;
				resultFile.flush();
#endif

				// t = t*t2
				mpz_mul(t, t, t2);
#ifdef LOG_ENABLE
				resultFile << "t = t * t2: " << mpz_get_str(NULL, 10, t) << endl;
				resultFile.flush();
#endif

				mpz_mul(X, X, t);
				//cout << "X: ";
				//mpz_out_str(stdout, 10, X);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "X: " << mpz_get_str(NULL, 10, X) << endl;
				resultFile.flush();
#endif

				// (ax + b)^2
				mpz_mul(tmp, tmp, tmp);
				//cout << "(ax + b)^2: ";
				//mpz_out_str(stdout, 10, tmp);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "(ax + b)^2: " << mpz_get_str(NULL, 10, tmp) << endl;
				resultFile.flush();
#endif

				// (a2x2 + b2)^2
				mpz_mul(tmp2, tmp2, tmp2);
#ifdef LOG_ENABLE
				resultFile << "(a2x2 + b2)^2: " << mpz_get_str(NULL, 10, tmp2) << endl;
				resultFile.flush();
#endif

				// u = (ax + b)^2 - N
				mpz_sub(tmp, tmp, data.N);
				mpz_set(u, tmp);
				//cout << "u = (ax + b)^2 - N: ";
				//mpz_out_str(stdout, 10, tmp);
				//cout << endl;
#ifdef LOG_ENABLE
				resultFile << "u = (ax + b)^2 - N: " << mpz_get_str(NULL, 10, u) << endl;
				resultFile.flush();
#endif

				// u2 = (a2x2 + b2)^2 - N
				mpz_sub(tmp2, tmp2, data.N);
				mpz_set(u2, tmp2);
#ifdef LOG_ENABLE
				resultFile << "u2 = (a2x2 + b2)^2 - N: " << mpz_get_str(NULL, 10, u2) << endl;
				resultFile.flush();
#endif

				// u = u * u2
				mpz_mul(u, u, u2);
#ifdef LOG_ENABLE
				resultFile << "u = u * u2: " << mpz_get_str(NULL, 10, u) << endl;
				resultFile.flush();
#endif

				mpz_mul(Y, Y, u);
#ifdef LOG_ENABLE
				resultFile << "Y: " << mpz_get_str(NULL, 10, Y) << endl << endl;
				resultFile.flush();
#endif
			}
		}

#ifdef LOG_ENABLE
		resultFile << "Going to square" << endl;
		resultFile.flush();
#endif
		mpz_sqrt(Y, Y);
		//cout << "sqrt(Y^2): ";
		//mpz_out_str(stdout, 10, Y);
		//cout << endl;
#ifdef LOG_ENABLE
		resultFile << "sqrt(Y): " << mpz_get_str(NULL, 10, Y) << endl;
		resultFile.flush();
#endif

		mpz_mod(X, X, data.N);
		//cout << "X mod n: ";
		//mpz_out_str(stdout, 10, X);
		//cout << endl;
#ifdef LOG_ENABLE
		resultFile << "X mod n: " << mpz_get_str(NULL, 10, X) << endl;
		resultFile.flush();
#endif

		mpz_mod(Y, Y, data.N);
		//cout << "Y mod n: ";
		//mpz_out_str(stdout, 10, Y);
		//cout << endl;
#ifdef LOG_ENABLE
		resultFile << "Y mod n: " << mpz_get_str(NULL, 10, Y) << endl;
		resultFile.flush();
#endif

		mpz_sub(tmp, X, Y);
		//cout << "X - Y: ";
		//mpz_out_str(stdout, 10, tmp);
		//cout << endl;
#ifdef LOG_ENABLE
		resultFile << "X - Y: " << mpz_get_str(NULL, 10, tmp) << endl;
		resultFile.flush();
#endif

		if (mpz_cmp_ui(tmp, 0) == 0)
			continue;

		// GCD(X - Y, N)
		Euklid_algorithm_BIG(tmp, tmpData, res);
		cout << "Divisor: ";
		mpz_out_str(stdout, 10, res);
		cout << endl;
#ifdef LOG_ENABLE
		resultFile << "Divisor: " << mpz_get_str(NULL, 10, res) << endl;
		resultFile << "==============================" << endl;
		resultFile.flush();

		resultFile.close();
#endif
		resCnt++;

		// Divisor found
		if (mpz_cmp_ui(res, 1) != 0 && mpz_cmp_ui(res, data.multiplier) != 0 && mpz_cmp(res, data.origN) != 0)
			break;
	}

	cout << "Real Divisor: ";
	if (data.multiplier > 1) {
		mpz_tdiv_qr_ui(tmp, tmpData, res, data.multiplier);
		if (mpz_cmp_ui(tmpData, 0) == 0)
			mpz_set(res, tmp);
	}
	mpz_out_str(stdout, 10, res);
	cout << endl;

	mpz_clear(X);
	mpz_clear(Y);
	mpz_clear(u);
	mpz_clear(u2);
	mpz_clear(t);
	mpz_clear(t2);
	mpz_clear(tmp);
	mpz_clear(tmp2);
	mpz_clear(tmpData);

	return;
}

void GetMultPrimes(vector<uint32_t> &multPrimes, string primeFileName)
{
	mpz_t prime;
	string line;
	uint32_t delta = 0;

	mpz_init(prime);

	ifstream primesFile(primeFileName);

	if (primesFile.is_open()) {
		for (uint32_t i = 0; i < 256; ++i) {
			getline(primesFile, line);
			delta = stoi(line);
			// Make a prime from the delta and the last prime
			mpz_add_ui(prime, prime, delta);

			multPrimes.push_back(mpz_get_ui(prime));
		}
		primesFile.close();
	}
	else
		cerr << "Unable to open " << primeFileName << endl;

	 // cleaning
	 mpz_clear(prime);

	 return;
}

// Function used by Knuth-Schroeppel Multiplier algorithm
double g(uint32_t p, uint32_t k, data_t &data, mpz_t nmodp, mpz_t kn, mpz_t knmodp)
{
	if (p == 2) {
		mpz_t knmod8;
		uint32_t knmod8_u;

		mpz_init(knmod8);

		mpz_mod_ui(knmod8, kn, 8);
		knmod8_u = mpz_get_ui(knmod8);
		mpz_clear(knmod8);

		if (knmod8_u == 1)
			return 2.0;
		else if (knmod8_u == 5)
			return 1.0;
		else if (knmod8_u == 3 || knmod8_u == 7)
			return 0.5;
		else
			return 0.0;
	}

	mpz_t mpzPrime;
	mpz_init_set_ui(mpzPrime, p);

	if (mpz_cmp_ui(knmodp, 0) == 0 || mpz_legendre(knmodp, mpzPrime) == 1) {
		if (mpz_cmp_ui(knmodp, 0) == 0) {
			mpz_clear(mpzPrime);
			return (1.0 / (double)(p - 1));
		}
		else {
			mpz_clear(mpzPrime);
			return (2.0 / (double)(p - 1));
		}
	}

	mpz_clear(mpzPrime);
	return 0.0;
}

// Knuth-Schroeppel Multiplier algorithm
// Algorithm was implemented as it's described in
// The Multiple Polynomial Quadratic Sieve by Robert D. Silverman
// But some details are not explained there, so the rest of algorithm
// is inspired by msieve.
uint32_t KnuthSchroeppelMultiplier(data_t &data, string primeFileName)
{
	const uint32_t p_cnt = 256;
	const uint8_t mult_list[] =
	{ 1,  2,  3,  5,  7, 11, 13, 17,
	 19, 23, 29, 31, 37, 41, 43, 47,
	 53, 59, 61, 67, 71, 73 };

	const uint32_t k_cnt = (sizeof(mult_list) / sizeof(uint8_t));

	uint32_t multiplier;
	double   best_res;

	vector<uint32_t> multPrimes;
	GetMultPrimes(multPrimes, primeFileName);

	// Scores
	double f_res[k_cnt];

	mpz_t nmodp;
	mpz_t kn;
	mpz_t knmodp;

	mpz_init(nmodp);
	mpz_init(kn);
	mpz_init(knmodp);

	for (uint32_t i = 0; i < k_cnt; ++i) {
		uint32_t k = mult_list[i];
		mpz_mul_ui(kn, data.N, k);
		double log_k = log(k) / 2;

		f_res[i] = -log_k;

		for (uint32_t l = 0; l < p_cnt; ++l) {
			uint32_t p = multPrimes.at(l);
			double log_p = log(p);

			mpz_mod_ui(nmodp, data.N, p);
			mpz_mod_ui(knmodp, kn, p);

			f_res[i] += g(p, k, data, nmodp, kn, knmodp)*log_p;
		}
	}

	multiplier = 1;
	best_res = 0.0;

	for (uint32_t i = 0; i < k_cnt; ++i) {
		if (f_res[i] > best_res) {
			best_res = f_res[i];
			multiplier = mult_list[i];
		}
	}

	// Cleaning
	multPrimes.clear();
	mpz_clear(nmodp);
	mpz_clear(kn);
	mpz_clear(knmodp);

	return multiplier;
}

/** Self-Initialization Quadratic Sieve
 *
 * @param primesFileName File where deltas of primes are stored
 * @param N              Number to factor
 * @param res            Factor of N
 * @return               Result is OK or not
 */
bool SIQS(string primesFileName, data_t &data, uint16_t digits, mpz_t res)
{ 
	uint32_t bits;
	uint32_t x = 0;

	mpz_t mpzM; // M at gmp form
	mpz_t idealDivisor;

	pol_t Q_ab; // Polynom

	bits = mpz_sizeinbase(data.N, 2);
	data.digits = digits;

	// Find optimal values
	while (prebuilt_params[x].bits < bits) ++x;

	if (digits < 30) {
		data.F = prebuilt_params[2].fb_size;
		data.M = prebuilt_params[2].sieve_size;
		data.largePrimeMultiplier = prebuilt_params[2].large_mult;

		data.divNum = prebuilt_params[2].div_num;
		data.k = 0;
	}
	else {
		data.F = prebuilt_params[x].fb_size;
		data.M = prebuilt_params[x].sieve_size;
		data.largePrimeMultiplier = prebuilt_params[x].large_mult;

		if (prebuilt_params[x].div_num) {
			data.divNum = prebuilt_params[x].div_num;
			data.k = data.divNum - 3;
		}
	}

	mpz_init(data.idealA);
	mpz_init_set(data.origN, data.N);
	mpz_init_set_ui(mpzM, data.M);
	mpz_init(idealDivisor);

	data.primesMPZ = (mpz_t *)malloc(sizeof(mpz_t)* data.F);
	data.quadraticResidue = (mpz_t *)malloc(sizeof(mpz_t)* data.F);

	// Get multiplayer to boost factorization
	data.multiplier = KnuthSchroeppelMultiplier(data, primesFileName);
	cout << "Multiplier k set to: " << data.multiplier << endl;
	mpz_mul_ui(data.N, data.N, data.multiplier);

	// Build factor base and get large primes for Single Large Prime Variation
	if (!GetPrimeDeltas(primesFileName, data)) {
		mpz_set(res, data.N);
		return OK;
	}

	// Init array for prime logarithms
	data.logs = (double *)malloc(sizeof(double) * (data.primes.at(data.primes.size() - 1) + 1));
	data.threshold = (double *)malloc(sizeof(double)* data.M);

	// Compute ideal a coeficient
	cout << "Ideal a is set to: ";
	mpz_mul_ui(data.idealA, data.N, 2);
	mpz_sqrt(data.idealA, data.idealA);
	mpz_tdiv_q(data.idealA, data.idealA, mpzM);
	mpz_out_str(stdout, 10, data.idealA);
	cout << endl;

	// Compute ideal divisor of a
	cout << "Ideal divisor of Ideal a: ";
	mpz_root(idealDivisor, data.idealA, data.divNum);
	mpz_out_str(stdout, 10, idealDivisor);
	cout << endl;

	if (data.divNum == 2) {
		// Set bounds
		data.leftBound = 25;
		data.rightBound = 100;

		// Get primes that will make Binary Search Tree that will contain
		// all tuples of these primes where the key is sum of logarithms
		// of primes that make a tuple
		for (uint32_t i = data.leftBound; i <= data.rightBound; ++i)
			data.treePrimes.push_back(data.primes.at(i));
	}
	else if (data.divNum == 3) {
		// Set bounds
		data.leftBound = 75;
		data.rightBound = 200;

		// Get primes that will make Binary Search Tree that will contain
		// all tuples of these primes where the key is sum of logarithms
		// of primes that make a tuple
		for (uint32_t i = data.leftBound; i <= data.rightBound; ++i)
			data.treePrimes.push_back(data.primes.at(i));
	}
	else {
		uint32_t idealPrime = mpz_get_ui(idealDivisor);
		
		// Find the closest prime to the idealDivisor
		uint32_t i;
		for (i = 0; i < data.primes.size(); ++i) {
			if (data.primes.at(i) > idealPrime)
				break;
		}
		--i;

		cout << "Ideal prime at index: " << i << endl;

		// Set bounds
		data.leftBound = i - 10;
		data.rightBound = i + 200;

		for (uint32_t i = 55; i <= data.rightBound; ++i) {
			// Get primes, that will be used by NEXKSB
			if ((i & 1) && (i >= data.leftBound))
				data.nexksbPrimes.push_back(data.primes.at(i));
			// Get primes that will make Binary Search Tree that will contain
			// all tuples of these primes where the key is sum of logarithms
			// of primes that make a tuple
			else
				data.treePrimes.push_back(data.primes.at(i));
		}
	}

	// Compute logarithm of ideal a
	data.logIdealA = log(mpz_get_d(data.idealA));

	// Compute logs for prime base
	ComputeLogs(data);

	// Compute how many relations we need to collect - 5% oversieve
	InitRelations(data);

	// Buil Binary Search Tree
	if (data.divNum == 2) {
		BSTInit(&(data.log2primes));

		Build2PrimesLogTree(data);
	}
	else {
		BSTInit(&(data.log3primes));

		Build3PrimesLogTree(data);
	}

#ifdef NO_PARALLEL
		InitPolynom(Q_ab, data);

		GeneratePolynom(data, Q_ab);

		Sieve(data, Q_ab);
#else
#pragma omp parallel shared(data) private(Q_ab)
	{
		InitPolynom(Q_ab, data);

		GeneratePolynom(data, Q_ab);

		Sieve(data, Q_ab);
	}
#endif

	// Linear algebra part
	uint64_t **matrixUI = NULL;
	uint32_t rows;
	uint32_t rowBlocks;
	uint32_t cols;
	uint32_t attempt = 0;
	vector<vector<uint32_t>> S;
	vector<relInfoUInt_t> singletons;
	mpz_set_ui(res, 1);

	while (!ProcessRelations(data, rows, cols, singletons)) {
		// Too few relations was gathered => bring singletons back and do more sieving
		cout << "Nasbirano prilis malo relaci :(" << endl;
		data.relations.insert(data.relations.begin(), singletons.begin(), singletons.end());
		data.relationsNeededPart = data.relations.size() + (data.relationsNeededPart / 8);
		cout << "Singletons num: " << singletons.size() << endl;
		singletons.clear();
#ifdef NO_PARALLEL
		InitPolynom(Q_ab, data);

		GeneratePolynom(data, Q_ab);

		Sieve(data, Q_ab);
#else
#pragma omp parallel shared(data) private(Q_ab)
		{
			InitPolynom(Q_ab, data);

			GeneratePolynom(data, Q_ab);

			Sieve(data, Q_ab);
		}
#endif
	}
	cout << "Singletons num: " << singletons.size() << endl;
	singletons.clear();

	// Create matrix for linear algebra
	matrixUI = InitMatrixUI(rows, rowBlocks, cols, data);

	// Fill it with relations
	FillMatrixUI(matrixUI, rows, rowBlocks, cols, data);

	// Do linear algebra
	FastGaussianUIParallel(matrixUI, cols, rows, rowBlocks, S, data);

	ComputeResult(S, data, res, attempt);

	// cleaning
	S.clear();
	free(data.logs);
	for (uint32_t i = 0; i < data.primes.size(); ++i) {
		mpz_clear(data.quadraticResidue[i]);
		mpz_clear(data.primesMPZ[i]);
	}
	free(data.quadraticResidue);
	free(data.primesMPZ);

	for (uint32_t i = 0; i < cols; ++i) {
		free(matrixUI[i]);
	}
	free(matrixUI);

	if (data.divNum == 2)
		BSTDispose(&(data.log2primes));
	else
		BSTDispose(&(data.log3primes));

	mpz_clear(data.idealA);
	mpz_clear(data.origN);
	mpz_clear(mpzM);
	mpz_clear(idealDivisor);

	return OK;
}
