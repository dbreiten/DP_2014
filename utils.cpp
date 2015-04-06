/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace celych cisel z pohledu lamani RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    utils.cpp                                                  **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

#include "utils.h"

using namespace std;

/** Implementation of great common divisor
*/
uint64_t Euklid_algorithm(uint64_t u, uint64_t w)
{
	uint64_t r;

	while (w != 0) {
		r = u % w;
		u = w;
		w = r;
	}

	return u;
}

/** Implementation of great common divisor using big numbers
*/
void Euklid_algorithm_BIG(mpz_t u, mpz_t w, mpz_t res)
{
	mpz_t r;

	mpz_init(r);

	while (mpz_cmp_ui(w, 0) != 0) {
		mpz_mod(r, u, w);
		mpz_set(u, w);
		mpz_set(w, r);
	}

	mpz_set(res, u);

	// cleaning
	mpz_clear(r);
}

/** This function represents polynom that is used by Pollard Rho Method
 *  Most common polynom is g(x) = x^2 + 1
 *
 *  Function is private
 *
 * @param x Input
 * @param N Z_N
 * @return  Result
 */
void Pollard_polynom(mpz_t x, mpz_t N) {

	// x^2
	mpz_mul(x, x, x);

	// +1
	mpz_add_ui(x, x, 1);

	// mod N
	mpz_mod(x,x,N);
}

/** Factorization method: Pollar Rho Method
 * This method is used for smaller numbers.
 *
 * @param N Number to factor
 * @return  Result
 */
void Pollard_rho_method_GMP(mpz_t res, mpz_t N) {
	mpz_t b_i;
	mpz_t b_2i;

#ifdef UTILS_DEBUG_PRINT_POLLARD_RHO
	cout << "Going to print Pollard Rho method debug msgs" << endl;
#endif

	// Initialize the number
	mpz_init(b_i);
	mpz_init(b_2i);
	mpz_set_ui(b_i, 2);
	mpz_set_ui(b_2i, 2);

	do {
		Pollard_polynom(b_i, N);
		Pollard_polynom(b_2i, N);
		Pollard_polynom(b_2i, N);

#ifdef UTILS_DEBUG_PRINT_POLLARD_RHO
		cout << " b_i = ";
		mpz_out_str(stdout, 10, b_i);
		cout << endl;

		cout << " b_2i = ";
		mpz_out_str(stdout, 10, b_2i);
		cout << endl;
#endif

		mpz_sub(res, b_i, b_2i);

#ifdef UTILS_DEBUG_PRINT_POLLARD_RHO
		cout << " b_i - b_2i = ";
		mpz_out_str(stdout, 10, res);
		cout << endl;
#endif

		mpz_abs(res,res);

#ifdef UTILS_DEBUG_PRINT_POLLARD_RHO
		cout << " abs = ";
		mpz_out_str(stdout, 10, res);
		cout << endl;
#endif

		mpz_gcd(res, res, N);

#ifdef UTILS_DEBUG_PRINT_POLLARD_RHO
		cout << " gcd(res, N) = ";
		mpz_out_str(stdout, 10, res);
		cout << endl;
#endif
	} while(mpz_cmp_ui(res,1) <= 0);

	mpz_clear(b_i);
	mpz_clear(b_2i);
}


uint32_t NEXTKSB(vector<uint64_t> &values, vector<uint64_t> &outIndexes, vector<uint64_t> &outValues, uint32_t n, uint32_t k, uint32_t &h, uint32_t &m)
{
#ifdef DEBUG_PRINT
	cout << "NEXTKSB()" << endl;
#endif

	// First entry
	if (outValues.size() == 0) {
		m = 0;
		h = k;

		for (uint32_t j = 1; j <= h; ++j) {
			outIndexes.push_back(m + j);
			outValues.push_back(values.at((m + j) - 1));
		}

		return OK;
	}

	// Later entries
	if (m < n - h)
		h = 0;

	h++;
	m = outIndexes.at((k + 1 - h) - 1);

	for (uint32_t j = 1; j <= h; ++j) {
		outIndexes.at((k + j - h) - 1) = m + j;
		outValues.at((k + j - h) - 1) = values.at((m + j) - 1);
	}

	if (outIndexes.at(0) == n - k + 1)
		return ERR;

	return OK;
}

void Build3PrimesLogTree(data_t &data)
{
	BSTInit(&(data.log3primes));

	uint64_t *content;
	double key;

	for (uint32_t i = 0; i < data.treePrimes.size() - 2; ++i) {
		for (uint32_t j = i + 1; j < data.treePrimes.size() - 1; ++j) {
			for (uint32_t k = j + 1; k < data.treePrimes.size(); ++k) {
				key = data.logs[data.treePrimes.at(i)] + data.logs[data.treePrimes.at(j)] + data.logs[data.treePrimes.at(k)];

				content = (uint64_t *)malloc(sizeof(uint64_t)* 3);
				content[0] = data.treePrimes.at(i);
				content[1] = data.treePrimes.at(j);
				content[2] = data.treePrimes.at(k);

				BSTInsert(&(data.log3primes), key, content);
			}
		}
	}

	cout << "Tree filled" << endl;

	return;
}

void Build2PrimesLogTree(data_t &data)
{
	BSTInit(&(data.log2primes));

	uint64_t *content;
	double key;

	for (uint32_t i = 0; i < data.treePrimes.size() - 2; ++i) {
		for (uint32_t j = i + 1; j < data.treePrimes.size() - 1; ++j) {
				key = data.logs[data.treePrimes.at(i)] + data.logs[data.treePrimes.at(j)];

				content = (uint64_t *)malloc(sizeof(uint64_t)* 2);
				content[0] = data.treePrimes.at(i);
				content[1] = data.treePrimes.at(j);

				BSTInsert(&(data.log2primes), key, content);
		}
	}

	cout << "Tree filled" << endl;

	return;
}

// find x^2 = q mod n
// return
// -1 q is quadratic non-residue mod n
//  1 q is quadratic residue mod n
//  0 q is congruent to 0 mod n
//
// Implemendation made by Ramón T. B. framontb at yahoo.es
int quadratic_residue(mpz_t x, mpz_t q, mpz_t n)
{
	int          leg;
	unsigned int mod4;
	mpz_t        tmp, ofac, nr, t, r, c, b;
	mp_bitcnt_t  twofac = 0, m, i, ix;

	mod4 = mpz_tstbit(n, 0);
	if (!mod4) // must be odd
		return 0;

	mod4 += 2 * mpz_tstbit(n, 1);

	leg = mpz_legendre(q, n);
	if (leg != 1)
		return leg;

	mpz_init_set(tmp, n);

	if (mod4 == 3) // directly, x = q^(n+1)/4 mod n
	{
		mpz_add_ui(tmp, tmp, 1UL);
		mpz_tdiv_q_2exp(tmp, tmp, 2);
		mpz_powm(x, q, tmp, n);
		mpz_clear(tmp);
	}
	else // Tonelli-Shanks
	{
		mpz_inits(ofac, t, r, c, b, NULL);

		// split n - 1 into odd number times power of 2 ofac*2^twofac
		mpz_sub_ui(tmp, tmp, 1UL);
		twofac = mpz_scan1(tmp, twofac); // largest power of 2 divisor
		if (twofac)
			mpz_tdiv_q_2exp(ofac, tmp, twofac); // shift right

		// look for non-residue
		mpz_init_set_ui(nr, 2UL);
		while (mpz_legendre(nr, n) != -1)
			mpz_add_ui(nr, nr, 1UL);

		mpz_powm(c, nr, ofac, n); // c = nr^ofac mod n

		mpz_add_ui(tmp, ofac, 1UL);
		mpz_tdiv_q_2exp(tmp, tmp, 1);
		mpz_powm(r, q, tmp, n); // r = q^(ofac+1)/2 mod n

		mpz_powm(t, q, ofac, n);
		mpz_mod(t, t, n); // t = q^ofac mod n

		if (mpz_cmp_ui(t, 1UL) != 0) // if t = 1 mod n we're done
		{
			m = twofac;
			do
			{
				i = 2;
				ix = 1;
				while (ix<m)
				{
					// find lowest 0 < ix < m | t^2^ix = 1 mod n
					mpz_powm_ui(tmp, t, i, n); // repeatedly square t
					if (mpz_cmp_ui(tmp, 1UL) == 0)
						break;
					i <<= 1; // i = 2, 4, 8, ...
					ix++; // ix is log2 i
				}
				mpz_powm_ui(b, c, 1 << (m - ix - 1), n); // b = c^2^(m-ix-1) mod n
				mpz_mul(r, r, b);
				mpz_mod(r, r, n); // r = r*b mod n
				mpz_mul(c, b, b);
				mpz_mod(c, c, n); // c = b^2 mod n
				mpz_mul(t, t, c);
				mpz_mod(t, t, n); // t = t b^2 mod n
				m = ix;
			} while (mpz_cmp_ui(t, 1UL) != 0); // while t mod n != 1
		}
		mpz_set(x, r);
		mpz_clears(tmp, ofac, nr, t, r, c, b, NULL);
	}

	return 1;
}