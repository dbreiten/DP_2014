/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace celych cisel z pohledu lamani RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    matrix.cpp                                                 **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

#include "matrix.h"

// OpenMP
#include <omp.h>

using namespace std;

void InitMatrix(bool ***matrix, uint32_t &rows, uint32_t &cols, data_t &data)
{
	if (*matrix == NULL) {
		// Matrix
		// +2 because +1 for -1 prime and +1 for dependency => we will get matrix A[n,m]
		//rows = data.primes.size() + 2;
		rows = data.relations.size();
		cols = data.primes.size() + 3;
		*matrix = (bool **)malloc(sizeof(bool *)*rows);
	}

	return;
}

void FillMatrixUI(uint64_t **matrix, uint32_t &rows, uint32_t &rowBlocks, uint32_t &cols, data_t &data)
{
	int32_t col;
	uint32_t accRowBlock = 0;
	uint32_t blockSize = sizeof(uint64_t)* 8;
	uint32_t blockIndex = 0;
	uint32_t realRows = blockSize * rowBlocks;
	uint32_t blockIndexRel;
	uint32_t blockOffsetRel;
	uint32_t accBlockRel;

	for (uint32_t row = 0; row < rows; ++row) {
#pragma omp parallel for private(blockIndexRel, blockOffsetRel, accBlockRel)
		for (col = 0; col < cols; ++col) {
			blockIndexRel = col % 32;
			blockOffsetRel = 0x80000000;
			blockOffsetRel >>= blockIndexRel;
			accBlockRel = col >> 5;

			matrix[col][accRowBlock] <<= 1;
			matrix[col][accRowBlock] |= (data.relations.at(row).relation[accBlockRel] & blockOffsetRel) ? 1 : 0;
		}
		blockIndex++;

		if (blockIndex == blockSize) {
			blockIndex = 0;
			accRowBlock++;
		}
	}

	uint32_t tmp = realRows - rows;
	for (uint32_t col = 0; col < cols; ++col) {
		matrix[col][accRowBlock] <<= tmp;
	}

#ifdef LOG_GAUSS_ENABLE
	ofstream matrixFile;

	matrixFile.open("matrixUI.txt");

	for (uint32_t i = 0; i < data.primes.size(); ++i) {
		if (data.primes.at(i) < 10) {
			matrixFile << "   " << data.primes.at(i) << " ";
		}
		else if (data.primes.at(i) < 100) {
			matrixFile << "  " << data.primes.at(i) << " ";
		}
		else if (data.primes.at(i) < 1000) {
			matrixFile << " " << data.primes.at(i) << " ";
		}
		else {
			matrixFile << data.primes.at(i) << " ";
		}
	}
	matrixFile << "  -1 ";
	matrixFile << "M ";
	matrixFile << "D" << endl;

	uint64_t offset = 0x8000000000000000;
	uint32_t rowBlock = 0;
	uint32_t rowBlockIndex = 0;
	bool val = 0;
	for (uint32_t i = 0; i < rows; ++i) {
		for (uint32_t j = 0; j < cols; ++j) {
			val = (matrix[j][rowBlock] & offset) ? true : false;
			if (j == cols - 1) {
				if (val == true)
					matrixFile << "+";
				else
					matrixFile << " ";
			}
			else if (j == cols - 2) {
				if (val == true)
					matrixFile << "* ";
				else
					matrixFile << "  ";
			}
			else
				matrixFile << "   " << val << " ";
		}
		offset >>= 1;
		rowBlockIndex++;
		if (rowBlockIndex >= 64) {
			offset = 0x8000000000000000;
			rowBlockIndex = 0;
			rowBlock++;
		}
		matrixFile << endl;
	}
	matrixFile.close();
#endif

	return;
}

void FillMatrix(uint64_t **matrix, uint32_t &rows, uint32_t &rowBlocks, uint32_t &cols, data_t &data)
{
	uint32_t accRowBlock = 0;
	uint32_t blockSize = sizeof(uint64_t)* 8;
	uint32_t blockIndex = 0;
	uint32_t realRows = blockSize * rowBlocks;

	for (uint32_t row = 0; row < rows; ++row) {
		for (uint32_t col = 0; col < cols; ++col) {
			matrix[col][accRowBlock] <<= 1;
			matrix[col][accRowBlock] |= data.relations.at(row).relation[col];
		}
		blockIndex++;

		if (blockIndex == blockSize) {
			blockIndex = 0;
			accRowBlock++;
		}
	}

	uint32_t tmp = realRows - rows;
	for (uint32_t col = 0; col < cols; ++col) {
		matrix[col][accRowBlock] <<= tmp;
	}

#ifdef LOG_GAUSS_ENABLE
	ofstream matrixFile;

	matrixFile.open("matrixUI.txt");

	for (uint32_t i = 0; i < data.primes.size(); ++i) {
		if (data.primes.at(i) < 10) {
			matrixFile << "   " << data.primes.at(i) << " ";
		}
		else if (data.primes.at(i) < 100) {
			matrixFile << "  " << data.primes.at(i) << " ";
		}
		else if (data.primes.at(i) < 1000) {
			matrixFile << " " << data.primes.at(i) << " ";
		}
		else {
			matrixFile << data.primes.at(i) << " ";
		}
	}
	matrixFile << "  -1 ";
	matrixFile << "M ";
	matrixFile << "D" << endl;

	uint64_t offset = 0x8000000000000000;
	uint32_t rowBlock = 0;
	uint32_t rowBlockIndex = 0;
	bool val = 0;
	for (uint32_t i = 0; i < rows; ++i) {
		for (uint32_t j = 0; j < cols; ++j) {
			val = (matrix[j][rowBlock] & offset) ? true : false;
			if (j == cols - 1) {
				if (val == true)
					matrixFile << "+";
				else
					matrixFile << " ";
			}
			else if (j == cols - 2) {
				if (val == true)
					matrixFile << "* ";
				else
					matrixFile << "  ";
			}
			else
				matrixFile << "   " << val << " ";
		}
		offset >>= 1;
		rowBlockIndex++;
		if (rowBlockIndex >= 64) {
			offset = 0x8000000000000000;
			rowBlockIndex = 0;
			rowBlock++;
		}
		matrixFile << endl;
	}
	matrixFile.close();
#endif

	return;
}

uint64_t **InitMatrixUI(uint32_t &rows, uint32_t &rowBlocks, uint32_t &cols, data_t &data)
{
	uint64_t **matrix;
	rows = data.relations.size();
	rowBlocks = ceil(rows / (double)(sizeof(uint64_t)* 8));
	cols = data.primes.size() + 3;
	matrix = (uint64_t **)malloc(sizeof(uint64_t *)*cols);

	for (uint32_t i = 0; i < cols; ++i) {
		matrix[i] = (uint64_t *)malloc(sizeof(uint64_t )*rowBlocks);
		for (uint32_t j = 0; j < rowBlocks; ++j) {
			matrix[i][j] = 0;
		}
	}

	return matrix;
}

void PrintMatrix(bool **matrix, uint32_t m, uint32_t n)
{
	uint32_t i;
	uint32_t j;

	for (i = 0; i < n; ++i) {
		for (j = 0; j <= m; ++j) {
			if (j == m) {
				if (matrix[i][j] == true)
					cout << "+" << endl;
				else
					cout << " " << endl;
			}
			else if (j == m - 1) {
				if (matrix[i][j] == true)
					cout << "* ";
				else
					cout << "  ";
			}
			else {
				if (matrix[i][j] == true)
					cout << "1 ";
				else
					cout << "0 ";
			}
		}
	}

	cout << "\n\n";

	return;
}

void PrintMatrixToFile(bool **matrix, uint32_t m, uint32_t n, data_t &data)
{
	uint32_t i;
	uint32_t j;
	ofstream matrixFile;

	matrixFile.open("matrix.txt");
	for (i = 0; i < data.primes.size(); ++i) {
		if (data.primes.at(i) < 10) {
			matrixFile << "   " << data.primes.at(i) << " ";
		}
		else if (data.primes.at(i) < 100) {
			matrixFile << "  " << data.primes.at(i) << " ";
		}
		else if (data.primes.at(i) < 1000) {
			matrixFile << " " << data.primes.at(i) << " ";
		}
		else {
			matrixFile << data.primes.at(i) << " ";
		}
	}
	matrixFile << "  -1 ";
	matrixFile << "M ";
	matrixFile << "D" << endl;

	for (i = 0; i < n; ++i) {
		for (j = 0; j <= m; ++j) {
			if (j == m) {
				if (matrix[i][j] == true)
					matrixFile << "+" << endl;
				else
					matrixFile << " " << endl;
			}
			else if (j == m - 1) {
				if (matrix[i][j] == true)
					matrixFile << "* ";
				else
					matrixFile << "  ";
			}
			else {
				if (matrix[i][j] == true)
					matrixFile << "   " << matrix[i][j] << " ";
				else
					matrixFile << "   " << matrix[i][j] << " ";
			}
		}
	}
	matrixFile.close();

	return;
}

void PrintSolvedMatrixToFile(bool **matrix, uint32_t m, uint32_t n, data_t &data)
{
	uint32_t i;
	uint32_t j;
	ofstream matrixFile;

	matrixFile.open("matrixSol.txt");
	for (i = 0; i < data.primes.size(); ++i) {
		if (data.primes.at(i) < 10) {
			matrixFile << "   " << data.primes.at(i) << " ";
		}
		else if (data.primes.at(i) < 100) {
			matrixFile << "  " << data.primes.at(i) << " ";
		}
		else if (data.primes.at(i) < 1000) {
			matrixFile << " " << data.primes.at(i) << " ";
		}
		else {
			matrixFile << data.primes.at(i) << " ";
		}
	}
	matrixFile << "  -1 ";
	matrixFile << "M ";
	matrixFile << "D" << endl;

	for (i = 0; i < n; ++i) {
		for (j = 0; j <= m; ++j) {
			if (j == m) {
				if (matrix[i][j] == true)
					matrixFile << "+" << endl;
				else
					matrixFile << " " << endl;
			}
			else if (j == m - 1) {
				if (matrix[i][j] == true)
					matrixFile << "* ";
				else
					matrixFile << "  ";
			}
			else {
				if (matrix[i][j] == true)
					matrixFile << "   " << matrix[i][j] << " ";
				else
					matrixFile << "   " << matrix[i][j] << " ";
			}
		}
	}
	matrixFile.close();

	return;
}

void DeleteDuplicantsUI(data_t &data)
{
	uint32_t sum;
	int32_t i;
	int32_t j;
	int32_t k;
	int32_t i_size = data.relations.size() - 1;
	int32_t j_size = data.relations.size();
	uint32_t cols = data.expVectorLengthAlloc;
	map<uint32_t, uint32_t> deleteIndexes;

	cout << "Deleting duplicants" << endl;

#pragma omp parallel for private(j,k,sum)
	for (i = 0; i < i_size; ++i) {
		for (j = i + 1; j < j_size; ++j) {
			sum = 0;
			auto data_relations_at_i = data.relations.at(i);
			auto data_relations_at_j = data.relations.at(j);
			for (k = 0; k < cols; ++k) {
				sum = (data_relations_at_i.relation[k] ^ data_relations_at_j.relation[k]);
				if (sum) break;
			}

			if (sum == 0) {
#pragma omp critical
				{
					deleteIndexes.insert(pair<uint32_t, uint32_t>(j, j));
				}
			}
		}
	}

	for (auto rit = deleteIndexes.crbegin(); rit != deleteIndexes.crend(); ++rit) {
		cout << "Deleting duplicant at: " << rit->first << endl;
		data.relations.erase(data.relations.begin() + rit->first);
	}

	return;
}

void DeleteDuplicants(data_t &data)
{
	uint32_t sum;
	int32_t i;
	int32_t j;
	int32_t k;
	int32_t i_size = data.relations.size() - 1;
	int32_t j_size = data.relations.size();
	uint32_t cols = data.primes.size() + 3;
	map<uint32_t, uint32_t> deleteIndexes;

	cout << "Deleting duplicants" << endl;

#pragma omp parallel for private(j,k,sum)
		for (i = 0; i < i_size; ++i) {
			for (j = i + 1; j < j_size; ++j) {
				sum = 0;
				auto data_relations_at_i = data.relations.at(i);
				auto data_relations_at_j = data.relations.at(j);
				for (k = 0; k < cols; ++k)
					sum += (data_relations_at_i.relation[k] ^ data_relations_at_j.relation[k]);

				if (sum == 0) {
#pragma omp critical
					{
						deleteIndexes.insert(pair<uint32_t, uint32_t>(j, j));
					}
				}
			}
		}

	for (auto rit = deleteIndexes.crbegin(); rit != deleteIndexes.crend(); ++rit) {
		cout << "Deleting duplicant at: " << rit->first << endl;
		data.relations.erase(data.relations.begin() + rit->first);
	}

	return;
}

void DeleteNullVectorsUI(data_t &data)
{
	int32_t i;
	int32_t j;
	uint32_t cols = data.expVectorLengthAlloc;
	vector<uint32_t> deleteIndexes;

//#pragma omp parallel for private(j)
	for (i = 0; i < data.relations.size(); ++i) {
		//cout << "OMP: " << omp_get_thread_num() << endl;
		for (j = 0; j < cols; ++j) {
			if (data.relations.at(i).relation[j])
				break;
		}

		if (j >= cols) {
//#pragma omp critical
//			{
				deleteIndexes.push_back(i);
//			}
		}
	}

	// Kdyz neco smazu, tak se samozrejme vektor zmensi a tak se musi prepocitavat i indexy
	// proto -i
	for (uint32_t i = 0; i < deleteIndexes.size(); ++i)
		data.relations.erase(data.relations.begin() + deleteIndexes.at(i) - i);

	return;
}

void DeleteNullVectors(data_t &data)
{
	uint32_t i;
	uint32_t j;
	uint32_t cols = data.primes.size() + 3;
	vector<uint32_t> deleteIndexes;

	for (i = 0; i < data.relations.size(); ++i) {
		for (j = 0; j < cols; ++j) {
			if (data.relations.at(i).relation[j])
				break;
		}

		if (j >= cols)
			deleteIndexes.push_back(i);
	}

	// Kdyz neco smazu, tak se samozrejme vektor zmensi a tak se musi prepocitavat i indexy
	// proto -i
	for (uint32_t i = 0; i < deleteIndexes.size(); ++i)
		data.relations.erase(data.relations.begin() + deleteIndexes.at(i) - i);

	return;
}

void GetRelationsFromPartialsUI(data_t &data)
{
	multimap<uint64_t, relInfoUInt_t>::iterator it = data.partialRelationsUInt.begin();
	multimap<uint64_t, relInfoUInt_t>::iterator it2 = data.partialRelationsUInt.begin();

	uint32_t index1 = 0;
	uint32_t index2 = 0;
	map<uint32_t, multimap<uint64_t, relInfoUInt_t>::iterator> deleteIndexes;

	cout << "Partials before: " << data.partialRelationsUInt.size() << endl;

	// So we can test two different items and because its multimap, it's legit
	it2++; index2++;

	while (it2 != data.partialRelationsUInt.end()) {
		if (it->first != it2->first) {
			// Primes don't match => Set to next pair
			it++; index1++;
			it2++; index2++;
			continue;
		}

		relInfoUInt_t relInfo;

		// Init space for exponent vector
		// +1 for -1 as a prime, + 2 is because flaging by Fast Gaussian method
		uint32_t *relation = (uint32_t *)malloc(sizeof(uint32_t)*data.expVectorLengthAlloc);

		for (uint32_t i = 0; i < data.expVectorLengthAlloc; ++i)
			relation[i] = it->second.relation[i] ^ it2->second.relation[i];

#ifdef LOG_ENABLE_SMOOTH
		ofstream relFile;

		string file2Name = "rel_comb" + to_string(omp_get_thread_num()) + ".txt";

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

		for (uint32_t j = 0; j < data.primes.size(); ++j) {
			relFile << "    ";
			if (relationBool[j])
				relFile << "1 ";
			else
				relFile << "0 ";
		}
		relFile << "    ";
		if (relationBool[data.primes.size()])
			relFile << "1 ";
		else
			relFile << "0 ";

		relFile << endl;
		relFile << "X1: " << it->second.x << "\tX2: " << it2->second.x << endl;
		relFile.close();
#endif

		// Save info about relation
		relInfo.x = it->second.x;
		relInfo.x2 = it2->second.x;
		relInfo.composed = true;
		relInfo.relation = relation;
		mpz_init_set(relInfo.a, it->second.a);
		mpz_init_set(relInfo.a2, it2->second.a);
		mpz_init_set(relInfo.b, it->second.b);
		mpz_init_set(relInfo.b2, it2->second.b);
		mpz_init_set(relInfo.c, it->second.c);
		mpz_init_set(relInfo.c2, it2->second.c);

#ifdef LOG_ENABLE_SMOOTH
		ofstream smoothFile;

		string fileName = "smooth_part_ui.txt";

		smoothFile.open(fileName, ios::out | ios::app);

		smoothFile << data.relationsNum << ":\t";
		smoothFile << "a: " << mpz_get_str(NULL, 10, relInfo.a);
		smoothFile << "\ta2: " << mpz_get_str(NULL, 10, relInfo.a2);
		smoothFile << "\tb: " << mpz_get_str(NULL, 10, relInfo.b);
		smoothFile << "\tb2: " << mpz_get_str(NULL, 10, relInfo.b2);
		smoothFile << "\tc: " << mpz_get_str(NULL, 10, relInfo.c);
		smoothFile << "\tc2: " << mpz_get_str(NULL, 10, relInfo.c2);
		smoothFile << "\tx: " << it->second.x;
		smoothFile << "\tx2: " << it2->second.x;
		smoothFile << "\t Divisors: ";
		blockIndex = 0;
		index = 0;
		tmp = 0;
		block = relation[blockIndex];
		for (uint32_t k = 0; k < 32; ++k) {
			tmp <<= 1;
			tmp |= (block & 1);
			block >>= 1;
		}
		for (uint32_t j = 0; j < data.primes.size(); ++j) {
			if (tmp & 1)
				smoothFile << data.primes.at(j) << " ";

			tmp >>= 1;
			++index;
			if (index == 32) {
				++blockIndex;
				index = 0;
				block = relation[blockIndex];
				tmp = 0;
				for (uint32_t k = 0; k < 32; ++k) {
					tmp <<= 1;
					tmp |= (block & 1);
					block >>= 1;
				}
			}
		}
		if (tmp)
			smoothFile << "-1";
		smoothFile << endl;

		smoothFile.close();
#endif

		data.relations.push_back(relInfo);

		deleteIndexes.insert(pair<uint32_t, multimap<uint64_t, relInfoUInt_t>::iterator>(index1, it));
		deleteIndexes.insert(pair<uint32_t, multimap<uint64_t, relInfoUInt_t>::iterator>(index2, it2));

		// Go behind tested
		it++; it++; index1++; index1++;
		it2++; it2++; index2++; index2++;
	}

	cout << "Partials to delete: " << deleteIndexes.size() << endl;
	for (auto rit = deleteIndexes.crbegin(); rit != deleteIndexes.crend(); ++rit) {
		data.partialRelationsUInt.erase(rit->second);
	}
	
	cout << "Partials after: " << data.partialRelationsUInt.size() << endl;

	return;
}

void DeleteSingletonsParallelUI(data_t &data, vector<relInfoUInt_t> &singletons)
{
	int32_t i;
	int32_t j;
	uint32_t lastIndex = 0;
	uint32_t primeCnt = 0;
	int32_t i_size = data.primes.size();
	int32_t j_size = data.relations.size();
	uint32_t blocksNum = data.expVectorLengthAlloc;
	uint32_t accBlock;
	uint32_t blockOffset;
	uint32_t blockIndex;
	map<uint32_t, uint32_t> deleteIndexes;

	cout << "Relation size with singletons: " << data.relations.size() << endl;

#pragma omp parallel for private(j, lastIndex, primeCnt, blockIndex, blockOffset, accBlock)
	for (i = 0; i < i_size; ++i) {
		blockIndex = i % 32;
		blockOffset = 0x80000000;
		blockOffset >>= blockIndex;
		accBlock = i >> 5;


		lastIndex = 0;
		primeCnt = 0;
		for (j = 0; j < j_size; ++j) {
			if (data.relations.at(j).relation[accBlock] & blockOffset) {
				primeCnt++;
				lastIndex = j;
			}
		}

		if (primeCnt == 1) {
#pragma omp critical
			{
				deleteIndexes.insert(pair<uint32_t, uint32_t>(lastIndex, lastIndex));
				singletons.push_back(data.relations.at(lastIndex));
			}
		}
	}

	for (auto rit = deleteIndexes.crbegin(); rit != deleteIndexes.crend(); ++rit) {
		cout << "Deleting singleton at: " << rit->first << endl;
		data.relations.erase(data.relations.begin() + rit->first);
	}

	cout << "Relation size without singletons: " << data.relations.size() << endl;
	cout << "Relations needed: " << data.primes.size() + 1 << endl;

	return;
}

/*
void DeleteSingletonsParallel(data_t &data, vector<relInfo_t> &singletons)
{
	int32_t i;
	int32_t j;
	uint32_t lastIndex = 0;
	uint32_t primeCnt = 0;
	int32_t i_size = data.primes.size();
	int32_t j_size = data.relations.size();
	map<uint32_t, uint32_t> deleteIndexes;

	cout << "Relation size with singletons: " << data.relations.size() << endl;

#pragma omp parallel for private(j, lastIndex, primeCnt)
	for (i = 0; i < i_size; ++i) {

		lastIndex = 0;
		primeCnt = 0;
		for (j = 0; j < j_size; ++j) {
			if (data.relations.at(j).relation[i]) {
				primeCnt++;
				lastIndex = j;
			}
		}

		if (primeCnt == 1) {
#pragma omp critical
			{
				deleteIndexes.insert(pair<uint32_t, uint32_t>(lastIndex, lastIndex));
			}
		}
	}

	for (auto rit = deleteIndexes.crbegin(); rit != deleteIndexes.crend(); ++rit) {
		cout << "Deleting singleton at: " << rit->first << endl;
		data.relations.erase(data.relations.begin() + rit->first);
	}

	cout << "Relation size without singletons: " << data.relations.size() << endl;
	cout << "Relations needed: " << data.primes.size() + 1 << endl;

	return;
}
*/

void DeleteSingletons(data_t &data)
{
	bool deleted = true;
	uint32_t lastIndex = 0;
	uint32_t primeCnt = 0;

	cout << "Relation size with singletons: " << data.relations.size() << endl;

	while (deleted) {
		deleted = false;

		for (uint32_t i = 0; i < data.primes.size(); ++i) {

			lastIndex = 0;
			primeCnt = 0;
			for (uint32_t j = 0; j < data.relations.size(); ++j) {
				if (data.relations.at(j).relation[i]) {
					primeCnt++;
					lastIndex = j;
				}
			}

			if (primeCnt == 1) {
				cout << "Singleton found at position: " << lastIndex << " Prime: " << data.primes.at(i) << endl;
				data.relations.erase(data.relations.begin() + lastIndex);
				deleted = true;
				break;
			}
		}
		
	}

	cout << "Relation size without singletons: " << data.relations.size() << endl;
	cout << "Relations needed: " << data.primes.size() + 1 << endl;
	//exit(0);
	return;
}

bool ProcessRelations(data_t &data, uint32_t &rows, uint32_t &cols, vector<relInfoUInt_t> &singletons)
{
	uint32_t i = 0;
	uint32_t j = 0;
	uint32_t last = 0;
	vector<uint32_t> emptyCols;
	map<uint32_t, bool> usedIndexes;

	GetRelationsFromPartialsUI(data);

	DeleteNullVectorsUI(data);
	
	DeleteDuplicantsUI(data);
	
	DeleteSingletonsParallelUI(data, singletons);

	// Nelze vytvorit matici, prilis malo relaci
	if (data.relations.size() <= data.primes.size())
		return false;

	return true;
}

void CheckVectors(data_t &data, vector<vector<uint32_t>> &S)
{
	bool *relation = (bool *)malloc(sizeof(bool)*(data.primes.size() + 1));
	for (uint32_t j = 0; j <= data.primes.size(); ++j)
		relation[j] = false;

	for (vector<vector<uint32_t>>::iterator it = S.begin(); it != S.end(); ++it) {
		for (uint32_t j = 0; j <= data.primes.size(); ++j)
			relation[j] = false;

		for (vector<uint32_t>::iterator it2 = (*it).begin(); it2 != (*it).end(); ++it2) {
			for (uint32_t i = 0; i <= data.primes.size(); ++i) {
				relation[i] ^= data.relations.at((*it2)).relation[i];
			}
		}

		uint32_t sum = 0;
		for (uint32_t i = 0; i <= data.primes.size(); ++i) {
			sum += relation[i];
		}

		if (sum == 0) {
			cout << "Vector OK" << endl;
		}
	}

	return;
}

/* Efective matrix is matrix[n,cols-2] because on index -2 is marking flag and on index -1 is dependency flag */
bool FastGaussianUIParallel(uint64_t **matrix, uint32_t cols, uint32_t n, uint32_t rowBlocks, vector<vector<uint32_t>> &S, data_t &data)
{
	uint32_t rowBlock;
	uint64_t rowPivot;
	uint64_t rowOffset;
	uint64_t tmp;
	uint64_t reversed;
	uint32_t pivotIndex;
	int32_t col;
	int32_t k;
	uint32_t x;

	uint32_t m = cols - 2;

#ifdef LOG_ENABLE_MATRIX
	PrintMatrixToFile(matrix, cols - 1, n, data);
#endif

	//PrintMatrix(matrix, cols - 1, n);
#ifdef LOG_GAUSS_ENABLE
	ofstream gaussFile;

	gaussFile.open("gaussUI.txt");
#endif
	for (col = 0; col < m; ++col) {
		// Find A_ij
		for (rowBlock = 0; rowBlock < rowBlocks; ++rowBlock) {
			if (matrix[col][rowBlock] & 0xffffffffffffffff)
				break;
		}

		// No A_ij was found
		if (rowBlock >= rowBlocks) {
			continue;
		}

		// Mark row i
		pivotIndex = (rowBlock + 1)*(sizeof(uint64_t)* 8);
		tmp = matrix[col][rowBlock];
		rowPivot = 1;
		tmp >>= 1;
		pivotIndex--;
		while (tmp) {
			rowPivot <<= 1;
			tmp >>= 1;
			pivotIndex--;
		}
		matrix[m][rowBlock] |= rowPivot;
#ifdef LOG_GAUSS_ENABLE
		gaussFile << "Pivot for column: " << col << " found at row: " << pivotIndex << endl;
#endif

#pragma omp parallel for private(x)
		for (k = 0; k < m; ++k) {
			if (k == col)
				continue;

			// Add column j to column k
			if (matrix[k][rowBlock] & rowPivot) {
#ifdef LOG_GAUSS_ENABLE
				gaussFile << "1 found at column: " << k << endl;
#endif
				for (x = 0; x < rowBlocks; ++x)
					matrix[k][x] ^= matrix[col][x];
			}
		}

		//PrintMatrix(matrix, cols - 1, n);
	}

#ifdef LOG_ENABLE_MATRIX
	PrintSolvedMatrixToFile(matrix, cols - 1, n, data);
#endif

	// Save S
#ifdef LOG_VECTOR_ENABLE
	ofstream vectorFile;

	vectorFile.open("vectors.txt");
#endif
	// Prochazi pres vsechny bloky a snazi se najit neoznaceny radek
	for (uint32_t rowBlock = 0; rowBlock < rowBlocks; ++rowBlock) {
		// Pokud nejaky z radku je neoznaceny, jinak vsechny budou oznacene a tim padem XOR povede k zisku 0, tedy false
		if (matrix[m][rowBlock] ^ 0xffffffffffffffff) {
			// pocatecni index
			uint32_t unmarkedRow = (rowBlock)*(sizeof(uint64_t)* 8);

			// Prohledani konkretniho bloku a nalezeni vsech neoznacenych radku
			tmp = matrix[m][rowBlock];
			reversed = 0;
			for (uint32_t i = 0; i < (sizeof(uint64_t)* 8) - 1; ++i) {
				reversed |= (tmp & 1);
				tmp >>= 1;
				reversed <<= 1;
			}
			reversed |= (tmp & 1);
			tmp = reversed;
			rowOffset = 0x8000000000000000;
			for (uint32_t i = 0; i < (sizeof(uint64_t)* 8); ++i) {
				// Bylo nasbirano dostatek moznosti
				if (S.size() >= 10)
					break;

				// Neoznaceny radek => nutne zpracovat
				if (!(tmp & 1)) {
#ifdef LOG_GAUSS_ENABLE
					gaussFile << "======================================" << endl;
					gaussFile << "NotMarked at: " << unmarkedRow << endl;
#endif
					vector<uint32_t> indexes;
					indexes.push_back(unmarkedRow);

#pragma omp parallel for private(reversed)
					for (col = 0; col < m; ++col) {
						// A[col][row] = 1 => nutne najit zavisle radky
						if (matrix[col][rowBlock] & rowOffset) {
#ifdef LOG_GAUSS_ENABLE
							gaussFile << "Column of dependency at: " << col << endl;
#endif
							// Projiti vsech radku
							bool found = false;
							for (uint32_t k = 0; k < rowBlocks; ++k) {
								// Mozna byl nalezen zavisly radek
								if (matrix[col][k] & 0xffffffffffffffff) {
									uint32_t markedRow = (k)*(sizeof(uint64_t)* 8);
									uint64_t markedRowOffset = 0x8000000000000000;
									uint64_t tmp2 = matrix[col][k];
									reversed = 0;
									for (uint32_t i = 0; i < (sizeof(uint64_t)* 8) - 1; ++i) {
										reversed |= tmp2 & 1;
										tmp2 >>= 1;
										reversed <<= 1;
									}
									reversed |= tmp2 & 1;
									tmp2 = reversed;
									for (uint32_t i = 0; i < (sizeof(uint64_t)* 8); ++i) {
										// Nalezen zavisly radek
										// Na radku se nachazi 1 && oznaceny neni shodny vybranym neoznacenym && je oznaceny
										if ((tmp2 & 1) && (markedRow != unmarkedRow) && (matrix[m][k] & markedRowOffset)) {
#ifdef LOG_GAUSS_ENABLE
											gaussFile << "Marked at: " << markedRow << endl;
#endif
#pragma omp critical
											{
												indexes.push_back(markedRow);
											}
											found = true;
											break;
										}
										markedRow++;
										markedRowOffset >>= 1;
										tmp2 >>= 1;
									}
									if (found)
										break;
								}
							}
						}
					}

					S.push_back(indexes);
					if (S.size() >= 10)
						break;
				}
				unmarkedRow++;
				tmp >>= 1;
				rowOffset >>= 1;
			}
			// Bylo nasbirano dostatek moznosti
			if (S.size() >= 10)
				break;
		}
	}

#ifdef LOG_GAUSS_ENABLE
	gaussFile.close();
#endif

	return true;
}

/* Efective matrix is matrix[n,cols-2] because on index -2 is marking flag and on index -1 is dependency flag */
bool FastGaussianUI(uint64_t **matrix, uint32_t cols, uint32_t n, uint32_t rowBlocks, vector<vector<uint32_t>> &S, data_t &data)
{
	//uint32_t row;
	uint32_t rowBlock;
	uint64_t rowPivot;
	uint64_t rowOffset;
	uint64_t tmp;
	uint64_t reversed;
	uint32_t pivotIndex;
	uint32_t col;
	uint32_t k;
	uint32_t x;

	uint32_t m = cols - 2;

#ifdef LOG_ENABLE_MATRIX
	PrintMatrixToFile(matrix, cols - 1, n, data);
#endif

	//PrintMatrix(matrix, cols - 1, n);
#ifdef LOG_GAUSS_ENABLE
	ofstream gaussFile;

	gaussFile.open("gaussUI.txt");
#endif
	for (col = 0; col < m; ++col) {
		// Find A_ij
		//for (row = 0; row < n; ++row) {
		//	if (matrix[i][j] == true)
		//		break;
		//}

		// Find A_ij
		for (rowBlock = 0; rowBlock < rowBlocks; ++rowBlock) {
			if (matrix[col][rowBlock] & 0xffffffffffffffff)
				break;
		}

		// No A_ij was found
		if (rowBlock >= rowBlocks) {
			continue;
		}

		// Mark row i
		//matrix[i][m] = true;

		// Mark row i
		//pivotIndex = (rowBlock) ? (rowBlock - 1)*(sizeof(uint64_t)* 8) : 0;
		pivotIndex = (rowBlock + 1)*(sizeof(uint64_t)* 8);
		tmp = matrix[col][rowBlock];
		rowPivot = 1;
		tmp >>= 1;
		pivotIndex--;
		while (tmp) {
			rowPivot <<= 1;
			tmp >>= 1;
			pivotIndex--;
		}
		matrix[m][rowBlock] |= rowPivot;
#ifdef LOG_GAUSS_ENABLE
		gaussFile << "Pivot for column: " << col << " found at row: " << pivotIndex << endl;
#endif
		//cout << "Pivot for column: " << col << " found at row: " << pivotIndex << endl;

		for (k = 0; k < m; ++k) {
			if (k == col)
				continue;

			// Add column j to column k
			//if (matrix[i][k] == 1) {
			//	for (x = 0; x < n; ++x)
			//		matrix[x][k] ^= matrix[x][j];
			//}

			// Add column j to column k
			if (matrix[k][rowBlock] & rowPivot) {
#ifdef LOG_GAUSS_ENABLE
				gaussFile << "1 found at column: " << k << endl;
#endif
				//cout << "1 found at column: " << k << endl;
				for (x = 0; x < rowBlocks; ++x)
					matrix[k][x] ^= matrix[col][x];
			}
		}

		//PrintMatrix(matrix, cols - 1, n);
	}

#ifdef LOG_ENABLE_MATRIX
	PrintSolvedMatrixToFile(matrix, cols - 1, n, data);
#endif

	// Save S
#ifdef LOG_VECTOR_ENABLE
	ofstream vectorFile;

	vectorFile.open("vectors.txt");
#endif
	// Prochazi pres vsechny bloky a snazi se najit neoznaceny radek
	for (uint32_t rowBlock = 0; rowBlock < rowBlocks; ++rowBlock) {
		// Pokud nejaky z radku je neoznaceny, jinak vsechny budou oznacene a tim padem XOR povede k zisku 0, tedy false
		if (matrix[m][rowBlock] ^ 0xffffffffffffffff) {
			// pocatecni index
			uint32_t unmarkedRow = (rowBlock)*(sizeof(uint64_t)* 8);

			// Prohledani konkretniho bloku a nalezeni vsech neoznacenych radku
			tmp = matrix[m][rowBlock];
			reversed = 0;
			for (uint32_t i = 0; i < (sizeof(uint64_t)* 8) - 1; ++i) {
				reversed |= (tmp & 1);
				tmp >>= 1;
				reversed <<= 1;
			}
			reversed |= (tmp & 1);
			tmp = reversed;
			rowOffset = 0x8000000000000000;
			for (uint32_t i = 0; i < (sizeof(uint64_t)* 8); ++i) {
				// Bylo nasbirano dostatek moznosti
				if (S.size() >= 10)
					break;

				// Neoznaceny radek => nutne zpracovat
				//if (!(tmp & 1)) {
				if (!(tmp & 1)) {
#ifdef LOG_GAUSS_ENABLE
					gaussFile << "======================================" << endl;
					gaussFile << "NotMarked at: " << unmarkedRow << endl;
#endif
					//cout << "======================================" << endl;
					//cout << "NotMarked at: " << unmarkedRow << endl;
					vector<uint32_t> indexes;
					indexes.push_back(unmarkedRow);

					for (col = 0; col < m; ++col) {
						// A[col][row] = 1 => nutne najit zavisle radky
						if (matrix[col][rowBlock] & rowOffset) {
#ifdef LOG_GAUSS_ENABLE
							gaussFile << "Column of dependency at: " << col << endl;
#endif
							//cout << "Column of dependency at: " << col << endl;
							// Projiti vsech radku
							bool found = false;
							for (uint32_t k = 0; k < rowBlocks; ++k) {
								// Mozna byl nalezen zavisly radek
								if (matrix[col][k] & 0xffffffffffffffff) {
									uint32_t markedRow = (k)*(sizeof(uint64_t)* 8);
									uint64_t markedRowOffset = 0x8000000000000000;
									uint64_t tmp2 = matrix[col][k];
									reversed = 0;
									for (uint32_t i = 0; i < (sizeof(uint64_t)* 8) - 1; ++i) {
										reversed |= tmp2 & 1;
										tmp2 >>= 1;
										reversed <<= 1;
									}
									reversed |= tmp2 & 1;
									tmp2 = reversed;
									for (uint32_t i = 0; i < (sizeof(uint64_t)* 8); ++i) {
										// Nalezen zavisly radek
										// Na radku se nachazi 1 && oznaceny neni shodny vybranym neoznacenym && je oznaceny
										if ((tmp2 & 1) && (markedRow != unmarkedRow) && (matrix[m][k] & markedRowOffset)) {
#ifdef LOG_GAUSS_ENABLE
											gaussFile << "Marked at: " << markedRow << endl;
#endif
											//cout << "Marked at: " << markedRow << endl;
											indexes.push_back(markedRow);
											found = true;
											break;
										}
										markedRow++;
										markedRowOffset >>= 1;
										tmp2 >>= 1;
									}
									if (found)
										break;
								}
							}
						}
					}

					S.push_back(indexes);
					if (S.size() >= 10)
						break;
				}
				unmarkedRow++;
				tmp >>= 1;
				rowOffset >>= 1;
			}
			// Bylo nasbirano dostatek moznosti
			if (S.size() >= 10)
				break;
		}
	}

#ifdef LOG_GAUSS_ENABLE
	gaussFile.close();
#endif

	return true;
}

/* Efective matrix is matrix[n,cols-2] because on index -2 is marking flag and on index -1 is dependency flag */
bool FastGaussian(bool **matrix, uint32_t cols, uint32_t n, vector<vector<uint32_t>> &S, data_t &data)
{
	int32_t i;
	int32_t j;
	int32_t k;
	int32_t x;

	uint32_t m = cols - 2;

#ifdef LOG_GAUSS_ENABLE_NO
	PrintMatrixToFile(matrix, cols - 1, n, data);
#endif

#ifdef LOG_GAUSS_ENABLE
	ofstream gaussFile;

	gaussFile.open("gaussSimple.txt");
#endif

	//PrintMatrix(matrix, cols - 1, n);

	for (j = 0; j < m; ++j) {
		// Find A_ij
		for (i = 0; i < n; ++i) {
			if (matrix[i][j] == true)
				break;
		}

		if (i >= n) {
			//PrintMatrixToFile(matrix, cols - 1, n, data);
			//return false;
			continue;
		}

		// Mark row i
		matrix[i][m] = true;
#ifdef LOG_GAUSS_ENABLE
		gaussFile << "Pivot for column: " << j << " found at row: " << i << endl;
#endif
		cout << "Pivot for column: " << j << " found at row: " << i << endl;


//#pragma omp parallel for private(x)
		for (k = 0; k < m; ++k) {
			if ((k != j) && (matrix[i][k] == 1)) {
#ifdef LOG_GAUSS_ENABLE
				gaussFile << "1 found at column: " << k << endl;
#endif
				cout << "1 found at column: " << k << endl;
				for (x = 0; x < n; ++x) matrix[x][k] ^= matrix[x][j];
			}
		}

		//PrintMatrix(matrix, cols - 1, n);
	}

#ifdef LOG_ENABLE
	PrintSolvedMatrixToFile(matrix, cols - 1, n, data);
#endif

	// Save S
#ifdef LOG_VECTOR_ENABLE
	ofstream vectorFile;

	vectorFile.open("vectors.txt");
#endif
	for (i = 0; i < n; ++i) {
		//cout << "row: " << i << " marked: " << matrix[i][m] << endl;
		// Unmarked => dependency finder
		if (matrix[i][m] == false) {
#ifdef LOG_GAUSS_ENABLE
			gaussFile << "======================================" << endl;
			gaussFile << "NotMarked at: " << i << endl;
#endif
			cout << "======================================" << endl;
			cout << "NotMarked at: " << i << endl;
			matrix[i][cols - 1] = true;
			vector<uint32_t> indexes;
			indexes.push_back(i);

#ifdef LOG_VECTOR_ENABLE
			if (i < 10)
				vectorFile << "  ";
			else if (i < 100)
				vectorFile << " ";
			vectorFile << i << " ";
			for (uint32_t x = 0; x <= m; ++x)
				vectorFile << matrix[i][x] << " ";
			vectorFile << endl;
#endif


			// Finding A_ij = 1
			for (j = 0; j < m; ++j) {
				if (matrix[i][j] == true) {
#ifdef LOG_GAUSS_ENABLE
					gaussFile << "Column of dependency at: " << j << endl;
#endif
					cout << "Column of dependency at: " << j << endl;
					// Projiti vsech radku, pokud je na danem radku v danem sloupci jedna, je zavisly
					for (k = 0; k < n; ++k) {
						if (k == i)
							continue;
						if (matrix[k][j] == true && matrix[k][m] == true) {
#ifdef LOG_GAUSS_ENABLE
							gaussFile << "Marked at: " << k << endl;
#endif
							cout << "Marked at: " << k << endl;
							indexes.push_back(k);
							matrix[k][cols - 1] = true;

#ifdef LOG_VECTOR_ENABLE
							if (k < 10)
								vectorFile << "  ";
							else if (k < 100)
								vectorFile << " ";
							vectorFile << k << " ";
							for (uint32_t x = 0; x <= m; ++x)
								vectorFile << matrix[k][x] << " ";
							vectorFile << endl;
#endif

							//PrintMatrix(matrix, cols - 1, n);
						}
					}
				}
			}
			S.push_back(indexes);
			// We need only like 10 combinations
			if(S.size() >= 10)
				break;
#ifdef LOG_VECTOR_ENABLE
			vectorFile << "\n\n";
#endif
		}
	}
#ifdef LOG_VECTOR_ENABLE
	vectorFile.close();
#endif

#ifdef LOG_GAUSS_ENABLE
	gaussFile.close();
#endif

	return true;
}


