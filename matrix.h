/*********************************************************************************
** Autor:            Dominik Breitenbacher, xbreit00@stud.fit.vutbr.cz          **
** Nazev projektu:   Paralelizace faktorizace celych cisel z pohledu lamani RSA **
** Vedouci projektu: Homoliak Ivan, Ing., UITS FIT VUT                          **
** Nazev souboru:    matrix.h                                                   **
** Posledni uprava:  2014-10-02                                                 **
**                                                                              **
** Popis:                                                                       **
**                                                                              **
*********************************************************************************/

#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>

// Big Numbers library
#if defined (_WIN32) || defined (_WIN64)
#include <mpir.h>
#else
#include <gmpxx.h>
#endif

// OpenMP
#include <omp.h>

#include "siqs.h"

using namespace std;

/******************************************************************************
** STRUCTS                                                                   **
******************************************************************************/
/******************************************************************************
** FUNCTIONS                                                                 **
******************************************************************************/

void InitMatrix(bool ***matrix, uint32_t &rows, uint32_t &cols, data_t &data);

uint64_t **InitMatrixUI(uint32_t &rows, uint32_t &rowBlocks, uint32_t &cols, data_t &data);

void FillMatrixUI(uint64_t **matrix, uint32_t &rows, uint32_t &rowBlocks, uint32_t &cols, data_t &data);

bool ProcessRelations(data_t &data, uint32_t &rows, uint32_t &cols, vector<relInfoUInt_t> &singletons);

bool FastGaussian(bool **matrix, uint32_t m, uint32_t n, vector<vector<uint32_t>> &S, data_t &data);

bool FastGaussianUI(uint64_t **matrix, uint32_t cols, uint32_t n, uint32_t rowBlocks, vector<vector<uint32_t>> &S, data_t &data);

bool FastGaussianUIParallel(uint64_t **matrix, uint32_t cols, uint32_t n, uint32_t rowBlocks, vector<vector<uint32_t>> &S, data_t &data);

#endif /* MATRIX_H_ */
