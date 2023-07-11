/*
 * SampleDecoder.h
 *
 * Any decoder must have the format below, i.e., implement the method decode(std::vector< double >&)
 * returning a double corresponding to the fitness of that vector. If parallel decoding is to be
 * used in the BRKGA framework, then the decode() method _must_ be thread-safe; the best way to
 * guarantee this is by adding 'const' to the end of decode() so that the property will be checked
 * at compile time.
 *
 * The chromosome inside the BRKGA framework can be changed if desired. To do so, just use the
 * first signature of decode() which allows for modification. Please use double values in the
 * interval [0,1) when updating, thus obeying the BRKGA guidelines.
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */

#ifndef SAMPLEDECODER_H
#define SAMPLEDECODER_H

#include <list>
#include <vector>
#include <algorithm>

using namespace std;

class SampleDecoder {
public:
	SampleDecoder();
    SampleDecoder(vector<int>, vector<int>, vector<int>, vector<vector<int>>, vector<int>, int, vector<double>, vector<double>, vector<double>);
	~SampleDecoder();

    vector<int> paresOD, nos, nosId, adj;
    vector<vector<int>> nosPares;
    vector<double> distancia, fo, fluxo;
    int nLinhas, nPares, busCapacity;
    double tempEspMax, tempEspRef, capRef, fMax, fMin, penality, a1, a2;

	double decode(const vector< double >& chromosome) const;

private:
};

#endif
