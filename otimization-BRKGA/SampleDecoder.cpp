/*
 * SampleDecoder.cpp
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */

#include <iostream>
#include "SampleDecoder.h"
#define MAX(x,y) ((x) > (y) ? (x) : (y))

using namespace std;

SampleDecoder::SampleDecoder() { };

SampleDecoder::SampleDecoder(vector<int> par, vector<int> _nos, vector<int> nosID, vector<vector<int>> _nosPares, vector<int> a, int nL, vector<double> d, vector<double> fl, vector<double> c)
{
	SampleDecoder::nLinhas = nL;
	SampleDecoder::nPares = par.size() / 2;
	SampleDecoder::paresOD = par;
	SampleDecoder::nos = _nos;
	SampleDecoder::nosId = nosID;
	SampleDecoder::nosPares = _nosPares;
	SampleDecoder::adj = a;
	SampleDecoder::distancia = d;
	SampleDecoder::fluxo = fl;
	SampleDecoder::busCapacity = 45;
	SampleDecoder::tempEspRef = c[0];
	SampleDecoder::tempEspMax = c[1];
	SampleDecoder::fMin = c[2];
	SampleDecoder::fMax = c[3];
	SampleDecoder::penality = c[4];
	SampleDecoder::capRef = c[5];
	SampleDecoder::a1 = c[6];
	SampleDecoder::a2 = c[7];
}

SampleDecoder::~SampleDecoder() { }

// Runs in \Theta(n \log n):
double SampleDecoder::decode(const vector< double >& chromosome) const {
	vector< pair< double, unsigned > > ranking(chromosome.size());

	int nNodes = nos.size();

	vector<double> frequencia(nLinhas);
	vector<double> capMax(nLinhas); // Max capacity of a bus line]
	vector<double> escolhaRotas(nPares * nLinhas, 0.);
	vector<double> perfilC(nLinhas, 0.);
	vector<double> pMax(nLinhas);
	double tempoEspera = 0.;
	double somatorioF, somatorioD = 0.;
	double diff = 0.; // Difference of values for penalities
	double pTE = 0.; // penality for unaccepted distance
	double pCap = 0.;  // penality for unaccepted capacity
	int posI = 0;

	// Frequency setting: [0;1] -> [fMin;fMax]
	for(int l = 0; l < nLinhas; l++)
	{
		frequencia[l] = chromosome[l]*(fMax - fMin) + fMin;
		capMax[l] = busCapacity * frequencia[l]; // Max capacity of a line
	}

	// Waiting Time
	for (int par = 0; par < nPares; par++)
	{
		somatorioF = 0.;

		for(int l = 0; l < nLinhas; l++) // Sum of frequencies from current OD-pair
			somatorioF += adj[par * nLinhas + l]*frequencia[l];

		// CONSTRAINT: Capacity of a bus line #1
		for(int l = 0; l < nLinhas; l++)// Construção da matriz Escolha de Rotas
			escolhaRotas[par * nLinhas + l] = fluxo[par] * frequencia[l] * adj[par * nLinhas + l] / somatorioF;

		// CONSTRAINT: Waiting Time
		tempoEspera += fluxo[par] / somatorioF;
	}
	tempoEspera *= 30.; // Waiting time is half of the headway

	if(tempEspRef < tempoEspera) // Waiting Time violation
		pTE += (tempoEspera - tempEspRef) / tempEspMax;


	// OBJECTIVE FUNCTION: Distance (cost) of the company
	for (int l = 0; l < nLinhas; l++)
	{
		somatorioD += frequencia[l] * distancia[l]; // Current 'spent' distance

		for(int j = 0; j < nosPares[0].size(); j++)
		{
			perfilC[l] += escolhaRotas[(0 + j) * nLinhas + l];// First node origin E[0][j][l]
		}
		pMax[l] = perfilC[l];
	}

	// CONSTRAINT: Capacity of a bus line #2
	for(int k = 1; k < nNodes - 1; k++)
	{
		int posK = nosId[k];
		for (int l = 0; l < nLinhas; l++)
		{
			for(int j = 0; j < nosPares[k].size(); j++) // [k][j] pairs
			{
				perfilC[l] += escolhaRotas[(posK + j) * nLinhas + l];// Sobe E[k][j]
			}
			for(int i = 0; i < k; i++) // [i][k] pairs
			{
				posI = nosId[i];

				for(int no = 0; no < nosPares[i].size(); no++)
				{
					if(nos[k] == paresOD[(posI + no) * 2 + 1])
					{
						perfilC[l] -= escolhaRotas[(posI + no)* nLinhas + l];// Desce E[i][k]
					}
				}
			}

			pMax[l] = MAX(pMax[l],perfilC[l]);
		}
	}

	// CONSTRAINT: Capacity of a bus line #3
	for(int l = 0; l < nLinhas; l++)
	{
		if(capMax[l] < pMax[l])
			pCap += (pMax[l] - capMax[l]) / (capRef * nLinhas);
	}
	/*if(pDist<0.00001)
	{
		cout << "pCap:\t" << pCap << "\tpDist:\t" << pDist << endl;


		cout << "frequencias:" << endl;
		for(int l=0;l<nLinhas;l++)
			cout << frequencia[l] << " ";

		cout <<"\npMax" << endl;
		for(int l=0;l<nLinhas;l++)
			cout << pMax[l] << " ";
		cout << endl;
	}*/
	// PENALITY
	diff = (a1*pTE + a2*pCap);// printar pDist e pCap 
//cout << diff << endl;
	if(diff > 0)
		somatorioD += penality * (1 + diff);

/*
cout << "Escolha  de rotas:" << endl;
	for(int par=0;par<nPares;par++)
	{
		for(int l=0;l<nLinhas;l++)
		{
			cout << escolhaRotas[par*nLinhas+l] <<" ";
		}
		cout << "\n";
	}*/

	return somatorioD;
}
