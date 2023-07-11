#include <iostream>
#include "SampleDecoder.h"
#include "MTRand.h"
#include "BRKGA.h"
#include <chrono>

#include <fstream>
#include <string>

#define MAX(x,y) ((x) > (y) ? (x) : (y))

using namespace std;

int main(int argc, char* argv[])
{
	//const unsigned n = 100;		// size of chromosomes
	const unsigned p = 1000;	// size of population
	const double pe = 0.3;		// fraction of population to be the elite-set
	const double pm = 0.15;		// fraction of population to be replaced by mutants
	const double rhoe = 0.75;	// probability that offspring inherit an allele from elite parent
	const unsigned K = 3;		// number of independent populations
	const unsigned MAXT = 2;	// number of threads for parallel decoding
	//-------------------------------------------------------------------------
	vector<int> linhas, paresOD, adj, nos, nosId;
	vector<vector<int>> nosPares;
	vector<double> distancia, fo, fluxo;
	int nLinhas, nParesOD, l, par, tMax, sd = 0;
	double distMax = 0., distRef = 0., fMax = 12, fMin = 1;
	string linha;
	double somatorioF, somatorioFo, tempEspRef = 0., tempEspMax = 0., a1 = 1.0, a2 = 1.;

	chrono::steady_clock sc;

	if(argc == 1) // Teste de argumento de entrada
		cout << "Ausência de endereço do arquivo!" << endl;
	if(argc > 2)
	{
		tMax = stoi(string(argv[2]));
		if(argc > 3)
		{
			sd = stoi(string(argv[3]));// seed

			if(argc == 6)
			{
				fMin = stoi(string(argv[4]));
				fMax = stoi(string(argv[5]));
			}
		}
		if(argc == 8)
		{
			a1 = stod(string(argv[6]));
			a2 = stoi(string(argv[8]));
		}
	}
	else
		tMax = 30;
	
	// Abertura de arquivo
	ifstream arquivo(argv[1]);

	if (arquivo.is_open())
	{
		while (getline(arquivo, linha, '\n'))
		{
			if (linha.compare("%LINHAS") == 0)
			{
				//Leitura do Número de Linhas
				arquivo >> nLinhas;

				// Realocação dos vetores [nLinhas]
				linhas.resize(nLinhas);
				distancia.resize(nLinhas);
				fo.resize(nLinhas);

				//Leitura das linhas
				for (l = 0; l < nLinhas; l++) 
					arquivo >> linhas[l];

				//Leitura das distancias
				for (l = 0; l < nLinhas; l++)
					arquivo >> distancia[l];	

				//Leitura das fo
				for (l = 0; l < nLinhas; l++)
				{
					arquivo >> fo[l];
					distRef += fo[l] * distancia[l];
					distMax += fMax * distancia[l];
				}
			}
			else if (linha.compare("%PAR_OD") == 0)
			{
				//Leitura do Número de Rotas
				arquivo >> nParesOD;

				// Realocação dos vetores [nParesOD]
				paresOD.resize(nParesOD * 2);
				fluxo.resize(nParesOD);
				adj.resize(nParesOD * nLinhas);

				for (par = 0; par < nParesOD; par++)
				{
					//Leitura do par OD
					arquivo >> paresOD[par * 2 + 0];
					arquivo >> paresOD[par * 2 + 1];

					//Leitura do fluxo
					arquivo >> fluxo[par];

					//Leitura da matriz de adjacência
					for (l = 0; l < nLinhas; l++)
						arquivo >> adj[par * nLinhas + l];
				}
			}
			else if (linha.compare("%END") == 0)
			{
				break;
			}
			else if (linha.compare("") == 0)
			{
				continue;
			}
			else
			{
				cout << "Erro no formato do arquivo!\n";
				exit(1);
			}
		}
		arquivo.close();
	}
	else
	{
		cout << "Não foi possível ler o arquivo!" << endl;
		exit(1);
	}
	//Size of chromosomes
	const unsigned n = nLinhas;		// size of chromosomes

	for (par = 0; par < nParesOD; par++)
	{
		somatorioF = 0.;
		somatorioFo = 0.;

		for (l = 0; l < nLinhas; l++)
		{
			// Sum of frequencies from current OD-pair
			somatorioF += adj[par*nLinhas + l]*fMin;// Min frequency sum (using fMin)
			somatorioFo += adj[par*nLinhas + l]*fo[l];
		}

		// Max Waiting time possible
		tempEspMax += fluxo[par]/somatorioF; // Waiting time is half of the headway
		tempEspRef += fluxo[par]/somatorioFo;
	}
	tempEspMax *= 30.; // Waiting time is half of the headway
	tempEspRef *= 30.;

	nos.push_back(paresOD[0 * 2 + 0]); // First node of the pair (O)
	nos.push_back(paresOD[0 * 2 + 1]); // First node of the pair (D)
	for(par = 1; par < nParesOD; par++)
	{
		int firstNode = 0, secondNode = 0;
		for(int no = 0; no < nos.size(); no++)
		{
			if(paresOD[par * 2 + 0] == nos[no])
				firstNode++;
			if(paresOD[par * 2 + 1] == nos[no])
				secondNode++;
		}
		if(firstNode == 0)
			nos.push_back(paresOD[par * 2 + 0]);
		if(secondNode == 0)
			nos.push_back(paresOD[par * 2 + 1]);
	}

	for(int no = 0; no < nos.size(); no++)
	{
		for(par = 0; par < nParesOD; par++)
		{
			if(nos[no] == paresOD[par * 2 + 0])
			{
				nosId.push_back(par);
				break;
			}
		}
	}
	for(int no = 0; no < nos.size() - 1; no++)
	{
		int idO = nosId[no];
		par = idO;
		vector<int> aux;

		while(paresOD[par * 2 + 0] == nos[no])
		{
			aux.push_back(paresOD[par * 2 + 1]);
			par++;
		}
		nosPares.push_back(aux);
	}
	//-------------------------------------------------------------------------------------
	vector<double> constantes(8);
	constantes[0] = tempEspRef;
	constantes[1] = tempEspMax;
	constantes[2] = fMin;
	constantes[3] = fMax;
	constantes[4] = distMax;
	constantes[5] = 45 * fMax;
	constantes[6] = a1;
	constantes[7] = a2;

	SampleDecoder decoder(paresOD, nos, nosId, nosPares, adj,nLinhas, distancia, fluxo, constantes);	// initialize the decoder
	
	const long unsigned rngSeed = sd;	// seed to the random number generator
	MTRand rng(rngSeed);				// initialize the random number generator
	
	// initialize the BRKGA-based heuristic
	BRKGA< SampleDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng);
	
	unsigned generation = 0;		// current generation
	const unsigned X_INTVL = 100;	// exchange best individuals at every 100 generations
	const unsigned X_NUMBER = 2;	// exchange top 2 best
	const unsigned MAX_GENS = 1000;	// run for 1000 gens
	
	double t = 0.;
	auto t1 = sc.now();
	do {
		algorithm.evolve();	// evolve the population for one generation
		
		if((++generation) % X_INTVL == 0)
		{
			algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
		}

		//cout << 100.*t/tMax << "%" << endl;

		auto t2 = sc.now();
		auto t_span = static_cast<chrono::duration<double>>(t2 - t1);
		t = t_span.count();
	} while (t < tMax);

	//----------------------------------------------------------------------------------
	// PRINTS
	cout <<  t << "\t" <<  generation << "\t" << algorithm.getBestFitness() << endl;
/*
	double d = 0;
	for(l = 0; l < nLinhas; l++)
	{
		double fL = algorithm.getBestChromosome()[l] * (fMax - fMin) + fMin;

		cout << linhas[l] << ":\t" << fL << endl;

		d += distancia[l] * fL;
	}

	cout << "\nTempo Espera Max:\t" << tempEspMax << endl;
	cout << "Tempo Espera Ref:\t" << tempEspRef << endl;
//	cout << "Cap Ref:\t" << constantes[5] << endl;
	cout << "Distância máxima:\t"  << distMax << endl;
	cout << "Distância Ref:\t" << distRef << "\nDistância calc:\t" << d << endl;
*/
	return 0;
}
