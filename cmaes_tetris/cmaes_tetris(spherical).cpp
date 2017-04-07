#define SCORELINE /// score = lines*10 (+placed block) switch
// #define TSLINE /// True Score Line: score = lines switch
// #define SCORECFF /// uses score coefficients

#define genpop 50
#define iterct 100
#define simct 16
#define simmove 1000000000

#define _USE_MATH_DEFINES

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <thread>
#include <random>
#include <chrono>
#include <omp.h>
#include <vector>
#include <iostream>
#include "TetrisGlobal.h"
#include "cmaes.h"

using namespace libcmaes;


std::mt19937 mt(0);

long long gaSimulate(const double *x);
bool isPossible(bool[27][10], st_tetro);
double calcMoveVal(bool[27][10], int, int, int, st_tetro, const double *x);
long long procMove(bool[27][10], int &, int &, int &, st_tetro);
double* toCartesian(const double *, const int);
FitFunc objFunc = [](const double *x, const int N)
{
	long long fit = 0;

	for (int i = 0; i < N; i++)
		printf("%.3lf ", x[i]);
	printf("\n");

	double *sph = toCartesian(x, N);

	#pragma omp parallel for schedule(dynamic) reduction(+:fit)
	for (int i = 0; i < simct; i++)
		fit += gaSimulate(sph);

	delete[] sph;

	printf("Simulation complete: fitness = %.2lf\n", (double)fit / simct);

	return -((double)fit / simct); // negation for maximizing value
};
ProgressFunc<CMAParameters<GenoPheno<pwqBoundStrategy>>,CMASolutions> progFunc = [](const CMAParameters<GenoPheno<pwqBoundStrategy>> &cmaparams, const CMASolutions &cmasols)
{
	if(cmasols.niter()==0)
		return 0;	
	
	Candidate bcand = cmasols.best_candidate();
	const double *x = cmaparams.get_gp().pheno(bcand.get_x_dvec()).data();

	printf("================================================================================\n");
	
 	FILE *fp = fopen("progress.txt", "at");

	fprintf(fp, "Generation %d\n", cmasols.niter());
	printf("Generation %d\n", cmasols.niter());

	fprintf(fp, "    Run Status: %d, Elapsed Time: %.2lf, ObjFunc Evals: %d\n", cmasols.run_status(), cmasols.elapsed_last_iter() / 1000.0, cmasols.fevals());
	printf("    Run Status: %d, Elapsed Time: %.2lf, ObjFunc Evals: %d\n", cmasols.run_status(), cmasols.elapsed_last_iter() / 1000.0, cmasols.fevals());

	double *sph = toCartesian(x, CFFCT - 1);
	
	fprintf(fp, "    Best Candidate: %.3lf", sph[0]);
	printf("    Best Candidate: %.3lf", sph[0]);
	for (int i = 1; i < CFFCT; i++)
	{
		fprintf(fp, ", %.3lf", sph[i]);
		printf(", %.3lf", sph[i]);
	}
	
	delete[] sph;

	fprintf(fp, " (%.2lf)\n", -bcand.get_fvalue());
	printf(" (%.2lf)\n", -bcand.get_fvalue());

	fclose(fp);
	
	printf("================================================================================\n");

	return 0;
};

int main(void)
{
	mt.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());

	int dim = CFFCT - 1; // sphere in 12-dimension, needs 11 angles
	double sigma = 1.0;
	double lbounds[CFFCT - 1], ubounds[CFFCT - 1];

	for(int i = 0; i < dim - 1; i++)
	{
		lbounds[i] = 0;
		ubounds[i] = M_PI;
	}
	lbounds[dim - 1] = 0;
	ubounds[dim - 1] = M_PI * 2;
	std::vector<double> x0(dim, M_PI_2);
	GenoPheno<pwqBoundStrategy> gp(lbounds, ubounds, dim);

	CMAParameters<GenoPheno<pwqBoundStrategy>> cmaparams(x0, sigma, genpop, std::chrono::high_resolution_clock::now().time_since_epoch().count(), gp);
	cmaparams.set_fplot("log.dat"); // log data for visualization (uses Python & Matplotlib)
	cmaparams.set_max_iter(iterct);
	cmaparams.set_elitism(1); // re-injects the best-ever seen solution
#ifndef SCORECFF
	cmaparams.set_fixed_p(7, 0.0);
	cmaparams.set_fixed_p(8, 0.0);
	cmaparams.set_fixed_p(9, 0.0);
#endif
	cmaparams.set_algo(aCMAES);
	cmaparams.set_mt_feval(false); // multithreading already implemented inside fitness function, turning this multithread off

	CMASolutions cmasols = cmaes<GenoPheno<pwqBoundStrategy>>(objFunc, cmaparams, progFunc);

	printf("Optimal Solution: ");
	cmasols.print(std::cout, 0, gp);
	printf("\n");
	printf("Expected Distance to Minimum: %.2lf\n", cmasols.edm());
	printf("Total Time Consumption: %lf seconds\n", cmasols.elapsed_time() / 1000.0);
	printf("Return Code: %d\n", cmasols.run_status());
	
	scanf("%*s"); // dummy input

	return cmasols.run_status();
}

long long gaSimulate(const double *x)
{
	bool board[27][10] = { 0 }; // board set to h=27, [0~21]: playfield, [22~26]: end zone
	st_tetro tester;
	long long score = 0, dsc;
	int b2bct = 0, b2btp = -1, comboct = 0;
	int tetromix[7] = { 0, 1, 2, 3, 4, 5, 6 };

	for (int i = 0; i < simmove; i++)
	{
		double optval = -987654321.0, mvval;
		st_tetro optmv;

		if (i % 7 == 0)
			std::shuffle(tetromix, tetromix + 7, mt);

		tester.type = type_tetro(tetromix[i % 7]); // randomly generate tetromino

		for (int k = 0; k <= 3; k++) // rotation state
		{
			tester.rot = k;
			for (int l = 0; l < 10; l++)
			{
				tester.crd.x = l;
				tester.crd.y = 24;

				if (!isPossible(board, tester)) // default fails (block goes over x-limit)
					continue;

				int lsc;
				for (lsc = 23; lsc >= 0; lsc--)
				{
					tester.crd.y = lsc;
					if (!isPossible(board, tester))
						break;
				}

				tester.crd.y = lsc + 1; // y = lsc fails (or is -1), then y = lsc + 1 must be a vaild move
				mvval = calcMoveVal(board, b2bct, b2btp, comboct, tester, x);
				if (mvval > optval)
				{
					optval = mvval;
					optmv = tester;
				}
			}
		}

		// place down block, update score (terminate if game end)
		if ((dsc = procMove(board, b2bct, b2btp, comboct, optmv)) == points[(int)type_score::gameover])
			return score + dsc;
		score += dsc;
	}

	return score /*+ points[(int)type_score::finish]*/;
}

bool isPossible(bool board[27][10], st_tetro move)
{
	int tx, ty;
	for (int i = 0; i < 4; i++)
	{
		tx = move.crd.x + tetsl[(int)move.type][move.rot][i].x - 2;
		ty = move.crd.y + 2 - tetsl[(int)move.type][move.rot][i].y;

		if (tx < 0 || tx >= 10 || ty < 0 || ty >= 27) // block goes over limit
			return false;

		if (board[ty][tx]) // block overlaps w/ onboard blocks
			return false;
	}

	return true;
}

double calcMoveVal(bool board[27][10], int b2bct, int b2btp, int comboct, st_tetro move, const double *x)
{
	// only valid moves are given for input

	double ret = 0;
	int listMaxH[10] = { 0 }, lcct = 0, pothole[10] = { 0 }, holect = 0, tmp, tx, ty, rt = 0, ct = 0;

	// place block
	for (int i = 0; i < 4; i++)
	{
		tx = move.crd.x + tetsl[(int)move.type][move.rot][i].x - 2;
		ty = move.crd.y + 2 - tetsl[(int)move.type][move.rot][i].y;
		board[ty][tx] = 1;
	}

	for (int i = 21; i >= 0; i--)
	{
		bool lineIsFull = true;

		for (int j = 0; j < 10; j++)
		{
			if (board[i][j])
			{
				if (listMaxH[j])
				{
					holect += pothole[j];
					pothole[j] = 0;
				}
				else
				{
					listMaxH[j] = i + 1;
					pothole[j] = 0;
				}
			}
			else
			{
				lineIsFull = false;
				pothole[j]++;

				if (j - 1 < 0 || board[i][j - 1])
					rt++;
				if (j + 1 >= 10 || board[i][j + 1])
					rt++;
				if (i - 1 < 0 || board[i - 1][j])
					ct++;
				if (i + 1 <= 21 && board[i + 1][j]) // ignoring topmost board limit
					ct++;
			}
		}

		if (lineIsFull)
			lcct++;
	}
	for (int i = 0; i < 10; i++)
		if (listMaxH[i])
			holect += pothole[i];

	tmp = 0;
	for (int i = 0; i < 10; i++)
	{
		ret += x[0] * listMaxH[i]; // smaxh, 0
		if (listMaxH[i] > tmp)
			tmp = listMaxH[i];
	}
	ret += x[1] * tmp; // maxh, 1

	ret += x[2] * lcct; // compl, 2

	ret += x[3] * holect; // hole, 3

	for (int i = 1; i < 10; i++)
		ret += x[4] * abs(listMaxH[i - 1] - listMaxH[i]); // bumpy, 4

	// remove block
	for (int i = 0; i < 4; i++)
	{
		tx = move.crd.x + tetsl[(int)move.type][move.rot][i].x - 2;
		ty = move.crd.y + 2 - tetsl[(int)move.type][move.rot][i].y;
		board[ty][tx] = 0;
	}

	for (int i = 0; i < 4; i++)
	{
		tx = move.crd.x + tetsl[(int)move.type][move.rot][i].x - 2;
		ty = move.crd.y + 2 - tetsl[(int)move.type][move.rot][i].y;

		if (tx - 1 < 0 || board[ty][tx - 1])
			ret += x[5]; // adjLR, 5
		if (tx + 1 >= 10 || board[ty][tx + 1])
			ret += x[5]; // adjLR, 5
		if (ty - 1 < 0 || board[ty - 1][tx])
			ret += x[6]; // adjD, 6
	}

	tmp = points[(int)type_score::place];
	if (lcct > 0)
	{
		tmp += points[(int)type_score(lcct - 1)];

		ret += x[7] * comboct; // rescombo, 7
		tmp += points[(int)type_score::combo] * comboct;
		if (b2btp == lcct)
		{
			ret += x[8] * b2bct * b2btp; // resb2b, 8 (takes b2b type into account)
			tmp += points[(int)type_score(b2btp - 1)] * b2bct / 2;
		}
	}

	bool breaker = false;
	for (int i = 22; i < 27; i++)
	{
		for (int j = 0; j < 10; j++)
		{
			if (board[i][j])
			{
				tmp = points[(int)type_score::gameover_penalty];
				breaker = true;
				break;
			}
		}

		if (breaker)
			break;
	}

	ret += x[9] * tmp; // dscore, 9

	ret += x[10] * rt; // rtans, 10

	ret += x[11] * ct; // ctrans, 11

	return ret;
}

long long procMove(bool board[27][10], int &b2bct, int &b2btp, int &comboct, st_tetro move)
{
	int tx, ty, score = 0, lclist[4], lcct = 0;

	// place block
	for (int i = 0; i < 4; i++)
	{
		tx = move.crd.x + tetsl[(int)move.type][move.rot][i].x - 2;
		ty = move.crd.y + 2 - tetsl[(int)move.type][move.rot][i].y;
		board[ty][tx] = 1;
	}

	for (int i = 0; i < 22; i++)
	{
		bool lineIsFull = true;

		for (int j = 0; j < 10; j++)
		{
			if (!board[i][j])
			{
				lineIsFull = false;
				break;
			}
		}

		if (lineIsFull)
			lclist[lcct++] = i;
	}

	// score retrieving phase
	score = points[(int)type_score::place];
	if (lcct > 0)
	{
		score += points[(int)type_score(lcct - 1)];

		score += points[(int)type_score::combo] * comboct++;

		if (b2btp == lcct)
			score += points[(int)type_score(b2btp - 1)] * b2bct++ / 2;
		else
		{
			b2btp = lcct;
			b2bct = 1;
		}
	}
	else
		comboct = b2bct = b2btp = 0;

#ifdef SCORELINE
	score = points[(int)type_score::place] + lcct*points[(int)type_score::one];
#endif
#ifdef TSLINE
	score = lcct;
#endif

	// gameover checker
	bool breaker = false;
	for (int i = 22; i < 27; i++)
	{
		for (int j = 0; j < 10; j++)
		{
			if (board[i][j])
			{
				score = points[(int)type_score::gameover];
				breaker = true;
				break;
			}
		}

		if (breaker)
			break;
	}

	// update board
	int curl = 0, curlc = 0;
	for (int i = 0; i < 22; i++)
	{
		if (curlc < lcct && i == lclist[curlc])
		{
			curlc++;
			continue;
		}
		else
		{
			for (int j = 0; j < 10; j++)
				board[curl][j] = board[i][j];
			curl++;
		}
	}
	for (int i = curl; i < 22; i++)
		for (int j = 0; j < 10; j++)
			board[i][j] = 0;

	return score;
}

double *toCartesian(const double *x, const int N)
{
	double *crd = new double[N + 1];
	double base = SPHRAD;

	for (int i = 0; i < N; i++)
	{
#ifndef SCORECFF
		if (i == 7 || i == 8 || i == 9)
		{
			crd[i] = 0;
			continue;
		}
#endif
		crd[i] = base*cos(x[i]);
		base *= sin(x[i]);
	}
	crd[N] = base;

	return crd; // always deallocate memory after use!
}