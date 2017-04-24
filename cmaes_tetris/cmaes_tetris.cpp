#define SCORELINE /// score = lines*10 (+placed block) switch
// #define TSLINE /// True Score Line: score = lines switch
// #define SCORECFF /// uses score coefficients

#define BOARDX 10
#define BOARDY 15

#define idx(x, y) ((y)*(1<<BOARDY)+(x)) // reversed for easier indexing
#define NOMINMAX

#define genpop 100
#define iterct 1000
#define simct 100
#define simmove LLONG_MAX
// #define restartct 100
#define SIGMA 1.5
#define algotype aCMAES
// #define ELITIST /// uses elitism

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <thread>
#include <vector>
#include <omp.h>
#include "TetrisGlobal.h"
#include "cmaes.h"

using namespace libcmaes;


std::mt19937 mt(0);
char* hashdat;

long long gaSimulate(const double *x);
inline bool isPossible(int[BOARDX], int[BOARDY + 5], st_tetro move);
double calcMoveVal(int[BOARDX], int[BOARDY + 5], int, int, int, st_tetro, const double *x);
long long procMove(int[BOARDX], int[BOARDY + 5], int &, int &, int &, st_tetro);
char* initHash();
inline int hashLineRemove(int, int, int);
FitFunc objFunc = [](const double *x, const int N)
{
	long long fit = 0;

	#pragma omp parallel for schedule(dynamic) reduction(+:fit)
	for (int i = 0; i < simct; i++)
		fit += gaSimulate(x);

	printf("Simulation complete: fitness = %.2lf\n", (double)fit / simct);

	return -((double)fit / simct); // negation for maximizing value
};
ProgressFunc<CMAParameters<GenoPheno<pwqBoundStrategy>>,CMASolutions> progFunc = [](const CMAParameters<GenoPheno<pwqBoundStrategy>> &cmaparams, const CMASolutions &cmasols)
{
	if(cmasols.niter()==0)
		return 0;
	
	CMASolutions tsols = cmasols; // temporary solution set to find best candidate
	tsols.sort_candidates();
	Candidate bcand = tsols.best_candidate();
	const double *x = cmaparams.get_gp().pheno(bcand.get_x_dvec()).data();

	printf("================================================================================\n");
	
 	FILE *fp = fopen("progress.txt", "at");

	fprintf(fp, "Generation %d\n", cmasols.niter());
	printf("Generation %d\n", cmasols.niter());

	fprintf(fp, "    Run Status: %d, Elapsed Time: %.2lf, ObjFunc Evals: %d\n", cmasols.run_status(), cmasols.elapsed_last_iter() / 1000.0, cmasols.fevals());
	printf("    Run Status: %d, Elapsed Time: %.2lf, ObjFunc Evals: %d\n", cmasols.run_status(), cmasols.elapsed_last_iter() / 1000.0, cmasols.fevals());
	
	fprintf(fp, "    Best Candidate: %.15lf", x[0]);
	printf("    Best Candidate: %.15lf", x[0]);
	for (int i = 1; i < CFFCT; i++)
	{
		fprintf(fp, ", %.15lf", x[i]);
		printf(", %.15lf", x[i]);
	}
	fprintf(fp, " (%.2lf)\n", -bcand.get_fvalue());
	printf(" (%.2lf)\n", -bcand.get_fvalue());

	fclose(fp);
	
	printf("================================================================================\n");

	return 0;
};

int main(void)
{
	printf("initHash() started.\n");
	hashdat = initHash();
	printf("initHash() completed. char* = %d\n", (int)hashdat);

	mt.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());

	int dim = CFFCT;
	double sigma = SIGMA;
	double lbounds[CFFCT], ubounds[CFFCT];

	for(int i = 0; i < dim; i++)
	{
		lbounds[i] = -5.0;
		ubounds[i] = 5.0;
	}
	std::vector<double> x0(dim, 0.0);
	GenoPheno<pwqBoundStrategy> gp(lbounds, ubounds, dim);

	CMAParameters<GenoPheno<pwqBoundStrategy>> cmaparams(x0, sigma, genpop, std::chrono::high_resolution_clock::now().time_since_epoch().count(), gp);
	cmaparams.set_fplot("log.dat"); // log data for visualization (uses Python & Matplotlib)
	cmaparams.set_max_iter(iterct);
#ifdef ELITIST
	cmaparams.set_elitism(1); // re-injects the best-ever seen solution
#endif
#ifndef SCORECFF
	cmaparams.set_fixed_p(7, 0.0);
	cmaparams.set_fixed_p(8, 0.0);
	cmaparams.set_fixed_p(9, 0.0);
#endif
	cmaparams.set_algo(algotype);
#ifdef restartct
	cmaparams.set_restarts(restartct);
#endif
	cmaparams.set_mt_feval(false);
	
	CMASolutions cmasols = cmaes<GenoPheno<pwqBoundStrategy>>(objFunc, cmaparams, progFunc);

	FILE *fp = fopen("progress.txt", "at");
	fprintf(fp, "\n");
	fprintf(fp, "Total Time Consumption: %lf seconds\n", cmasols.elapsed_time() / 1000.0);
	fprintf(fp, "Return Code: %d\n", cmasols.run_status());
	fclose(fp);
	
	printf("Total Time Consumption: %lf seconds\n", cmasols.elapsed_time() / 1000.0);
	printf("Return Code: %d\n", cmasols.run_status());
	
	scanf("%*s"); // dummy input

	return cmasols.run_status();
}

long long gaSimulate(const double *x)
{
	int cboard[BOARDX] = { 0 }, rboard[BOARDY + 5] = { 0 };
	st_tetro tester;
	long long score = 0, dsc;
	int b2bct = 0, b2btp = -1, comboct = 0;
	int tetromix[7] = { 0, 1, 2, 3, 4, 5, 6 };

	for (long long i = 0; i < simmove; i++)
	{
		double optval = -987654321.0, mvval;
		st_tetro optmv;

		if (i % 7 == 0)
			std::shuffle(tetromix, tetromix + 7, mt);

		tester.type = type_tetro(tetromix[i % 7]); // randomly generate tetromino

		for (int k = 0; k <= 3; k++) // rotation state
		{
			tester.rot = k;
			for (int l = 0; l < BOARDX; l++)
			{
				tester.crd.x = l;
				tester.crd.y = BOARDY + 2;

				if (!isPossible(cboard, rboard, tester)) // default fails (block goes over x-limit)
					continue;

				tester.crd.y = 0;
				for (int m = 0; m < 4; m++)
					tester.crd.y = std::max(tester.crd.y, hashdat[idx(cboard[tester.crd.x + tetsl[(int)tester.type][tester.rot][m].x - 2], 0)] + tetsl[(int)tester.type][tester.rot][m].y - 2);

				mvval = calcMoveVal(cboard, rboard, b2bct, b2btp, comboct, tester, x);
				if (mvval > optval)
				{
					optval = mvval;
					optmv = tester;
				}
			}
		}

		// place down block, update score (terminate if game end)
		if ((dsc = procMove(cboard, rboard, b2bct, b2btp, comboct, optmv)) == points[(int)type_score::gameover])
			return score + dsc;
		score += dsc;
	}

	return score /*+ points[(int)type_score::finish]*/;
}

inline bool isPossible(int cboard[BOARDX], int rboard[BOARDY + 5], st_tetro move)
{
	int tx, ty;
	for (int i = 0; i < 4; i++)
	{
		tx = move.crd.x + tetsl[(int)move.type][move.rot][i].x - 2;
		ty = move.crd.y + 2 - tetsl[(int)move.type][move.rot][i].y;

		if (tx < 0 || tx >= BOARDX || ty < 0 || ty >= BOARDY + 5) // block goes over limit
			return false;

		if ((cboard[tx] & (1 << ty)) || (rboard[ty] & (1 << tx))) // block overlaps w/ onboard blocks
			return false;
	}

	return true;
}

double calcMoveVal(int cboard[BOARDX], int rboard[BOARDY + 5], int b2bct, int b2btp, int comboct, st_tetro move, const double *x)
{
	// only valid moves are given for input

	double ret = 0;
	int lcct = 0, holect = 0, tmp, tx, ty;

	// place block
	for (int i = 0; i < 4; i++)
	{
		tx = move.crd.x + tetsl[(int)move.type][move.rot][i].x - 2;
		ty = move.crd.y + 2 - tetsl[(int)move.type][move.rot][i].y;
		cboard[tx] ^= (1 << ty);
		rboard[ty] ^= (1 << tx);

		if (rboard[ty] == (1 << BOARDX) - 1)
			lcct++;
	}

	for (int i = BOARDY; i < BOARDY + 5; i++)
	{
		if (rboard[i] != 0)
		{
			for (int j = 0; j < 4; j++)
			{
				tx = move.crd.x + tetsl[(int)move.type][move.rot][j].x - 2;
				ty = move.crd.y + 2 - tetsl[(int)move.type][move.rot][j].y;
				cboard[tx] ^= (1 << ty);
				rboard[ty] ^= (1 << tx);
			}

			return points[(int)type_score::gameover_penalty];
		}
	}

	tmp = 0;
	for (int i = 0; i < BOARDX; i++)
	{
		ret += x[0] * hashdat[idx(cboard[i], 0)]; // smaxh, 0
		
		ret += x[3] * hashdat[idx(cboard[i], 1)]; // hole, 3
		
		ret += x[11] * hashdat[idx(cboard[i], 2)]; // ctrans, 11

		if (hashdat[idx(cboard[i], 0)] > tmp) // finding maximum height
			tmp = hashdat[idx(cboard[i], 0)];

		if (i > 0)
			ret += x[4] * abs(hashdat[idx(cboard[i - 1], 0)] - hashdat[idx(cboard[i], 0)]); // bumpy, 4
	}

	ret += x[1] * tmp; // maxh, 1

	ret += x[2] * lcct; // compl, 2

	for (int i = 0; i < tmp; i++)
		ret += x[10] * hashdat[idx(rboard[i] + (1 << BOARDX), 2)]; // rtrans, 10
	ret += x[10] * (BOARDY - tmp) * 2; // rtrans compensation for uncounted rows, for compatibility w/ previous definition of rtrans

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

	// remove block
	for (int i = 0; i < 4; i++)
	{
		tx = move.crd.x + tetsl[(int)move.type][move.rot][i].x - 2;
		ty = move.crd.y + 2 - tetsl[(int)move.type][move.rot][i].y;
		cboard[tx] ^= (1 << ty);
		rboard[ty] ^= (1 << tx);
	}

	ret += x[9] * tmp; // dscore, 9

	for (int i = 0; i < 4; i++)
	{
		tx = move.crd.x + tetsl[(int)move.type][move.rot][i].x - 2;
		ty = move.crd.y + 2 - tetsl[(int)move.type][move.rot][i].y;

		if (tx - 1 < 0 || (rboard[ty] & (1 << (tx - 1))))
			ret += x[5]; // adjLR, 5
		if (tx + 1 >= BOARDX || (rboard[ty] & (1 << (tx + 1))))
			ret += x[5]; // adjLR, 5
		if (ty - 1 < 0 || (cboard[tx] & (1 << (ty - 1))))
			ret += x[6]; // adjD, 6
	}

	return ret;
}

long long procMove(int cboard[BOARDX], int rboard[BOARDY + 5], int &b2bct, int &b2btp, int &comboct, st_tetro move)
{
	int tx, ty, score = 0, lclist[4], lcct = 0, maxh = 0;

	// place block
	for (int i = 0; i < 4; i++)
	{
		tx = move.crd.x + tetsl[(int)move.type][move.rot][i].x - 2;
		ty = move.crd.y + 2 - tetsl[(int)move.type][move.rot][i].y;
		cboard[tx] ^= (1 << ty);
		rboard[ty] ^= (1 << tx);

		if (rboard[ty] == (1 << BOARDX) - 1)
			lclist[lcct++] = ty;
	}

	// gameover checker
	for (int i = BOARDY; i < BOARDY + 5; i++)
		if (rboard[i] != 0)
			return points[(int)type_score::gameover];

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
#elif defined(TSLINE)
	score = lcct;
#endif

	// update board
	if (lcct > 0)
	{
		int curl = 0, curlc = 0, rc = 0, maxh = 0;
		bool tmp[4] = { 0 };

		std::sort(lclist, lclist + lcct);
		curl = lclist[0];

		for (int i = 0; i < lcct; i++)
			tmp[lclist[i] - lclist[0]] = 1;

		for (int i = 0; i < 4; i++)
			rc += (1 << i)*tmp[i];
		rc = (rc - 1) / 2;

		for (int i = 0; i < BOARDX; i++)
		{
			maxh = std::max(maxh, (int)hashdat[idx(cboard[i], 0)]);
			cboard[i] = hashLineRemove(cboard[i], lclist[0] + 1, rc);
		}

		for (int i = lclist[0]; i < maxh; i++)
		{
			if (curlc < lcct && i == lclist[curlc])
			{
				curlc++;
				continue;
			}
			else
				rboard[curl++] = rboard[i];
		}
		for (int i = curl; i < maxh; i++)
			rboard[i] = 0;
	}

	return score;
}

char* initHash()
{
	char *dat = new char[(1 << BOARDY) * 3];

	int log2_cmp = 1, log2_val = 0;
	dat[idx(0, 0)] = 0, dat[idx(0, 1)] = 0, dat[idx(0, 2)] = 1;
	for (int i = 1; i < (1 << BOARDY); i++)
	{
		// height, 0
		dat[idx(i, 0)] = (log2_cmp == i ? (log2_cmp <<= 1, ++log2_val) : log2_val);

		int hct = 0, ctr = 0;
		for (int j = 0; j < dat[idx(i, 0)]; j++)
		{
			if (!(i&(1 << j))) // if empty
			{
				hct++;
				ctr += ((int)(j - 1 < 0 || (i&(1 << (j - 1)))) + (int)(j + 1 <= BOARDY - 1 && (i&(1 << (j + 1)))));
			}
		}

		// hole, 1
		dat[idx(i, 1)] = hct;

		// ctrans, 2 (use this for rtrans as rtrans[n] = ctrans[n+(1<<BOARDX)] for row hash n)
		dat[idx(i, 2)] = ctr;
	}

	return dat;
}

inline int hashLineRemove(int hv, int lo, int rc) // lo = n when lowest is n-th digit
{
	switch (rc) // bit-twiddling
	{
		case 0:
			return (hv&((1 << (lo - 1)) - 1)) + ((hv&((1 << BOARDY) - (1 << lo))) >> 1);
		case 1:
			return (hv&((1 << (lo - 1)) - 1)) + ((hv&((1 << BOARDY) - (1 << (lo + 1)))) >> 2);
		case 2:
			return (hv&((1 << (lo - 1)) - 1)) + ((hv&(1 << lo)) >> 1) + ((hv&((1 << BOARDY) - (1 << (lo + 2)))) >> 2);
		case 3:
			return (hv&((1 << (lo - 1)) - 1)) + ((hv&((1 << BOARDY) - (1 << (lo + 2)))) >> 3);
		case 4:
			return (hv&((1 << (lo - 1)) - 1)) + ((hv&((1 << lo) + (1 << (lo + 1)))) >> 1) + ((hv&((1 << BOARDY) - (1 << (lo + 3)))) >> 2);
		case 5:
			return (hv&((1 << (lo - 1)) - 1)) + ((hv&(1 << (lo + 1))) >> 2) + ((hv&((1 << BOARDY) - (1 << (lo + 3)))) >> 3);
		case 6:
			return (hv&((1 << (lo - 1)) - 1)) + ((hv&(1 << lo)) >> 1) + ((hv&((1 << BOARDY) - (1 << (lo + 3)))) >> 3);
		case 7:
			return (hv&((1 << (lo - 1)) - 1)) + ((hv&((1 << BOARDY) - (1 << (lo + 3)))) >> 4);
		default:
			return -1; // error case
	}
}