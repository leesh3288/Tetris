#define DBG /// debug switch
// #define DBG_MBM /// debug move-by-move
// #define SCORELINE /// score = lines*10 (+placed block) switch
#define TSLINE /// True Score Line: score = lines switch
// #define SCORECFF /// uses score coefficients
// #define OPTBEST /// output best entity of each generation to optimal.txt
#define CHECKPT /// checkpoint switch

#define BOARDX 10
#define BOARDY 20

#define idx(x, y) ((y)*(1<<BOARDY)+(x)) // reversed for easier indexing
#define NOMINMAX

#include <cstdio>
#include <cmath>
#include <climits>
#include <thread>
#include <atomic>
#include <random>
#include <chrono>
#include <algorithm>
#include <omp.h>
#include "TetrisGlobal.h"
#include "rlutil.h" // cross-platform getch(), gotoxy(), setColor(), etc.


std::atomic<int> fct;
char* hashdat;

st_func gaSolver(int, int, int, int, double, double);
void gaInitPool(st_func *, int, std::mt19937 &);
void gaFitnessCalc(st_func *, int, int, int, std::mt19937 &);
void gaEvolve(st_func *&, int, int, int, double, double, std::mt19937 &);
long long gaSimulate(st_func *, int, long long, std::mt19937 &);
void gaCrossover(st_func *, int, int, int, int, double, double, std::mt19937 &);
void gaNormalize(st_func *, int);
inline bool isPossible(int[BOARDX], int[BOARDY + 5], st_tetro move);
double calcMoveVal(int[BOARDX], int[BOARDY + 5], int, int, int, st_tetro, st_func &);
long long procMove(int[BOARDX], int[BOARDY + 5], int &, int &, int &, st_tetro);
char* initHash();
inline int hashLineRemove(int, int, int);

int main(void)
{

	printf("initHash() started.\n");
	hashdat = initHash();
	printf("initHash() completed.\n", (int)hashdat);

#ifndef DBG
	st_func optimal = gaSolver(50, 50, 20, 1'000'000'000, 0.1, 1.5);
#endif

#ifdef DBG
#ifdef DBG_MBM
	rlutil::saveDefaultColor();
#endif

	/// DEBUG
	st_func *t = new st_func[1];
	double dat[CFFCT] = { -1.227216147591351, -0.181782005756474, 3.243910450182081, -3.744504339902469, -0.879149441484679, -0.519702872040023, 0.090352749470464, 0.000000000000000, 0.000000000000000, 0.000000000000000, -2.238796977210718, -4.666693888260132 };
	for (int i = 0; i < CFFCT; i++)
		t[0].cff[i] = dat[i];
	// gaNormalize(t, 0);
	std::random_device rd;
	std::mt19937 mt(rd());
	
	FILE *fs = fopen("testing.txt", "at");
	fprintf(fs, "Parameters: (");
	for(int i=0; i<CFFCT - 1; i++)
		fprintf(fs, "%.15lf, ", t[0].cff[i]);
	fprintf(fs, "%.15lf)\n", t[0].cff[CFFCT - 1]);
	fprintf(fs, "Map Size: (%d, %d)\n", BOARDX, BOARDY);
	fprintf(fs, "Test Mode:");

#ifdef SCORELINE
	fprintf(fs, " SCORELINE");
#endif
#ifdef TSLINE
	fprintf(fs, " TSLINE");
#endif
#ifdef SCORECFF
	fprintf(fs, " SCORECFF");
#endif
	fprintf(fs, "\n\n");
	fclose(fs);
	fs = fopen("testing.csv", "at");
	fprintf(fs, "Test#,Score,AvgScore\n");
	fclose(fs);

	long long int totline = 0;
	int counts = 0;
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < 100; i++)
	{
		long long int curline;
		curline = gaSimulate(t, 0, LLONG_MAX, mt);

		#pragma omp critical
		{
			totline += curline;
			counts++;
			printf("Test %d / %d: Average Score %lld\n", counts, 100, totline / counts);
			FILE *fp = fopen("testing.txt", "at");
			fprintf(fp, "Test %d / %d: Average Score %lld\n", counts, 100, totline / counts);
			fclose(fp);
			fp = fopen("testing.csv", "at");
			fprintf(fp, "%d,%lld,%lld\n", counts, curline, totline / counts);
			fclose(fp);
		}
	}
	///

	// gaSimulate(t, 0, 1'000'000'000, mt);
	
	delete[] t;
	///
#endif

	return 0;
}

// population | evolve # | simulation # | max mvmt in simulation | mutation chance | max mutation deviation
st_func gaSolver(int genpop, int evolvect, int simct, int simmove, double mutch, double mutdev)
{
	st_func *entity = new st_func[genpop], opt;

	std::random_device rd;
	std::mt19937 mt(rd());

	gaInitPool(entity, genpop, mt);

	for (int i = 0; i < evolvect; i++)
	{
		auto timebase = std::chrono::high_resolution_clock::now();
		
		printf("==============================================================================\n");
		printf("Generation #%d fitness calculation started\n", i);
		printf("==============================================================================\n");
		
		gaFitnessCalc(entity, genpop, simct, simmove, mt);
		
		long long Midx = 0, avg = 0;
		for (int j = 1; j < genpop; j++)
		{
			// printf("Entity #%d: fitness = %lld\n", i, entity[i].fitness); /// DEBUG
			if (entity[j].fitness > entity[Midx].fitness)
				Midx = j;
			avg += entity[j].fitness;
		}

		FILE *fp = fopen("log.txt", "at");
		fprintf(fp, "%d,%lld,%lld\n", i, entity[Midx].fitness, avg / genpop);
		fclose(fp);

		printf("==============================================================================\n");
		printf("Generation #%d fitness calculation finished with: \n", i);
		printf(" >> optimal fitness = %lld, average fitness = %lld\n", entity[Midx].fitness, avg / genpop);
		printf(" >> optimal: %.2lf", entity[Midx].cff[0]);
		for (int j = 1; j < CFFCT; j++)
			printf(", %.2lf", entity[Midx].cff[j]);
		printf("\n");
		printf("==============================================================================\n");

#ifdef OPTBEST
		fp = fopen("optimal.txt", "at");
		fprintf(fp, "(opt) %lld: ", entity[Midx].fitness);
		fprintf(fp, "%.2lf", entity[Midx].cff[0]);
		for (int j = 1; j < CFFCT; j++)
			fprintf(fp, ", %.2lf", entity[Midx].cff[j]);
		fprintf(fp, "\n");
		fclose(fp);
#endif
		
		printf("==============================================================================\n");
		printf("Generation #%d -> #%d evolution started\n", i, i + 1);
		printf("==============================================================================\n");

		gaEvolve(entity, genpop, simct, simmove, mutch, mutdev, mt);
		
		auto timelen = std::chrono::high_resolution_clock::now() - timebase;
		printf("==============================================================================\n");
		printf("Generation #%d -> #%d completed (Total %lld.%09llds)\n", i, i + 1, timelen.count() / 1'000'000'000LL, timelen.count() % 1'000'000'000LL); /// DEBUG
		printf("==============================================================================\n");

		// rlutil::getkey(); /// DEBUG
	}
	gaFitnessCalc(entity, genpop, simct, simmove, mt);

	long long mfit = entity[0].fitness;
	int midx = 0;

	for (int i = 1; i < genpop; i++)
	{
		if (entity[i].fitness > mfit)
		{
			mfit = entity[i].fitness;
			midx = i;
		}
	}
	opt = entity[midx];

	delete[] entity;

	return opt;
}

void gaInitPool(st_func *entity, int genpop, std::mt19937 &mt)
{
	std::uniform_real_distribution<double> rnd(-50.0, 50.0);

	for (int i = 0; i < genpop; i++)
	{
		for (int j = 0; j < CFFCT; j++)
		{
			entity[i].cff[j] = rnd(mt);
#ifndef SCORECFF
			if (j == 7 || j == 8 || j == 9)
				entity[i].cff[j] = 0.0;
#endif
		}
		gaNormalize(entity, i);
	}

	return;
}

void gaFitnessCalc(st_func *entity, int genpop, int simct, int simmove, std::mt19937 &mt)
{
	for (int i = 0; i < genpop; i++)
	{
		if (entity[i].fitness != 0) // fitness already calculated
			continue;

		fct = 0;
		long long fit = 0;
		#pragma omp parallel for schedule(dynamic) reduction(+:fit)
		for (int j = 0; j < simct; j++)
			fit += gaSimulate(entity, i, simmove, mt); // random engine works fine multithreaded
		entity[i].fitness = fit / simct;

		if (fct >= simct * 4 / 5)
		{
			FILE *fp = fopen("optimal.txt", "at");
			fprintf(fp, "%d / %d (%lld): %.2lf", fct.load(), simct, entity[i].fitness, entity[i].cff[0]);
			for (int k = 1; k < CFFCT; k++)
				fprintf(fp, ", %.2lf", entity[i].cff[k]);
			fprintf(fp, "\n");
			fclose(fp);
			printf("Entity succeeded %d times, wrote results to optimal.txt\n", fct.load());
		}

		printf("Entity #%d simulation complete: fitness = %lld\n", i, entity[i].fitness); /// DEBUG
	}

	return;
}

void gaEvolve(st_func *&entity, int genpop, int simct, int simmove, double mutch, double mutdev, std::mt19937 &mt)
{
	/*
	1. randomly shuffle
	2. take two, choose the fitter, add to next gen (total ~50% entity filled)
	3. from the parents, randomly choose 2, produce offspring (fill leftover 50%~ with offsprings)
	*/
	std::shuffle(entity, entity + genpop, mt);

	st_func *nextgen = new st_func[genpop];

	int par = 0;
	for (int i = 1; i < genpop; i += 2)
	{
		if (entity[i - 1].fitness > entity[i].fitness)
			nextgen[par++] = entity[i - 1];
		else
			nextgen[par++] = entity[i];
	}
	delete[] entity;

	// nextgen[0 ~ par-1]: parents
	std::uniform_int_distribution<int> r1(0, par - 2), r2(0, par - 1);
	int p1, p2;
	for (int i = par; i < genpop; i++)
	{
		/*
		select 2 unique random numbers from [0, par - 1] in uniform distribution
		1. randomly select 1 number from [0, par - 2]
		2. randomly select 1 number from [0, par - 1]. if overlaps, set it to par - 1.
		*/
		p1 = r1(mt);
		p2 = r2(mt);
		if (p1 == p2)
			p2 = par - 1;

		gaCrossover(nextgen, i, p1, p2, 4, mutch, mutdev, mt);
	}

	entity = nextgen;

	return;
}

long long gaSimulate(st_func *entity, int cursim, long long simmove, std::mt19937 &mt)
{
	// single simulation (with evaluation function entity[cursim] & maximum move # simmove)

	int cboard[BOARDX] = { 0 }, rboard[BOARDY + 5] = { 0 };
#ifdef DBG
	int bcboard[BOARDX] = { 0 }, brboard[BOARDY + 5] = { 0 }; /// DEBUG
#endif
	st_tetro tester;
	long long score = 0, dsc;
	int b2bct = 0, b2btp = -1, comboct = 0;
	int tetromix[7] = { 0, 1, 2, 3, 4, 5, 6 };

	for (long long i = 0; i < simmove; i++)
	{
#ifdef CHECKPT
		if (i % 100'000 == 0)
			printf("Checkpoint: %lld, Current Score: %lld\n", i, score);
#endif
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

				mvval = calcMoveVal(cboard, rboard, b2bct, b2btp, comboct, tester, entity[cursim]);
				if (mvval > optval)
				{
					optval = mvval;
					optmv = tester;
				}
			}
		}

		// place down block, update score (terminate if game end)
		if ((dsc = procMove(cboard, rboard, b2bct, b2btp, comboct, optmv)) == points[(int)type_score::gameover])
		{
#ifdef DBG
			printf("Total Score: %lld\n", score);
			printf("Total Moves: %lld\n", i);

			/// TESTING
			FILE *fp = fopen("testing.txt", "at");
			fprintf(fp, "Total Score: %lld\n", score);
			fprintf(fp, "Total Moves: %lld\n", i);
			fclose(fp);
			///
#endif
			return score + dsc;
		}
		score += dsc;

#ifdef DBG_MBM
		/// DEBUG
		rlutil::cls();

		for (int i = BOARDY - 1; i >= 0; i--)
		{
			for (int j = 0; j < BOARDX; j++)
			{
				if (brboard[i] & (1 << j))
					printf("бс");
				else
				{
					bool bP = false;
					for (int k = 0; k < 4; k++)
						if (j == optmv.crd.x + tetsl[(int)optmv.type][optmv.rot][k].x - 2 && i == optmv.crd.y + 2 - tetsl[(int)optmv.type][optmv.rot][k].y)
							bP = true;
					if (bP)
					{
						rlutil::setColor(12);
						printf("бс");
						rlutil::resetColor();
					}
					else
						printf("бр");
				}
			}
			printf("\n");
		}
		printf("\n");

		for (int i = BOARDY - 1; i >= 0; i--)
		{
			for (int j = 0; j < BOARDX; j++)
			{
				if (rboard[i] & (1 << j))
					printf("бс");
				else
					printf("бр");
			}
			brboard[i] = rboard[i];
			printf("\n");
		}
		printf("Score: %lld\n", score);
		printf("Moves: %lld\n", i + 1);

		// std::this_thread::sleep_for(std::chrono::milliseconds(50));
		rlutil::getkey();
		///
#endif
	}

	// printf("Entity #%d successfully finishes all %d moves\n", cursim, simmove); /// DEBUG
	fct++;

#ifdef DBG
	printf("Total Score: %lld\n", score);
	printf("Total Moves: %lld\n", simmove);

	/// TESTING
	FILE *fp = fopen("testing.txt", "at");
	fprintf(fp, "Total Score: %lld\n", score);
	fprintf(fp, "Total Moves: %lld\n", simmove);
	fclose(fp);
	///
#endif

	return score /*+ points[(int)type_score::finish]*/;
}

void gaCrossover(st_func *entity, int child, int parent1, int parent2, int cotype, double mutch, double mutdev, std::mt19937 &mt)
{
	switch (cotype)
	{
		case 1: // Average Crossover
		{
			for (int i = 0; i < CFFCT; i++)
				entity[child].cff[i] = entity[parent1].cff[i] + entity[parent2].cff[i];
		}
			break;
		case 2: // Weighted Average Crossover
		{
			for (int i = 0; i < CFFCT; i++)
				entity[child].cff[i] = entity[parent1].cff[i] * entity[parent1].fitness + entity[parent2].cff[i] * entity[parent2].fitness;
		}
			break;
		case 3: // Equal Probability Random Crossover
		{
			std::bernoulli_distribution bd(0.5);
			for (int i = 0; i < CFFCT; i++)
			{
				if (bd(mt))
					entity[child].cff[i] = entity[parent1].cff[i];
				else
					entity[child].cff[i] = entity[parent2].cff[i];
			}
		}
			break;
		case 4: // Weighted Probability Random Crossover
		{
			std::bernoulli_distribution bd((double)entity[parent1].fitness / (entity[parent1].fitness + entity[parent2].fitness));
			for (int i = 0; i < CFFCT; i++)
			{
				if (bd(mt))
					entity[child].cff[i] = entity[parent1].cff[i];
				else
					entity[child].cff[i] = entity[parent2].cff[i];
			}
		}
			break;
		default:
			break;
	}

	std::bernoulli_distribution mutate(mutch);
	std::uniform_real_distribution<double> delta(-1.0, 1.0);
	for (int i = 0; i < CFFCT; i++)
		if (mutate(mt))
			entity[child].cff[i] *= (1 + delta(mt)*mutdev);
	
	gaNormalize(entity, child);

	return;
}

void gaNormalize(st_func *entity, int target)
{
	double rlen = 0;

	for (int i = 0; i < CFFCT; i++)
		rlen += entity[target].cff[i] * entity[target].cff[i];

	rlen = sqrt(rlen / SQRAD);
	for (int i = 0; i < CFFCT; i++)
		entity[target].cff[i] /= rlen;

	return;
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

double calcMoveVal(int cboard[BOARDX], int rboard[BOARDY + 5], int b2bct, int b2btp, int comboct, st_tetro move, st_func &valfunc)
{
	// only valid moves are given for input

	double ret = 0;
	int lcct = 0, holect = 0, tmp, tx, ty;

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
	for (int i = 0; i < 10; i++)
	{
		ret += valfunc.cff[0] * hashdat[idx(cboard[i], 0)]; // smaxh, 0
		
		ret += valfunc.cff[3] * hashdat[idx(cboard[i], 1)]; // hole, 3
		
		ret += valfunc.cff[11] * hashdat[idx(cboard[i], 2)]; // ctrans, 11

		if (hashdat[idx(cboard[i], 0)] > tmp) // finding maximum height
			tmp = hashdat[idx(cboard[i], 0)];

		if (i > 0)
			ret += valfunc.cff[4] * abs(hashdat[idx(cboard[i - 1], 0)] - hashdat[idx(cboard[i], 0)]); // bumpy, 4
	}

	ret += valfunc.cff[1] * tmp; // maxh, 1

	ret += valfunc.cff[2] * lcct; // compl, 2

	for (int i = 0; i < tmp; i++)
		ret += valfunc.cff[10] * hashdat[idx(rboard[i] + (1 << BOARDX), 2)]; // rtrans, 10
	ret += valfunc.cff[10] * (BOARDY - tmp) * 2; // rtrans compensation for uncounted rows

	tmp = points[(int)type_score::place];
	if (lcct > 0)
	{
		tmp += points[(int)type_score(lcct - 1)];

		ret += valfunc.cff[7] * comboct; // rescombo, 7
		tmp += points[(int)type_score::combo] * comboct;
		if (b2btp == lcct)
		{
			ret += valfunc.cff[8] * b2bct * b2btp; // resb2b, 8 (takes b2b type into account)
			tmp += points[(int)type_score(b2btp - 1)] * b2bct / 2;
		}
	}

	ret += valfunc.cff[9] * tmp; // dscore, 9

	// remove block
	for (int i = 0; i < 4; i++)
	{
		tx = move.crd.x + tetsl[(int)move.type][move.rot][i].x - 2;
		ty = move.crd.y + 2 - tetsl[(int)move.type][move.rot][i].y;
		cboard[tx] ^= (1 << ty);
		rboard[ty] ^= (1 << tx);
	}

	for (int i = 0; i < 4; i++)
	{
		tx = move.crd.x + tetsl[(int)move.type][move.rot][i].x - 2;
		ty = move.crd.y + 2 - tetsl[(int)move.type][move.rot][i].y;

		if (tx - 1 < 0 || (rboard[ty] & (1 << (tx - 1))))
			ret += valfunc.cff[5]; // adjLR, 5
		if (tx + 1 >= BOARDX || (rboard[ty] & (1 << (tx + 1))))
			ret += valfunc.cff[5]; // adjLR, 5
		if (ty - 1 < 0 || (cboard[tx] & (1 << (ty - 1))))
			ret += valfunc.cff[6]; // adjD, 6
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