#define DBG /// debug switch
#define DBG_MBM /// debug move-by-move
// #define SCORELINE /// score = lines*10 (+placed block) switch
#define TSLINE /// True Score Line: score = lines switch
#define SCORECFF /// uses score coefficients
// #define OPTBEST /// output best entity of each generation to optimal.txt

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <thread>
#include <atomic>
#include <random>
#include <chrono>
#include <omp.h>
#include "TetrisGlobal.h"
#include "rlutil.h" // cross-platform getch(), gotoxy(), setColor(), etc.


std::atomic<int> fct;

st_func gaSolver(int, int, int, int, double, double);
void gaInitPool(st_func *, int, std::mt19937 &);
void gaFitnessCalc(st_func *, int, int, int, std::mt19937 &);
void gaEvolve(st_func *&, int, int, int, double, double, std::mt19937 &);
long long gaSimulate(st_func *, int, int, std::mt19937 &);
void gaCrossover(st_func *, int, int, int, int, double, double, std::mt19937 &);
void gaNormalize(st_func *, int);
bool isPossible(bool[27][10], st_tetro);
double calcMoveVal(bool[27][10], int, int, int, st_tetro, st_func &);
long long procMove(bool[27][10], int &, int &, int &, st_tetro);
/*
void tPlayerAIMain();
void tPlayerMain();
void tPlayerActionWatch(const bool &, std::atomic<bool> &, std::atomic<int> &);
void clearScreen();
void drawScreen();
*/

int main(void)
{
	/*
	printf("t to trigger automatic timer\n");
	printf("r to reset timer\n");
	printf("c to reset timer counts\n");
	printf("q to trigger quit command (tPAW -> tPM -> tPAW)\n");

	std::thread tP1_PM(tPlayerMain); // thread parent->child relationship: tP1_PM -> tP1_PAW
	
	tP1_PM.join();
	*/

#ifndef DBG
	st_func optimal = gaSolver(50, 50, 20, 1'000'000'000, 0.1, 1.5);
#endif

#ifdef DBG
	/// DEBUG
	st_func *t = new st_func[1];
	double dat[CFFCT] = { -40.64, -24.96, 37.15, -45.96, -2.24, 44.61, 25.35, 34.15, -20.70, -1.02 };
	for (int i = 0; i < CFFCT; i++)
		t[0].cff[i] = dat[i];
	gaNormalize(t, 0);
	std::random_device rd;
	std::mt19937 mt(rd());
	
	/*
	/// TESTING
	long long int totline = 0;
	for (int i = 0; i < 100; i++)
	{
		totline += gaSimulate(t, 0, 1'000'000'000, mt);
		printf("Test %d / %d: Average Score %lld\n", i + 1, 100, totline / (i + 1));
		FILE *fp = fopen("testing.txt", "at");
		fprintf(fp, "Test %d / %d: Average Score %lld\n", i + 1, 100, totline / (i + 1));
		fclose(fp);
	}
	///
	*/
	gaSimulate(t, 0, 1'000'000'000, mt);
	
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
			fit += gaSimulate(entity, i, simmove, mt); // random engine "seems" to work fine multithreaded
		entity[i].fitness = fit / simct;

		if (fct >= simct * 4 / 5)
		{
			FILE *fp = fopen("optimal.txt", "at");
			fprintf(fp, "%d / %d (%lld): %.2lf", fct, simct, entity[i].fitness, entity[i].cff[0]);
			for (int k = 1; k < CFFCT; k++)
				fprintf(fp, ", %.2lf", entity[i].cff[k]);
			fprintf(fp, "\n");
			fclose(fp);
			printf("Entity succeeded %d times, wrote results to optimal.txt\n", fct);
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

long long gaSimulate(st_func *entity, int cursim, int simmove, std::mt19937 &mt)
{
	// single simulation (with evaluation function entity[cursim] & maximum move # simmove)

	bool board[27][10] = { 0 }; // board set to h=27, [0~21]: playfield, [22~26]: end zone
#ifdef DBG
	bool backboard[27][10] = { 0 }; /// DEBUG
#endif
	st_tetro tester;
	long long score = 0, dsc;
	int b2bct = 0, b2btp = -1, comboct = 0;
	int tetromix[7] = { 0, 1, 2, 3, 4, 5, 6 };

	for (int i = 0; i < simmove; i++)
	{
#ifdef DBG
#ifndef DBG_MBM
		if (i % 100'000 == 0)
			printf("Checkpoint: %d, Current Score: %lld\n", i, score);
#endif
#endif
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
				mvval = calcMoveVal(board, b2bct, b2btp, comboct, tester, entity[cursim]);
				if (mvval > optval)
				{
					optval = mvval;
					optmv = tester;
				}
			}
		}

		// place down block, update score (terminate if game end)
		if ((dsc = procMove(board, b2bct, b2btp, comboct, optmv)) == points[(int)type_score::gameover])
		{
#ifdef DBG
			printf("Total Score: %lld\n", score);
			printf("Total Moves: %d\n", i);

			/// TESTING
			FILE *fp = fopen("testing.txt", "at");
			fprintf(fp, "Total Score: %lld\n", score);
			fprintf(fp, "Total Moves: %d\n", i);
			fclose(fp);
			///
#endif
			return score + dsc;
		}
		score += dsc;

#ifdef DBG_MBM
		/// DEBUG
		rlutil::cls();

		for (int i = 21; i >= 0; i--)
		{
			for (int j = 0; j < 10; j++)
			{
				if (backboard[i][j])
					printf("бс");
				else
					printf("бр");
			}
			printf("\n");
		}
		printf("\n");

		for (int i = 21; i >= 0; i--)
		{
			for (int j = 0; j < 10; j++)
			{
				if (board[i][j])
					printf("бс");
				else
					printf("бр");

				backboard[i][j] = board[i][j];
			}
			printf("\n");
		}
		printf("Score: %lld\n", score);
		printf("Moves: %d\n", i + 1);

		// std::this_thread::sleep_for(std::chrono::milliseconds(50));
		rlutil::getkey();
		///
#endif
	}

	// printf("Entity #%d successfully finishes all %d moves\n", cursim, simmove); /// DEBUG
	fct++;

#ifdef DBG
	printf("Total Score: %lld\n", score);
	printf("Total Moves: %d\n", simmove);

	/// TESTING
	FILE *fp = fopen("testing.txt", "at");
	fprintf(fp, "Total Score: %lld\n", score);
	fprintf(fp, "Total Moves: %d\n", simmove);
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

double calcMoveVal(bool board[27][10], int b2bct, int b2btp, int comboct, st_tetro move, st_func &valfunc)
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
		ret += valfunc.cff[0] * listMaxH[i]; // smaxh, 0
		if (listMaxH[i] > tmp)
			tmp = listMaxH[i];
	}
	ret += valfunc.cff[1] * tmp; // maxh, 1

	ret += valfunc.cff[2] * lcct; // compl, 2

	ret += valfunc.cff[3] * holect; // hole, 3

	for (int i = 1; i < 10; i++)
		ret += valfunc.cff[4] * abs(listMaxH[i - 1] - listMaxH[i]); // bumpy, 4

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
			ret += valfunc.cff[5]; // adjLR, 5
		if (tx + 1 >= 10 || board[ty][tx + 1])
			ret += valfunc.cff[5]; // adjLR, 5
		if (ty - 1 < 0 || board[ty - 1][tx])
			ret += valfunc.cff[6]; // adjD, 6
	}

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

	ret += valfunc.cff[9] * tmp; // dscore, 9

	ret += valfunc.cff[10] * rt; // rtans, 10

	ret += valfunc.cff[11] * ct; // ctrans, 11

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


/*
void tPlayerAIMain()
{
}

void tPlayerMain() // busy-waiting for block drop / action update from tPlayerActionWatch()
{
	bool killPAW = false;
	std::atomic<bool> actflag = false;
	std::atomic<int> actdat;

	bool autotimer = false;
	int timercount = 0;

	std::thread tP1_PAW(tPlayerActionWatch, std::ref(killPAW), std::ref(actflag), std::ref(actdat));

	auto timebase = std::chrono::high_resolution_clock::now();

	while (!killPAW)
	{
		if (autotimer && (std::chrono::high_resolution_clock::now() - timebase).count() > 1000000000LL * 1)
		{
			printf("automatic 1s timer triggered from tPM\n");
			timebase = std::chrono::high_resolution_clock::now();
			printf("timebase reset\n");

			timercount++;
			printf("timercount incremented to %d\n", timercount);
			
			if (timercount == 10)
			{
				printf("automatic timer count reached max limit (10)\n");
				printf("killPAW set to true\n");
				killPAW = true;
			}
		}

		if (actflag) // input received, need processing
		{
			switch (actdat)
			{
				case 'q':
				{
					killPAW = true;
					printf("killPAW set to true\n");
					break;
				}
				case 'r':
				{
					timebase = std::chrono::high_resolution_clock::now();
					printf("timer reset\n");
					break;
				}
				case 't':
				{
					if (autotimer)
					{
						autotimer = false;
						printf("autotimer turned off\n");
					}
					else
					{
						autotimer = true;
						printf("autotimer turned on\n");
						timebase = std::chrono::high_resolution_clock::now();
						printf("timebase reset\n");
					}
					break;
				}
				case 'c':
				{
					timercount = 0;
					printf("timercount reset to 0\n");
					break;
				}
				default:
				{
					printf("command '%c'(%d) not defined, input ignored from tPM\n", (char)actdat, (int)actdat);
					break;
				}
			}

			actflag = false; // time delay can possibly lead to spammed keys being ignored
		}
	}

	printf("waiting for PAW thread to return\n");

	if(tP1_PAW.joinable())
		tP1_PAW.join();

	printf("tPAW joined, tPM now returning\n");

	return;
}

void tPlayerActionWatch(const bool &killsw, std::atomic<bool> &actionflag, std::atomic<int> &actiondat) // busy-waiting for input
{
	while (!killsw)
	{
		actiondat = getch();
		actionflag = true;
	}

	printf("tPAW now returning\n");

	return;
}
*/
/*
void clearScreen()
{
	rlutil::cls();
}

void drawScreen() // needs to be completely mutually exclusive(mutex) - input ignorance/queueing, atomic<bool> for function run check (redundant?)
{
	
}
*/


/*
Controls: 
LArr/RArr for movement
DArr for soft drop
SpaceBar for hard drop
X for right rotation
Z for left rotation
Shift for hold
P for pause
*/

/*
Offset Data / Tetromino Shape viewing code

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 5; k++)
				printf("(%d, %d) ", offset[i][j][k].x, offset[i][j][k].y);
			printf("\n");
		}
		printf("\n");
	}

	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				for (int l = 0; l < 5; l++)
				{
					if (j == 2 && l == 2)
						printf("в├");
					else if (tetroshape[i][k].tetrogrid[j][l])
						printf("бс");
					else
						printf("  ");
				}

				printf("    ");
			}
			printf("\n");
		}
		printf("\n");
	}
*/