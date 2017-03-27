#include <cstdio>
#include <algorithm>
#include <thread>
#include <atomic>
#include <random>
#include <chrono>
#include "TetrisGlobal.h"
#include "rlutil.h" // cross-platform getch(), gotoxy(), setColor(), etc.



st_func gaSolver(int, int, int, int, double, double, double);
void gaInitPool(st_func *, int, double, std::mt19937 &);
void gaFitnessCalc(st_func *, int, int, int, std::mt19937 &);
void gaEvolve(st_func *&, int, int, int, double, double, std::mt19937 &);
double gaSimulate(st_func *, int, int, std::mt19937 &);
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

	st_func optimal = gaSolver(1000, 1000, 500, 5000, 10.0, 50.0, 50.0);

	return 0;
}

// population | evolve # | simulation # | max mvmt in simulation | mutation chance (%) | max mutation deviation (%) | max abs init val
st_func gaSolver(int genpop, int evolvect, int simct, int simmove, double mutch, double mutdev, double absiv)	
{
	st_func *entity = new st_func[genpop], opt;

	std::random_device rd;
	std::mt19937 mt(rd());

	gaInitPool(entity, genpop, absiv, mt);

	for (int i = 0; i < evolvect; i++)
	{
		gaFitnessCalc(entity, genpop, simct, simmove, mt);
		gaEvolve(entity, genpop, simct, simmove, mutch, mutdev, mt);
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

void gaInitPool(st_func *entity, int genpop, double absiv, std::mt19937 &mt)
{
	std::uniform_real_distribution<double> rnd(-absiv, absiv);

	for (int i = 0; i < genpop; i++)
		for (int j = 0; j < 10; j++)
			entity[i].cff[j] = rnd(mt);

	return;
}

void gaFitnessCalc(st_func *entity, int genpop, int simct, int simmove, std::mt19937 &mt)
{
	for (int i = 0; i < genpop; i++)
	{
		if (entity[i].fitness != 0.0) // fitness already calculated
			continue;

		for (int j = 0; j < simct; j++)
			entity[i].fitness += gaSimulate(entity, j, simmove, mt);
		entity[i].fitness /= simct;
	}

	return;
}

void gaEvolve(st_func *&entity, int genpop, int simct, int simmove, double mutch, double mutdev, std::mt19937 &mt)
{
	/*
	1. randomly shuffle
	2. take two, choose the fitter, add to next gen (total 50% entity filled)
	3. from the parents, randomly choose 2, produce offspring (fill leftover 50% with offsprings)
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
	for (int i = par; i < genpop; i++)
	{
		/*
		select 2 unique random numbers from [0, par - 1]
		1. randomly select 1 number from [0, par - 2]
		2. randomly select 1 number from [0, par - 1]. if overlaps, set it to par - 1.
		*/
	}

	entity = nextgen;

	return;
}

double gaSimulate(st_func *entity, int cursim, int simmove, std::mt19937 &mt)
{
	// single simulation (with evaluation function entity[cursim] & maximum move # simmove)
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