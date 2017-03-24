#include <cstdio>
#include <thread>
#include <atomic>
#include <chrono>
#include <array>
#include "TetrisGlobal.h"
#include "rlutil.h" // cross-platform getch(), gotoxy(), setColor(), etc.


std::atomic<int> tetGrid[10][22]; // use bitmask to separate data values

void tPlayerAIMain();
void tPlayerMain();
void tPlayerActionWatch(const bool &, std::atomic<bool> &, std::atomic<int> &);
// void clearScreen();
// void drawScreen();

int main(void)
{
	printf("t to trigger automatic timer\n");
	printf("r to reset timer\n");
	printf("c to reset timer counts\n");
	printf("q to trigger quit command (tPAW -> tPM -> tPAW)\n");

	std::thread tP1_PM(tPlayerMain); // thread parent->child relationship: tP1_PM -> tP1_PAW
	
	tP1_PM.join();

	return 0;
}

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