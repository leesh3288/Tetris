#define CFFCT 12
#define SPHRAD 100.0 // radius of (CFFCT)-dimensional sphere in solution space

enum class type_tetro
{
	undef = -1, I, O, J, L, S, T, Z
};

enum class type_score
{
	one, two, three, four, combo, gameover, place, finish, gameover_penalty
	// b2b bonus = (score of performed sequential move)*(performed times(-1))*0.5
	// combo bonus = (value of combo)*(performed times(-1))
};

struct st_crd
{
	int x, y;
};

struct st_tetro
{
	type_tetro type;
	int rot;
	st_crd crd;
};

struct st_grid
{
	bool tetg[5][5];
};

struct st_gbdat
{
	int color;
	bool exists, isVirtual;
};

/*
Coefficients for value function
1. sum of total maximum height (smaxh, 0)
2. maximum height (maxh, 1)
3. completed lines (compl, 2)
4. holes (hole, 3)
5. sum of height differences between adjacent columns (bumpy, 4)
6. adjacent sides (L/R) (adjLR, 5)
7. adjacent sides (D) (adjD, 6)
8. resultant combo counts (rescombo, 7)
9. resultant b2b counts (resb2b, 8)
10. score change (dscore, 9)
11. row transition (rtrans, 10)
12. column transition (ctrans, 11)
*/
struct st_func
{
	double cff[CFFCT];
	long long fitness = 0;

	bool operator<(const st_func &cmp) const
	{
		return fitness < cmp.fitness;
	}
};

extern const st_crd offset[3][4][5];
extern const st_grid tets[7][4];
extern const st_crd tetsl[7][4][4];
extern const int points[9];