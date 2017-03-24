enum class type_tetro
{
	undef = -1, I, O, J, L, S, T, Z
};

enum class type_score
{
	one, two, three, four, combo
	// b2b bonus = (score of performed sequential move)*(performed times)*0.5
	// combo bonus = (value of combo)*(performed times)
};

struct st_crd
{
	int x, y;

	st_crd(int x = 0, int y = 0) : x(x), y(y) {}
};

struct st_tetro
{
	type_tetro type;
	int rot;
	st_crd crd;

	st_tetro(type_tetro type = type_tetro::undef, int rot = 0, st_crd crd = st_crd()) : type(type), rot(rot), crd(crd) {}
};

struct st_grid
{
	bool tetrogrid[5][5];
};

struct st_gbdat
{
	int color;
	bool exists, isVirtual;
};

/*
Coefficients for value function
1. sum of total maximum height (smaxh)
2. maximum height (maxh)
3. completed lines (compl)
4. holes (hole)
5. sum of height differences between adjacent columns (bumpy)
6. adjacent sides (L/R) (adjLR)
7. adjacent sides (D) (adjD)
8. previous combo counts (prevcombo)
9. previous b2b counts (prevb2b)
10. score change (dscore)
*/
struct st_func
{
	double smaxh, maxh, compl, hole, bumpy, adjLR, adjD, prevcombo, prevb2b, dscore;
};

extern const st_crd offset[3][4][5];
extern const st_grid tetroshape[7][4];
extern const int points[5];