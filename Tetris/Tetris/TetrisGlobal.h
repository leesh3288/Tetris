enum class type_tetro
{
	undef = -1, I, O, J, L, S, T, Z
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

extern const st_crd offset[3][4][5];
extern const st_grid tetroshape[7][4];