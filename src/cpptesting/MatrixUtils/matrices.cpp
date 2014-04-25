#include <vector>
using std::vector;

#define HEIGHT 5
#define WIDTH 3

int main() {
	vector<vector<double> > array2D;

	// Set up sizes. (HEIGHT X WIDTH)
	array2D.resize(HEIGHT);
	for (int i = 0; i < HEIGHT; ++i) {
		array2D[i].resize(WIDTH);
	}

	// Put some values in
	array2D[1][2] = 6.0;
	array2D[3][1] = 5.5;

	return 0;
}
