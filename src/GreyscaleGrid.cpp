#include <stdlib.h>

class GreyscaleGrid {
  public:
    GreyscaleGrid(R2Image& img);
    double get(int i, int j);
    void set(int i, int j, double val);

  private:
    int width;
    int height;
    int size;
    double* grid;
};

GreyscaleGrid::
GreyscaleGrid(R2Image& img) {
  width = img.Width();
  height = img.Height();
  size = width * height;
  grid = (double*) malloc(sizeof(double)*size);
  assert(grid);

  double* comps;
  for(int i = 0; i < width; i++) {
    for(int j = 0; j < height; j++) {
      double sum = 0;
      comps = img[i][j].Components();
      for(int k = 0; k < 3; k++) {
        sum += comps[i];
      }
      grid[i*height+ j] = sum/3;
    }
  }
}

double GreyscaleGrid::
get(int i, int j) {
  return grid[i*height + j];
}

void GreyscaleGrid::
set(int i, int j, double val) {
  grid[i*height + j] = val;
}
