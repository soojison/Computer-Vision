// Source file for image class



// Include files

#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include "GreyscaleGrid.cpp"

#include <vector>
#include <algorithm>
#include <queue>
#include <vector>
#include <unordered_map>
#include <limits>
#include <iostream>
#include <random>



////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////


R2Image::
R2Image(void)
  : pixels(NULL),
    npixels(0),
    width(0),
    height(0)
{
}



R2Image::
R2Image(const char *filename)
  : pixels(NULL),
    npixels(0),
    width(0),
    height(0)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height)
  : pixels(NULL),
    npixels(width * height),
    width(width),
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
  : pixels(NULL),
    npixels(width * height),
    width(width),
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width),
    height(image.height)

{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}


void R2Image::
svdTest(void)
{
	// fit a 2D conic to five points
	R2Point p1(1.2,3.5);
	R2Point p2(2.1,2.2);
	R2Point p3(0.2,1.6);
	R2Point p4(0.0,0.5);
	R2Point p5(-0.2,4.2);

	// build the 5x6 matrix of equations
	double** linEquations = dmatrix(1,5,1,6);

	linEquations[1][1] = p1[0]*p1[0];
	linEquations[1][2] = p1[0]*p1[1];
	linEquations[1][3] = p1[1]*p1[1];
	linEquations[1][4] = p1[0];
	linEquations[1][5] = p1[1];
	linEquations[1][6] = 1.0;

	linEquations[2][1] = p2[0]*p2[0];
	linEquations[2][2] = p2[0]*p2[1];
	linEquations[2][3] = p2[1]*p2[1];
	linEquations[2][4] = p2[0];
	linEquations[2][5] = p2[1];
	linEquations[2][6] = 1.0;

	linEquations[3][1] = p3[0]*p3[0];
	linEquations[3][2] = p3[0]*p3[1];
	linEquations[3][3] = p3[1]*p3[1];
	linEquations[3][4] = p3[0];
	linEquations[3][5] = p3[1];
	linEquations[3][6] = 1.0;

	linEquations[4][1] = p4[0]*p4[0];
	linEquations[4][2] = p4[0]*p4[1];
	linEquations[4][3] = p4[1]*p4[1];
	linEquations[4][4] = p4[0];
	linEquations[4][5] = p4[1];
	linEquations[4][6] = 1.0;

	linEquations[5][1] = p5[0]*p5[0];
	linEquations[5][2] = p5[0]*p5[1];
	linEquations[5][3] = p5[1]*p5[1];
	linEquations[5][4] = p5[0];
	linEquations[5][5] = p5[1];
	linEquations[5][6] = 1.0;

	printf("\n Fitting a conic to five points:\n");
	printf("Point #1: %f,%f\n",p1[0],p1[1]);
	printf("Point #2: %f,%f\n",p2[0],p2[1]);
	printf("Point #3: %f,%f\n",p3[0],p3[1]);
	printf("Point #4: %f,%f\n",p4[0],p4[1]);
	printf("Point #5: %f,%f\n",p5[0],p5[1]);

	// compute the SVD
	double singularValues[7]; // 1..6
	double** nullspaceMatrix = dmatrix(1,6,1,6);
	svdcmp(linEquations, 5, 6, singularValues, nullspaceMatrix);

	// get the result
	printf("\n Singular values: %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6]);

	// find the smallest singular value:
	int smallestIndex = 1;
	for(int i=2;i<7;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

	// solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
	printf("Conic coefficients: %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex]);

	// make sure the solution is correct:
	printf("Equation #1 result: %f\n",	p1[0]*p1[0]*nullspaceMatrix[1][smallestIndex] +
										p1[0]*p1[1]*nullspaceMatrix[2][smallestIndex] +
										p1[1]*p1[1]*nullspaceMatrix[3][smallestIndex] +
										p1[0]*nullspaceMatrix[4][smallestIndex] +
										p1[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #2 result: %f\n",	p2[0]*p2[0]*nullspaceMatrix[1][smallestIndex] +
										p2[0]*p2[1]*nullspaceMatrix[2][smallestIndex] +
										p2[1]*p2[1]*nullspaceMatrix[3][smallestIndex] +
										p2[0]*nullspaceMatrix[4][smallestIndex] +
										p2[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #3 result: %f\n",	p3[0]*p3[0]*nullspaceMatrix[1][smallestIndex] +
										p3[0]*p3[1]*nullspaceMatrix[2][smallestIndex] +
										p3[1]*p3[1]*nullspaceMatrix[3][smallestIndex] +
										p3[0]*nullspaceMatrix[4][smallestIndex] +
										p3[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #4 result: %f\n",	p4[0]*p4[0]*nullspaceMatrix[1][smallestIndex] +
										p4[0]*p4[1]*nullspaceMatrix[2][smallestIndex] +
										p4[1]*p4[1]*nullspaceMatrix[3][smallestIndex] +
										p4[0]*nullspaceMatrix[4][smallestIndex] +
										p4[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #5 result: %f\n",	p5[0]*p5[0]*nullspaceMatrix[1][smallestIndex] +
										p5[0]*p5[1]*nullspaceMatrix[2][smallestIndex] +
										p5[1]*p5[1]*nullspaceMatrix[3][smallestIndex] +
										p5[0]*nullspaceMatrix[4][smallestIndex] +
										p5[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	R2Point test_point(0.34,-2.8);

	printf("A point off the conic: %f\n",	test_point[0]*test_point[0]*nullspaceMatrix[1][smallestIndex] +
											test_point[0]*test_point[1]*nullspaceMatrix[2][smallestIndex] +
											test_point[1]*test_point[1]*nullspaceMatrix[3][smallestIndex] +
											test_point[0]*nullspaceMatrix[4][smallestIndex] +
											test_point[1]*nullspaceMatrix[5][smallestIndex] +
											nullspaceMatrix[6][smallestIndex]);

	return;
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
  // Brighten the image by multiplying each pixel component by the factor.
  // This is implemented for you as an example of how to access and set pixels
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) *= factor;
      Pixel(i,j).Clamp();
    }
  }
}

/*
 * input: int x, -1 <= x <= 1
 *        int y, -1 <= y <= 1
 * output: weights, which are calculated as such:
 *  {-1, 0, 1}       {-1 x 1, 0 x 1, 1 x 1} x = -1
 *  {-2, 0, 2}  ==>  {-1 x 2, 0 x 2, 1 x 2} x = 0
 *  {-1, 0, 1}       {-1 x 1, 0 x 1, 1 x 1} x = 1
 *                    y = -1  y = 0, y = 1
 */
int SobelYKernel(int x, int y) {
  if(x == 0) {
    return y * 2;
  } else {
    return y * 1;
  }
}

void R2Image::
SobelY(void)
{
  // alloc image
  R2Image oldImg(*this);

  // half pixel to be added after computation
  R2Pixel halfPix(0.5,0.5,0.5,1);

  // for all the pixels in the image
  for(int x = 1; x < width - 3; x++) {
    for(int y = 1; y < height -3; y++) {
      // the pixel has to be a blank slate before we put the new value
      Pixel(x,y).Reset(0,0,0,1);
      // find the weight thru kernel
      for(int i = -1; i < 2; i++) {
        for(int j = -1; j < 2; j++) {
          // new image's pixel is calculated using the kernel
          Pixel(x,y) += oldImg.Pixel(x+i,y+j) * SobelYKernel(i,j);
          Pixel(i,j) += halfPix; // add half pix
        }
      }
      //Pixel(x,y).Clamp();
    }
  }

}

/*
 * input: int x, -1 <= x <= 1
 *        int y, -1 <= y <= 1
 * output: weights, which are calculated as such:
 *  {-1, -2, +1}       {-1 x 1, -1 x 2, -1 x 1} x = -1
 *  { 0,  0,  0}  ==>  { 0 x 1,  0 x 2,  0 x 1} x = 0
 *  {+1, +2, +1}       { 1 x 1,  1 x 2,  1 x 1} x = 1
 *                      y = -1   y = 0,  y = 1
 */
int SobelXKernel(int x, int y) {
  if(y == 0) {
    return x * 2;
  } else {
    return x * 1;
  }
}

void R2Image::
SobelX(void)
{
  // alloc image
  R2Image oldImg(*this);

  // half pixel to be added after computation
  R2Pixel halfPix(0.5,0.5,0.5,1);

  // for all the pixels in the image
  for(int x = 1; x < width - 3; x++) {
    for(int y = 1; y < height -3; y++) {
      // the pixel has to be a blank slate before we put the new value
      Pixel(x,y).Reset(0,0,0,1);
      // find the weight thru kernel
      for(int i = -1; i < 2; i++) {
        for(int j = -1; j < 2; j++) {
          // new image's pixel is calculated using the kernel
          Pixel(x,y) += oldImg.Pixel(x+i,y+j) * SobelXKernel(i,j);
          Pixel(i,j) += halfPix; // add half pix
        }
      }
      //Pixel(x,y).Clamp();
    }
  }

}

void R2Image::
LoG(void)
{
  // Apply the LoG oprator to the image

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "LoG() not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{
  // Changes the saturation of an image
  // Find a formula that changes the saturation without affecting the image brightness

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "ChangeSaturation(%g) not implemented\n", factor);
}

// Linear filtering ////////////////////////////////////////////////

/*
 * Gaussian distribution in 1-D form:
 * G(x) = (1 / sqrt(2 * pi) * sigma) * e^(-x^2 / 2 * sigma^2)
 */
double Gaussian(double sigma, int i) {
  double x = i - (3 * sigma);
  // exp(double x) returns the base-e exponential function of x
  double numerator = exp(-(x * x) / (2 * sigma * sigma));
  // sqrt(x) = x^(1/2) = x^(0.5)
  double denominator = pow(2 * M_PI, 0.5) * sigma;
  return numerator / denominator;
}

/*
 * input: int x, -1 <= x <= 1
 *        int y, -1 <= y <= 1
 * output: weights, which are calculated as such:
 *  {     0, -1/4,     0}
 *  { -1/4,     2,  -1/4}
 *  {     0, -1/4,     0}
 */
 double SharpenKernel(int x, int y) {
   if(x == 0 && y == 0) {
     return 2.0;
   } else if(x == 0 || y == 0) {
     return -0.25;
   } else {
     return 0;
   }
 }

void R2Image::
Sharpen()
{
  // alloc image
  R2Image oldImg(*this);

  // half pixel to be added after computation
  R2Pixel halfPix(0.5,0.5,0.5,1);

  // for all the pixels in the image
  for(int x = 1; x < width - 3; x++) {
    for(int y = 1; y < height -3; y++) {
      // the pixel has to be a blank slate before we put the new value
      Pixel(x,y).Reset(0,0,0,1);
      // find the weight thru kernel
      for(int i = -1; i < 2; i++) {
        for(int j = -1; j < 2; j++) {
          // new image's pixel is calculated using the kernel
          Pixel(x,y) += oldImg.Pixel(x+i,y+j) * SharpenKernel(i,j);
          Pixel(i,j) += halfPix; // add half pix
        }
      }
      Pixel(x,y).Clamp();
    }
  }

}

/*
 * returns a value that is within the bound
 * if the value is less than the minimum, it returns the minimum
 * if the value is greater than the maximum, it returns the maximum
 * if the value is within the bounds, it returns the value
 */
int ValWithinBound(int val, int min, int max) {
  if(val < min) {
    return min;
  } else if (val > max) {
    return max;
  } else {
    return val;
  }
}


void R2Image::
Blur(double sigma)
{
  // Initial calculations
  int sigmaInt = (int) sigma;
  int kernelSize = (6 * sigmaInt) + 1;
  // Better way of doing: double weights[kernelSize];
  std::vector<double> weights;
  weights.resize(kernelSize);
  double sumOfWeights = 0;

  // Create the kernel with appropriate weights
  for (int i = 0; i < kernelSize; i++) {
    double weight = Gaussian(sigma, i);
    weights[i] = weight;
    sumOfWeights += weight;
  }

  // Normalize kernel
  for (int i = 0; i < kernelSize; i++) {
    weights[i] /= sumOfWeights;
  }

  // Create a blank temp image
  R2Image tempImg(width, height);

  // First pass in the x direction
  for(int x = 0; x < width; x++) {
    for(int y = 0; y < height; y++) {
      R2Pixel pix;
      for(int lx = -3 * sigmaInt; lx <= 3 * sigmaInt; lx++) {
        // Border processing: if x+lx is less than 0, then use 0.
        //                    if greater than width-1, use width-1
        //                    otherwise, use the value of x+lx
        int val = ValWithinBound(x+lx, 0, width-1);
        pix += Pixel(val, y) * weights[lx + (3 * sigmaInt)];
      }
    tempImg.SetPixel(x, y, pix);
    }
  }
  // second pass in the y direction
  for(int x = 0; x < width; x++) {
    for(int y = 0; y < height; y++) {
      R2Pixel pix;
      for(int ly = -3 * sigmaInt; ly <= 3 * sigmaInt; ly++) {
        // Same border processing stuff
        int val = ValWithinBound(y+ly, 0, height-1);
        pix += tempImg.Pixel(x, val) * weights[ly + (3 * sigmaInt)];
      }
      SetPixel(x, y, pix);
    }
  }
}

void R2Image::
HighPassSharpen(double sigma, double contrast)
{
  R2Image original(*this);

  Blur(sigma);

  R2Image highPass(width, height);
  for(int i = 0; i < width; i++) {
    for(int j = 0; j < height; j++) {
      highPass.Pixel(i,j) = original.Pixel(i,j) - Pixel(i,j);
    }
  }

  R2Image final(width, height);
  for(int i = 0; i < width; i++) {
    for(int j = 0; j < height; j++) {
      final.Pixel(i,j) = highPass.Pixel(i,j) * contrast + Pixel(i,j);
      final.Pixel(i,j).Clamp();
    }
  }

  (*this) = final;
}

R2Image generateHarrisImage(R2Image* orig, double sigma) {
  printf("Getting Harris image... ");
  R2Image Ix2(*orig);
  R2Image Iy2(*orig);
  int width = orig->Width();
  int height = orig->Height();
  R2Image IxIy(width, height);
  R2Image harrisImg(width, height);

  Ix2.SobelX();
  Iy2.SobelY();

  R2Pixel ix;
  R2Pixel iy;
  for(int i = 0; i < width; i++) {
    for(int j = 0; j < height; j++) {
      ix = Ix2[i][j];
      iy = Iy2[i][j];
      IxIy[i][j] = ix * iy;
      Ix2[i][j] = ix * ix;
      Iy2[i][j] = iy * iy;
    }
  }

  Ix2.Blur(sigma);
  Iy2.Blur(sigma);
  IxIy.Blur(sigma);

  R2Pixel gray(0.5, 0.5, 0.5, 1);
  R2Pixel tmp;
  for(int i = 0; i < width; i++) {
    for(int j = 0; j < height; j++) {
      tmp = Ix2[i][j] * Iy2[i][j]
                             - IxIy[i][j] * IxIy[i][j]
                             - (0.04 * (Ix2[i][j] + Iy2[i][j])
                                     * (Ix2[i][j] + Iy2[i][j]))
                             + gray ;
      tmp.Clamp();
      harrisImg.SetPixel(i,j,tmp);
    }
  }
  printf("Done\n");
  return harrisImg;
}

typedef struct Point {
  int x;
  int y;
  float RGBsum;
} Point;

// if p1 comes before p2, returns true
struct PointComparator {
  bool operator()(Point a, Point b) {
    return a.RGBsum < b.RGBsum;
  }
};



void getFeaturePoints(R2Image* harris, std::vector<Point> &featurePoints, int numFeaturePoints, int borderX, int borderY ) {
  int width = harris->Width();
  int height = harris->Height();
  const int minDistance = 20;

  // initializing the data structure that tells us whether
  // the point is a valid feature point or not
  bool allowed[width][height];
  for(int i = 0; i < width; i++) {
    for(int j = 0; j < height; j++) {
      allowed[i][j] = true;
    }
  }

  // second argument vector is there bc
  // the C++ STL PQ is a container adapter
  std::priority_queue<Point, std::vector<Point>, PointComparator> pq;

  // fill out the PQ with the pixels' info
  // since we are using a PQ, the points should be ordered
  // based on their priority -- intensity of the pixel
  for(int i = borderX; i < width-borderX; i++) {
    for(int j = borderY; j < height - borderY; j++) {
      R2Pixel cur = harris->Pixel(i,j);
      Point p;
      p.x = i;
      p.y = j;
      p.RGBsum = cur.Red() + cur.Green() + cur.Blue();
      pq.push(p);
    }
  }

  int numPointsSoFar = 0;
  printf("Fetching feature points... ");
  while(numPointsSoFar < numFeaturePoints) {
    // get the point with highest priority so far
    Point p = pq.top();
    // check whether the point is at least 10px away from another FP
    if(allowed[p.x][p.y]) {
      featurePoints[numPointsSoFar] = p;
      // invalidate the points within 10 pixels distance
      for(int i = -1 * minDistance; i <= minDistance; i++) {
        for(int j = -1 * minDistance; j <= minDistance; j++) {
          int xi = ValWithinBound(p.x+i, 0, width-1);
          int yj = ValWithinBound(p.y+j, 0, height-1);
          allowed[xi][yj] = false;
        }
      }
      numPointsSoFar++;
    }
    pq.pop();
  }
  printf("Done\n");
}

void MarkPoints(R2Image &img, Point p, R2Pixel color) {

  const int size = 3;
  for(int i = -1 * size; i <= size; i++) {
    for(int j = -1 * size; j <= size; j++) {
      int x = ValWithinBound(p.x + i, 0, img.Width());
      int y = ValWithinBound(p.y + j, 0, img.Height());
      img.SetPixel(x,y,color);
    }
  }
}

void R2Image::line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
  if(x0>x1)
  {
    int x=y1;
    y1=y0;
    y0=x;

    x=x1;
    x1=x0;
    x0=x;
  }
  int deltax = x1 - x0;
  int deltay = y1 - y0;
  float error = 0;
  float deltaerr = 0.0;
  if(deltax!=0) deltaerr =fabs(float(float(deltay) / deltax));    // Assume deltax != 0 (line is not vertical),
  // note that this division needs to be done in a way that preserves the fractional part
  int y = y0;
  for(int x=x0;x<=x1;x++)
  {
    Pixel(x,y).Reset(r,g,b,1.0);
    error = error + deltaerr;
    if(error>=0.5)
    {
      if(deltay>0) y = y + 1;
      else y = y - 1;

      error = error - 1.0;
    }
  }
  if(x0>3 && x0<width-3 && y0>3 && y0<height-3)
  {
    for(int x=x0-3;x<=x0+3;x++)
    {
      for(int y=y0-3;y<=y0+3;y++)
      {
        Pixel(x,y).Reset(r,g,b,1.0);
      }
    }
  }
}


void R2Image::
Harris(double sigma)
{
  R2Image harris = generateHarrisImage(this, sigma);
  int width = (this)->Width();
  int height = (this)->Height();
  int windowX = (int) (0.2f * width);
  int windowY = (int) (0.2f * height);

  const int numFeaturePoints = 150;
  std::vector<Point> featurePoints(numFeaturePoints);
  //without Border
  getFeaturePoints(&harris, featurePoints, numFeaturePoints, 0, 0);
  //with Border
  //getFeaturePoints(&harris, featurePoints, numFeaturePoints, windowX, windowY);
  for(int i = 0; i < numFeaturePoints; i++) {
    MarkPoints(*this, featurePoints[i], R2Pixel(1, 0, 1, 0.5));
  }
  //(*this) = harris;
}

double GetSSDOf(R2Image* I0, R2Image* I1, Point feature,
  int pointx, int pointy, int radius) {
  double* comp0;
  double* comp1;
  double diff = 0;
  for(int i = -1 * radius; i <= radius; i++) {
    for(int j = -1 * radius; j <= radius; j++) {
      comp0 = I0->Pixel(feature.x+i, feature.y+j).Components();
      comp1 = I1->Pixel(pointx + i,pointy + j).Components();
      for(int k = 0; k < 3; k++) {
        diff += (comp0[k] - comp1[k]) * (comp0[k] - comp1[k]);
      }
    }
  }
  return diff;
}

void track(R2Image * featureImage, R2Image * compareImage, int numFeaturePoints,
    std::vector<Point> &features, std::unordered_map<int,Point> &trackedFeatures) {
  int width = featureImage->Width();
  int height = featureImage->Height();
  // 20% image size search window
  int windowX = (int) (0.2f * width);
  int windowY = (int) (0.2f * height);

  int sigma = 2;
  int radius = 6 * sigma + 1;

  R2Image harris = generateHarrisImage(featureImage, sigma);

  //without Border
  getFeaturePoints(&harris, features, numFeaturePoints, 0, 0);
  //with Border
  //getFeaturePoints(&harris, features, numFeaturePoints, windowX, windowY);

  double DOUBLE_MAX = std::numeric_limits<double>::max();
  Point bestSoFar;

  for(int i = 0; i < numFeaturePoints; i++) {
    printf("Tracking point %d of 150", i+1);
    Point curFeature = features[i];
    int startX = std::max(radius, curFeature.x - windowX/2);
    int startY = std::max(radius, curFeature.y - windowY/2);
    int endX = std::min(width-radius, curFeature.x + windowX/2);
    int endY = std::min(height-radius, curFeature.y + windowY/2);
    bestSoFar.RGBsum = DOUBLE_MAX;
    double ssd;

    for(int j = startX; j < endX; j++) {
      for(int k = startY; k < endY; k++) {
        ssd = GetSSDOf(featureImage, compareImage, curFeature, j, k, radius);
        if(ssd < bestSoFar.RGBsum) {
          bestSoFar.RGBsum = ssd;
          bestSoFar.x = j;
          bestSoFar.y = k;
        }
      }
    }
    if(bestSoFar.RGBsum != DOUBLE_MAX) {
      trackedFeatures[i] = bestSoFar;
    }
    std::cout<<'\r';
    std::cout.flush();
  }
  printf("Tracking point 150 of 150... Done\n");
  return;
}

void R2Image::
trackWithRANSAC(R2Image *other) {

  int numFeaturePoints = 150;
  std::vector<Point> features(numFeaturePoints);
  std::unordered_map<int, Point> trackedFeatures(numFeaturePoints);
  track(this, other, numFeaturePoints, features, trackedFeatures);

  int numTrials = 100;
  double acceptThreshold = 4;
  int numInliers = 4;

  int bestNumMatches = -1;
  double bestAvgDeltaX = 0;
  double bestAvgDeltaY = 0;

  printf("Running RANSAC... ");
  for(int trial = 0; trial < numTrials; trial++) {
    int inlier[numInliers];
    double avgDeltaX = 0;
    double avgDeltaY = 0;
    for(int i = 0; i < numInliers; i++) {
      inlier[i] = std::rand() % numFeaturePoints;
      avgDeltaX += trackedFeatures[inlier[i]].x - features[inlier[i]].x;
      avgDeltaY += trackedFeatures[inlier[i]].y - features[inlier[i]].y;
    }
    avgDeltaX /= numInliers;
    avgDeltaY /= numInliers;

    int numMatches = 0;
    for(int i = 0; i < numFeaturePoints; i++) {
      if(trackedFeatures.count(inlier[i]) > 0) {
        double deltaX =  trackedFeatures[inlier[i]].x - features[inlier[i]].x;

        double deltaY =  trackedFeatures[inlier[i]].y - features[inlier[i]].y;
        if(std::abs(deltaX - avgDeltaX) <= acceptThreshold
            && std::abs(deltaY - avgDeltaY) <= acceptThreshold) {
          numMatches++;
        }
      }
    }
    if (numMatches > bestNumMatches) {
      bestNumMatches = numMatches;
      bestAvgDeltaX = avgDeltaX;
      bestAvgDeltaY = avgDeltaY;
    }
  }
  printf("Done\n");

  for(int i = 0; i < numFeaturePoints; i++) {
    if(trackedFeatures.count(i) > 0) {
      double deltaX = trackedFeatures[i].x - features[i].x;
      double deltaY = trackedFeatures[i].y - features[i].y;
      if(std::abs(deltaX - bestAvgDeltaX) <= acceptThreshold
          && std::abs(deltaY - bestAvgDeltaY) <= acceptThreshold) {
        line(features[i].x, trackedFeatures[i].x,
            features[i].y, trackedFeatures[i].y,
            0, 1, 0);
      } else {
         line(features[i].x, trackedFeatures[i].x,
            features[i].y, trackedFeatures[i].y,
            1, 0, 0);
      }
    }
  }
}

void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
  const int numFeaturePoints = 150;
  std::vector<Point> features(numFeaturePoints);
  std::unordered_map<int, Point> trackedFeatures(numFeaturePoints);
  track(this, otherImage, numFeaturePoints, features, trackedFeatures);
  *this = *otherImage;

  for(int i = 0; i < numFeaturePoints; i++) {
    if(trackedFeatures.count(i) > 0) {
      line(features[i].x, trackedFeatures[i].x,
        features[i].y, trackedFeatures[i].y,
        0, 1, 0);
      //MarkPoints(*this, features[i], R2Pixel(1, 0, 1, 1));
    }
  }
}

void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching homography, and blend the transformed "otherImage" into this image with a 50% opacity.
	fprintf(stderr, "fit other image using a homography and blend imageB over imageA\n");
	return;
}

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp);
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);

  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);

  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);

  // Check info header
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }

  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }

  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
	// Read pixel values
	int red, green, blue;
	if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
	  fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
	  fclose(fp);
	  return 0;
	}

	// Assign pixel values
	double r = (double) red / max_value;
	double g = (double) green / max_value;
	double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }

    // Close file
    fclose(fp);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
  extern "C" {
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}




int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width; 	/* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 95, TRUE);
  jpeg_start_compress(&cinfo, TRUE);

  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}
