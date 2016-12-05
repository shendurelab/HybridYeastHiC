#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

int main(int argc, char** argv) {

  bool time_run = false;
  double contact_dist;

  clock_t t1, t2;
  t1 = clock();

  // arg 1: bins file
  ifstream bins_f(argv[1]);
  string line;
  vector<int> mapping;
  while (getline(bins_f,line)) {
    stringstream lstream(line);
    int binn;
    lstream >> binn;
    mapping.push_back(binn);
  }
  bins_f.close();
  int nbin = mapping.back() + 1;
  
  // arg 2: max distance for contact
  double maxdist;
  stringstream dist_arg(argv[2]);
  dist_arg >> maxdist;

  // arg 3: minimum average value across row
  stringstream minrowavg_arg(argv[3]);
  double minrowavg;
  minrowavg_arg >> minrowavg;
  long minrowsum = (long) (minrowavg * nbin);
  
  // arg 4: interactions file
  ifstream interactions_f(argv[4]);

  // arg 5: output row sum file
  ofstream outrowsumf(argv[5]);

  // arg 6: output matrix
  ofstream outmatrf(argv[6]);
  
  // initialize matrices
  long* outmatr = new long[nbin*nbin];
  double* normmatr = new double[nbin*nbin];
  for (long i = 0; i < nbin*nbin; i++) {
    outmatr[i]=0;
    normmatr[i]=nan("");
  }

  // fill matrix from interactions file
  while (getline(interactions_f,line)) {
    stringstream lstream(line);
    int b1;
    int b2;
    double dist;
    lstream >> b1;
    lstream >> b2;
    lstream >> dist;
    // exclude interactions in same bin or not within max interaction distance
    if ((dist < maxdist) && (mapping[b1] != mapping[b2])) {
      outmatr[mapping[b1]*nbin + mapping[b2]]++;
      outmatr[mapping[b2]*nbin + mapping[b1]]++;
    }
  }
  interactions_f.close();

  // calculate row sums
  long* rowsums = new long[nbin];
  for (int i = 0; i < nbin; i++) {
    rowsums[i] = 0;
    for (int j = 0; j < nbin; j++) {
      rowsums[i] += outmatr[nbin*i + j];
    }
    outrowsumf << i << '\t' << rowsums[i] << '\n';
  }

  // calculate total sum, excluding rows that don't meet minimum
  long allsum = 0;
  for (int i = 0; i < nbin; i++) {
    if (rowsums[i] >= minrowsum) {
      allsum += rowsums[i];
    }
  }

  // calculate coverage-normalized matrix
  for (int i = 0; i < nbin; i++) {
    if (rowsums[i] >= minrowsum) {
      for (int j = 0; j < nbin; j++) {
        if (rowsums[j] >= minrowsum) {
          normmatr[nbin*i + j] = ((double)outmatr[nbin*i + j]*allsum)/((double)rowsums[i]*rowsums[j]);
        }
      }
    }
  }

  // print matrix, one line per entry of matrix
  for (int i = 0; i < nbin; i++) {
    for (int j = 0; j < nbin; j++) {
      outmatrf << i << '\t' << j << '\t' << outmatr[nbin*i + j] << '\t' << normmatr[nbin*i + j] << '\n';
    }
  }

  delete[] outmatr;
  delete[] normmatr;
  delete[] rowsums;
  
  t2 = clock();
  if (time_run) { 
    float diff = (float) t2 - (float) t1;
    float seconds = diff / CLOCKS_PER_SEC;
    cout << "Time to parse: " << seconds << '\n';
  }
}
