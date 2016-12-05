#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <map>

using namespace std;

int main(int argc, char** argv) {

  // arg 1: bin annotations
  ifstream binannot(argv[1]);

  // arg 2: homologous bins
  ifstream hombins(argv[2]);
  
  // arg 3: maxdist
  stringstream maxdarg(argv[3]);
  int maxd;
  maxdarg >> maxd;

  // arg 4: new homologous bins
  ofstream newhom(argv[4]);

  // arg 5: new homologous bins + neighbors
  ofstream newneighb(argv[5]);

  // parse bin annotations 
  string line;
  vector<string> chrs;
  string chr;
  int bin;
  int nbin = 0;
  while (getline(binannot,line)) {
    stringstream ss(line);
    ss >> bin;
    ss >> chr;
    chrs.push_back(chr);
    if (bin > nbin) nbin = bin;
  }
  binannot.close();
  nbin++;

  // initialize matrices
  int* outmatr = new int[nbin*nbin];
  int* neighbors = new int[nbin*nbin];
  int* newneighbors = new int[nbin*nbin];
  for (long i = 0; i < nbin*nbin; i++) {
    outmatr[i]=0;
    neighbors[i]=0;
    newneighbors[i]=0;
  }
  
  // parse homologous bins
  int bin1, bin2;
  while (getline(hombins,line)) {
    stringstream ss(line);
    ss >> bin1;
    ss >> bin2;
    outmatr[nbin*bin1+bin2] = 1;
  }

  // count neighbors (including self) for each bin
  for (int i = 0; i < nbin; i++) {
    for (int j = 0; j < nbin; j++) {
      for (int a = -maxd; a <= maxd; a++) {
        for (int b = -maxd; b <= maxd; b++) {
          if (((i+a) >= 0 && (i+a) < nbin) && ((j+b) >= 0 && (j+b) < nbin)) {
            if (outmatr[nbin*(i+a) + (j+b)]) neighbors[nbin*i + j]++;
          }
        }
      }
    }
  }

  // remove homologous bins with fewer than 2 other neighbors in maxd radius
  for (int i = 0; i < nbin; i++) {
    for (int j = 0; j < nbin; j++) {
      if (neighbors[nbin*i + j] < 3) outmatr[nbin*i + j] = 0;
    }
  }

  for (int i = 0; i < nbin; i++) {
    for (int j = 0; j < nbin; j++) {
      if (outmatr[nbin*i + j]) {
        newhom << i << '\t' << j << '\n';
        for (int a = -maxd; a <= maxd; a++) {
          for (int b = -maxd; b <= maxd; b++) {
            if (((i+a) >= 0 && (i+a) < nbin) && ((j+b) >= 0 && (j+b) < nbin)) {
              if ((chrs[i].compare(chrs[i+a])==0) && (chrs[j].compare(chrs[j+b])==0)) {
                newneighbors[nbin*(i+a) + (j+b)] = 1;
              }
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < nbin; i++) {
    for (int j = 0; j < nbin; j++) {
      if (newneighbors[nbin*i + j]) {
        newneighb << i << '\t' << j << '\n';
      }
    }
  }

  delete[] outmatr;
  delete[] neighbors;

  newhom.close();
  newneighb.close();
}
