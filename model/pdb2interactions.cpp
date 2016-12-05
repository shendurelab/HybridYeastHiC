#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

// data structure for 3D coordinates
struct coord {
  double x;
  double y;
  double z;
  string chr;
  coord(double a, double b, double c, string chrom) {
    x = a;
    y = b;
    z = c;
    chr = chrom;
  }
};

// distance function for two 3D coordinates
double dist(coord a, coord b) {
  return pow(pow(a.x-b.x,2) + pow(a.y-b.y,2) + pow(a.z-b.z,2),0.5);
}

int main(int argc, char** argv) {

  bool time_run = false;
  double contact_dist = 45;

  clock_t t1, t2;
  t1 = clock();

  stringstream arg(argv[1]);
  arg >> contact_dist;

  // read PDB file
  ifstream pdb(argv[2]);
  string line;
  vector<coord> coords;
  while (getline(pdb,line)) {
    stringstream lstream(line);
    double x, y, z;
    string chrom;
    lstream.ignore(30);
    lstream >> x;
    lstream >> y;
    lstream >> z;
    lstream >> chrom;
    coords.push_back(coord(x,y,z,chrom));
  }
  pdb.close();

  // check for interactions
  for (int i = 0; i < coords.size(); i++) {
    for (int j = i+1; j < coords.size(); j++) {
      double d = dist(coords[i],coords[j]);
      if (d < contact_dist) {
        if ((coords[i].chr.compare(coords[j].chr)!=0) || (j - i > 2)) {
          printf("%d\t%d\t%.3f\n",i,j,d);
        }
      }
    }
  }
  
  t2 = clock();
  if (time_run) { 
    float diff = (float) t2 - (float) t1;
    float seconds = diff / CLOCKS_PER_SEC;
    cout << "Time to parse: " << seconds << '\n';
  }

  // print matrix
/*  cout << "seqpos";
  for (int i = 0; i < coords.size(); i++)
    cout << '\t' << i;
  cout << '\n';
  for (int i = 0; i < coords.size(); i++) {
    cout << i;
    for (int j = 0; j < coords.size(); j++) {
      cout << '\t' << matrix[i][j];
    }
    cout << '\n';
  }
*/
 

}
