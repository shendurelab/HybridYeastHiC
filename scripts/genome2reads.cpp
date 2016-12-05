/* genome2reads.cpp
 * read in genome, every WINSIZE bases take READLEN bases
 * USAGE: genome2reads GENOME WINSIZE READLEN
 * 
 * Seungsoo Kim
 */
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char** argv) {

  int readlen;
  int winsize;
  if (argc != 3) {
    cout << "USAGE: genome2reads GENOME WINSIZE READLEN\n";
    return 1;
  }

  // arg 1: genome fasta
  ifstream genome(argv[1]);
  string line;

  // arg 2: window size
  stringstream arg1(argv[2]);
  arg1 >> winsize;

  // arg 3: read length
  stringstream arg2(argv[3]);
  arg2 >> readlen;

  string chr;
  string seq;
  while (getline(genome,line)) {
    if (line[0]=='>') {
      if (seq.length() != 0) {
        for (int p = 0; p < seq.length()-readlen; p += winsize) {
          cout << '>' << chr << '-' << p + 1 << '\n';
          cout << seq.substr(p,readlen) << '\n';
        }
      }
      chr = line.substr(1);
      seq.clear();
    }
    else {
      seq += line;
    }
  }
  for (int p = 0; p < seq.length()-readlen; p += winsize) {
    cout << '>' << chr << '-' << p + 1 << '\n';
    cout << seq.substr(p,readlen) << '\n';
  }

  genome.close();
}
