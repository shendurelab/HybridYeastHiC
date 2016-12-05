/* digest_genome.cpp
 * in silico restriction digestion of genome
 * 
 * Seungsoo Kim
 * USAGE: digest_genome GENOME_FILE RESTRICTION_SITE RESTRICTION_OFFSET
 */
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

int main(int argc, char** argv) {

  // check number of arguments
  if (argc < 3) {
    cout << "USAGE: digest_genome GENOME_FILE RESTRICTION_SITE RESTRICTION_OFFSET\n";
  }

  // arg 1: genome fasta
  ifstream genome_f(argv[1]);
  string line;
  string prevline;
  string chrom;
  long line_pos;

  // arg 2: restriction site
  string resite;
  stringstream ss(argv[2]);
  ss >> resite;

  // arg 3: restriction offset (distance between start of binding site and cut site)
  stringstream parseoffset(argv[3]);
  int offset;
  parseoffset >> offset;

  vector<string> chroms;
  vector<vector<long> > resites;
  vector<long> chrom_length;
  while (getline(genome_f,line)) {
    // new chromosome
    if (line[0]=='>') {
      // extract chromosome name
      chrom = line.substr(1);
      // add to list of chromosomes
      chroms.push_back(chrom);
      // start counting length of chromosome
      chrom_length.push_back(0);
      // start list of restriction sites
      resites.push_back(vector<long>(1,0));
      // reset counter variables
      line_pos = 0;
      prevline = "";
    }
    // continuing on same chromosome
    else {
      // add to chromosome length
      chrom_length.back() += line.length();
      // if not first real line of sequence, buffer sequence with previous line
      if (line_pos > 0)
        line = prevline.substr(prevline.length()-resite.length()) + line;
      size_t pos = 0;
      // repeat for all restriction sites:
      do {
        // find restriction site
        pos = line.find(resite,pos+1);
        if (pos != string::npos) {
          resites.back().push_back(line_pos + pos + offset);
        }
      } while (pos != string::npos);
      line_pos += line.length() - resite.length();
      prevline = line;
    }
  }
  genome_f.close();

  // print restriction fragments  
  for (int i = 0; i < chroms.size(); i++) {
    for (int j = 0; j < resites[i].size() - 1; j++) {
      cout << chroms[i] << '\t' << resites[i][j] << '\t' << resites[i][j+1] << '\t' << chroms[i] << '_' << j + 1 << "\t0\t+\n";
    }
    cout << chroms[i] << '\t' << resites[i].back() << '\t' << chrom_length[i] << '\t' << chroms[i] << '_' << resites[i].size() << "\t0\t+\n";
  }
}
