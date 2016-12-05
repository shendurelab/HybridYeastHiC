/* process_pairs_dual.cpp
 * 
 * Seungsoo Kim
 */
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <utility>
#include <unordered_set>


template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

using namespace std;

// function to find end position of read alignment from CIGAR string
int readend(string cig) {
  stringstream ss(cig);
  int pos = 0;
  int nts;
  char type;
  while (ss >> nts, ss >> type, !ss.eof()) {
    if (type=='M' || type=='D') {
      pos += nts;
    }
  }
  return pos;
}

// data type for read pairs
struct readpair {
  bool strand1;
  string chrom1;
  long coord1;
  bool strand2;
  string chrom2;
  long coord2;
  bool operator==(const readpair & rp) const {
    return (rp.chrom1.compare(chrom1) == 0 &&
      rp.chrom2.compare(chrom2) == 0 &&
      rp.strand1 == strand1 &&
      rp.strand2 == strand2 &&
      rp.coord1 == coord1 &&
      rp.coord2 == coord2);
  }
};

namespace std {
  template <> struct hash<readpair> {
    inline size_t operator()(const readpair & a) const {
      size_t seed = 0;
      hash_combine(seed, a.chrom1);
      hash_combine(seed, a.chrom2);
      hash_combine(seed, a.strand1);
      hash_combine(seed, a.strand2);
      hash_combine(seed, a.coord1);
      hash_combine(seed, a.coord2);
      return seed;
    }
  };
}

int main(int argc, char** argv) {
  // parameters

  if (argc < 9) {
    cout << "USAGE: ./process_pairs_dual RESTRICTION_SITE GENOME_DIGEST MAPQ READ1A_FILE READ1B_FILE READ2A_FILE READ2B_FILE READPAIR_OUT STAT_OUT\n";
    return 1;
  }

  // build map (binary search tree) of restriction fragments, with a vector for each chromosome
  map<string, vector<pair<long,long> > > restr_map;

  // arg 1: restriction site
  string resite;
  stringstream ss(argv[1]);
  ss >> resite;

  // arg 2: restriction digest file
  ifstream genome_digest(argv[2]);
  
  // arg 3: mapq
  int mapq_filt;
  stringstream ss2(argv[3]);
  ss2 >> mapq_filt;

  // arg 4: shotgun threshold
  int shot_filt;
  stringstream ss3(argv[4]);
  ss3 >> shot_filt;

  // arg 5: read 1a file
  ifstream read1afile(argv[5]);
  // arg 6: read 1b file
  ifstream read1bfile(argv[6]);
  // arg 7: read 2a file
  ifstream read2afile(argv[7]);
  // arg 8: read 2b file
  ifstream read2bfile(argv[8]);
  // arg 9: output file for read pairs
  ofstream readout(argv[9]);
  // arg 10: output file for stats
  ofstream statout(argv[10]);

  string line;
  string chr;
  long pos1, pos2;
  long pos1a, pos1b, pos2a, pos2b;

  // build restriction map
  while (getline(genome_digest,line)) {
    // parse line of restriction digest
    stringstream ss(line);
    ss >> chr;
    ss >> pos1;
    ss >> pos2;
    // if new chr, add new vector
    if (restr_map.find(chr) != restr_map.end()) {
      restr_map.insert(make_pair(chr,vector<pair<long,long> >()));
    }
    // add new element to vector for each restriction fragment
    restr_map[chr].push_back(make_pair(pos1+1,pos2));
  }
  genome_digest.close();

  string ignore_entry;
  string rname1a, rname1b, rname2a, rname2b, rname1, rname2;
  int dir1a, dir1b, dir2a, dir2b, dir1, dir2;
  string chr1a, chr1b, chr2a, chr2b, chr1, chr2;
  int mapq1a, mapq1b, mapq2a, mapq2b, mapq1, mapq2;
  string cig1a, cig1b, cig2a, cig2b, cig1, cig2;
  string seq1a, seq1b, seq2a, seq2b, seq1, seq2;
  string qual1a, qual1b, qual2a, qual2b, qual1, qual2;
  int as1a, as1b, as2a, as2b;
  int xm1a, xm1b, xm2a, xm2b;
  int xo1a, xo1b, xo2a, xo2b;
  
  // data structure to store all read pairs, to check for duplicates
  unordered_set<readpair> readpairs;

  // counters for read pair categories
  long fail_mapq, dup, intra, inter, shot, same_re, same_dir, self_lig = 0;

  size_t aspos, xmpos, xopos;

  while (getline(read1afile,line)) {
    // parse line from read 1a file
    stringstream ss1a(line);
    ss1a >> rname1a;
    ss1a >> dir1a;
    ss1a >> chr1a;
    ss1a >> pos1a;
    ss1a >> mapq1a;
    ss1a >> cig1a;
    ss1a >> ignore_entry;
    ss1a >> ignore_entry;
    ss1a >> ignore_entry;
    ss1a >> seq1a;
    ss1a >> qual1a;
    aspos = line.find("AS:i:");
    if (aspos != string::npos) {
      stringstream ss(line.substr(aspos+5));
      ss >> as1a;
    }
    xmpos = line.find("XM:i:");
    if (xmpos != string::npos) {
      stringstream ss(line.substr(xmpos+5));
      ss >> xm1a;
    }
    xopos = line.find("XO:i:");
    if (xopos != string::npos) {
      stringstream ss(line.substr(xopos+5));
      ss >> xo1a;
    }

    // parse line from read 1b file
    getline(read1bfile,line);
    stringstream ss1b(line);
    ss1b >> rname1b;
    ss1b >> dir1b;
    ss1b >> chr1b;
    ss1b >> pos1b;
    ss1b >> mapq1b;
    ss1b >> cig1b;
    ss1b >> ignore_entry;
    ss1b >> ignore_entry;
    ss1b >> ignore_entry;
    ss1b >> seq1b;
    ss1b >> qual1b;
    aspos = line.find("AS:i:");
    if (aspos != string::npos) {
      stringstream ss(line.substr(aspos+5));
      ss >> as1b;
    }
    xmpos = line.find("XM:i:");
    if (xmpos != string::npos) {
      stringstream ss(line.substr(xmpos+5));
      ss >> xm1b;
    }
    xopos = line.find("XO:i:");
    if (xopos != string::npos) {
      stringstream ss(line.substr(xopos+5));
      ss >> xo1b;
    }

    // parse line from read 2a file
    getline(read2afile,line);
    stringstream ss2a(line);
    ss2a >> rname2a;
    ss2a >> dir2a;
    ss2a >> chr2a;
    ss2a >> pos2a;
    ss2a >> mapq2a;
    ss2a >> cig2a;
    ss2a >> ignore_entry;
    ss2a >> ignore_entry;
    ss2a >> ignore_entry;
    ss2a >> seq2a;
    ss2a >> qual2a;
    aspos = line.find("AS:i:");
    if (aspos != string::npos) {
      stringstream ss(line.substr(aspos+5));
      ss >> as2a;
    }
    xmpos = line.find("XM:i:");
    if (xmpos != string::npos) {
      stringstream ss(line.substr(xmpos+5));
      ss >> xm2a;
    }
    xopos = line.find("XO:i:");
    if (xopos != string::npos) {
      stringstream ss(line.substr(xopos+5));
      ss >> xo2a;
    }

    // parse line from read 2b file
    getline(read2bfile,line);
    stringstream ss2b(line);
    ss2b >> rname2b;
    ss2b >> dir2b;
    ss2b >> chr2b;
    ss2b >> pos2b;
    ss2b >> mapq2b;
    ss2b >> cig2b;
    ss2b >> ignore_entry;
    ss2b >> ignore_entry;
    ss2b >> ignore_entry;
    ss2b >> seq2b;
    ss2b >> qual2b;
    aspos = line.find("AS:i:");
    if (aspos != string::npos) {
      stringstream ss(line.substr(aspos+5));
      ss >> as2b;
    }
    xmpos = line.find("XM:i:");
    if (xmpos != string::npos) {
      stringstream ss(line.substr(xmpos+5));
      ss >> xm2b;
    }
    xopos = line.find("XO:i:");
    if (xopos != string::npos) {
      stringstream ss(line.substr(xopos+5));
      ss >> xo2b;
    }

    // for each read, must map perfectly to one genome and at least 1 SNP and 2 SNP + indel to other genome
    if ((mapq1a < 30) || (mapq1b < 30) || (mapq2a < 30) || (mapq2b < 30)) {
      fail_mapq++;
      continue;
    }
    if (!(((as1a == 0) && (xm1b >= 1) && (xm1b + xo1b >= 2)) || ((as1b == 0) && (xm1a >= 1) && (xm1a + xo1a >= 2)))) {
      fail_mapq++;
      continue;
    }
    if (!(((as2a == 0) && (xm2b >= 1) && (xm2b + xo2b >= 2)) || ((as2b == 0) && (xm2a >= 1) && (xm2a + xo2a >= 2)))) {
      fail_mapq++;
      continue;
    }

    if (xm1a == 0) {
      rname1 = rname1a;
      dir1 = dir1a;
      chr1 = chr1a;
      pos1 = pos1a;
      mapq1 = mapq1a;
      cig1 = cig1a;
    }
    else {
      rname1 = rname1b;
      dir1 = dir1b;
      chr1 = chr1b;
      pos1 = pos1b;
      mapq1 = mapq1b;
      cig1 = cig1b;
    }
    if (xm2a == 0) {
      rname2 = rname2a;
      dir2 = dir2a;
      chr2 = chr2a;
      pos2 = pos2a;
      mapq2 = mapq2a;
      cig2 = cig2a;
    }
    else {
      rname2 = rname2b;
      dir2 = dir2b;
      chr2 = chr2b;
      pos2 = pos2b;
      mapq2 = mapq2b;
      cig2 = cig2b;
    }

    // find end positions of read alignments
    long pos1end = pos1 - 1 + readend(cig1);
    long pos2end = pos2 - 1 + readend(cig2);

    // determine end of read alignment further from end of fragment that was ligated
    long map1end = (dir1 == 0) ? pos1 : pos1end;
    long map2end = (dir2 == 0) ? pos2 : pos2end;

    // if not a duplicate, insert read pair
    readpair newpair;
    newpair.strand1 = (dir1 == 0);
    newpair.chrom1 = chr1;
    newpair.coord1 = map1end;
    newpair.strand2 = (dir2 == 0);
    newpair.chrom2 = chr2;
    newpair.coord2 = map2end;

    if (readpairs.count(newpair)==0) {
      readpairs.insert(newpair);
    }
    else {
      dup++;
      continue;
    }

    // find restriction fragment ends
    long re1start, re1end, re2start, re2end;

    // binary search for restriction fragment 1 start
    int upper = restr_map[chr1].size();
    int lower = 0;
    while (upper != lower) {
      int test = (upper+lower)/2;
      if (restr_map[chr1][test].first <= pos1) {
        lower = test;
      }
      if (restr_map[chr1][test].second >= pos1) {
        upper = test;
      }
    }
    re1start = restr_map[chr1][lower].first;

    // search for restriction fragment 1 end starting from fragment 1 start
    while (restr_map[chr1][upper].second + resite.length() < pos1end) {
      upper++;
    }
    re1end = restr_map[chr1][upper].second + resite.length();

    // binary search for restriction fragment 2 start
    upper = restr_map[chr2].size();
    lower = 0;
    while (upper != lower) {
      int test = (upper+lower)/2;
      if (restr_map[chr2][test].first <= pos2) {
        lower = test;
      }
      if (restr_map[chr2][test].second >= pos2) {
        upper = test;
      }
    }
    re2start = restr_map[chr2][lower].first;
    // search for restriction fragment 2 end starting from fragment 2 start
    while (restr_map[chr2][upper].second + resite.length() < pos2end) {
      upper++;
    }
    re2end = restr_map[chr2][upper].second + resite.length();

    // determine length of fragment
    int fraglen = ((dir1 == 0) ? (re1end - pos1 + 1) : (pos1end - re1start + 1)) + ((dir2 == 0) ? (re2end - pos2 + 1) : (pos2end - re2start + 1));
    
    string pairtype;
    // if on same chromosome
    if (chr1.compare(chr2)==0) {
      // if facing toward and insert size less than threshold
      if ((dir1 != dir2) && (((dir1==0) && (pos1 < pos2end) && (pos2end-pos1+1 < shot_filt)) || ((dir2==0) && (pos2 < pos1end) && (pos1end-pos2+1 < shot_filt)))) {
        pairtype = "invalid_shotgun";
        fraglen = (dir1 == 0) ? (pos2end-pos1+1) : (pos1end-pos2+1);
        shot++;
      }
      // if both reads are in same RE frag
      else if ((re1start == re2start) || (re1end == re2end)) {
        // if facing same way
        if (dir1 == dir2) {
          pairtype = "invalid_samedir";
          same_dir++;
        }
        else {
          // if facing toward
          if (((dir1 == 0) && (pos1 < pos2)) || ((dir2 == 0) && (pos2 < pos1))) {
            pairtype = "invalid_sameRE";
            fraglen = (dir1 == 0) ? (pos2end-pos1+1) : (pos1end-pos2+1);
            same_re++;
          }
          // if facing away
          else {
            pairtype = "invalid_selflig";
            self_lig++;
          }
        }
      }
      else {
        pairtype = "intrachromosomal";
        intra++;
      }
    }
    // if on different chromosomes
    else {
      pairtype = "interchromosomal";
      inter++;
    }
    // print read pair info
    readout << chr1 << '\t' << pos1 << '\t' << pos1end << '\t' << map1end << '\t' << ((dir1 == 0) ? re1end : re1start) << '\t' << dir1 << '\t';
    readout << chr2 << '\t' << pos2 << '\t' << pos2end << '\t' << map2end << '\t' << ((dir2 == 0) ? re2end : re2start) << '\t' << dir2 << '\t';
    readout << fraglen << '\t' << pairtype << '\n';
  }
  readout.close();

  // print read pair statistics
  statout << "low_mapq\t" << fail_mapq << '\n';
  statout << "duplicates\t" << dup << '\n';
  statout << "invalid_shotgun\t" << shot << '\n';
  statout << "invalid_samedir\t" << same_dir << '\n';
  statout << "invalid_sameRE\t" << same_re << '\n';
  statout << "invalid_selflig\t" << self_lig << '\n';
  statout << "intrachromosomal\t" << intra << '\n';
  statout << "interchromosomal\t" << inter << '\n';
  statout.close();
}
