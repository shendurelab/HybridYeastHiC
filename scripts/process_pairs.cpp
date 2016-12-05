/* process_pairs.cpp
 * takes mapped reads, deduplicates them, and assigns to restriction fragments, and identifies valid ligations
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

// hash function for deduplicating read pairs
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

  if (argc < 8) {
    cout << "USAGE: ./process_pairs RESTRICTION_SITE GENOME_DIGEST MAPQ SHOTGUN_THR READ1_FILE READ2_FILE READPAIR_OUT STAT_OUT\n";
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
  
  // arg 3: mapq filter
  int mapq_filt;
  stringstream ss2(argv[3]);
  ss2 >> mapq_filt;

  // arg 4: shotgun threshold
  int shot_filt;
  stringstream ss3(argv[4]);
  ss3 >> shot_filt;

  // arg 5: read 1 file
  ifstream read1file(argv[5]);
  // arg 6: read 2 file
  ifstream read2file(argv[6]);
  // arg 7: output file for read pairs
  ofstream readout(argv[7]);
  // arg 8: output file for stats
  ofstream statout(argv[8]);

  // alignment score filter
  int as_filt = 10;

  string line;
  string chr;
  long pos1, pos2;

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

  string rname1, rname2;
  int dir1, dir2;
  string chr1, chr2;
  int mapq1, mapq2;
  string cig1, cig2;
 
  // data structure to store all read pairs, to check for duplicates
  unordered_set<readpair> readpairs;

  // counters for read pair categories
  long fail_mapq, shot, dup, intra, inter, same_re, same_dir, self_lig = 0;

  int as1, xs1, as2, xs2;
  size_t aspos, xspos;

  while (getline(read1file,line)) {
    // parse line from read 1 file
    stringstream ss(line);
    ss >> rname1;
    ss >> dir1;
    ss >> chr1;
    ss >> pos1;
    ss >> mapq1;
    ss >> cig1;

    as1 = 0;
    aspos = line.find("AS:i:");
    if (aspos != string::npos) {
      stringstream ss(line.substr(aspos+5));
      ss >> as1;
    }
    xs1 = -40;
    xspos = line.find("XS:i:");
    if (xspos != string::npos) {
      stringstream ss(line.substr(xspos+5));
      ss >> xs1;
    }
    // parse line from read 2 file
    getline(read2file,line);
    stringstream ss2(line);
    ss2 >> rname2;
    ss2 >> dir2;
    ss2 >> chr2;
    ss2 >> pos2;
    ss2 >> mapq2;
    ss2 >> cig2;
    
    as2 = 0;
    aspos = line.find("AS:i:");
    if (aspos != string::npos) {
      stringstream ss2(line.substr(aspos+5));
      ss2 >> as2;
    }
    xs2 = -40;
    xspos = line.find("XS:i:");
    if (xspos != string::npos) {
      stringstream ss2(line.substr(xspos+5));
      ss2 >> xs2;
    }
  
    // if mapq < threshold or as - xs < threshold for either read, ignore read pair
    if ((mapq1 < mapq_filt) || (as1-xs1 < as_filt) || (mapq2 < mapq_filt) || (as2-xs2 < as_filt)) {
 //   if ((mapq1 < mapq_filt) || (mapq2 < mapq_filt)) {
      fail_mapq++;
      continue;
    }

    // find end positions of read alignments
    long pos1end = pos1 - 1 + readend(cig1);
    long pos2end = pos2 - 1 + readend(cig2);

    // determine end of read alignment further from end of fragment that was ligated
    long map1end = (dir1 == 0) ? pos1 : pos1end;
    long map2end = (dir2 == 0) ? pos2 : pos2end;

    // if not a duplicate, insert read pair into data structure
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
    long fraglen = ((dir1 == 0) ? (re1end - pos1 + 1) : (pos1end - re1start + 1)) + ((dir2 == 0) ? (re2end - pos2 + 1) : (pos2end - re2start + 1));
    

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
