# include <iostream>
# include <cmath>
# include <cstdlib>
# include <ctime>
# include <limits>
# include <vector>
# include "count_min_sketch.hpp"
using namespace std;

std::vector<CountMinSketch> feature_counters;
int current_counter_index;
/**
   Class definition for CountMinSketch.
   public operations:
   // overloaded updates
   void update(uint64_t item, int c);
   void update(const char *item, int c);
   // overloaded estimates
   uint32_t estimate(uint32_t item);
   uint32_t estimate(const char *item);
**/


// CountMinSketch constructor
// ep -> error 0.01 < ep < 1 (the smaller the better)
// gamma -> probability for error (the smaller the better) 0 < gamm < 1
CountMinSketch::CountMinSketch(float ep, float gamm) {
  if (!(0.009 <= ep && ep < 1)) {
    cout << "eps must be in this range: [0.01, 1)" << endl;
    exit(EXIT_FAILURE);
  } else if (!(0 < gamm && gamm < 1)) {
    cout << "gamma must be in this range: (0,1)" << endl;
    exit(EXIT_FAILURE);
  }
  eps = ep;
  gamma = gamm;
  w = ceil(exp(1)/eps);
  d = ceil(log(1/gamma));
  total = 0;
  // initialize counter array of arrays, C
  C = new uint32_t *[d];
  uint32_t i, j;
  for (i = 0; i < d; i++) {
    C[i] = new uint32_t[w];
    for (j = 0; j < w; j++) {
      C[i][j] = 0;
    }
  }
  // initialize d pairwise independent hashes
  srand(time(NULL));
  hashes = new uint64_t* [d];
  for (i = 0; i < d; i++) {
    hashes[i] = new uint64_t[2];
    genajbj(hashes, i);
  }
}

// CountMinSkectch destructor
CountMinSketch::~CountMinSketch() {
  // free array of counters, C
  unsigned int i;
  for (i = 0; i < d; i++) {
    delete[] C[i];
  }
  delete[] C;

  // free array of hash values
  for (i = 0; i < d; i++) {
    delete[] hashes[i];
  }
  delete[] hashes;
}

// CountMinSketch totalcount returns the
// total count of all items in the sketch
uint32_t CountMinSketch::totalcount() {
  return total;
}

// countMinSketch update item count (int)
void CountMinSketch::update(uint64_t item, int c) {
  total = total + c;
  uint32_t hashval = 0;
  uint32_t count = estimate(item);
  for (uint32_t j = 0; j < d; j++) {
    hashval = (hashes[j][0]*item+hashes[j][1])%LONG_PRIME%w;
    if (C[j][hashval] < count + c)  // conservative update
        C[j][hashval] = count + c;
  }
}

// countMinSketch update item count (string)
void CountMinSketch::update(const char *str, int c) {
  uint64_t hashval = hashstr(str);
  update(hashval, c);
}

// CountMinSketch estimate item count (int)
uint32_t CountMinSketch::estimate(uint64_t item) {
  uint32_t minval = numeric_limits<uint32_t>::max();
  uint64_t hashval = 0;
  for (uint32_t j = 0; j < d; j++) {
    hashval = (hashes[j][0]*item+hashes[j][1])%LONG_PRIME%w;
    minval = std::min(minval, C[j][hashval]);
  }
  return minval;
}

// CountMinSketch estimate item count (string)
uint32_t CountMinSketch::estimate(const char *str) {
  uint64_t hashval = hashstr(str);
  return estimate(hashval);
}

// generates aj,bj from field Z_p for use in hashing
void CountMinSketch::genajbj(uint64_t** hashes, int i) {
  uint64_t s = uint64_t(float(rand())*LONG_PRIME/RAND_MAX + 1);
#ifdef DEBUG
  cout << "set hashes[" << i << "][0] to " << s << endl;
#endif
  hashes[i][0] = s;
  s = uint64_t(float(rand())*LONG_PRIME/RAND_MAX + 1);
#ifdef DEBUG
  cout << "set hashes[" << i << "][1] to " << s << endl;
#endif
  hashes[i][1] = s;
}

// generates a hash value for a sting
// same as djb2 hash function
uint64_t CountMinSketch::hashstr(const char *str) {
  uint64_t hash = 5381;
  int c;
  while (c = *str++) {
    hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
  }
  return hash;
}

void CountMinSketch::reset() {
  for (uint32_t i = 0; i < d; i++) {
    for (uint32_t j = 0; j < w; j++) {
      C[i][j] = 0;
    }
  }
}

