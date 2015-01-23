#ifndef LOCALALIGNMENT_H_
#define LOCALALIGNMENT_H_

#include "sdk.hpp"
#include "bio_util.hpp"
#include "evalue.hpp"
#include <time.h>

#include <vector>
using std::vector;

namespace local_alignment {

//#define I (i - 1)
//#define J (j - 1)
//#define P (p - 1)
//#define Q (q - 1)

#define DIAG ('a')
#define UP ('b')
#define LEFT ('c')
#define STOPorSTART ('d')

class LocalAlignment {
 public:
  LocalAlignment(const Evalue* evalue, const uint32_t& max_rows,
                 const uint32_t& max_cols);
  ~LocalAlignment();

  bool RunLocalAlignment(const char * U, const char * V, const uint32_t & Ul,
                         const uint32_t & Vl, M8Results& res);
 private:
  char * rU;
  char * rV;
  char * midline;
  int rs;

  clock_t start_t;
  uint64_t sum_time;

  vector<vector<int> > s;
  vector<vector<char> > l;

  vector<vector<int> > g;
  vector<vector<int> > h;

  uint32_t n;
  uint32_t m;

  int alignScore;
  int nIdentity;
  double e_value;

  int gapopen;
  int gapextension;

  const Evalue* evalue;

  uint32_t max_rows;
  uint32_t max_cols;

  void DisplayAlignment();
  void stringReverse(char * str, const int & n);
  int MaxOfFour(const int & s1, const int & s2, const int & s3, const int & s4);
  char Direction(const int & s1, const int & s2, const int & s3,
                 const int & s4);
};

}  // namespace local_alignment

#endif /* LOCALALIGNMENT_H_ */
