#ifndef GLOCALALIGNMENT_H_
#define GLOCALALIGNMENT_H_

#include "util/bio_util.hpp"

namespace global_alignment {

class GlobalAlignment {
 public:
  GlobalAlignment();
  ~GlobalAlignment();
  int RunGlobalAlignment(const char * U, const char * V, const uint32_t & n,
                         const uint32_t & m);

 private:

  vector<vector<int> > s;
  vector<vector<int> > g;
  vector<vector<int> > h;

  int gapopen;
  int gapextension;

};

}  // namespace global_alignment

#endif /* GLOCALALIGNMENT_H_ */
