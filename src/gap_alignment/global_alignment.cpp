#include "global_alignment.hpp"
#include "option.hpp"

#include <vector>
#include <limits.h>
#include "bio_util.hpp"

namespace global_alignment {
#define MAX_ALIGNMENT_LEN 6000

GlobalAlignment::GlobalAlignment() {
  Option::GetOption("-gopen", gapopen, -11);
  Option::GetOption("-gext", gapextension, -1);

  s.resize(MAX_ALIGNMENT_LEN);
  g.resize(MAX_ALIGNMENT_LEN);
  h.resize(MAX_ALIGNMENT_LEN);
  for (uint32_t i = 0; i < MAX_ALIGNMENT_LEN; i++) {
    s[i].resize(MAX_ALIGNMENT_LEN);
    g[i].resize(MAX_ALIGNMENT_LEN);
    h[i].resize(MAX_ALIGNMENT_LEN);
  }
}

GlobalAlignment::~GlobalAlignment() {

}

int GlobalAlignment::RunGlobalAlignment(const char * U, const char * V,
                                        const uint32_t & n,
                                        const uint32_t & m) {
#ifdef TESTCODE
  for (uint32_t i = 0; i < n; i++) {
    cout << U[i];
  }
  cout << endl;
  for (uint32_t i = 0; i < m; i++) {
    cout << V[i];
  }
  cout << endl;
  //start_t_ = clock();
#endif
  if (n >= MAX_ALIGNMENT_LEN || m >= MAX_ALIGNMENT_LEN) {
    ERROR_INFO("The length of the two sequences is too long.");
    return 0;
  }

  s[0][0] = 0;
  g[0][0] = INT_MIN + 10000;
  h[0][0] = INT_MIN + 10000;
  for (uint32_t i = 1; i <= m; i++) {
    s[0][i] = gapopen_ + i * gapextension;
    g[0][i] = INT_MIN + 10000;
    h[0][i] = INT_MIN + 10000;
  }
  for (uint32_t i = 1; i <= n; i++) {
    s[i][0] = gapopen_ + i * gapextension;
    g[i][0] = INT_MIN + 10000;
    h[i][0] = INT_MIN + 10000;
  }

  int stmp;
  for (uint32_t i = 1; i <= n; i++) {
    for (uint32_t j = 1; j <= m; j++) {
      g[i][j] = max(s[i][j - 1] + gapopen_ + gapextension,
                    g[i][j - 1] + gapextension);
      h[i][j] = max(s[i - 1][j] + gapopen_ + gapextension,
                    h[i - 1][j] + gapextension);

      stmp = s[i - 1][j - 1]
          + BLOSUM62[base[U[i - 1] - 'A']][base[V[j - 1] - 'A']];
      s[i][j] = max(max(g[i][j], h[i][j]), stmp);
    }
  }

  return s[n][m];
}

}  // namespace local_alignment
