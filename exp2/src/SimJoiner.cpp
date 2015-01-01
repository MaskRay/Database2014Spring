#include "SimJoiner.h"

//#define BRUTE_FORCE
#define DIVIDE_SKIP
//#define FILTER_VERIFY

// static members {{{1

int CharHistogram::charmap[256];

// common {{{1

u ipow(u a, u n)
{
  u s = 1;
  for (; n; n >>= 1) {
    if (n & 1)
      s *= a;
    a *= a;
  }
  return s;
}

template<typename UINT>
int levenshtein_bit_vector(const char *a, int m, const char *b, int n)
{
  int r = m;
  UINT vp = (UINT)-1;
  UINT vn = 0;
  UINT hp = 0;
  UINT hn = 0;
  int bl = min<int>(n, sizeof(UINT) * 8);
  for (int j = 0; j < m; j++) {
    UINT pm = 0;
    for (int i = bl-1; i >= 0; i--) {
      pm <<= 1;
      pm |= b[i] == a[j];
    }
    UINT d = (((pm & vp) + vp)) ^ vp | pm | vn;
    hp = (vn | ~(d | vp));
    UINT hpw = (hp << 1) | 1;
    hn = d & vp;
    vp = (hn << 1) | ~(d | hpw);
    vn = d & hpw;
  }
  for (int i = 0; i < bl; i++) {
    r += (vp & 1) - (vn & 1);
    vp >>= 1;
    vn >>= 1;
  }
  return r;
}

// @param th > 0
int levenshteinDistance(const char *a, int m, const char *b, int n, int th)
{
  static int dp[2][N];
  //while (m > 0 && n > 0 && a[0] == b[0])
  //a++, b++, m--, n--;
  //while (m > 0 && n > 0 && a[m-1] == b[n-1])
  //m--, n--;
  if (m < n)
    swap(a, b), swap(m, n);
  if (m-n > th)
    return th+1;
  if (n <= 64)
    return levenshtein_bit_vector<uint64_t>(a, m, b, n);
  if (n <= 128)
    return levenshtein_bit_vector<unsigned __int128>(a, m, b, n);
#if 0
  REP(j, n+1)
    dp[0][j] = j;
  REP1(i, m) {
    dp[i&1][0] = i;
    REP1(j, n)
      dp[i&1][j] = a[i-1] == b[j-1] ? dp[i-1&1][j-1] : min(min(dp[i-1&1][j], dp[i&1][j-1]), dp[i-1&1][j-1]) + 1;
  }
#else
  REP(j, min(n, th) + 1)
    dp[0][j] = j;
  REP1(i, m) {
    int l = max(i-th, 0), h = min(i+th, n);
    if (! l)
      dp[i&1][l++] = i;
    FOR(j, l, h+1)
      dp[i&1][j] = a[i-1] == b[j-1] ? dp[i-1&1][j-1] : min(min(j-(i-1) > th ? th : dp[i-1&1][j], i-(j-1) > th ? th : dp[i&1][j-1]), dp[i-1&1][j-1]) + 1;
  }
#endif
  return dp[m&1][n];
}

JaccardJoinResult mkJaccard(int id1, int id2, ft s)
{
  JaccardJoinResult res;
  res.id1 = id1;
  res.id2 = id2;
  res.s = s;
  return res;
}

EDJoinResult mkED(u id1, u id2, u s)
{
  EDJoinResult res;
  res.id1 = id1;
  res.id2 = id2;
  res.s = s;
  return res;
}

static int intersect(Counter &c, const int *a, const int *aa, const int *b, const int *bb)
{
  int r = 0;
  c.reset();
  while (a < aa)
    c[*a++]++;
  for (; b < bb; b++) {
    int &d = c[*b];
    if (d > 0) {
      d--;
      r++;
    }
  }
  return r;
}

bool SimJoiner::cmp_gram(int a, int b)
{
  if (df[a] != df[b]) return df[a] < df[b];
  return a < b;
}

bool SimJoiner::partition(int *x, int m, int w, int lo, int hi, int &lx, int &rx, int &diff)
{
  int saved_hi = hi;
  while (lo < hi) {
    int mi = (lo+hi) >> 1;
    if (cmp_gram(x[mi], w)) lo = mi+1;
    else hi = mi;
  }
  lx = lo;
  if (lo < saved_hi && x[lo] == w) {
    diff = 0;
    rx = m-(lo+1);
  } else {
    diff = 1;
    rx = m-lo;
  }
  return true;
}

int SimJoiner::suffix_filter(int *x, int m, int *y, int n, int hamming_ub, int dep)
{
  if (! m || ! n || dep >= SUFFIX_FILTER_DEPTH) return abs(m-n);
  int mid = (n-1)/2, w = y[mid];
  int diff, lx, rx, ly, ry;
  int o = (hamming_ub-abs(m-n))>>1;
  partition(y, n, w, 0, mid+1, ly, ry, diff);
  if (! partition(x, m, w, max(ly-o-(m<n?abs(m-n):0), 0), min(ly+o+(m<n?0:abs(m-n))+1, m), lx, rx, diff))
    return hamming_ub+1;
  int h = abs(lx-ly)+diff+abs(rx-ry);
  if (h > hamming_ub) return h;
  int lh = suffix_filter(x, lx, y, ly, hamming_ub-abs(rx-ry)-diff, dep+1);
  h = lh + diff + abs(rx-ry);
  if (h > hamming_ub) return h;
  return lh + diff + suffix_filter(x+m-rx, rx, y+n-ry, ry, hamming_ub-lh-diff, dep+1);
}

void SimJoiner::ppjoin(VS R, VS S, QueryType qt, ft th, int delta, vector<JaccardJoinResult> *res_jaccard, vector<EDJoinResult> *res_ed)
{
  s2g = new IdAllocator<NB>(q, CAP);
  prefg2sid.clear();
  gss.clear();
  df.clear();
  lenL.clear();

  // extract q-grams
  REP(i, S.size()) {
    string &y = S[i];
    gss.eb();
    VI &gs = gss.back();
    if (y.size() >= lenL.size())
      lenL.resize(y.size()+1);
    lenL[y.size()].pb(i);
    B_RABIN_KARP(y.c_str(), y.size(), j)
      int ngid = s2g->size, gid = s2g->find(h, &y[j]);
      if (gid == ngid) {
        df.pb(0);
        prefg2sid.eb();
      }
      df[gid]++;
      gs.pb(gid);
    E_RABIN_KARP(y.c_str(), y.size(), j)
  }

  // initialize frequency histogram
  freqs.assign(s2g->size, 0);
  tick.assign(s2g->size, false);
  flip = false;
  Counter counter(s2g->size);
  Counter A(S.size());

  // construct inverted lists with heuristic
  REP(i, gss.size()) {
    VI &gs = gss[i];
    int ly = gs.size();
    int prefy = min(ly, qt == JACCARD ? int(ly-ceil(th*ly-EPS)+1) : q*delta+1);
    sort(ALL(gs), [&](int a, int b) {
      if (df[a] != df[b]) return df[a] < df[b];
      return a < b;
    });
    REP(j, prefy)
      prefg2sid[gs[j]].pb(mp(i, j));
  }

  // probe
  REP(i, R.size()) {
    VI gs;
    string &x = R[i];
    B_RABIN_KARP(x.c_str(), x.size(), j)
      int gid = s2g->get(h, &x[j]);
      if (gid != -1)
        gs.pb(gid);
    E_RABIN_KARP(x.c_str(), x.size(), j)
    sort(ALL(gs), [&](int a, int b) {
      if (df[a] != df[b]) return df[a] < df[b];
      return a < b;
    });

    int lx = gs.size();
    int prefx = min(lx, qt == JACCARD ? int(lx-ceil(th*lx-EPS)+1) : q*delta+1);
    A.reset();
    REP(j, prefx) {
      int w = gs[j];
      if (j > 0 && w == gs[j-1]) continue; // one q-gram in y can only be mapped once
      REP(k, prefg2sid[w].size()) {
        int sid = prefg2sid[w][k].fi, pref_pos = prefg2sid[w][k].se;
        int ly = gss[sid].size();
        // multiple q-grams may be mapped by one in x
        if (lx >= ly*th && ly >= lx*th) {
          int alpha = qt == JACCARD ? int(ceil((lx+ly)*th/(1.0+th)-EPS)) : max(lx, ly)-delta*q,
              ub = 1+min(lx-j-1, ly-pref_pos-1);
          int &v = A[sid];
          if (v >= 0) {
            if (v+ub >= alpha) {
              VI &gs2 = gss[sid];
              int hamming_ub = max(lx+ly-2*int(ceil((lx+ly)*th/(1.0+th)-EPS))-(j+pref_pos), 0);
              if (v == 0 && lx+ly > 1000 && suffix_filter(&gs[0]+j+1, gs.size()-j-1, &gs2[0]+pref_pos+1, gs2.size()-pref_pos-1, hamming_ub, 0) > hamming_ub)
                v = -1;
              else
                v++;
            } else
              v = -1;
          }
        }
      }
    }

    sort(A.b, A.b+A.num);
    REP(j, A.num) {
      int sid = A.b[j], o = A.c[sid];
      if (o >= 0) {
        VI &gs2 = gss[sid];
        int ly = gs2.size();
        int prefy = min(ly, qt == JACCARD ? int(ly-ceil(th*ly-EPS)+1) : q*delta+1);
        int wx = gs[prefx-1], wy = gs2[prefy-1];
        int alpha = qt == JACCARD ? int(ceil((lx+ly)*th/(1.0+th)-EPS)) : max(lx, ly)-delta*q;
        int ub = cmp_gram(wx, wy) ? min(o+lx-prefx, ly) : min(o+ly-prefy, lx);
        if (ub >= alpha) {
          o = intersect(counter, &gs[0], &*gs.end(), &gs2[0], &*gs2.end());
          if (o >= alpha) {
            if (qt == JACCARD) {
              ft jaccard = ft(o) / (gs.size()+gs2.size()-o);
              if (jaccard >= th)
                res_jaccard->pb(mkJaccard(i, sid, jaccard));
            } else {
              int ed = levenshteinDistance(x.c_str(), x.size(), S[sid].c_str(), S[sid].size(), delta);
              if (ed <= delta)
                res_ed->pb(mkED(i, sid, ed));
            }
          }
        }
      }
    }

    if (qt == ED) {
      int ub = min(min(int(x.size()+delta), q*delta)+1, int(lenL.size()));
      FOR(l, max(int(x.size()-delta), 0), ub)
        for (auto sid: lenL[l]) {
          int ed = levenshteinDistance(x.c_str(), x.size(), S[sid].c_str(), S[sid].size(), delta);
          if (ed <= delta)
            res_ed->pb(mkED(i, sid, ed));
        }
    }
  }

  if (qt == ED) {
    sort(ALL(*res_ed), [](const EDJoinResult &a, const EDJoinResult &b) {
         if (a.id1 != b.id1) return a.id1 < b.id1;
         return a.id2 < b.id2;
         });
    res_ed->erase(unique(ALL(*res_ed), [](const EDJoinResult &a, const EDJoinResult &b) {
                         return a.id1 == b.id1 && a.id2 == b.id2;
                         }), res_ed->end());
  }
  delete s2g;
}

static pair<VS, VS> read_file(const char *file1, const char *file2)
{
  char buf[N];
  pair<VS, VS> r;
  FILE *h = fopen(file1, "r");
  while (fgets(buf, sizeof buf, h)) {
    buf[strlen(buf)-1] = '\0';
    r.fi.pb(buf);
  }
  fclose(h);
  h = fopen(file2, "r");
  while (fgets(buf, sizeof buf, h)) {
    buf[strlen(buf)-1] = '\0';
    r.se.pb(buf);
  }
  fclose(h);
  return r;
}

int SimJoiner::joinJaccard(const char *file1, const char *file2, unsigned q, double threshold, vector<JaccardJoinResult> &res)
{
  this->q = q;
  q_pow_p1 = ipow(BASE, q-1);
  pair<VS, VS> RS = read_file(file1, file2);
  ppjoin(RS.fi, RS.se, JACCARD, threshold, -1, &res, nullptr);
  return SUCCESS;
}

int SimJoiner::joinED(const char *file1, const char *file2, unsigned q, unsigned threshold, vector<EDJoinResult> &res)
{
  this->q = q;
  q_pow_p1 = ipow(BASE, q-1);
  pair<VS, VS> RS = read_file(file1, file2);
  ppjoin(RS.fi, RS.se, ED, -1.0, threshold, nullptr, &res);
  return SUCCESS;
}
