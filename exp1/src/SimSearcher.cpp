#include "SimSearcher.h"

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
u levenshteinDistance(const char *a, const char *b, int th)
{
  static int dp[2][N];
  int m = (int)strlen(a), n = (int)strlen(b);
  //while (m > 0 && n > 0 && a[0] == b[0])
  //a++, b++, m--, n--;
  //while (m > 0 && n > 0 && a[m-1] == b[n-1])
  //m--, n--;
  if (m < n)
    swap(a, b), swap(m, n);
  if (m-n > th)
    return th;
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

void tokenize(const char *s, vector<string> &res)
{
  res.clear();
  for (u i = 0; *s; ) {
    for (i = 0; s[i] && s[i] != ' '; i++);
    res.PB(string(s, i));
    s += i;
    if (*s == ' ') s++;
  }
}

// QGram {{{1

QGram::QGram(u q) : q(q), pp(ipow(BASE, q-1)), allo(IdAllocator<NB>(q, CAP))
{}

u QGram::bf_overlap(const char *a, u n, const char *b, u m)
{
  flip = ! flip;
  u h = 0, pp = 1;
  REP(j, q-1) {
    h = h * BASE + a[j];
    pp *= BASE;
  }
  FOR(j, q-1, n) {
    h = h * BASE + a[j];
    int id = allo.get(h, a+(j-q+1));
    if (id != -1)
      freq(id)++;
    h -= pp * a[j-q+1];
  }

  u overlap = 0;
  h = 0;
  REP(j, q-1)
    h = h * BASE + b[j];
  FOR(j, q-1, m) {
    h = h * BASE + b[j];
    int id = allo.get(h, b+(j-q+1));
    if (id != -1 && freq(id) > 0) {
      freq(id)--;
      overlap++;
    }
    h -= pp * b[j-q+1];
  }
  return overlap;
}

u QGram::bf_overlap2(const char *b, u m)
{
  u overlap = 0, h = 0;
  REP(j, q-1)
    h = h * BASE + b[j];
  FOR(j, q-1, m) {
    h = h * BASE + b[j];
    int id = allo.get(h, b+(j-q+1));
    if (id != -1 && freq(id)-- > 0)
      overlap++;
    h -= pp * b[j-q+1];
  }

  h = 0;
  REP(j, q-1)
    h = h * BASE + b[j];
  FOR(j, q-1, m) {
    h = h * BASE + b[j];
    int id = allo.get(h, b+(j-q+1));
    if (id != -1)
      freq(id)++;
    h -= pp * b[j-q+1];
  }
  return overlap;
}

// BruteForce {{{1

BruteForce::BruteForce(const char *filename, u q) : QGram(q)
{
  char line[N];
  FILE *fp = fopen(filename, "r");
  strs.clear();
  while (fgets(line, sizeof line, fp)) {
    u n = strlen(line);
    line[n-1] = '\0';
    strs.PB(line);
    u h = 0;
    REP(j, q-1)
      h = h * BASE + line[j];
    FOR(j, q-1, n) {
      h = h * BASE + line[j];
      allo.find(h, &strs.back()[j-q+1]);
      h -= pp * line[j-q+1];
    }
  }
  fclose(fp);
}

int BruteForce::searchJaccard(const char *query, ft threshold, vector<pair<u, ft> > &res)
{
  res.clear();

  u n = strlen(query);
  REP(i, strs.size()) {
    u m = strs[i].size();
    if (n >= q && m >= q) {
      u m = strs[i].size(), overlap = bf_overlap(query, n, strs[i].c_str(), m);
      ft th = ft(overlap) / (n-q+1 + m-q+1 - overlap);
      if (th >= threshold)
        res.PB(MP(i, th));
    } else if (threshold == .0)
      res.PB(MP(i, .0));
  }
  return 0;
}

int BruteForce::searchED(const char *query, u threshold, VPUU &res)
{
  res.clear();
  REP(i, strs.size()) {
    u d = levenshteinDistance(query, strs[i].c_str(), threshold+1);
    if (d <= threshold)
      res.push_back(MP(i, d));
  }
  return 0;
}

// Tournament {{{1

template<typename Derived>
Tournament<Derived>::Tournament(const char *filename, u q) : QGram(q)
{
  char line[N];
  FILE *fp = fopen(filename, "r");
  VVI idx;
  gs_min = N;
  while (fgets(line, sizeof line, fp)) {
    u n = strlen(line)-1;
    line[n] = '\0';
    strs.PB(line);
    u h = 0;
    REP(i, q-1)
      h = h * BASE + line[i];
    FOR(i, q-1, n) {
      h = h * BASE + line[i];
      u id = allo.find(h, &strs.back()[i-q+1]);
      mapping[id].PB(strs.size()-1);
      h -= pp * line[i-q+1];
    }
    if (n >= q)
      gs_min = min(gs_min, n-q+1);
  }
  for (auto &i: mapping)
    i.second.PB(INT_MAX);
  frequency.assign(strs.size(), 0);
  tick.assign(strs.size(), false);
  flip = false;
  fclose(fp);
}

struct Compare
{
  bool operator()(const pair<u*,u*> &a, const pair<u*,u*> &b) const {
    return *a.fi < *b.fi;
  }
};

// @param overlap > 0
// @return len >= q
template<typename Derived>
void Tournament<Derived>::search(const char *query, u n, u overlap, VPUU &res)
{
  u h = 0;
  Heap<pair<u*,u*>, Compare> L;
  REP(i, q-1)
    h = h * BASE + query[i];
  FOR(i, q-1, n) {
    h = h * BASE + query[i];
    int id = allo.get(h, query+(i-q+1));
    if (id != -1 && mapping.count(id))
      L.PB(MP(&mapping[id][0], &mapping[id][0]+mapping[id].size()));
    h -= pp * query[i-q+1];
  }

  if (L.size() < overlap)
    return;

  flip = ! flip;
  //VU freq(strs.size(), 0);
  L.make();
  for (auto x: L.a)
    freq(*x.fi)++;
  while (*L[0].fi < INT_MAX) {
    if (freq(*L[0].fi) >= overlap)
      res.PB(MP(*L[0].fi, freq(*L[0].fi)));
    int t = *L[0].fi;
    while (*L[0].fi == t) {
      int old = *L[0].fi++;
      if (old < *L[0].fi) {
        freq(old)--;
        if (*L[0].fi < INT_MAX)
          freq(*L[0].fi)++;
        L.down(0);
      }
    }
  }
}

template<typename Derived>
int Tournament<Derived>::searchJaccard(const char *query, ft threshold, vector<pair<u, ft> > &res)
{
  u n = strlen(query), gq = n >= q ? n-q+1 : 0,
    overlap = (u)ceil(max(threshold*gq, (gq + gs_min)*threshold/(1+threshold)) - eps);
  VPUU r;
  res.clear();
  if (overlap > 0) {
    static_cast<Derived*>(this)->search(query, n, overlap, r);
    filter_jaccard(query, n, threshold, r, res);
  } else {
    flip = ! flip;
    u h = 0;
    REP(j, q-1)
      h = h * BASE + query[j];
    FOR(j, q-1, n) {
      h = h * BASE + query[j];
      int id = allo.get(h, query+(j-q+1));
      if (id != -1)
        freq(id)++;
      h -= pp * query[j-q+1];
    }

    REP(i, strs.size())
      if (strs.size() >= q) {
        //u m = strs[i].size(), overlap = bf_overlap2(strs[i].c_str(), m);
        u m = strs[i].size(), overlap = bf_overlap(query, n, strs[i].c_str(), m);
        ft th = ft(overlap) / (n-q+1 + m-q+1 - overlap);
        if (th >= threshold)
          res.PB(MP(i, th));
      } else if (threshold == .0)
        res.PB(MP(i, .0));
  }
  return 0;
}

template<typename Derived>
int Tournament<Derived>::searchED(const char *query, u threshold, VPUU &res)
{
  u n = strlen(query), gq = n-q+1,
    overlap = gq >= threshold*q ? gq-threshold*q : 0;
  VPUU r;
  res.clear();
  if (overlap > 0) {
    static_cast<Derived*>(this)->search(query, n, overlap, r);
    filter_edit_distance(query, n, threshold, r, res);
  } else
    REP(i, strs.size()) {
      u d = levenshteinDistance(query, strs[i].c_str(), threshold+1);
      if (d <= threshold)
        res.PB(MP(i, d));
    }
  return 0;
}

template<typename Derived>
void Tournament<Derived>::filter_jaccard(const char *query, u n, ft threshold, const VPUU &cand, vector<pair<u, ft>> &res)
{
  for (auto i: cand) {
    ft th = ft(i.se) / (n-q+1 + strs[i.fi].size()-q+1 - i.se);
    if (th >= threshold)
      res.PB(MP(i.fi, th));
  }
}

template<typename Derived>
void Tournament<Derived>::filter_edit_distance(const char *query, u n, u threshold, const VPUU &cand, VPUU &res)
{
  for (auto i: cand) {
    u d = levenshteinDistance(query, strs[i.fi].c_str(), threshold+1);
    if (d <= threshold)
      res.PB(MP(i.fi, d));
  }
}

// MergeSkip

MergeSkip::MergeSkip(const char *filename, u q) : Tournament<MergeSkip>(filename, q)
{}

// @param overlap > 0
// @return len >= q
void MergeSkip::search(const char *query, u n, u overlap, VPUU &res)
{
  u h = 0;
  Heap<pair<u*,u*>, Compare> L;
  REP(i, q-1)
    h = h * BASE + query[i];
  FOR(i, q-1, n) {
    h = h * BASE + query[i];
    int id = allo.get(h, query+(i-q+1));
    if (id != -1 && mapping.count(id))
      L.PB(MP(&mapping[id][0], &mapping[id][0]+mapping[id].size()));
    h -= pp * query[i-q+1];
  }

  if (L.size() < overlap)
    return;

  flip = ! flip;
  //VU freq(strs.size(), 0);
  L.make();
  for (auto x: L.a)
    freq(*x.fi)++;
  while (*L[0].fi < INT_MAX) {
    if (freq(*L[0].fi) >= overlap)
      res.PB(MP(*L[0].fi, freq(*L[0].fi)));

//    FOR(i, L.size()-(overlap-1), L.size()) {
    //  L.emplace();
    //}

    int t = *L[0].fi;
    while (*L[0].fi == t) {
      int old = *L[0].fi++;
      if (old < *L[0].fi) {
        freq(old)--;
        if (*L[0].fi < INT_MAX) {
          freq(*L[0].fi)++;
          L.down(0);
        } else
          L.pop();
      }
    }
    if (L.size() < overlap) break;

    REP(i, overlap-1) {
      if (*L[0].fi < INT_MAX)
        freq(*L[0].fi)--;
      L.pop_to_back();
    }

    t = *L[0].fi;
    FOR(i, L.size()-(overlap-1), L.size()) {
      u *l = L[i].fi, *h = L[i].se;
      while (l < h) {
        u *x = l + ((h-l)>>1);
        if (*x < t) l = x+1;
        else h = x;
      }
      L[i].fi = l;
      //while (*L[i].fi < t)
        //L[i].fi++;
      if (*L[i].fi < INT_MAX)
        freq(*L[i].fi)++;
      L.emplace();
    }
  }
}

// DivideSkip {{{1
DivideSkip::DivideSkip(const char *filename, u q) : Tournament<DivideSkip>(filename, q)
{}

void DivideSkip::search(const char *query, u n, u overlap, VPUU &res)
{
  u h = 0;
  u maxlen = 0;
  vector<pair<u*,u*> > L;
  Heap<pair<u*,u*>, Compare> H;
  REP(i, q-1)
    h = h * BASE + query[i];
  FOR(i, q-1, n) {
    h = h * BASE + query[i];
    int id = allo.get(h, query+(i-q+1));
    if (id != -1 && mapping.count(id)) {
      L.PB(MP(&mapping[id][0], &mapping[id][0]+mapping[id].size()));
      maxlen = max(maxlen, (u)mapping[id].size());
    }
    h -= pp * query[i-q+1];
  }

  if (L.size() < overlap)
    return;

  const ft MU = 0.007;

  u nlong = min(u(overlap/(MU*log(ft(maxlen))+1)), overlap-1);
  u nshort = L.size() - nlong;
  assert(overlap > nlong);

  // sort invertex indices by length
  sort(ALL(L), [](const pair<u*,u*> &a, const pair<u*,u*> &b) {
       return a.se-a.fi < b.se-b.fi;
       });
  REP(i, nshort)
    H.PB(L[i]);
  H.make();

  u min_freq_in_short = overlap-nlong;
  flip = ! flip;
  for (auto x: H.a)
    freq(*x.fi)++;
  while (*H[0].fi < INT_MAX) {
    u cnt = freq(*H[0].fi);
    if (cnt >= min_freq_in_short) {
      FOR(i, nshort, L.size()) {
        u *pos = lower_bound(L[i].fi, L[i].se, *H[0].fi);
        if (*pos == *H[0].fi)
          cnt++;
      }
      if (cnt >= overlap)
        res.PB(MP(*H[0].fi, cnt));
    }

    int t = *H[0].fi;
    while (*H[0].fi == t) {
      int old = *H[0].fi++;
      if (old < *H[0].fi) {
        freq(old)--;
        if (*H[0].fi < INT_MAX)
          freq(*H[0].fi)++;
        H.down(0);
      }
    }
    //if (H.size() < overlap) break;

    REP(i, min_freq_in_short-1) {
      if (*H[0].fi < INT_MAX)
        freq(*H[0].fi)--;
      H.pop_to_back();
    }

    t = *H[0].fi;
    FOR(i, H.size()-(min_freq_in_short-1), H.size()) {
      H[i].fi = lower_bound(H[i].fi, H[i].se, t);
      //while (*L[i].fi < t)
        //L[i].fi++;
      if (*H[i].fi < INT_MAX)
        freq(*H[i].fi)++;
      H.emplace();
    }
  }
}

// SimSearcher {{{1

SimSearcher::SimSearcher()
{
}

SimSearcher::~SimSearcher()
{
  if (impl)
    delete impl;
}

int SimSearcher::createIndex(const char *filename, u q)
{
  if (impl)
    delete impl;
#ifdef BRUTE_FORCE
  impl = new BruteForce(filename, q);
#elif defined(DIVIDE_SKIP)
  impl = new DivideSkip(filename, q);
#elif defined(MERGE_SKIP)
  impl = new MergeSkip(filename, q);
#else
  impl = new Tournament<Tournament>(filename, q);
#endif
  return 0;
}

int SimSearcher::searchJaccard(const char *query, double threshold, vector<pair<u, ft> > &result)
{
  return impl->searchJaccard(query, threshold, result);
}

int SimSearcher::searchED(const char *query, u threshold, VPUU &result)
{
  return impl->searchED(query, threshold, result);
}
