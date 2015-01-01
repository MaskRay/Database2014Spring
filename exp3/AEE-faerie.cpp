#include "AEE.h"
using namespace std;

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

tuple<int, char *, size_t> read_file(const char *filename)
{
  int fd = open(filename, O_RDONLY);
  if (fd < 0)
    errx(1, "failed to open %s: %s", filename, strerror(errno));
  struct stat s;
  int status = fstat(fd, &s);
  if (status < 0)
    errx(1, "failed to stat %s: %s", filename, strerror(errno));
  char *mapped = (char *)mmap(0, s.st_size, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);
  return tuple<int, char *, size_t>(fd, mapped, s.st_size);
}

void AEE::close_file()
{
  if (get<0>(tup) != -1) {
    munmap(get<1>(tup), get<2>(tup));
    close(get<0>(tup));
    get<0>(tup) = 0;
    get<1>(tup) = nullptr;
    get<2>(tup) = 0;
  }
}

int AEE::createIndex(const char *filename, unsigned q)
{
  ngid = 0;
  maxL = 0;
  this->q = q;
  q_pow_p1 = ipow(BASE, q-1);
  e.clear();
  eL.clear();
  gid2e.clear();
  if (s2gid)
    delete s2gid;
  s2gid = new IdAllocator<NB>(q, CAP);

  close_file();

  tup = read_file(filename);
  char *mapped = get<1>(tup);

  int file_size = int(get<2>(tup)), eid = 0;
  for (int i = 0, j = 0; j < file_size; i = j, eid++) {
    for (; j < file_size && mapped[j] != '\n'; j++);
    B_RABIN_KARP((mapped+i), j-i, k)
      int ngid = s2gid->size;
      int gid = s2gid->find(h, mapped+i+k);
      if (gid == ngid) {
        gid2e.emplace_back();
        gid2e[gid].pb(eid);
      } else {
        if (gid2e[gid].size() && gid2e[gid].back() != eid)
          gid2e[gid].pb(eid);
      }
    E_RABIN_KARP((mapped+i), j-i, k)
    e.eb(mapped+i, j-i);
    eL.pb(j-i);
    maxL = max(maxL, j-i);
    if (j < file_size && mapped[j] == '\n')
      j++;
  }
  for (auto &e: gid2e)
    e.pb(INT_MAX);

  freqs.assign(s2gid->size, 0);
  ticks.assign(s2gid->size, 0);
  flip = false;
	return SUCCESS;
}

struct Cmp
{
  bool operator()(const pair<int *, int> &a, const pair<int *, int> &b) const {
    if (*a.fi != *b.fi) return *a.fi < *b.fi;
    return a.se < b.se;
  }
};

int AEE::overlap(const char *a, int m, const char *b, int n)
{
  flip = ! flip;
  B_RABIN_KARP(a, m, i)
    int id = s2gid->get(h, &a[i]);
    if (id != -1)
      freq(id)++;
  E_RABIN_KARP(a, n, i)

  int o = 0;
  B_RABIN_KARP(b, n, i)
    int id = s2gid->get(h, &b[i]);
    if (id != -1 && freq(id) > 0) {
      freq(id)--;
      o++;
    }
  E_RABIN_KARP(b, m, i)
  return o;
}

#if 0
int AEE::overlap(const char *a, int m, const char *b, int n)
{
  unordered_map<string, int> freq;
  int o = 0;
  REP(i, m-q+1) {
    string s(a+i, q);
    freq[s]++;
  }
  REP(i, n-q+1) {
    string s(b+i, q);
    auto it = freq.find(s);
    if (it != freq.end()) {
      o++;
      if (! --it->second)
        freq.erase(it);
    }
  }
  return o;
}
#endif

void AEE::faerie(const char *doc, ft threshold, bool flag_ed, vector<EDExtractResult> *result_ed, vector<JaccardExtractResult> *result_jaccard)
{
  VI empty(1, INT_MAX);
  Heap<pair<int *, int>, Cmp> L;
  if (flag_ed)
    result_ed->clear();
  else
    result_jaccard->clear();

  // make heap
  int n = int(strlen(doc)), ng = n-q+1, th = int(threshold);
  B_RABIN_KARP(doc, n, i)
    int gid = s2gid->get(h, doc+i);
    if (gid == -1)
      L.eb(&empty[0], i);
    else
      L.eb(&gid2e[gid][0], i);
  E_RABIN_KARP(doc, n, i)
  L.make();

  VI p;
  int eeid = *L[0].fi;
  do {
    int pos = L[0].se, eid = *L[0].fi++;
    L.down(0);
    if (eid == eeid)
      p.pb(pos);
    else {
      int Be, Te, Tl;
      if (flag_ed) {
        Be = max(eL[eeid]-q+1-th, 1);
        Te = eL[eeid]-q+1+th;
        Tl = max(eL[eeid]-q+1-th*q, 1);
      } else {
        Be = max(int(ceil(eL[eeid]-q+1)*threshold), 1);
        Te = int(floor((eL[eeid]-q+1)/threshold));
        Tl = Be;
      }
      if (p.size() >= Tl)
        for (int i = 0; i <= int(p.size())-Tl; ) { // i <= |p|-Tl
          int j = i+Tl-1; // j < |p|
          if (p[j]-p[i]+1 <= Te) {
            int lo = j, hi = min(i+Te-1, int(p.size())-1);
            while (lo < hi) { // [lo, hi]
              int mi = (lo+hi+1)>>1;
              if (p[mi]-p[i]+1 <= Te) lo = mi;
              else hi = mi-1;
            }
            if (lo == hi && p[lo]-p[i]+1 <= Te) {
              int bgn_lo = max(p[j]-Te+1, i ? p[i-1]+1 : 0), end_hi = min(p[i]+Te-1, lo == p.size()-1 ? ng-1 : p[lo+1]-1);
              if (flag_ed) {
                EDExtractResult res;
                res.id = eeid;
                FORC(bgn, bgn_lo, p[i]) {
                  res.pos = bgn;
                  int min_l = max(p[j], bgn+Be-1)-bgn+q, max_l = end_hi-bgn+q;
                  const char *a = e[eeid].c_str(), *b = doc+bgn;

                  static int dp[2][N];
                  int m = eL[eeid];
                  REP(j, min(max_l, th) + 1)
                    dp[0][j] = j;
                  REP1(i, m) {
                    int l = max(i-th, 0), h = min(i+th, max_l);
                    if (! l)
                      dp[i&1][l++] = i;
                    FOR(j, l, h+1)
                      dp[i&1][j] = a[i-1] == b[j-1] ? dp[i-1&1][j-1] : min(min(j-(i-1) > th ? th : dp[i-1&1][j], i-(j-1) > th ? th : dp[i&1][j-1]), dp[i-1&1][j-1]) + 1;
                  }
                  FORC(l, min_l, max_l) {
                    if (abs(l-eL[eeid]) <= th && dp[eL[eeid]&1][l] <= th) {
                      res.len = l;
                      res.sim = dp[eL[eeid]&1][l];
                      result_ed->pb(res);
                    }
                  }
                  //FORC(end, max(p[j], bgn+Be-1), end_hi) {
                    //res.len = end-bgn+q;
                    //res.sim = levenshteinDistance(e[eeid].c_str(), eL[eeid], doc+bgn, res.len, th+1);
                    //if (res.sim <= th)
                      //result_ed->pb(res);
                  //}
                }
              } else {
                JaccardExtractResult res;
                res.id = eeid;
                FORC(bgn, bgn_lo, p[i]) {
                  //int overlap = j-i, pp = j, end = p[j];
                  res.pos = bgn;
                  //for (; end < bgn+Be-1; end++)
                    //if (pp < p.size() && end == p[pp]) {
                      //overlap++;
                      //pp++;
                    //}
                  FORC(end, max(p[j], bgn+Be-1), end_hi) {
                    //if (pp < p.size() && end == p[pp]) {
                      //overlap++;
                      //pp++;
                    //}
                    res.len = end-bgn+q;
                    int o = overlap(e[eeid].c_str(), eL[eeid], doc+bgn, res.len);
                    res.sim = ft(o) / (res.len-q+1 + eL[eeid]-q+1 - o);
                    if (res.sim >= threshold)
                      result_jaccard->pb(res);
                  }
                }
              }
            }
            i++;
          } else {
            do {
              j = i+Tl-1; // j < |p|
              int lo = i, hi = j+1;
              while (lo < hi) { // [lo,hi)
                int mi = (lo+hi)>>1;
                if (p[j]+(mi-i)-p[mi]+1 > Te) lo = mi+1; // imply p[mi+j-i]-p[mi]+1 > Te
                else hi = mi;
              }
              i = lo;
              j = i+Tl-1;
            } while (i <= int(p.size())-Tl && p[j]-p[i]+1 > Te);
          }
        }
      eeid = eid;
      p.clear();
      p.pb(pos);
    }
  } while (eeid < INT_MAX);
}

int AEE::aeeED(const char *doc, unsigned threshold, vector<EDExtractResult> &result)
{
  faerie(doc, ft(threshold), true, &result, nullptr);
  return 0;
}

int AEE::aeeJaccard(const char *doc, double threshold, vector<JaccardExtractResult> &result)
{
  faerie(doc, threshold, false, nullptr, &result);
  return 0;
}
