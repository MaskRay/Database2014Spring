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
  this->q = q;
  q_pow_p1 = ipow(BASE, q-1);
  e.clear();

  close_file();

  tup = read_file(filename);
  char *mapped = get<1>(tup);

  int file_size = int(get<2>(tup)), eid = 0;
  for (int i = 0, j = 0; j < file_size; i = j, eid++) {
    for (; j < file_size && mapped[j] != '\n'; j++);
    e.eb(mapped+i, j-i);
    re.eb(mapped+i, j-i);
    reverse(ALL(re.back()));
    if (j < file_size && mapped[j] == '\n')
      j++;
  }
	return SUCCESS;
}

struct Cmp
{
  bool operator()(const pair<int *, int> &a, const pair<int *, int> &b) const {
    if (*a.fi != *b.fi) return *a.fi < *b.fi;
    return a.se < b.se;
  }
};

void AEE::trie_traverse(TrieNode *p, function<void(vector<pair<PII,int>>&)> fn)
{
  if (! p) return;
  if (p->sub)
    fn(*p->sub);
  REP(i, ALPHABET)
    if (p->ch[i] && p->ch[i]->dep == p->dep+1)
      trie_traverse(p->ch[i], fn);
}

void AEE::trie_t(TrieNode *p, int d)
{
  if (! p) return;
  REP(i,d)printf("  ");
  printf("%p(", p);
  REP(i, ALPHABET)
    if (p->ch[i] && p->ch[i]->dep == p->dep+1)
      printf("%p,",p->ch[i]);
  puts(")");
  REP(i, ALPHABET)
    if (p->ch[i] && p->ch[i]->dep == p->dep+1)
      trie_t(p->ch[i], d+1);
}

int AEE::aeeED(const char *doc, unsigned threshold, vector<EDExtractResult> &result)
{
  int tau = threshold;
  if (tau != last_tau) {
    int tag = 0;
    if (last_tau != -1)
      delete trie_root;
    trie_root = new TrieNode;
    trie_root->dep = 0;
    REP(i, e.size()) {
      int n = e[i].size(), base_len = n/(tau+1), larger = n%(tau+1), pos = 0;
      REP(j, tau+1) {
        int pos2 = pos+base_len+(j<larger);
        TrieNode *p = trie_root;
        FOR(g, pos, pos2) {
          u8 c = e[i][g];
          if (c == ' ') c = 36;
          else if (unsigned(c-'0') < 10) c = c-'0'+26;
          else c -= 'a';
          TrieNode *&q = p->ch[c];
          if (! q) {
            q = new TrieNode;
            q->dep = g-pos+1;
          }
          p = q;
        }
        if (! p->sub)
          p->sub = new vector<pair<PII, int>>;
        p->sub->eb(mp(i, pos2), tag);
        p->sub->eb(mp(i, ~ pos), tag++);
        pos = pos2;
      }
    }

    // aho-corasick automaton
    queue<TrieNode*> q;
    REP(i, ALPHABET)
      if (trie_root->ch[i])
        trie_root->ch[i]->pi = trie_root, q.push(trie_root->ch[i]);
      else
        trie_root->ch[i] = trie_root;
    while (! q.empty()) {
      TrieNode *p = q.front();
      q.pop();
      REP(i, ALPHABET) {
        if (p->ch[i])
          p->ch[i]->pi = p->pi->ch[i], q.push(p->ch[i]);
        else
          p->ch[i] = p->pi->ch[i];
      }
    }

    //trie_t(trie_root, 0);

    trie_traverse(trie_root, [&](vector<pair<PII,int>> &a) {
                  if (a.empty()) return;
                  // pre (<0) | post (>=0)
                  sort(ALL(a), [](const pair<PII,int> &b, const pair<PII,int> &c) {
                       return b.fi.se < c.fi.se;
                       });
                  // sort pre
                  sort(a.begin(), a.begin()+a.size()/2, [&](const pair<PII,int> &b, const pair<PII,int> &c) {
                       return strcmp(&re[b.fi.fi][re[b.fi.fi].size()-~b.fi.se], &re[c.fi.fi][re[c.fi.fi].size()-~c.fi.se]) < 0;
                       });
                  // sort post
                  sort(a.begin()+a.size()/2, a.end(), [&](const pair<PII,int> &b, const pair<PII,int> &c) {
                       return strcmp(&e[b.fi.fi][b.fi.se], &e[c.fi.fi][c.fi.se]) < 0;
                       });
                  });
    last_tau = tau;
  }
  result.clear();

  static int dp[N][N];
  int m = strlen(doc);
  char *rdoc = strdup(doc);
  reverse(rdoc, rdoc+m);
  TrieNode *q = trie_root, *p;
  for (int i, ii = 0, j = 0; j < m; ) {
    u8 c = doc[j++];
    if (c == ' ') c = 36;
    else if (unsigned(c-'0') < 10) c = c-'0'+26;
    else c -= 'a';
    ii += q->dep+1-q->ch[c]->dep;
    q = q->ch[c];
    for (i = ii, p = q; p; i += p->pi ? p->dep-p->pi->dep : 0, p = p->pi)
    if (p->sub) {
      int eeid = -1, eer = 0;
      const char *aa;

      unordered_map<int, vector<VI>> mapping;
      FOR(l, 0, tau+1)
        dp[0][l] = l;
      REP(eidx, p->sub->size()) {
        bool left = eidx < p->sub->size()/2;
        int eid = p->sub->at(eidx).fi.fi, pos = p->sub->at(eidx).fi.se,
            tag = p->sub->at(eidx).se;
        int er;
        int lcp = 0;
        const char *a, *b;
        if (left) {
          er = ~pos;
          a = &re[eid][re[eid].size()-er];
          b = rdoc+m-i;
        } else {
          a = &e[eid][pos];
          b = doc+j;
          er = e[eid].size()-pos;
        }
        if (eidx == p->sub->size()/2)
          eeid = -1;
        if (eeid != -1) {
          int maxl = min(er, eer);
          while (lcp < maxl && a[lcp] == aa[lcp])
            lcp++;
        }
        eeid = eid;
        eer = 0;
        aa = a;
        for (int k = lcp+1; k <= er; k++) {
          bool stop = true;
          int l = max(k-tau, 0), h = min(k+tau, left ? i : m-j);
          if (! l)
            dp[k][l++] = k;
          FOR(o, l, h+1) {
            dp[k][o] = a[k-1] == b[o-1] ? dp[k-1][o-1] : min(min(o-(k-1) > tau ? tau : dp[k-1][o], k-(o-1) > tau ? tau : dp[k][o-1]), dp[k-1][o-1]) + 1;
            if (dp[k][o] <= tau)
              stop = false;
          }
          eer = k;
          if (stop) goto L1;
        }
        {
          int l = max(er-tau, 0), h = min(er+tau, left ? i : m-j);
          vector<VI> B(tau+1);
          if (left) {
            FOR(o, l, h+1)
              if (dp[er][o] <= tau && i-o >= 0)
                B[dp[er][o]].pb(i-o);
            mapping[tag] = B;
          } else {
            if (! mapping.count(tag)) continue;
            vector<VI> &BB = mapping[tag];
            FOR(o, l, h+1)
              if (dp[er][o] <= tau && j+o <= m)
                B[dp[er][o]].pb(j+o);
            REP(u, tau+1)
              REP(v, tau-u+1)
              for (auto uu: BB[u])
                for (auto vv: B[v]) {
                  EDExtractResult res;
                  res.id = eid;
                  res.pos = uu;
                  res.len = vv-uu;
                  res.sim = u+v;
                  result.pb(res);
                }
          }
        }
L1:;
      }
    }
  }
  sort(ALL(result), [](const EDExtractResult &a, const EDExtractResult &b) {
       if (a.id != b.id) return a.id < b.id;
       if (a.pos != b.pos) return a.pos < b.pos;
       return a.len < b.len;
       });
  result.erase(unique(ALL(result), [](const EDExtractResult &a, const EDExtractResult &b) {
         return a.id == b.id && a.pos == b.pos && a.len == b.len;
         }), result.end());
  free(rdoc);
  return 0;
}

int AEE::aeeJaccard(const char *doc, double threshold, vector<JaccardExtractResult> &result)
{
  return 0;
}
