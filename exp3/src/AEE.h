#pragma once
#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <err.h>
#include <fcntl.h>
#include <functional>
#include <map>
#include <queue>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <type_traits>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;

typedef double ft;
typedef unsigned u;
typedef uint8_t u8;
typedef uint32_t u32;
typedef vector<int> VI;
typedef vector<VI> VVI;
typedef pair<int, int> PII;
typedef vector<PII> VPII;
typedef vector<u> VU;
typedef vector<VU> VVU;
typedef pair<u, u> PUU;
typedef vector<PUU> VPUU;
typedef vector<string> VS;
#define FOR(i, a, b) for (remove_cv<decltype(b)>::type i = (a); i < (b); i++)
#define FORC(i, a, b) for (decltype(b) i = (a); i <= (b); i++)
#define REP(i, n) FOR(i, 0, n)
#define REP1(i, n) for (decltype(n) i = 1; i <= (n); i++)
#define mp make_pair
#define pb push_back
#define eb emplace_back
#define ALL(x) (x).begin(), (x).end()
#define fi first
#define se second

const int NB = 50000017;
const int CAP = 20000003;
const int BASE = 31;
const int SUCCESS = 0;
const int FAILURE = 1;
const u N = 2048+2;
const int ALPHABET = 26+10+1;
const ft eps = 1e-8;

#define B_RABIN_KARP(str, n, i)    \
  {                                \
    u h = 0;                       \
    REP(i, q-1)                    \
      h = h * BASE + str[i];       \
    if (n >= q-1)                  \
      REP(i, n-(q-1)) {            \
        h = h * BASE + str[i+q-1];

#define E_RABIN_KARP(str, n, i)    \
        h -= q_pow_p1 * str[i];    \
      }                            \
  }

template<u NB>
struct IdAllocator
{
  u q, cap, size = 0;
  int *b;
  struct Node { int next; u k; const char *kk; } *a;

  IdAllocator(u q, u cap) : q(q), cap(cap), b(new int[NB]), a((Node *)malloc(sizeof(Node)*cap)) {
    fill_n(b, NB, -1);
  }

  ~IdAllocator() {
    delete[] b;
    free(a);
  }

  u hash(u k) const {
    return k % NB;
  }

  int get(u k, const char *kk) {
    u h = hash(k);
    int i;
    for (i = b[h]; i != -1; i = a[i].next)
      if (a[i].k == k && ! memcmp(a[i].kk, kk, q))
        return i;
    return -1;
  }

  int find(u k, const char *kk) {
    u h = hash(k);
    int i;
    for (i = b[h]; i != -1; i = a[i].next)
      if (a[i].k == k && ! memcmp(a[i].kk, kk, q))
        return i;

    if (size == cap) {
      cap <<= 1;
      a = (Node *)realloc(a, sizeof(Node)*cap);
    }

    i = size++;
    a[i].k = k;
    a[i].kk = kk;
    a[i].next = b[h];
    b[h] = i;
    return i;
  }
};

template<typename T, typename Cmp>
struct Heap
{
  void down(u x) {
    T key = a[x];
    for (u y; y = 2*x+1, y < n; x = y) {
      if (y+1 < n && cmp(a[y+1], a[y])) y++;
      if (! cmp(a[y], key)) break;
      a[x] = a[y];
    }
    a[x] = key;
  }

  void up(u x) {
    T key = a[x];
    for (; x && Cmp()(key, a[(x-1)/2]); x = (x-1)/2)
      a[x] = a[(x-1)/2];
    a[x] = key;
  }

  void make() {
    for (u i = n/2; i > 0; )
      down(--i);
  }

  template<typename... Args>
  void emplace_back(Args&&... args) {
    a.eb(forward<Args>(args)...);
    n++;
  }

  void push_back(const T &v) {
    a.pb(v);
    n++;
  }

  u size() const {
    return a.size();
  }

  bool empty() const {
    return a.empty();
  }

  void pop() {
    a[0] = a.back();
    a.pop_back();
    n--;
    down(0);
  }

  void pop_to_back() {
    swap(a[0], a[--n]);
    down(0);
  }

  void emplace() {
    n++;
    up(n-1);
  }

  T &operator[](u x) {
    return a[x];
  }

  vector<T> a;
  Cmp cmp;
  u n = 0;
};

tuple<int, char *, size_t> read_file(const char *filename);


template <typename _IDType, typename _PosType, typename _LenType, typename _SimType>
struct ExtractResult
{
	_IDType id;
	_PosType pos;
	_LenType len;
	_SimType sim;
};

typedef ExtractResult<unsigned, unsigned, unsigned, unsigned> EDExtractResult;
typedef ExtractResult<unsigned, unsigned, unsigned, double> JaccardExtractResult;

class AEE {
public:
  AEE() : tup(-1, nullptr, 0) {}
  ~AEE() { close_file(); }

  int createIndex(const char *filename, unsigned q);
  int aeeJaccard(const char *doc, double threshold, std::vector<JaccardExtractResult> &result);
  int aeeED(const char *doc, unsigned threshold, std::vector<EDExtractResult> &result);
  int overlap(const char *a, int m, const char *b, int n);

  void faerie(const char *doc, ft threshold, bool flag_ed, vector<EDExtractResult> *, vector<JaccardExtractResult> *);
  void close_file();

  int q, q_pow_p1, ngid;
  VS e;
  VS re;
  tuple<int, char *, size_t> tup;

  VI freqs;
  vector<bool> ticks;
  bool flip;
  int &freq(int id) {
    if (ticks[id] != flip) {
      ticks[id] = flip;
      freqs[id] = 0;
    }
    return freqs[id];
  }

  int last_tau = -1;
  struct TrieNode {
    TrieNode *ch[ALPHABET], *pi;
    vector<pair<PII, int>> *sub;
    int dep;
    TrieNode() {
      fill_n(ch, ALPHABET, nullptr);
      sub = nullptr;
    }
    ~TrieNode() {
      REP(i, ALPHABET)
        delete ch[i];
    }
  } *trie_root = nullptr;
  void trie_traverse(TrieNode *p, function<void(vector<pair<PII,int>>&)> fn);
  void trie_t(TrieNode *p, int d);
};
