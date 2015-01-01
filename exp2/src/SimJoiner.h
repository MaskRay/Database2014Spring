#pragma once
#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <climits>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <x86intrin.h>
using namespace std;

typedef __m128i mmi;
typedef uint8_t u8;
typedef double ft;
typedef unsigned u;
typedef vector<int> VI;
typedef vector<VI> VVI;
typedef pair<int, int> PII;
typedef vector<PII> VPII;
typedef vector<u> VU;
typedef vector<VU> VVU;
typedef pair<u, u> PUU;
typedef vector<PUU> VPUU;
typedef vector<string> VS;
#define FOR(i, a, b) for (remove_reference<remove_cv<decltype(b)>::type>::type i = (a); i < (b); i++)
#define REP(i, n) FOR(i, 0, n)
#define REP1(i, n) for (decltype(n) i = 1; i <= (n); i++)
#define mp make_pair
#define pb push_back
#define eb emplace_back
#define ALL(x) (x).begin(), (x).end()
#define fi first
#define se second

#define B_READ(file)                       \
  {                                        \
    char line[N];                          \
    FILE *fp = fopen(file, "r");           \
    while (fgets(line, sizeof line, fp)) { \
      u n = strlen(line)-1;                \
      line[n] = 0;

#define E_READ   \
    }            \
    fclose(fp);  \
  }

#define B_RABIN_KARP(str, n, i)      \
  {                                  \
    u h = 0;                         \
    REP(i, q-1)                      \
      h = h * BASE + (u8)str[i];     \
    REP(i, int(n)-(q-1)) {           \
      h = h * BASE + (u8)str[i+q-1];

#define E_RABIN_KARP(str, n, i)      \
      h -= q_pow_p1 * (u8)str[i];    \
    }                                \
  }

template <typename _IDType, typename _SimType>
struct JoinResult
{
	_IDType id1;
	_IDType id2;
	_SimType s;
};

typedef JoinResult<unsigned, double> JaccardJoinResult;
typedef JoinResult<unsigned, unsigned> EDJoinResult;

const int SUCCESS = 0;
const int FAILURE = 1;

const u N = 2048+2;
const int NB = 50000017;
const int CAP = 20000003;
const int BASE = 31;
const ft EPS = 1e-8;
const int SUFFIX_FILTER_DEPTH = 1;

struct CharHistogram
{
  union { u8 f[32]; mmi p[2]; };
  static int charmap[256];

  CharHistogram(const char *s) {
    p[0] = p[1] = _mm_setzero_si128();
    for (; *s; s++)
      f[charmap[*(u8 *)s]]++;
  }

  u delta(const CharHistogram &o) const {
    mmi a = _mm_add_epi64(_mm_sad_epu8(p[0], o.p[0]), _mm_sad_epu8(p[1], o.p[1]));
    //printf("%08x %08x %08x %08x\n", *(int*)&a, *((int*)&a+1), *((int*)&a+2), *((int*)&a+3));
    return _mm_cvtsi128_si32(_mm_srli_si128(a, 16)) + _mm_cvtsi128_si32(a);
  }

  static void initialize() {
    REP(i, 256)
      charmap[i] = 31;
    for (int i = 'a'; i <= 'z'; i++) {
      charmap[i] = i-'a';
      charmap[toupper(i)] = i-'a';
    }
    charmap[' '] = 26;
    charmap[','] = 27;
    charmap['.'] = 28;
    charmap[':'] = 29;
    charmap['?'] = 30;
  }
};

struct Counter
{
  Counter(int n) : n(n), a((int*)calloc(n, sizeof(int))), b((int*)malloc(n*sizeof(int))), c((int*)malloc(n*sizeof(int))), num(0) {}
  ~Counter() {
    free(a);
    free(b);
    free(c);
  }
  void reset() { num = 0; }
  int &operator[](int x) {
    if (a[x] < num && b[a[x]] == x)
      return c[x];
    a[x] = num;
    b[num++] = x;
    return c[x] = 0;
  }
  int n, num, *a, *b, *c;
};

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

template<u NB>
struct StrAllocator
{
  u q, cap, size = 0;
  int *b;
  struct Node { int next; const char *k; u len; } *a;

  StrAllocator(u cap) : cap(cap), b(new int[NB]), a((Node *)malloc(sizeof(Node)*cap)) {
    fill_n(b, NB, -1);
  }

  ~StrAllocator() {
    delete[] b;
    free(a);
  }

  u hash(const char *k, u len) const {
    u h = 0;
    REP(i, len)
      h = h*BASE+((unsigned char *)k)[i];
    return h%NB;
  }

  int get(const char *k, u len) {
    u h = hash(k, len);
    int i;
    for (i = b[h]; i != -1; i = a[i].next)
      if (! memcmp(a[i].k, k, len))
        return i;
    return -1;
  }

  int find(const char *k, u len) {
    u h = hash(k);
    int i;
    for (i = b[h]; i != -1; i = a[i].next)
      if (! memcmp(a[i].k, k, len))
        return i;

    if (size == cap) {
      cap <<= 1;
      a = (Node *)realloc(a, sizeof(Node)*cap);
    }

    i = size++;
    a[i].k = k;
    a[i].len = len;
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

enum QueryType { JACCARD, ED };

// original {{{1

class SimJoiner
{
public:
  SimJoiner() {}
  ~SimJoiner() {}
  void ppjoin(VS R, VS S, QueryType qt, ft th, int delta, vector<JaccardJoinResult> *res_jaccard, vector<EDJoinResult> *res_ed);
  int overlap(const char *a, u n, const char *b, u m);
  inline bool cmp_gram(int a, int b);
  inline bool partition(int *x, int m, int w, int lo, int hi, int &lx, int &rx, int &diff);
  int suffix_filter(int *x, int m, int *y, int n, int hamming_ub, int dep);

  int q, q_pow_p1;
  IdAllocator<NB> *s2g;
  vector<VPII> prefg2sid;
  VI df;
  vector<VI> gss;
  vector<VI> lenL;

  VI freqs;
  vector<bool> tick;
  bool flip;
  int &freq(int id) {
    if (tick[id] != flip) {
      tick[id] = flip;
      freqs[id] = 0;
    }
    return freqs[id];
  }

  int joinJaccard(const char *filename1, const char *filename2, unsigned q, double threshold, std::vector<JaccardJoinResult> &result);
  int joinED(const char *filename1, const char *filename2, unsigned q, unsigned threshold, std::vector<EDJoinResult> &result);
};
