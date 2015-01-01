#pragma once

#if 0
#define DEBUG
#else
#endif
#define MERGE_SKIP 1
#define DIVIDE_SKIP 1

#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <utility>
#include <vector>
using namespace std;

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
#define FOR(i, a, b) for (decltype(b) i = (a); i < (b); i++)
#define REP(i, n) FOR(i, 0, n)
#define REP1(i, n) for (decltype(n) i = 1; i <= (n); i++)
#define MP make_pair
#define PB push_back
#define ALL(x) (x).begin(), (x).end()
#define fi first
#define se second

const int SUCCESS = 0;
const int FAILURE = 1;
const int NB = 50000017;
const int CAP = 20000003;
const int BASE = 31;

const u N = 2048+2;
const ft eps = 1e-8;

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

  void push_back(const T &v) {
    a.PB(v);
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

struct QGram
{
  QGram(u q);
  IdAllocator<NB> allo;
  u q, pp;
  vector<string> strs;

  VU frequency;
  vector<bool> tick;
  bool flip;
  u &freq(u id) {
    if (tick[id] != flip) {
      tick[id] = flip;
      frequency[id] = 0;
    }
    return frequency[id];
  }

  u bf_overlap(const char *a, u n, const char *b, u m);
  u bf_overlap2(const char *b, u m);
};

struct BruteForce : QGram
{
public:
  BruteForce(const char *filename, u q);
  int searchJaccard(const char *query, ft threshold, vector<pair<u, ft> > &res);
  int searchED(const char *query, u threshold, VPUU &res);
};

template<typename Derived>
struct Tournament : QGram
{
  Tournament(const char *filename, u q);
  void search(const char *query, u n, u overlap, VPUU &res);
  int searchJaccard(const char *query, ft threshold, vector<pair<u, ft> > &res);
  int searchED(const char *query, u threshold, VPUU &res);

  u nrows, gs_min;
  map<u, VU> mapping;

  void filter_jaccard(const char *query, u n, ft threshold, const VPUU &cand, vector<pair<u, ft>> &res);
  void filter_edit_distance(const char *query, u n, u threshold, const VPUU &cand, VPUU &res);
};

struct MergeSkip : Tournament<MergeSkip>
{
  MergeSkip(const char *filename, u q);
  void search(const char *query, u n, u overlap, VPUU &res);
};

struct DivideSkip : Tournament<DivideSkip>
{
  DivideSkip(const char *filename, u q);
  void search(const char *query, u n, u overlap, VPUU &res);
};

class SimSearcher
{
public:
  SimSearcher();
  ~SimSearcher();

  int createIndex(const char *filename, u q);
  int searchJaccard(const char *query, ft threshold, vector<pair<u, double> > &result);
  int searchED(const char *query, u threshold, VPUU &result);

protected:
  vector<string> strs;
#ifdef BRUTE_FORCE
  BruteForce *impl = NULL;
#elif defined(DIVIDE_SKIP)
  DivideSkip *impl = NULL;
#elif defined(MERGE_SKIP)
  MergeSkip *impl = NULL;
#else
  Tournament *impl = NULL;
#endif
};

u levenshteinDistance(const char *a, const char *b, int th);
