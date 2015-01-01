#include <algorithm>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <err.h>
#include <errno.h>
#include <functional>
#include <inttypes.h>
#include <jsoncpp/json/json.h>
#include <map>
#include <netinet/in.h>
#include <queue>
#include <set>
#include <stack>
#include <stdint.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <unordered_set>
#include <vector>
using namespace std;

#define REP(i, n) for (int i = 0; i < (n); i++)
#define REP1(i, n) for (int i = 1; i <= (n); i++)
#define FOR(i, a, b) for (int i = (a); i < (b); i++)
#define ROF(i, a, b) for (int i = (b); --i >= (a); )
#define mp make_pair
#define fi first
#define se second
#define pb push_back
#define eb emplace_back
#define ALL(x) (x).begin(), (x).end()

typedef double ft;
typedef vector<int> VI;
typedef pair<int, int> PII;

const int PORT = 1337;
const int D = 2;
const int MAX_MATCH = 10;

int RI()
{
  int x;
  scanf("%d", &x);
  return x;
}

int RI(int &x)
{
  return scanf("%d", &x);
}

struct V {
  int id;
  ft a[D];
  ft &operator[](int axis) {
    return a[axis];
  }
  ft operator[](int axis) const {
    return a[axis];
  }
  ft get_min(int axis) const {
    return a[axis];
  }
  ft get_max(int axis) const {
    return a[axis];
  }
  ft dist2(const V &o) const {
    ft r = 0;
    REP(i, D)
      r += (a[i]-o.a[i])*(a[i]-o.a[i]);
    return r;
  }
};

struct MBR {
  vector<ft> min, max;
  bool valid() const {
    REP(i, D)
      if (min[i] > max[i])
        return false;
    return true;
  }
  ft get_min(int axis) const {
    return min[axis];
  }
  ft get_max(int axis) const {
    return max[axis];
  }
  bool intersect(const V &o) const {
    REP(i, D)
      if (o[i] < min[i] || max[i] < o[i])
        return false;
    return true;
  }
  bool intersect(const MBR &o) const {
    REP(i, D)
      if (o.max[i] < min[i] || max[i] < o.min[i])
        return false;
    return true;
  }
  ft dist2(const V &o) const {
    ft r = 0;
    REP(i, D) {
      if (o[i] < min[i])
        r += (min[i]-o[i])*(min[i]-o[i]);
      else if (max[i] < o[i])
        r += (o[i]-max[i])*(o[i]-max[i]);
    }
    return r;
  }
};

struct Node {
  Node(bool leaf) : leaf(leaf) {}
  bool leaf;
  MBR mbr;
  ft get_min(int axis) {
    return mbr.get_min(axis);
  }
  ft get_max(int axis) {
    return mbr.get_max(axis);
  }
};
struct Leaf : Node {
  Leaf() : Node(true) {}
  vector<V> a;
};
struct Inner : Node {
  Inner() : Node(false) {}
  vector<Node *> a;
};

template<typename T>
struct NodeCmpGen {
  const vector<T> &arr;
  NodeCmpGen(const vector<T> &arr) : arr(arr) {}
  function<bool(int, int)> get_min_cmp(int axis) {
    return [=](int a, int b) -> bool {
      return arr[a].get_min(axis) < arr[b].get_min(axis);
    };
  }
  function<bool(int, int)> get_max_cmp(int axis) {
    return [=](int a, int b) -> bool {
      return arr[a].get_max(axis) > arr[b].get_max(axis);
    };
  }
};
template<>
struct NodeCmpGen<Node *> {
  const vector<Node *> &arr;
  NodeCmpGen(const vector<Node *> &arr) : arr(arr) {}
  function<bool(int, int)> get_min_cmp(int axis) {
    return [=](int a, int b) -> bool {
      return arr[a]->get_min(axis) < arr[b]->get_min(axis);
    };
  }
  function<bool(int, int)> get_max_cmp(int axis) {
    return [=](int a, int b) -> bool {
      return arr[a]->get_max(axis) > arr[b]->get_max(axis);
    };
  }
};

const int B = 10;

template<typename T, typename L>
struct Builder {
  L build(const vector<T> &);
};
template<>
struct Builder<V, Leaf *> {
  Leaf *build(const vector<V> &data) {
    auto r = new Leaf;
    r->a = data;
    r->mbr.min = vector<ft>(D);
    r->mbr.max = vector<ft>(D);
    REP(i, D) {
      r->mbr.min[i] = DBL_MAX;
      for (auto &x: data)
        r->mbr.min[i] = min(r->mbr.min[i], x[i]);
      r->mbr.max[i] = -DBL_MAX;
      for (auto &x: data)
        r->mbr.max[i] = max(r->mbr.max[i], x[i]);
    }
    return r;
  }
};
template<>
struct Builder<Node *, Inner *> {
  Inner *build(const vector<Node *> &data) {
    auto r = new Inner;
    r->a = data;
    r->mbr.min = vector<ft>(D);
    r->mbr.max = vector<ft>(D);
    REP(i, D) {
      r->mbr.min[i] = DBL_MAX;
      for (auto x: data)
        r->mbr.min[i] = min(r->mbr.min[i], x->get_min(i));
      r->mbr.max[i] = -DBL_MAX;
      for (auto &x: data)
        r->mbr.max[i] = max(r->mbr.max[i], x->get_max(i));
    }
    return r;
  }
};

void traverse(Node *x, int dep)
{
  if (!x) return;
  if (x->leaf) {
    Leaf *y = (Leaf *)x;
    MBR &mbr = y->mbr;
    REP(i,dep)printf("  ");
    printf("leaf (%lg,%lg) (%lg,%lg)\n", mbr.min[0], mbr.min[1], mbr.max[0], mbr.max[1]);
    for (auto &p: ((Leaf*)x)->a) {
      REP(i,dep+1)printf("  ");
      printf("%lg,%lg\n", p[0], p[1]);
    }
  } else {
    Inner *y = (Inner *)x;
    MBR &mbr = y->mbr;
    REP(i,dep)printf("  ");
    printf("inner (%lg,%lg) (%lg,%lg)\n", mbr.min[0], mbr.min[1], mbr.max[0], mbr.max[1]);
    for (auto &p: ((Inner*)x)->a)
      traverse(p, dep+1);
  }
}

class PRTree
{
public:
  PRTree() {
  }
  ~PRTree() {
  }
  Node *root;

  void load(vector<V> &data) {
    vector<Node *> nodes, nodes2;
    build_leaves<V, Leaf *, Node *>(data, nodes);
    //for (auto &x: nodes) {
    //  printf("Node\n");
    //  for (auto &y: ((Leaf*)x)->a)
    //    printf("  %lg %lg\n", y[0], y[1]);
    //}
    while (nodes.size() > B) {
      nodes2.clear();
      build_leaves<Node *, Inner *, Node *>(nodes, nodes2);
      nodes.swap(nodes2);
    }
    if (nodes.size() == 0)
      root = Builder<V, Leaf *>().build(data);
    else if (nodes.size() == 1)
      root = nodes[0];
    else
      root = Builder<Node *, Inner *>().build(nodes);
  }

  static void find(Node *rt, const MBR &mbr, VI &res) {
    if (rt->leaf) {
      Leaf *x = (Leaf *)rt;
      for (auto &p: x->a)
        if (mbr.intersect(p))
          res.pb(p.id);
    } else {
      Inner *x = (Inner *)rt;
      for (auto p: x->a)
        if (mbr.intersect(p->mbr))
          find(p, mbr, res);
    }
  }

  void find(MBR mbr, VI &res) {
    res.clear();
    find(root, mbr, res);
  }

  struct Entity {
    bool obj;
    union { V v; Node *n; };
  };
  typedef pair<ft, Entity> QElem;

  friend bool operator>(const QElem &a, const QElem &b) {
    return a.fi > b.fi;
  }

  void k_nearest_neighbor(const V &q, int k, VI &res, function<bool(int)> filter) {
    res.clear();
    if (! root->mbr.valid()) return;
    priority_queue<QElem, vector<QElem>, greater<QElem>> pq;
    Entity ent;
    ent.obj = false;
    ent.n = root;
    pq.push(mp(root->mbr.dist2(q), ent));
    while (! pq.empty() && res.size() < k) {
      ent = pq.top().se;
      pq.pop();
      if (ent.obj) {
        if (filter(ent.v.id))
          res.pb(ent.v.id);
      } else if (ent.n->leaf)
        for (auto &p: ((Leaf*)ent.n)->a) {
          Entity ent2;
          ent2.obj = true;
          ent2.v = p;
          pq.push(mp(p.dist2(q), ent2));
        }
      else
        for (auto p: ((Inner*)ent.n)->a) {
          Entity ent2;
          ent2.obj = false;
          ent2.n = p;
          pq.push(mp(p->mbr.dist2(q), ent2));
        }
    }
  }

  void k_nearest_neighbor(const V &q, int k, VI &res) {
    k_nearest_neighbor(q, k, res, [](int) { return true; });
  }

  struct Slice {
    int dep, id, num;
    VI pos;
    Slice(int dep, int id, int num, VI pos) : dep(dep), id(id), num(num), pos(pos) {}
  };

  template<typename T, typename L, typename N>
  void build_leaves(vector<T> &data, vector<N> &leaves) {
    NodeCmpGen<T> cmpgen(data);
    Builder<T, L> builder;

    vector<VI> ids(D*2, VI(data.size()));
    vector<function<bool(int, int)>> fns;
    REP(i, D*2)
      iota(ALL(ids[i]), 0);
    REP(i, D)
      fns.push_back(cmpgen.get_min_cmp(i));
    REP(i, D)
      fns.push_back(cmpgen.get_max_cmp(i));
    REP(i, D*2)
      sort(ALL(ids[i]), fns[i]);

    map<int, int> owner;
    REP(i, data.size())
      owner[i] = 0;
    stack<Slice> slice_stack;
    slice_stack.push(Slice(0, 0, data.size(), VI(D*2, 0)));
    while (! slice_stack.empty()) {
      Slice slice = slice_stack.top();
      slice_stack.pop();
      REP(i, D*2) {
        int t = min(slice.num, B);
        if (! t) break;
        vector<T> node_data;
        REP(j, t) {
          while (slice.pos[i] < data.size() && owner[ids[i][slice.pos[i]]] != slice.id)
            slice.pos[i]++;
          node_data.pb(data[ids[i][slice.pos[i]]]);
          owner[ids[i][slice.pos[i]++]] = -1;
        }
        leaves.pb(builder.build(node_data));
        slice.num -= t;
      }
      if (slice.num > 0) {
        int axis = slice.dep%(D*2);
        int num0 = slice.num/2, num1 = slice.num-num0;
        int k = slice.pos[axis];
        REP(j, num0) {
          while (k < data.size() && owner[ids[axis][k]] != slice.id)
            k++;
          owner[ids[axis][k]] = slice.id*2+1;
        }
        slice_stack.push(Slice(slice.dep+1, slice.id*2+1, num0, slice.pos));
        slice.pos[axis] = k;
        REP(j, num1) {
          while (k < data.size() && owner[ids[axis][k]] != slice.id)
            k++;
          owner[ids[axis][k]] = slice.id*2+2;
        }
        slice_stack.push(Slice(slice.dep+1, slice.id*2+2, num1, slice.pos));
      }
    }
  }
};

struct Place
{
  string addr, name;
  int pcode;
  ft lat, lng;
};

namespace KoAluru
{
  bool *t;
  int *b;

  template<typename T>
  void bucket(T a[], int n, int k, bool end)
  {
    fill_n(b, k, 0);
    REP(i, n) b[a[i]]++;
    if (end)
      FOR(i, 1, k) b[i] += b[i-1];
    else {
      int s = 0;
      REP(i, k)
        s += b[i], b[i] = s-b[i];
    }
  }

  template<typename T>
  void plus_to_minus(T a[], int sa[], int n, int k)
  {
    bucket(a, n, k, false);
    sa[b[a[n-1]]++] = n-1;
    REP(i, n-1) {
      int j = sa[i]-1;
      if (j >= 0 && ! t[j])
        sa[b[a[j]]++] = j;
    }
  }

  template<typename T>
  void minus_to_plus(T a[], int sa[], int n, int k)
  {
    bucket(a, n, k, true);
    ROF(i, 0, n) {
      int j = sa[i]-1;
      if (j >= 0 && t[j])
        sa[--b[a[j]]] = j;
    }
  }

  template<typename T>
  void ka(T a[], int sa[], int n, int k)
  {
    t[n-1] = false;
    ROF(i, 0, n-1)
      t[i] = a[i] < a[i+1] || a[i] == a[i+1] && t[i+1];
    bool minor = 2 * count(t, t+n, false) > n;

    bucket(a, n, k, minor);
    fill_n(sa, n, -1);
    if (minor) {
      REP(i, n)
        if (t[i])
          sa[--b[a[i]]] = i;
      plus_to_minus(a, sa, n, k);
      minus_to_plus(a, sa, n, k);
    } else {
      sa[b[a[n-1]]++] = n-1;
      REP(i, n-1)
        if (! t[i])
          sa[b[a[i]]++] = i;
      minus_to_plus(a, sa, n, k);
      plus_to_minus(a, sa, n, k);
    }

    int last = -1, name = 0, nn = count(t, t+n, minor);
    int *sa2, *pi;
    if (minor)
      sa2 = sa, pi = sa+n-nn;
    else
      sa2 = sa+n-nn, pi = sa;
    fill_n(b, n, -1);
    REP(i, n)
      if (sa[i] >= 0 && minor == t[sa[i]]) {
        bool diff = last == -1;
        int p = sa[i];
        if (! diff)
          REP(j, n) {
            if (last+j >= n || p+j >= n || a[last+j] != a[p+j] || t[last+j] != t[p+j]) {
              diff = true;
              break;
            } else if (j > 0 && (minor == t[last+j] || minor == t[p+j]))
              break;
          }
        if (diff) {
          name++;
          last = p;
        }
        b[p] = name-1;
      }
    nn = 0;
    REP(i, n)
      if (b[i] >= 0)
        pi[nn++] = b[i];

    if (name < nn)
      ka(pi, sa2, nn, name);
    else
      REP(i, nn)
        sa2[pi[i]] = i;

    ROF(i, 0, nn)
      t[i] = a[i] < a[i+1] || a[i] == a[i+1] && t[i+1];

    nn = 0;
    bucket(a, n, k, minor);
    if (minor) {
      REP(i, n)
        if (minor == t[i])
          pi[nn++] = i;
      REP(i, nn)
        sa[i] = pi[sa2[i]];
      ROF(i, 0, nn) {
        int j = sa[i];
        sa[i] = -1;
        sa[--b[a[j]]] = j;
      }
    } else {
      REP(i, n)
        if (minor == t[i])
          pi[nn++] = i;
      ROF(i, 0, nn)
        sa[n-nn+i] = pi[sa2[i]];
      REP(i, nn) {
        int j = sa[n-nn+i];
        sa[n-nn+i] = -1;
        sa[b[a[j]]++] = j;
      }
    }
    if (minor)
      plus_to_minus(a, sa, n, k);
    else
      minus_to_plus(a, sa, n, k);
  }

  template<typename T>
  void main(T a[], int sa[], int b[], int n, int k)
  {
    if (n > 0) {
      KoAluru::b = b;
      t = new bool[n];
      ka(a, sa, n, k);
      delete[] t;
    }
  }

  template<typename T>
  void calc_rank_lcp(T a[], int sa[], int n, int rank[], int lcp[])
  {
    REP(i, n)
      rank[sa[i]] = i;
    int k = 0;
    lcp[0] = 0;
    FOR(i, 0, n)
      if (rank[i]) {
        for (int j = sa[rank[i]-1]; i+k < n && j+k < n && a[i+k] == a[j+k]; k++);
        lcp[rank[i]] = k;
        k && k--;
      }
  }

  void calc_child(int lcp[], int n, int child[]) {
    stack<int> st;
    st.push(0);
    int last = -1;
    fill_n(child, n, -1);
    FOR(i, 1, n) {
      while (lcp[i] < lcp[st.top()]) {
        last = st.top();
        st.pop();
        if (lcp[i] <= lcp[st.top()] && lcp[st.top()] != lcp[last])
          child[st.top()] = last;
      }
      if (last != -1) {
        child[i-1] = last;
        last = -1;
      }
      st.push(i);
    }
    while (0 < lcp[st.top()]) {
      last = st.top();
      st.pop();
      if (0 <= lcp[st.top()] && lcp[st.top()] != lcp[last])
        child[st.top()] = last;
    }

    while (! st.empty())
      st.pop();
    st.push(0);
    FOR(i, 1, n) {
      while (lcp[i] < lcp[st.top()])
        st.pop();
      if (lcp[i] == lcp[st.top()]) {
        child[st.top()] = i;
        st.pop();
      }
      st.push(i);
    }
  }

  int get_lcp(int lcp[], int child[], int i, int j)
  {
    if (i == j-1) return lcp[j];
    int k = child[j-1]; // up[j]
    if (i < k && k <= j)
      return lcp[k];
    return child[i] != -1 ? lcp[child[i]] : -1; // down[j]
  }

  void get_child_intervals(int lcp[], int child[], int l, int h)
  {
    printf("(%d %d)\n", l, h);
    if (l >= h-1) return;
    int i = l < child[h-1] && child[h-1] < h ? child[h-1] : child[l];
    get_child_intervals(lcp, child, l, i);
    for (; child[i] > i && lcp[child[i]] == lcp[i]; i = child[i]) // next[i]
      get_child_intervals(lcp, child, i, child[i]);
    get_child_intervals(lcp, child, i, h);
  }

  template<typename T>
  PII get_interval(T a[], int sa[], int lcp[], int child[], int n, int l, int h, int d, T c)
  {
    int i = l == 0 && h == n ? child[0] : l < child[h-1] && child[h-1] < h ? child[h-1] : child[l];
    if (sa[l]+d < n && a[sa[l]+d] == c)
      return mp(l, i);
    for (; child[i] > i && lcp[child[i]] == lcp[i]; i = child[i]) { // next[i]
      if (a[sa[i]+d] == c)
        return mp(i, child[i]);
    }
    if (a[sa[i]+d] == c)
      return mp(i, h);
    return mp(-1,-1);
  }

  void top_down_traversal(int lcp[], int child[], int n)
  {
    get_child_intervals(lcp, child, 0, n);
  }

  template<typename T>
  PII search(T a[], int sa[], int lcp[], int child[], int n, T s[], int m)
  {
    if (m == 0)
      return mp(0, n);
    PII sub = get_interval(a, sa, lcp, child, n, 0, n, 0, s[0]);
    int l = sub.fi, h = sub.se, i = 1;
    bool found = true;
    while (found && i < m) {
      if (l < h-1) {
        int j = min(get_lcp(lcp, child, l, h), m);
        FOR(k, i, j) {
          if (a[sa[l]+k] != s[k]) {
            found = false;
            break;
          }
        }
        i = j;
        if (i < m) {
          sub = get_interval(a, sa, lcp, child, n, l, h, i, s[i]);
          l = sub.fi;
          h = sub.se;
        }
      } else {
        FOR(k, i, m) {
          if (a[sa[l]+k] != s[k]) {
            found = false;
            break;
          }
        }
        i = m;
      }
    }
    return found ? mp(l, h) : mp(-1, -1);
  }
};

struct Timer
{
  const char *s;
  clock_t bgn;
  Timer(const char *s) : s(s), bgn(clock()) {}
  ~Timer() { printf("%s: %lf\n", s, ft(clock()-bgn)/CLOCKS_PER_SEC); }
};

int main()
{
  int yes = 1, sockfd = socket(AF_INET, SOCK_STREAM, 0);
  if (sockfd == -1)
    errx(2, "failed to socket: %s", strerror(errno));
  if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof yes) == -1)
    errx(2, "failed to setsocket: %s", strerror(errno));
  struct sockaddr_in sin;
  sin.sin_family = AF_INET;
  sin.sin_addr.s_addr = INADDR_ANY;
  sin.sin_port = htons(PORT);
  if (bind(sockfd, (const struct sockaddr *)&sin, sizeof sin) == -1)
    errx(2, "failed to bind: %s", strerror(errno));
  if (listen(sockfd, 1) == -1)
    errx(2, "failed to listen: %s", strerror(errno));

  vector<Place> places;
  vector<V> vs;
  {
    Timer timer("read file");
    FILE *f = fopen("zipcode.json", "r");
    char buf[1024];
    while (fgets(buf, sizeof buf, f)) {
      Json::Value root;
      Json::Reader reader;
      bool flag = reader.parse(buf, root);
      if (! flag) {
        puts(buf);
        return 1;
      }

      Place p;
      p.addr = root["addr"].asString();
      p.name = root["name"].asString();
      p.pcode = root["pcode"].asInt();
      p.lat = root["latlng"][0].asDouble();
      p.lng = root["latlng"][1].asDouble();

      V v;
      v[0] = p.lat;
      v[1] = p.lng;
      v.id = vs.size();
      places.pb(p);
      vs.pb(v);
      if (vs.size() % 1000 == 0)
        printf("loaded %zd\n", vs.size());
    }
    fclose(f);
  }

  int len = 0;
  for (auto &p: places)
    len += p.name.size()+1;
  auto names = new unsigned char[len];
  len = 0;
  for (auto &p: places) {
    strcpy((char *)names+len, p.name.c_str());
    len += p.name.size()+1;
  }
  int *sa = new int[len], *id = new int[len], *rnk = new int[len], *lcp = new int[len], *child = new int[len];
  {
    Timer timer("Ko-Aluru suffix array");
    KoAluru::main(names, sa, rnk, len, 256);
    KoAluru::calc_rank_lcp(names, sa, len, rnk, lcp);
    KoAluru::calc_child(lcp, len, child);
    len = 0;
    REP(i, places.size()) {
      FOR(j, len, len+places[i].name.size()+1)
        id[rnk[j]] = i;
      len += places[i].name.size()+1;
    }
  }

  PRTree prtree;
  {
    Timer timer("Priority R-tree");
    prtree.load(vs);
  }
  Json::FastWriter writer;
  for(;;) {
    struct sockaddr_storage ss;
    socklen_t sslen = sizeof ss;
    int clifd = accept(sockfd, (struct sockaddr *)&ss, &sslen);
    if (clifd == -1) break;
    pid_t pid = fork();
    if (pid == -1) break;
    if (pid == 0) {
      char buf[1024], query[1024];
      VI res;
      V q;
      bool first = true;
      int qt;
      alarm(60);
      FILE *f = fdopen(clifd, "r");
      if (! fgets(buf, sizeof buf, f) || sscanf(buf, "%d", &qt) != 1)
        goto err;
      if (qt == 0) {
        query[0] = '\0';
        if (! fgets(buf, sizeof buf, f) || sscanf(buf, "%lf%lf %[^\n]", &q[0], &q[1], query) < 2)
          goto err;
        prtree.k_nearest_neighbor(q, MAX_MATCH, res, [&](int id) {
                                  return strstr(places[id].name.c_str(), query) != NULL;
                                  });
      } else {
        if (! fgets(buf, sizeof buf, f))
          goto err;
        int m = int(strlen(buf));
        if (m)
          buf[--m] = '\0';
        PII r = KoAluru::search(names, sa, lcp, child, len, (unsigned char *)buf, m);
        unordered_set<int> used;
        FOR(i, r.fi, r.se) {
          used.insert(id[i]);
          if (used.size() == MAX_MATCH) break;
        }
        for (auto i: used)
          res.pb(i);
      }
      dprintf(clifd, "[");
      for (auto i: res) {
        Json::Value json;
        json["lat"] = places[i].lat;
        json["lng"] = places[i].lng;
        json["addr"] = places[i].addr;
        json["name"] = places[i].name;
        json["pcode"] = places[i].pcode;
        string r = writer.write(json);
        dprintf(clifd, "%s%s", first ? "" : ",", r.c_str());
        first = false;
      }
      dprintf(clifd, "]");
err:
      fclose(f);
    }
    close(clifd);
  }

  //for (auto i: res)
  //  printf("%lf %lf: %s\n", places[i].lat, places[i].lng, places[i].addr.c_str());
}
