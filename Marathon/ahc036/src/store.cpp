#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
typedef pair<ll, ll> pll;
typedef pair<pll, ll> ppll;
typedef vector<ppll> vpp;
typedef vector<pll> vp;
typedef vector<vp> vvp;
typedef vector<char> vc;
typedef vector<int> vi;
typedef vector<string> vs;
typedef vector<vs> vvs;
#define rep(i, l, n) for (ll i = (ll)(l); i < (ll)(n); i++)
#define repd(i, n, l) for (ll i = (ll)(n); i > (ll)(l); i--)
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
template <class T, class S>
inline bool chmax(T &a, const S &b) {
  return (a < b ? a = b, 1 : 0);
}
template <class T, class S>
inline bool chmin(T &a, const S &b) {
  return (a > b ? a = b, 1 : 0);
}
#define INF INT_MAX
#define INFF (9223372036854775800)
#define TIME_LIMIT (3.5)
#define def (20101)
// #define MOD (1000000007)
// #define MAX (2147483647)
#define MAX (1073741823)
#define MOD (998244353)

ll N, M, T, LA, LB;
vll u, v;
vll t;
vll x, y;
vvll g, gMat;

struct UnionFind {
  // 根のとき-1,iの親はpar[i]
  vi par;
  vll sum;
  UnionFind(ll N) {
    par.resize(N);
    sum.assign(N, 1);
    rep(i, 0, N) par[i] = i;
  }
  // cの根っこを返す
  ll root(ll c) {
    if (par[c] == c) return c;
    ll a = root(par[c]);
    par[c] = a;
    return a;
  }
  // aをbに結合(違う時)
  bool unite(ll a, ll b) {
    ll ra = root(a);
    ll rb = root(b);
    if (ra == rb) return 0;
    par[ra] = rb;
    sum[rb] += sum[ra];
    return 1;
  }
  bool same(ll a, ll b) {
    if (root(a) == root(b)) return 1;
    return 0;
  }
};

class Event {
 public:
  char type;
  ll l, PA, PB;
  ll v;

  void setSignal(ll l, ll PA, ll PB) {
    this->type = 's';
    cin >> this->l >> this->PA >> this->PB;
  }

  void setMove(vll v) {
    this->type = 'm';
    cin >> this->v;
  }

  void printSignal() {
    cerr << this->type << " " << this->l << " ";
    cerr << this->PA << " " << this->PB << endl;
  }

  void printMove() {
    cerr << this->type << " " << this->v << endl;
  }
};

struct PressG {
  vvll G;
  map<ll, vll> clusterToVertices;

  // 圧縮後の頂点c1,c2が接続しているとき、圧縮前の接続する頂点を返す。
  pll edgeComposeVertice(ll c1, ll c2) {
    pll ret = {-1, -1};
    for (ll v1 : clusterToVertices[c1]) {
      for (ll v2 : clusterToVertices[c2]) {
        if (gMat[v1][v2] == 1) {
          ret = {v1, v2};
          return ret;
        }
      }
    }
    return ret;
  }

  void constructPressG(vll &A) {
    vll ok;
    rep(i, 0, A.size() - LB + 1) {
      map<ll, ll> ma;
      rep(j, i, i + LB) ma[A[j]] = 1;
      queue<ll> qu;
      qu.push(A[i]);
      vll seen(N, 0);
      seen[A[i]] = 1;
      while (!qu.empty()) {
        ll cu = qu.front();
        qu.pop();
        for (ll to : g[cu]) {
          if (ma.count(to) == 0) continue;
          if (seen[to] == 1) continue;
          qu.push(to);
          seen[to] = 1;
        }
      }
      bool flag = 1;
      rep(j, i, i + LB) if (seen[A[j]] == 0) flag = 0;
      if (flag) ok.push_back(i);
    }

    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    std::shuffle(ok.begin(), ok.end(), engine);

    cerr << "ok.size() : " << ok.size() << endl;
    G.resize((ll)min((ll)ok.size(), (ll)100));
    rep(i, 0, min((ll)ok.size(), (ll)100)) {
      rep(j, ok[i], ok[i] + LB) {
        clusterToVertices[i].push_back(A[j]);
        cerr << A[j] << " ";
      }
      cerr << endl;
    }
    rep(i, 0, G.size()) {
      rep(j, i + 1, G.size()) {
        pll _ = edgeComposeVertice(i, j);
        if (_.first != -1) {
          G[i].push_back(j);
          G[j].push_back(i);
        }
      }
    }
  }

};

struct Solver {
  void dfs(ll v, vvll &tree, vll &visited) {
    visited[v] = 1;
    for (ll c : g[v]) {
      if (visited[c] == 1) continue;
      visited[c] = 1;
      tree[v].push_back(c);
      tree[c].push_back(v);
      dfs(c, tree, visited);
    }
  }

  void constructTree(vvll &tree) {
    vll visited(N, 0);
    dfs(0, tree, visited);
  }

  void constructADfs(ll v, vll &A, vvll &tree, vll &visited) {
    // cerr << v << " ";
    visited[v] = 1;
    A.push_back(v);
    for (ll c : g[v]) {
      if (visited[c] == 1) continue;
      visited[c] = 1;
      constructADfs(c, A, tree, visited);
    }
  }

  void initConstructA(vll &A, vvll &tree) {
    vll visited(N, 0);
    constructADfs(0, A, tree, visited);
  }

  void _initAns(vll &A, vector<Event> &ans) {
    vvll gMat(N, vll(N, 0));
    rep(i, 0, N) {
      for (ll c : g[i]) {
        gMat[i][c] = 1, gMat[c][i] = 1;
      }
    }
    vvll dp(1e4, vll(N, -2));
    dp[0][0] = -1;
    vvp track(1e4, vp(N, {-1, -1}));
    rep(i, 0, 1e5 - 1) {
      rep(j, 0, N) {
        if (dp[i][j] == -2) continue;
        if (j % 10 == 0) cerr << "i,j" << i << " " << j;
        // 配るDB
        // dp[i][j] := 信号操作i回終了後、頂点jにいて、t[]の最大インデックス
        // 頂点jにいて、A[k],...,A[k+LB-1]をBに適用したとき、dp[i+1][連結する頂点]に分配する。
        rep(k, 0, A.size()) {
          UnionFind uf(LB + 1);
          // {j, A[k],...,A[k+LB-1]}
          vll sub(LB + 1, j);
          rep(l, 1, LB + 1) sub[l] = A[k + l - 1];
          map<ll, ll> ma;
          rep(l, 0, LB + 1) ma[sub[l]] = l;
          rep(_i, 0, LB + 1) {
            rep(_j, _i + 1, LB + 1) {
              uf.unite(ma[sub[_i]], ma[sub[_j]]);
            }
          }
          vll isUnite(LB, 0);
          rep(l, 0, LB) {
            if (uf.same(0, l + 1)) isUnite[l] = 1;
          }
          // t[]のインデックスがどんだけ進めるか
          ll proc = dp[i][j];
          while (1) {
            ll c = t[proc + 1];
            // t[]の次の値が、B[]に入っているかつ頂点jと繋がっている
            if (0 < ma.count(c) && uf.same(0, ma[c])) {
              proc++;
              if (proc == T - 1) break;
            } else break;
          }
          rep(l, 0, LB + 1) {
            if (dp[i + 1][sub[l]] < proc) {
              dp[i + 1][sub[l]] = proc;
              track[i + 1][sub[l]] = {j, k};
            }
          }

          // cerr << "i, j, k, dp[i][];" << i << " " << j << " " << k << " " << dp[i][j] << endl;
          // rep(_j, 0, N) cerr << dp[i + 1][_j] << " ";
          // cerr << endl;
        }
      }

      int flag = 0;
      rep(j, 0, N) if (dp[i + 1][j] == T - 1) flag = 1;
      cerr << "i" << i << " ";

      if (flag) cerr << i << endl;
      if (flag) break;
    }
  }

  void initAns(vll &A, vector<Event> &ans) {}

  void solve(vector<Event> &ans) {
    vvll tree(N);
    constructTree(tree);

    vll A;
    initConstructA(A, tree);
    cerr << "A[] : ";
    rep(i, 0, LA) cerr << A[i] << " ";
    cerr << endl;

    PressG pressG;
    pressG.constructPressG(A);

    // initAns(A, ans);
  }
};

void print_ans(vector<Event> &ans) {
  rep(i, 0, ans.size()) {
    if (ans[i].type == 's') ans[i].printSignal();
    else ans[i].printMove();
  }
}

int main() {
  cin >> N >> M >> T >> LA >> LB;

  u.resize(M);
  v.resize(M);
  rep(i, 0, M) cin >> u[i] >> v[i];
  t.resize(T);
  rep(i, 0, T) cin >> t[i];
  x.resize(N);
  y.resize(N);
  rep(i, 0, N) cin >> x[i] >> y[i];
  g.resize(N);
  rep(i, 0, M) {
    g[u[i]].push_back(v[i]);
    g[v[i]].push_back(u[i]);
  }
  gMat.resize(N);
  rep(i, 0, N) gMat[i].assign(N, 0);
  rep(i, 0, N) {
    for (ll cu : g[i]) {
      gMat[i][cu] = 1;
      gMat[cu][i] = 1;
    }
  }

  vector<Event> ans;
  Solver solver;
  solver.solve(ans);

  print_ans(ans);
}

// 8/24(土)
// constructPressG -> initConstrcutAの順に処理するため一旦保存