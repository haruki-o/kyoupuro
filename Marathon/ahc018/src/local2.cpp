#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef pair<ll, ll> P;
typedef pair<P, ll> Pll;
typedef pair<ll, P> llP;
typedef tuple<ll, ll, ll> T;
typedef tuple<ll, ll, ll, ll> F;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vvi> vvvi;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
typedef vector<vvvll> vvvvll;
typedef vector<vvvvll> vvvvvll;
typedef vector<vector<double>> vvd;
typedef vector<double> vd;
typedef vector<P> vP;
typedef vector<vP> vvP;
typedef vector<T> vT;
typedef vector<vT> vvT;
typedef vector<F> vF;
typedef vector<char> vc;
typedef vector<vc> vvc;
typedef vector<string> vs;
typedef vector<vs> vvs;
#define rep(i, l, n) for (ll i = (ll)(l); i < (ll)(n); i++)
#define repd(i, n, l) for (ll i = (ll)(n); i > (ll)(l); i--)
#define Copy(from, to) copy(from.begin(), from.end(), to.begin())
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
#define Unique(a) sort(a.begin(), a.end())
#define dPQ priority_queue<ll, vll, greater<ll>>
#define PQ priority_queue<ll, vll>
template <class T, class S> inline bool chmax(T &a, const S &b) {
  return (a < b ? a = b, 1 : 0);
}
template <class T, class S> inline bool chmin(T &a, const S &b) {
  return (a > b ? a = b, 1 : 0);
}
#define INF INT_MAX
#define INFF (9223372036854775800)
#define TIME_LIMIT (5.96)
#define def (10101010)
// #define MOD (1000000007)
#define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision(10) << << endl;
// 偏角ソートはlong ddouble!
// auto ite = s.lower_bound("B");
// cout << __builtin_popcount(N) << endl;

// 下,右,上,左
vll dy = {1, 0, -1, 0};
vll dx = {0, 1, 0, -1};
ll N, W, K, C;
vll a, b, c, d;
vvll S;
ll score = 0;
// 四分位数は[100,500,1000]
vll predict_boder = {100, 200, 200, 300};
vvll breaked_rock(200, vll(200, 0));

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
    if (par[c] == c)
      return c;
    ll a = root(par[c]);
    par[c] = a;
    return a;
  }
  // aをbに結合(違う時)
  bool unite(ll a, ll b) {
    ll ra = root(a);
    ll rb = root(b);
    if (ra == rb)
      return 0;
    par[ra] = rb;
    sum[rb] += sum[ra];
    return 1;
  }
  bool same(ll a, ll b) {
    if (root(a) == root(b))
      return 1;
    return 0;
  }
};

ll local_predict_S(ll _y, ll _x) {
  ll ret = 0;
  for (ll border : predict_boder) {
    if (breaked_rock[_y][_x] == 1)
      return breaked_rock[_y][_x];
    cout << _y << " " << _x << " " << border << endl;
    cout.flush();
    ret += border;
    score += border + C;

    ll r;
    S[_y][_x] -= border;
    if (S[_y][_x] <= 0) {
      r = 1;
    } else
      r = 0;

    if (r == 1) {
      breaked_rock[_y][_x] = ret;
      return ret;
    }
  }

  return ret;
}

ll predict_S(ll _y, ll _x) {
  ll ret = 0;
  for (ll border : predict_boder) {
    if (breaked_rock[_y][_x] != 0)
      return breaked_rock[_y][_x];
    cout << _y << " " << _x << " " << border << endl;
    cout.flush();
    ret += border;
    score += border + C;

    ll r;
    cin >> r;
    if (r == -1 || r == 2)
      return -1;

    if (r == 1) {
      breaked_rock[_y][_x] = ret;
      return ret;
    }
  }

  return ret;
}

struct Graph {
  ll N;
  vvll g;
  vll y, x;
  vll houses;
  vll sw;
  map<ll, ll> ma_h;
  map<ll, ll> ma_s;

  void only_Graph() {
    N = W + K;
    g.resize(N);

    for (ll _a : a)
      y.push_back(_a);
    for (ll _b : b)
      x.push_back(_b);
    for (ll _c : c)
      y.push_back(_c);
    for (ll _d : d)
      x.push_back(_d);

    rep(i, 0, N) {
      rep(j, 0, W) {
        if (a[j] == y[i] && b[j] == x[i]) {
          sw.push_back(i);
          ma_s[j] = i;
        }
      }
      rep(j, 0, K) {
        if (c[j] == y[i] && d[j] == x[i]) {
          houses.push_back(i);
          ma_h[j] = i;
        }
      }
    }

    rep(i, 0, N) { rep(j, 0, N) if (i != j) g[i].push_back(j); }
  }

  void lattice_Graph() {
    set<ll> se_y, se_x;
    rep(i, 0, W) se_y.insert(a[i]);
    rep(i, 0, W) se_x.insert(b[i]);
    rep(i, 0, K) se_y.insert(c[i]);
    rep(i, 0, K) se_x.insert(d[i]);
    N = (ll)se_y.size() * (ll)se_x.size();
    g.resize(N);

    for (auto cu_y : se_y) {
      for (auto cu_x : se_x) {
        y.push_back(cu_y);
        x.push_back(cu_x);
      }
    }

    rep(i, 0, N) {
      rep(j, 0, W) {
        if (a[j] == y[i] && b[j] == x[i]) {
          sw.push_back(i);
          ma_s[j] = i;
        }
      }
      rep(j, 0, K) {
        if (c[j] == y[i] && d[j] == x[i]) {
          houses.push_back(i);
          ma_h[j] = i;
        }
      }
    }

    rep(i, 0, N) {
      auto ite = se_y.find(y[i]);
      if (ite != se_y.begin()) {
        g[i].push_back(i - (ll)se_x.size());
      }
      ite++;
      if (ite != se_y.end())
        g[i].push_back(i + (ll)se_x.size());

      ite = se_x.find(x[i]);
      if (ite != se_x.begin()) {
        g[i].push_back(i - 1);
      }
      ite++;
      if (ite != se_x.end())
        g[i].push_back(i + 1);
    }
  }

  void dijkstra(ll sta_sw, vll _houses, vP &breack_rocks) {
    vll dp(N, -1);
    vi sub(N, 0);
    sub[ma_s[sta_sw]] = 1;
    for (ll house : _houses)
      sub[ma_h[house]] = 1;
    vvll tree_G(N);
    priority_queue<T, vT, greater<T>> pq;
    pq.push({(ll)0, ma_s[sta_sw], -1});
    while (!pq.empty()) {
      auto [cu_cost, cu, origin] = pq.top();
      pq.pop();
      if (dp[cu] != -1)
        continue;
      dp[cu] = cu_cost;
      if (origin != -1) {
        tree_G[cu].push_back(origin);
        tree_G[origin].push_back(cu);
      }
      for (ll to : g[cu]) {
        if (dp[to] != -1)
          continue;
        if (sub[to] == 0)
          continue;
        ll dist = abs(y[cu] - y[to]) + abs(x[cu] - x[to]);
        pq.push({dist, to, cu});
      }
    }

    map<P, int> breaked;
    vi seen(N, 0);
    queue<ll> qu;
    qu.push(ma_s[sta_sw]);
    seen[ma_s[sta_sw]] = 1;
    breack_rocks.push_back({y[ma_s[sta_sw]], x[ma_s[sta_sw]]});
    breaked[{y[ma_s[sta_sw]], x[ma_s[sta_sw]]}] = 1;
    while (!qu.empty()) {
      ll cu = qu.front();
      qu.pop();
      for (ll to : tree_G[cu]) {
        if (seen[to] == 1)
          continue;
        seen[to] = 1;
        qu.push(to);
        ll cu_y = y[cu];
        ll cu_x = x[cu];
        ll _dy = y[to] - cu_y;
        ll _dx = x[to] - cu_x;

        ll pred1 = 0, pred2 = 0;
        if (_dx != 0)
          pred1 += local_predict_S(cu_y, cu_x + _dx / 2);
        if (_dy != 0)
          pred1 += local_predict_S(cu_y + _dy / 2, cu_x + _dx);
        if (_dy != 0)
          pred2 += local_predict_S(cu_y + _dy / 2, cu_x);
        if (_dx != 0)
          pred2 += local_predict_S(cu_y + _dy, cu_x + _dx / 2);

        // if (_dx != 0)
        //   pred1 += predict_S(cu_y, cu_x + _dx / 2);
        // if (_dy != 0)
        //   pred1 += predict_S(cu_y + _dy / 2, cu_x + _dx);
        // if (_dy != 0)
        //   pred2 += predict_S(cu_y + _dy / 2, cu_x);
        // if (_dx != 0)
        //   pred2 += predict_S(cu_y + _dy, cu_x + _dx / 2);
        if (pred1 > pred2) {
          rep(i, 0, (ll)abs(_dy)) {
            cu_y += _dy / abs(_dy);
            if (!breaked[{cu_y, cu_x}])
              breack_rocks.push_back({cu_y, cu_x});
            breaked[{cu_y, cu_x}] = 1;
          }
          rep(i, 0, (ll)abs(_dx)) {
            cu_x += _dx / abs(_dx);
            if (!breaked[{cu_y, cu_x}])
              breack_rocks.push_back({cu_y, cu_x});
            breaked[{cu_y, cu_x}] = 1;
          }
        } else {
          rep(i, 0, (ll)abs(_dx)) {
            cu_x += _dx / abs(_dx);
            if (!breaked[{cu_y, cu_x}])
              breack_rocks.push_back({cu_y, cu_x});
            breaked[{cu_y, cu_x}] = 1;
          }
          rep(i, 0, (ll)abs(_dy)) {
            cu_y += _dy / abs(_dy);
            if (!breaked[{cu_y, cu_x}])
              breack_rocks.push_back({cu_y, cu_x});
            breaked[{cu_y, cu_x}] = 1;
          }
        }
      }
    }
  }

  void prim(ll sta_sw, vll _houses, vP &breack_rocks) {
    // cout << "sw : " << ma_s[sta_sw] << endl;
    // for (ll house : _houses)
    //   cout << ma_h[house] << " ";
    // cout << endl;
    vvll tree_G(N);
    UnionFind uf(N);
    vT at;
    rep(i, 0, N) {
      for (ll cu : g[i]) {
        if (cu < i)
          continue;
        ll dist = abs(y[i] - y[cu]) + abs(x[i] - x[cu]);
        at.push_back({dist, i, cu});
      }
    }
    Sort(at);
    for (auto cu : at) {
      auto [cost, u, v] = cu;
      if (uf.unite(u, v)) {
        tree_G[u].push_back(v);
        tree_G[v].push_back(u);
      }
    }
    // cout << N << " ";
    // ll _sum = 0;
    // rep(j, 0, N) _sum += (ll)tree_G[j].size();
    // _sum /= 2;
    // cout << _sum << endl;

    // 無駄な枝をそぐ
    queue<ll> qu;
    rep(i, 0, N) {
      if ((ll)tree_G[i].size() == 1) {
        int flag = 1;
        for (ll house : _houses) {
          if (ma_h[house] == i)
            flag = 0;
        }
        if (ma_s[sta_sw] == i)
          flag = 0;
        if (!flag)
          continue;
        qu.push(i);
      }
    }
    while (!qu.empty()) {
      ll cu = qu.front();
      qu.pop();
      // cout << cu << " ";
      for (ll del : tree_G[cu]) {
        for (auto ite = tree_G[del].begin(); ite != tree_G[del].end(); ite++) {
          if (*ite == cu) {
            tree_G[del].erase(ite);
            break;
          }
        }
        if ((ll)tree_G[del].size() == 1) {
          int flag = 1;
          for (ll house : _houses) {
            if (ma_h[house] == del)
              flag = 0;
          }
          if (ma_s[sta_sw] == del)
            flag = 0;
          if (!flag)
            continue;
          qu.push(del);
        }
      }
      tree_G[cu].erase(tree_G[cu].begin());
    }

    // cout << N << " ";
    // ll _sum = 0;
    // rep(j, 0, N) _sum += (ll)tree_G[j].size();
    // _sum /= 2;
    // cout << _sum << endl;

    // 出力
    map<P, int> breaked;
    vi seen(N, 0);
    qu.push(ma_s[sta_sw]);
    seen[ma_s[sta_sw]] = 1;
    breack_rocks.push_back({y[ma_s[sta_sw]], x[ma_s[sta_sw]]});
    breaked[{y[ma_s[sta_sw]], x[ma_s[sta_sw]]}] = 1;
    while (!qu.empty()) {
      ll cu = qu.front();
      qu.pop();
      for (ll to : tree_G[cu]) {
        if (seen[to] == 1)
          continue;
        seen[to] = 1;
        qu.push(to);
        ll cu_y = y[cu];
        ll cu_x = x[cu];
        ll _dy = y[to] - cu_y;
        ll _dx = x[to] - cu_x;
        rep(i, 0, (ll)abs(_dy)) {
          cu_y += _dy / abs(_dy);
          if (!breaked[{cu_y, cu_x}])
            breack_rocks.push_back({cu_y, cu_x});
          breaked[{cu_y, cu_x}] = 1;
        }
        rep(i, 0, (ll)abs(_dx)) {
          cu_x += _dx / abs(_dx);
          if (!breaked[{cu_y, cu_x}])
            breack_rocks.push_back({cu_y, cu_x});
          breaked[{cu_y, cu_x}] = 1;
        }
      }
    }
  }
};

void local_solve_mining(vP &break_rocks) {
  ll set_power = 50;
  ll all_mine_sum = 0;
  for (auto rock : break_rocks) {
    if (breaked_rock[rock.first][rock.second] != 0)
      continue;
    breaked_rock[rock.first][rock.second] = 1;
    ll mine_sum = 0;
    while (1) {
      cout << rock.first << " " << rock.second << " " << set_power << endl;
      cout.flush();
      score += set_power + C;
      all_mine_sum++;
      ll r;
      S[rock.first][rock.second] -= set_power;
      if (S[rock.first][rock.second] <= 0) {
        r = 1;
      } else
        r = 0;

      if (r == 0) {
        mine_sum++;
      }
      if (r == 1) {
        if (mine_sum == 0)
          set_power /= 3;
        else
          set_power *= (mine_sum);
        chmax(set_power, (ll)10);
        chmin(set_power, (ll)500);
        break;
      }
    }
  }
}

void solve_mining(vP &break_rocks) {
  ll set_power = 50;
  for (auto rock : break_rocks) {
    if (breaked_rock[rock.first][rock.second] != 0)
      continue;
    breaked_rock[rock.first][rock.second] = 1;
    ll mine_sum = 0;
    while (1) {
      cout << rock.first << " " << rock.second << " " << set_power << endl;
      cout.flush();

      ll r;
      cin >> r;

      if (r == -1 || r == 2)
        return;

      if (r == 0) {
        mine_sum++;
      }
      if (r == 1) {
        if (mine_sum == 0)
          set_power /= 3;
        else
          set_power *= (mine_sum);
        chmax(set_power, (ll)10);
        chmin(set_power, (ll)500);
        break;
      }
    }
  }
}

int main() {
  cin >> N >> W >> K >> C;

  S.resize(N);
  rep(i, 0, N) {
    S[i].resize(N);
    rep(j, 0, N) cin >> S[i][j];
  }

  a.resize(W);
  b.resize(W);
  c.resize(K);
  d.resize(K);
  rep(i, 0, W) cin >> a[i] >> b[i];
  rep(i, 0, K) cin >> c[i] >> d[i];

  vvll house_W(W);
  vi selected_house(K, 0);
  rep(i, 0, K) {
    P _best = {0, 0};
    ll best_dist = 1e10;
    // house jを選ぶ
    rep(j, 0, K) {
      if (selected_house[j] == 1)
        continue;
      ll mi = 1e10;
      ll best_W = 0;
      // 水源kに所属している家・水源kからの距離を調べる
      rep(k, 0, W) {
        // 水源kから
        ll _dist = abs(a[k] - c[j]) + abs(b[k] - d[j]);
        if (chmin(mi, _dist))
          best_W = k;
        // 水源kに所属する家
        for (ll house : house_W[k]) {
          _dist = abs(c[house] - c[j]) + abs(d[house] - d[j]);
          if (chmin(mi, _dist))
            best_W = k;
        }
      }
      if (chmin(best_dist, mi))
        _best = {j, best_W};
    }
    house_W[_best.second].push_back(_best.first);
    selected_house[_best.first] = 1;
  }

  Graph G;
  G.only_Graph();
  // G.lattice_Graph();
  vP break_rocks;
  // rep(i, 0, W) G.prim(i, house_W[i], break_rocks);
  rep(i, 0, W) G.dijkstra(i, house_W[i], break_rocks);

  local_solve_mining(break_rocks);
  // solve_mining(break_rocks);

  // cerr << score << endl;
  cerr << W << " " << K << " " << C << endl;
}

// local2.cpp
// 縦と横のどちらを先に選ぶか
// local1.cppと比較して1.08倍