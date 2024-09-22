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
#define INF INT_MAX
#define INFF (9223372036854775800)
#define TIME_LIMIT (5.95)
#define def (10101010)
// #define MOD (1000000007)
#define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision(10) << << endl;
// 偏角ソートはlong ddouble!
// auto ite = s.lower_bound("B");
// cout << __builtin_popcount(N) << endl;

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

ll calc_init_all_cost(vvll &G, map<P, ll> &G_cost) {
  ll sum = 0;
  ll N = G.size();
  // 始点iのダイクストラ
  rep(i, 0, N) {
    vll dp(N, -1);
    priority_queue<P, vP, greater<P>> pq;
    pq.push({0, i});
    while (!pq.empty()) {
      auto [cu_cost, cu] = pq.top();
      pq.pop();
      if (dp[cu] != -1)
        continue;
      dp[cu] = cu_cost;
      // if (i < cu)
      sum += cu_cost;
      for (ll to : G[cu]) {
        ll l = min(cu, to);
        ll r = max(cu, to);
        if (dp[to] != -1)
          continue;
        pq.push({cu_cost + G_cost[{l, r}], to});
      }
    }
  }
  return sum;
}

ll calc_diff_cost(vvll &G, map<P, ll> &G_cost, vll &r, map<P, ll> &r_idx, vll f,
                  ll all_d) {
  ll D = (ll)f.size();
  ll N = (ll)G.size();
  // d日は使わないとき
  rep(d, 0, D) {
    auto startClock = system_clock::now();
    // 始点iのダイクストラ
    rep(i, 0, N) {
      vll dp(N, -1);
      priority_queue<P, vP, greater<P>> pq;
      pq.push({0, i});
      while (!pq.empty()) {
        auto [cu_cost, cu] = pq.top();
        pq.pop();
        if (dp[cu] != -1)
          continue;
        dp[cu] = cu_cost;
        for (ll to : G[cu]) {
          ll l = min(cu, to);
          ll r = max(cu, to);
          if (r_idx[{l, r}] == d)
            continue;
          if (dp[to] != -1)
            continue;
          pq.push({cu_cost + G_cost[{l, r}], to});
        }
      }
      rep(j, 0, N) f[d] += dp[j] == -1 ? 1e9 : dp[j];
      // rep(j, 0, N) if (dp[j] == -1) cerr << "a";
    }
  }
  // rep(i, 0, D) cerr << "day : " << i << "  f_k : " << (f[i] - all_d) / (N *
  // (N - 1)) << endl;
  ll ret = 0;
  rep(i, 0, D) ret += (f[i] - all_d) / (N * (N - 1));
  ret = 1e3 * ret / D;
  return ret;
}

void init_r(vll &r, ll M, ll D, ll K) {
  ll id = M / D + 3;
  r.assign(M, D - 1);
  rep(i, 0, D - 1) {
    rep(j, 0, id) { r[min(i * id + j, M - 1)] = i; }
  }
}

void modify(vll &r, map<P, ll> &r_idx, vll &u, vll &v, vvll &G, ll D,
            vll &fix) {
  ll M = (ll)r.size();
  vll fix_idx;
  vll shuffle_r;
  rep(i, 0, M) {
    if (fix[r[i]] == 0) {
      fix_idx.push_back(i);
      shuffle_r.push_back(r[i]);
    }
  }
  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());
  std::shuffle(shuffle_r.begin(), shuffle_r.end(), engine);
  rep(i, 0, (ll)shuffle_r.size()) { r[fix_idx[i]] = shuffle_r[i]; }
  rep(i, 0, M) r_idx[{u[i], v[i]}] = r[i];
}

bool is_connect(vvll &G, vll &r, map<P, ll> &r_idx, ll D) {
  ll N = (ll)G.size();
  vll dp(N, -1);
  dp[0] = 0;
  queue<ll> qu;
  qu.push(0);
  ll sum = 0;
  while (!qu.empty()) {
    ll cu = qu.front();
    qu.pop();
    for (ll to : G[cu]) {
      ll mi = min(to, cu);
      ll ma = max(to, cu);
      if (dp[to] != -1)
        continue;
      if (r_idx[{mi, ma}] == D)
        continue;
      qu.push(to);
      dp[to] = 0;
      sum++;
    }
  }
  if (N - 1 == sum)
    return true;

  return false;
}

void modify_swap(vll &r, map<P, ll> &r_idx, vll &u, vll &v, ll &M1, ll &M2,
                 ll &D1, ll &D2, ll flag) {
  ll M = (ll)r.size();
  if (flag == 1) {
    while (1) {
      M1 = rand() % M;
      M2 = (M1 + rand() % M) % M;
      D1 = r_idx[{u[M1], v[M1]}];
      D2 = r_idx[{u[M2], v[M2]}];
      if (D1 != D2)
        break;
    }
    swap(r[M1], r[M2]);
    r_idx[{u[M1], v[M1]}] = r[M1];
    r_idx[{u[M2], v[M2]}] = r[M2];
  } else {
    swap(r[M1], r[M2]);
    r_idx[{u[M1], v[M1]}] = r[M1];
    r_idx[{u[M2], v[M2]}] = r[M2];
  }
}

// 始点sta[] のday Dのコスト
ll calc_d_cost(vvll &G, map<P, ll> &G_cost, vll &r, map<P, ll> &r_idx, ll D,
               vP sta) {
  ll N = (ll)G.size();
  ll ret = 0;
  priority_queue<P, vP, greater<P>> pq;
  for (auto cu : sta) {
    ll i = cu.first;
    ll p_sum = cu.second;
    vll dp(N, -1);
    pq.push({0, i});
    while (!pq.empty()) {
      auto [cu_cost, cu] = pq.top();
      pq.pop();
      if (dp[cu] != -1)
        continue;
      dp[cu] = cu_cost;
      ret += cu_cost * p_sum;
      for (ll to : G[cu]) {
        ll l = min(cu, to);
        ll r = max(cu, to);
        if (r_idx[{l, r}] == D)
          continue;
        if (dp[to] != -1)
          continue;
        pq.push({cu_cost + G_cost[{l, r}], to});
      }
    }
  }

  return ret;
}

void select_start_point(vP &start_idx, vll &x, vll &y) {
  ll N = (ll)x.size();
  ll cen0_idx = 0;
  ll cen1_idx = 0;
  ll cen2_idx = 0;
  ll mi0_center = 1e9;
  ll mi1_center = 1e9;
  ll mi2_center = 1e9;
  rep(i, 0, N) {
    if (abs(500 - x[i]) + abs(800 - y[i]) < mi0_center) {
      mi0_center = abs(500 - x[i]) + abs(800 - y[i]);
      cen0_idx = i;
    }
    if (abs(250 - x[i]) + abs(350 - y[i]) < mi1_center) {
      mi1_center = abs(250 - x[i]) + abs(350 - y[i]);
      cen1_idx = i;
    }
    if (abs(750 - x[i]) + abs(350 - y[i]) < mi2_center) {
      mi2_center = abs(750 - x[i]) + abs(350 - y[i]);
      cen2_idx = i;
    }
  }
  ll sum0 = 0;
  ll sum1 = 0;
  ll sum2 = 0;
  rep(i, 0, N) {
    ll d = abs(x[cen0_idx] - x[i]) * abs(x[cen0_idx] - x[i]);
    d += abs(x[cen0_idx] - y[i]) * abs(x[cen0_idx] - y[i]);
    if (d < (ll)200000)
      sum0++;

    d = abs(x[cen1_idx] - x[i]) * abs(x[cen1_idx] - x[i]);
    d += abs(x[cen1_idx] - y[i]) * abs(x[cen1_idx] - y[i]);
    if (d < (ll)200000)
      sum1++;
    d = abs(x[cen2_idx] - x[i]) * abs(x[cen2_idx] - x[i]);
    d += abs(x[cen2_idx] - y[i]) * abs(x[cen2_idx] - y[i]);
    if (d < (ll)200000)
      sum2++;
  }
  start_idx.push_back({cen0_idx, sum0});
  start_idx.push_back({cen1_idx, sum1});
  start_idx.push_back({cen2_idx, sum2});
  for (auto cu : start_idx)
    cerr << "{" << x[cu.first] << " " << y[cu.first] << "}, " << cu.second
         << endl;
}

int main() {
  srand((unsigned int)time(NULL));
  auto startClock = system_clock::now();
  ll N, M, D, K;
  cin >> N >> M >> D >> K;
  vvll G(N);
  vll u(M), v(M), w(M);
  map<P, ll> G_cost;
  rep(i, 0, M) {
    cin >> u[i] >> v[i] >> w[i];
    u[i]--;
    v[i]--;
    w[i]--;
    G[u[i]].push_back(v[i]);
    G[v[i]].push_back(u[i]);
    G_cost[{u[i], v[i]}] = w[i];
  }
  vll x(N), y(N);
  rep(i, 0, N) cin >> x[i] >> y[i];

  vP start_idx;
  select_start_point(start_idx, x, y);

  vll r(M);
  map<P, ll> r_idx;
  vll f(D, 0);

  init_r(r, M, D, K);
  rep(i, 0, M) r_idx[{u[i], v[i]}] = r[i];

  ll turn = 0;
  ll update_turn = 0;
  vll fix(D, 0);
  ll ju = 0;
  // 連結になるように振り分ける
  while (1) {
    turn++;
    const double time =
        duration_cast<microseconds>(system_clock::now() - startClock).count() *
        1e-6;
    if (time > TIME_LIMIT)
      break;
    vll new_r(M);
    new_r = r;
    map<P, ll> new_r_idx;
    rep(i, 0, M) r_idx[{u[i], v[i]}] = r[i];
    modify(new_r, new_r_idx, u, v, G, D, fix);
    vll new_fix(D, 0);
    rep(i, 0, D) new_fix[i] = is_connect(G, new_r, new_r_idx, i);

    ll fix_sum = 0;
    ll new_fix_sum = 0;
    for (ll cu : fix)
      fix_sum += cu;
    for (ll cu : new_fix)
      new_fix_sum += cu;
    if (new_fix_sum == D - 1)
      continue;
    if (fix_sum < new_fix_sum) {
      r = new_r;
      r_idx = new_r_idx;
      fix = new_fix;
      update_turn = turn;
    }

    if (new_fix_sum == D) {
      ju = 1;
      break;
    }
  }

  if (ju == 0) {
    rep(i, 0, M) cout << r[i] + 1 << " ";
    return 0;
  }

  // cerr
  //     << duration_cast<microseconds>(system_clock::now() -
  //     startClock).count() *
  //            1e-6
  //     << " " << turn << endl;

  // ll d = calc_init_all_cost(G, G_cost);
  // cerr << calc_diff_cost(G, G_cost, r, r_idx, f, d) << " " << D << endl;

  // 2辺のみ交換して、連結かどうか判断し、頂点1からの距離の最小化を目指す
  // startClock = system_clock::now();
  vll d_cost(D, 0);
  rep(i, 0, D) d_cost[i] = calc_d_cost(G, G_cost, r, r_idx, i, start_idx);
  turn = 0;
  while (1) {
    turn++;
    const double time =
        duration_cast<microseconds>(system_clock::now() - startClock).count() *
        1e-6;
    if (time > TIME_LIMIT)
      break;
    // D1とD2の辺が交換される
    ll D1, D2;
    ll M1, M2;
    modify_swap(r, r_idx, u, v, M1, M2, D1, D2, 1);

    // 交換しても連結してるか
    bool connect_flag = is_connect(G, r, r_idx, D1);
    connect_flag &= is_connect(G, r, r_idx, D2);

    if (!connect_flag) {
      modify_swap(r, r_idx, u, v, M1, M2, D1, D2, 0);
      continue;
    }

    ll update_cost = 0;
    ll D1_cost = calc_d_cost(G, G_cost, r, r_idx, D1, start_idx);
    ll D2_cost = calc_d_cost(G, G_cost, r, r_idx, D2, start_idx);
    update_cost += d_cost[D1] - D1_cost;
    update_cost += d_cost[D2] - D2_cost;

    // if (update_cost)
    //   cerr << update_cost << endl;
    if (0 < update_cost) {
      d_cost[D1] = D1_cost;
      d_cost[D2] = D2_cost;
    } else {
      modify_swap(r, r_idx, u, v, M1, M2, D1, D2, 0);
    }
    // break;
  }

  rep(i, 0, M) cout << r[i] + 1 << " ";

  cerr
      << duration_cast<microseconds>(system_clock::now() - startClock).count() *
             1e-6
      << " " << turn << endl;
  // ll d = calc_init_all_cost(G, G_cost);
  // cerr << calc_diff_cost(G, G_cost, r, r_idx, f, d) << " " << D << endl;
}

// 3角形の3点の周辺の個数をかけました