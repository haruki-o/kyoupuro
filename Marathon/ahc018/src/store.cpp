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
vll predict_border;
map<P, int> breaked;
// 予測できない時
ll upper = 1000;
vll stat_P(4, 0);
vll stat_P_border = {0, 100, 500, 1000, 5000};

struct PredictBan {
  vvll ban;
  vll se_y;
  vll se_x;

  // se_y,se_xをもとに更新
  void update_ban() {
    // 格子状
    rep(i, 0, (ll)se_x.size()) {
      rep(j, 0, (ll)se_y.size() - 1) {
        ll val1 = ban[se_y[j]][se_x[i]];
        ll val2 = ban[se_y[j + 1]][se_x[i]];
        long double log_val1 = log10(val1);
        long double log_val2 = log10(val2);
        long double unit = abs(log_val2 - log_val1) / (se_y[j + 1] - se_y[j]);
        rep(k, se_y[j] + 1, se_y[j + 1]) {
          long double _val;
          if (log_val1 < log_val2)
            _val = pow(10, log_val1 + unit * (k - se_y[j]));
          else
            _val = pow(10, log_val1 - unit * (k - se_y[j]));
          ban[k][se_x[i]] = _val;
        }
      }
    }
    rep(i, 0, (ll)se_y.size()) {
      rep(j, 0, (ll)se_x.size() - 1) {
        ll val1 = ban[se_y[i]][se_x[j]];
        ll val2 = ban[se_y[i]][se_x[j + 1]];
        long double log_val1 = log10(val1);
        long double log_val2 = log10(val2);
        long double unit = abs(log_val2 - log_val1) / (se_x[j + 1] - se_x[j]);
        rep(k, se_x[j] + 1, se_x[j + 1]) {
          long double _val;
          if (log_val1 < log_val2)
            _val = pow(10, log_val1 + unit * (k - se_x[j]));
          else
            _val = pow(10, log_val1 - unit * (k - se_x[j]));
          ban[se_y[i]][k] = _val;
        }
      }
    }
    // ここまで格子状
    rep(i, 0, 200) {
      rep(j, 0, (ll)se_y.size() - 1) {
        ll val1 = ban[se_y[j]][i];
        ll val2 = ban[se_y[j + 1]][i];
        long double log_val1 = log10(val1);
        long double log_val2 = log10(val2);
        long double unit = abs(log_val2 - log_val1) / (se_y[j + 1] - se_y[j]);
        rep(k, se_y[j] + 1, se_y[j + 1]) {
          long double _val;
          if (log_val1 < log_val2)
            _val = pow(10, log_val1 + unit * (k - se_y[j]));
          else
            _val = pow(10, log_val1 - unit * (k - se_y[j]));
          ban[k][i] = _val;
        }
      }
    }
  }

  // flag := 0 ほんばん, 1 local
  PredictBan(ll flag) {
    if (C <= 2) {
      rep(i, 0, 20) predict_border.push_back(25);
    } else if (C <= 8) {
      rep(i, 0, 10) predict_border.push_back(50);
    } else
      predict_border = {50, 100, 200, 200};

    ban.resize(200);
    rep(i, 0, 200) ban[i].assign(200, 0);
    rep(_y, 0, 200) if (_y % 30 == 0 || _y == 199) se_y.push_back(_y);
    rep(_x, 0, 200) if (_x % 30 == 0 || _x == 199) se_x.push_back(_x);
    map<P, int> skip;
    ll mi_y = 1e10;
    ll mi_x = 1e10;
    ll ma_y = 0;
    ll ma_x = 0;
    for (ll _y : a) {
      chmin(mi_y, _y);
      chmax(ma_y, _y);
    }
    for (ll _y : c) {
      chmin(mi_y, _y);
      chmax(ma_y, _y);
    }
    for (ll _x : b) {
      chmin(mi_x, _x);
      chmax(ma_x, _x);
    }
    for (ll _x : d) {
      chmin(mi_x, _x);
      chmax(ma_x, _x);
    }
    mi_y -= 20;
    mi_x -= 20;
    ma_y += 20;
    ma_x += 20;
    for (auto _y : se_y) {
      for (auto _x : se_x) {
        if (_y < mi_y || ma_y < _y)
          skip[{_y, _x}] = 1;
        if (_x < mi_x || ma_x < _x)
          skip[{_y, _x}] = 1;
      }
    }
    for (auto _y : se_y) {
      for (auto _x : se_x) {
        if (skip.count({_y, _x})) {
          ban[_y][_x] = upper;
          // continue;
        }
        ll ret = 0;
        for (ll border : predict_border) {
          cout << _y << " " << _x << " " << border << endl;
          cout.flush();
          ret += border;
          score += border + C;
          ll r;
          if (flag == 0) {
            cin >> r;

            if (r == -1 || r == 2)
              return;
          }
          if (flag == 1) {
            S[_y][_x] -= border;
            if (S[_y][_x] <= 0) {
              r = 1;
            } else
              r = 0;
          }

          if (r == 1) {
            breaked[{_y, _x}] = 1;
            ban[_y][_x] = ret;
            break;
          }
          ban[_y][_x] = upper;
        }
      }
    }
    update_ban();
  }

  // se_y,se_xと新たな点のbanを更新
  void update_se(vP &ava) {}
};

void local_unit_mining(ll _y, ll _x, PredictBan &p_ban) {
  if (breaked.count({_y, _x}))
    return;
  ll set_power = min((ll)500, p_ban.ban[_y][_x]);
  if (C <= 8)
    set_power /= 10;
  else if (C <= 16)
    set_power /= 5;
  else if (C <= 64)
    set_power /= 2;
  chmax(set_power, (ll)10);
  ll mine_sum = 1;
  rep(i, 0, (ll)stat_P_border.size() - 1) {
    if (stat_P_border[i] <= S[_y][_x] && S[_y][_x < stat_P_border[i + 1]])
      stat_P[i]++;
  }
  while (1) {
    cout << _y << " " << _x << " " << set_power << endl;
    cout.flush();
    score += set_power + C;
    ll r;
    S[_y][_x] -= set_power;
    if (S[_y][_x] <= 0) {
      r = 1;
    } else
      r = 0;

    if (r == 0) {
      mine_sum++;
    }
    if (r == 1) {
      p_ban.ban[_y][_x] = mine_sum * set_power;
      breaked[{_y, _x}] = 1;
      break;
    }
  }
}

void unit_mining(ll _y, ll _x, PredictBan &p_ban) {
  if (breaked.count({_y, _x}))
    return;
  ll set_power = min((ll)500, p_ban.ban[_y][_x]);
  if (C <= 8)
    set_power /= 10;
  else if (C <= 16)
    set_power /= 5;
  else if (C <= 64)
    set_power /= 2;
  chmax(set_power, (ll)10);
  ll mine_sum = 1;
  while (1) {
    cout << _y << " " << _x << " " << set_power << endl;
    cout.flush();
    score += set_power + C;
    ll r;
    cin >> r;

    if (r == -1 || r == 2)
      return;

    if (r == 0) {
      mine_sum++;
    }
    if (r == 1) {
      p_ban.ban[_y][_x] = mine_sum * set_power;
      breaked[{_y, _x}] = 1;
      break;
    }
  }
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

  void calc_dist_dijkstra(ll cu, vll &dist, PredictBan &p_ban) {
    vvll dp(200, vll(200, -1));
    priority_queue<T, vT, greater<T>> pq;
    pq.push({0, y[cu] * 200 + x[cu], -1});
    while (!pq.empty()) {
      auto [cu_cost, cu, origin] = pq.top();
      pq.pop();
      ll cu_y = cu / 200;
      ll cu_x = cu % 200;
      if (dp[cu_y][cu_x] != -1)
        continue;
      dp[cu_y][cu_x] = cu_cost;

      rep(dir, 0, 4) {
        ll to_y = cu_y + dy[dir];
        ll to_x = cu_x + dx[dir];
        if (to_y < 0 || to_y >= 200 || to_x < 0 || to_x >= 200)
          continue;
        if (dp[to_y][to_x] != -1)
          continue;
        ll at_dist = cu_cost;
        if (!breaked.count({to_y, to_x}))
          at_dist += p_ban.ban[to_y][to_x];
        pq.push({at_dist, to_y * 200 + to_x, cu});
      }
    }
    for (ll house : houses) {
      dist[house] = dp[y[house]][x[house]];
    }
  }

  void break_load(ll cu, ll to, PredictBan &p_ban) {
    vvll dp(200, vll(200, -1));
    vvll dir(200, vll(200, -1));
    priority_queue<T, vT, greater<T>> pq;
    pq.push({0, y[cu] * 200 + x[cu], -1});
    while (!pq.empty()) {
      auto [cu_cost, cu, origin] = pq.top();
      pq.pop();
      ll cu_y = cu / 200;
      ll cu_x = cu % 200;
      if (dp[cu_y][cu_x] != -1)
        continue;
      dp[cu_y][cu_x] = cu_cost;
      dir[cu_y][cu_x] = origin;

      rep(_dir, 0, 4) {
        ll to_y = cu_y + dy[_dir];
        ll to_x = cu_x + dx[_dir];
        if (to_y < 0 || to_y >= 200 || to_x < 0 || to_x >= 200)
          continue;
        if (dp[to_y][to_x] != -1)
          continue;
        ll at_dist = cu_cost;
        if (!breaked.count({to_y, to_x}))
          at_dist += p_ban.ban[to_y][to_x];
        pq.push({at_dist, to_y * 200 + to_x, cu});
      }
    }

    ll at_y = y[to];
    ll at_x = x[to];
    vP _mining;
    while (1) {
      local_unit_mining(at_y, at_x, p_ban);
      // unit_mining(at_y, at_x, p_ban);
      _mining.push_back({at_y, at_x});
      ll _to = dir[at_y][at_x];
      if (_to == -1)
        break;
      at_y = _to / 200;
      at_x = _to % 200;
    }
    p_ban.update_se(_mining);
  }

  void dijkstra(PredictBan &p_ban) {
    vll dp(N, -1);
    priority_queue<T, vT, greater<T>> pq;
    for (ll _sw : sw)
      pq.push({(ll)0, _sw, -1});

    while (!pq.empty()) {
      auto [cu_cost, cu, origin] = pq.top();
      pq.pop();
      if (dp[cu] != -1)
        continue;
      dp[cu] = cu_cost;
      if (origin != -1) {
        break_load(cu, origin, p_ban);
      }

      vll _dist(N, 0);
      calc_dist_dijkstra(cu, _dist, p_ban);
      for (ll to : houses) {
        if (dp[to] != -1)
          continue;
        pq.push({_dist[to], to, cu});
      }
    }
  }
};

int main() {
  auto startClock = system_clock::now();
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

  PredictBan p_ban(1);
  Graph G;
  G.only_Graph();

  G.dijkstra(p_ban);

  // cerr << score << endl;

  const double time =
      duration_cast<microseconds>(system_clock::now() - startClock).count() *
      1e-6;
  cerr << score << " ";
  // cerr << time << endl;
  for(ll _i : stat_P){
    cerr << _i << " ";
  }
  cerr << endl;

  // fstream output_fstream;
  // output_fstream.open("./output2.txt", std::ios_base::out);
  // rep(i, 0, 200) {
  //   rep(j, 0, 200) output_fstream << p_ban.ban[i][j] << " ";
  //   output_fstream << endl;
  // }
  // output_fstream.close();
}

// p_banをy = 10^xとしてる()