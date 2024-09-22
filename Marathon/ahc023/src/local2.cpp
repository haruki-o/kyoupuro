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
typedef vector<pll> vp;
typedef vector<vp> vvp;
#define rep(i, l, n) for (ll i = (ll)(l); i < (ll)(n); i++)
#define repd(i, n, l) for (ll i = (ll)(n); i > (ll)(l); i--)
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
template <class T, class S> inline bool chmax(T &a, const S &b) {
  return (a < b ? a = b, 1 : 0);
}
template <class T, class S> inline bool chmin(T &a, const S &b) {
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

double stat_time = 0.0;
double stat_loop = 0;
ll stat_ans = 0;

ll T, H, W, i0;
vvll h, v;
ll K;
vll S, D;

vvll g;
ll sg;

vll dy = {0, 1, 0, -1}, dx = {1, 0, -1, 0};

double ela_times(system_clock::time_point &clock) {
  return duration_cast<microseconds>(system_clock::now() - clock).count() *
         1e-6;
}

bool is_out(ll i, ll j) {
  if (i < 0 || i >= H || j < 0 || j >= W) return true;
  return false;
}

bool can_through(ll cu, ll to) {
  if (cu % W != to % W) {
    ll mi = min(cu % W, to % W);
    if (v[cu / W][mi] == 0) return true;
  }
  if (cu / W != to / W) {
    ll mi = min(cu / W, to / W);
    if (h[mi][cu % W] == 0) return true;
  }
  return false;
}

struct Work {
  ll k, i, j, s;

  Work(ll k, ll i, ll j, ll s) : k(k), i(i), j(j), s(s) {}
};

struct Pos {
  ll i, j;

  Pos() {}
  Pos(ll i, ll j) : i(i), j(j) {}
};

void read_input() {
  cin >> T >> H >> W >> i0;
  h.resize(H - 1), v.resize(H);
  rep(i, 0, H - 1) h[i].resize(W);
  rep(i, 0, H) v[i].resize(W - 1);
  rep(i, 0, H - 1) {
    string s;
    cin >> s;
    rep(j, 0, W) h[i][j] = s[j] - '0';
  }
  rep(i, 0, H) {
    string s;
    cin >> s;
    rep(j, 0, W - 1) v[i][j] = s[j] - '0';
  }
  cin >> K;
  S.resize(K), D.resize(K);
  rep(i, 0, K) cin >> S[i] >> D[i];
}

struct Graph {
  vvll g;
  const ll sg = i0 * W;
  vll dist;

  Graph() {
    g.resize(H * W);
    rep(i, 0, H) {
      rep(j, 0, W) {
        if (i != H - 1) {
          if (h[i][j] == 0) g[i * W + j].push_back((i + 1) * W + j);
        }
        if (j != W - 1) {
          if (v[i][j] == 0) g[i * W + j].push_back(i * W + j + 1);
        }
      }
    }
    rep(i, 0, H * W) {
      for (ll j : g[i])
        if (i < j) g[j].push_back(i);
    }

    random_device seed_gen;
    rep(i, 0, H * W) {
      std::mt19937 engine(seed_gen());
      std::shuffle(g[i].begin(), g[i].end(), engine);
    }
  }

  ll bfs(ll s_bfs) {
    dist.assign(H * W, -1);
    queue<ll> qu;
    qu.push(s_bfs);
    dist[qu.front()] = 0;
    while (!qu.empty()) {
      ll cu = qu.front();
      qu.pop();
      for (ll to : g[cu]) {
        if (dist[cu] != -1) continue;
        qu.push(to);
        dist[to] = dist[cu] + 1;
      }
    }

    ll ma = 0;
    rep(i, 0, H * W) if (dist[ma] < dist[i]) ma = i;
    return ma;
  }
};

// 一直線だった時, 採用数が最大化します.
void adopt_plan(vll plan, vll &_adopt_plan, ll max_period) {
  // plan[i] := plantのidx
  ll N = (ll)plan.size();
  vp period(N);
  rep(i, 0, N) period[i] = {D[plan[i]], i};
  Sort(period);

  ll sum = 0;
  ll pre = -1;
  rep(i, 0, N) {
    if (S[plan[period[i].second]] <= pre) continue;

    sum++;
    _adopt_plan[period[i].second] = 1;
    if (sum == max_period) {
      sum = 0;
      pre = period[i].first;
    }
  }
}

struct Stat {
  void shrink() {
    vvll g(H, vll(W, -1));
    rep(i, 0, H) {
      rep(j, 0, W) {
        rep(dir, 0, 4) if (is_out(i + dy[dir], j + dx[dir])) g[i][j] = 0;
        if (i != H - 1) {
          if (h[i][j] == 1) g[i][j] = g[i + 1][j] = 0;
        }
        if (j != W - 1) {
          if (v[i][j] == 1) g[i][j] = g[i][j + 1] = 0;
        }
      }
    }
    rep(color, 0, 50) {
      rep(i, 0, H) {
        rep(j, 0, W) {
          if (g[i][j] != -1) continue;
          rep(dir, 0, 4) {
            ll toi = i + dy[dir];
            ll toj = j + dx[dir];
            if (is_out(toi, toj)) continue;
            // 横
            if (dir % 2 == 0) {
              ll mi_j = min(j, toj);
              ll ma_j = max(j, toj);
              if (v[i][mi_j] == 1) continue;
            }
            // 縦
            if (dir % 2 == 1) {
              ll mi_i = min(i, toi);
              ll ma_i = max(i, toi);
              if (h[mi_i][j] == 1) continue;
            }

            if (g[toi][toj] == color) g[i][j] = color + 1;
          }
        }
      }
    }
    rep(i, 0, H) {
      rep(j, 0, W) cout << g[i][j] << " ";
      cout << endl;
    }
  }

  void range_in_num() {
    vvll sum(T, vll(T, 0));
    rep(i, 0, K) sum[S[i] - 1][D[i] - 1]++;
    rep(i, 0, T) {
      rep(j, 0, T) cout << sum[i][j] << " ";
      cout << endl;
    }
  }

  void side_color(vvll &side) {
    vvll g(H, vll(W, 0));
    ll N = side.size();
    vll shu(N);
    rep(i, 0, N) shu[i] = i;

    random_device seed_gen;
    mt19937 engine(seed_gen());
    shuffle(shu.begin(), shu.end(), engine);

    rep(i, 0, N) {
      ll color = shu[i];
      rep(j, 0, (ll)side[i].size()) {
        g[side[i][j] / W][side[i][j] % W] = color + 1;
      }
    }
    rep(i, 0, H) {
      rep(j, 0, W) cout << g[i][j] << " ";
      cout << endl;
    }
  }

  void side_num(vvll &side) {
    ll N = side.size();
    vp num(N);
    rep(i, 0, N) num[i] = {side.size(), i};
    gSort(num);

    rep(i, 0, N) cerr << num[i].first << " ";
    cerr << endl;
  }

  void side_occupy(vvll &side, vector<Work> &ans) {
    ll N = side.size();
    vll all(N);
    rep(i, 0, N) all[i] = side[i].size() * T;

    vll sum(N, 0);
    vvll g(H, vll(W, -1));
    rep(i, 0, N) {
      for (auto cu : side[i])
        g[cu / W][cu % W] = i;
    }

    for (auto cu : ans) {
      sum[g[cu.i][cu.j]] += D[cu.k] - S[cu.k] + 1;
    }

    vector<pair<double, ll>> stat(N);
    rep(i, 0, N) stat[i] = {sum[i] / (double)all[i], side[i].size()};
    gSort(stat);
    rep(i, 0, N) {
      cerr << setprecision(2) << stat[i].first;
      cerr << " " << stat[i].second << endl;
    }
  }
};

struct SA_side {
  double time;

  SA_side(double time) : time(time) {}

  void init(Graph &G, vll &road, vvll &side) {
    vll seen(H * W, 0);
    for (ll cu : road)
      seen[cu] = 1;

    queue<pll> qu;
    ll idx = 0;
    rep(i, 0, (ll)road.size()) {
      for (ll to : G.g[road[i]]) {
        if (seen[to] == 1) continue;
        seen[to] = 1;
        qu.push({idx, to});
        idx++;
      }
    }
    vvll _side(idx);

    queue<pll> qu2;
    while (!qu.empty()) {
      auto [cu_idx, cu] = qu.front();
      qu.pop();
      _side[cu_idx].push_back(cu);
      qu2.push({cu_idx, cu});
    }
    while (!qu2.empty()) {
      auto cu = qu2.front();
      qu2.pop();
      qu.push(cu);
    }

    while (!qu.empty()) {
      auto [cu_idx, cu] = qu.front();
      qu.pop();
      for (ll to : G.g[cu]) {
        if (seen[to] == 1) continue;
        qu.push({cu_idx, to});
        _side[cu_idx].push_back(to);
        seen[to] = 1;
      }
    }

    rep(i, 0, idx) {
      if (!_side[i].empty()) side.push_back(_side[i]);
    }
  }

  void dfs(ll c, ll p, vll &ord, vll &low, vll &g, vll &is_joint, ll &cnt) {
    ord[c] = cnt;
    rep(dir, 0, 4) {
      ll toi = c / W + dy[dir];
      ll toj = c % W + dx[dir];
      if (is_out(toi, toj)) continue;
      ll to = toi * W + toj;
      if (!can_through(c, to)) continue;
      if (p == to) continue;
      if (g[to] == 0) continue;
      if (ord[to] != -1) {
        low[c] = min(low[c], ord[to]);
        continue;
      }
      cnt++;
      dfs(to, c, ord, low, g, is_joint, cnt);

      low[c] = min(low[c], low[to]);
      if (ord[c] <= low[to]) is_joint[c] = 1;
    }
  }

  bool modify(vvll &side) {
    ll N = side.size();
    vp num(N);
    rep(i, 0, N) num[i] = {side[i].size(), i};
    gSort(num);
    ll idx = num[rand() % min((ll)5, N)].second;

    vll g(H * W, 0);
    for (ll cu : side[idx])
      g[cu] = 1;

    vll is_joint(H * W, 0);
    is_joint[side[idx][0]] = 1;
    vll ord(H * W, -1);
    vll low(H * W, 1e5);
    ll cnt = 0;
    dfs(side[idx][0], -1, ord, low, g, is_joint, cnt);

    vll color(H * W, -1);
    rep(i, 0, N) {
      for (ll cu : side[i])
        color[cu] = i;
    }

    ll era_cu = -1;
    ll era_to = -1;
    ll best = 1e5;
    for (ll cu : side[idx]) {
      if (is_joint[cu]) continue;
      rep(dir, 0, 4) {
        ll toi = cu / W + dy[dir];
        ll toj = cu % W + dx[dir];
        if (is_out(toi, toj)) continue;
        ll to = toi * W + toj;
        if (color[to] == -1 || color[to] == color[cu]) continue;
        if (!can_through(cu, to)) continue;
        if (chmin(best, (ll)side[color[to]].size())) {
          era_cu = cu;
          era_to = to;
        }
      }
    }

    if (era_cu == -1) return false;
    auto ite = side[idx].begin();
    while (ite != side[idx].end()) {
      if (*ite == era_cu) ite = side[idx].erase(ite);
      else ite++;
    }
    side[color[era_to]].push_back(era_cu);

    return true;
  }

  void adjust(vvll &side) {
    ll N = side.size();
    vll g(H * W, -1);
    rep(i, 0, N) {
      for (ll cu : side[i])
        g[cu] = i;
    }

    vll seen(H * W, 0);
    rep(i, 0, N) {
      queue<ll> qu;
      qu.push(side[i][0]);
      side[i].clear();
      seen[qu.front()] = 1;
      side[i].push_back(qu.front());

      while (!qu.empty()) {
        ll cu = qu.front();
        qu.pop();
        rep(dir, 0, 4) {
          ll toi = cu / W + dy[dir];
          ll toj = cu % W + dx[dir];
          if (is_out(toi, toj)) continue;
          ll to = toi * W + toj;
          if (seen[to] == 1) continue;
          if (g[cu] != g[to]) continue;
          if (!can_through(cu, to)) continue;
          qu.push(to);
          side[i].push_back(to);
          seen[to] = 1;
        }
      }
    }
  }

  void solve(Graph &G, vll &road, vvll &side) {
    auto start_time = system_clock::now();

    init(G, road, side);

    while (1) {
      if (time < ela_times(start_time)) break;

      // if (modify(side)) break;
      modify(side);
    }
    adjust(side);
  }
};

struct SA_plant {
  ll N;
  vvll side;
  double time;

  SA_plant(vvll &side, double time) : side(side), time(time) {
    N = (ll)side.size();
  }

  void init(vvll &period) {
    rep(i, 0, N) {
      period[i].push_back(1);
      period[i].push_back(T + 1);

      vll shu;
      rep(j, 3, T - 1) shu.push_back(j);
      random_device seed_gen;
      mt19937 engine(seed_gen());
      shuffle(shu.begin(), shu.end(), engine);
      rep(j, 0, rand() % 4 + 5) period[i].push_back(shu[j]);
      Sort(period[i]);
    }
  }

  void modify(vvll &period) {
    ll idx_i = rand() % N;
    ll idx_j = rand() % period[idx_i].size();
    map<ll, ll> ma;
    for (ll cu : period[idx_i])
      ma[cu] = 1;
    ll val = period[idx_i][idx_j];
    if (rand() % 2) {
      if (val != T + 1 && !ma.count(val + 1)) period[idx_i][idx_j]++;
    } else {
      if (val != 1 && !ma.count(val - 1)) period[idx_i][idx_j]--;
    }
  }

  ll calc_score(vvll &period) {
    vvll num(T + 1, vll(T + 1, 0));
    rep(i, 0, K) num[S[i]][D[i]]++;

    ll ret = 0;
    rep(i, 0, N) {
      rep(j, 0, (ll)period[i].size() - 1) {
        // [l, r)
        ll l = period[i][j];
        ll r = period[i][j + 1];

        vp range;
        rep(_num, 0, 5) {
          rep(_l, 0, 5) {
            rep(_r, 0, 5) {
              if (r - 1 - _r <= l + _l) break;
              if (_num == _l + _r) range.push_back({l + _l, r - 1 - _r});
            }
          }
        }

        ll rest = side[i].size();
        for (auto [_l, _r] : range) {
          ll diff = (ll)side[i].size();
          diff = min(diff, num[_l][_r]);
          diff = min(diff, rest);

          rest -= diff;
          num[_l][_r] -= diff;
          ret += (_r - _l + 1) * diff;
          if (rest == 0) break;
        }
        if (rest) {
          cout << i << " " << side[i].size() << " " << rest << " ";
          cout << "{" << l << " " << r << "}" << endl;
        }
      }
    }
    return ret;
  }

  void decide_ans(vvll &period, vector<Work> &ans) {
    vvvll num(T + 1, vvll(T + 1));
    rep(i, 0, K) num[S[i]][D[i]].push_back(i);

    rep(i, 0, N) {
      rep(j, 0, (ll)period[i].size() - 1) {
        // [l, r)
        ll l = period[i][j];
        ll r = period[i][j + 1];

        vp range;
        rep(_num, 0, 5) {
          rep(_l, 0, 5) {
            rep(_r, 0, 5) {
              if (r - 1 - _r <= l + _l) break;
              if (_num == _l + _r) range.push_back({l + _l, r - 1 - _r});
            }
          }
        }
        vp _range;
        for (auto cu : range)
          _range.push_back({cu.second, cu.first});
        Sort(_range);
        range.clear();
        for (auto cu : _range)
          range.push_back({cu.second, cu.first});

        ll rest = side[i].size();
        ll _idx = 0;
        for (auto [_l, _r] : range) {
          ll diff = (ll)side[i].size();
          diff = min(diff, (ll)num[_l][_r].size());
          diff = min(diff, rest);

          rest -= diff;
          auto ite = num[_l][_r].begin();
          ll _rest = diff;
          while (ite != num[_l][_r].end()) {
            if (_rest == 0) break;
            ans.push_back(Work(*ite, side[i][_idx] / W, side[i][_idx] % W, l));
            _idx++;
            ite = num[_l][_r].erase(ite);
            _rest--;
          }
          if (rest == 0) break;
        }
      }
    }
  }

  void solve(vector<Work> &ans) {
    auto start_time = system_clock::now();

    // [period[i][l], period[i][l + 1]) :=
    // 範囲の中にside[i].size()個以下を詰め込む
    vvll period(N);
    init(period);
    ll best_score = calc_score(period);
    double s_temp = 30000, e_temp = 10000;

    while (1) {
      if (time < ela_times(start_time)) break;

      vvll _period(N);
      _period = period;
      modify(_period);
      stat_loop++;
      ll score = calc_score(_period);
      // 温度関数
      double temp = s_temp + (e_temp - s_temp) * ela_times(start_time) / time;
      // 遷移確率関数(最大化の場合)
      double prob = exp((score - best_score) / temp);
      // if (prob > (rand() % INF) / (double)INF) { // 確率probで遷移する
      if (best_score < score) {
        period = _period;
        best_score = score;
      }
    }

    decide_ans(period, ans);
  }
};

struct Solver {
  Solver() {}

  void construct_road(Graph &G, vll &road) {
    road.push_back(G.sg);
    while (1) {
      vll dist(H * W, -1);
      vll pre(H * W, -1);
      queue<ll> qu;
      for (ll cu : road) {
        qu.push(cu);
        dist[cu] = 0;
      }

      while (!qu.empty()) {
        ll cu = qu.front();
        qu.pop();
        for (ll to : G.g[cu]) {
          if (dist[to] != -1) continue;
          qu.push(to);
          dist[to] = dist[cu] + 1;
          pre[to] = cu;
        }
      }

      ll eg = 0;
      rep(i, 0, H * W) if (dist[eg] < dist[i]) eg = i;

      if (dist[eg] <= 10) break;

      ll sum = 0;
      while (eg != -1) {
        sum++;
        if (5 < sum && eg != G.sg) road.push_back(eg);
        eg = pre[eg];
      }
    }
  }

  void solve(vector<Work> &ans) {
    Graph g;
    vll road;
    construct_road(g, road);
    vvll side;
    SA_side sa_s(0.5);
    sa_s.solve(g, road, side);

    SA_plant sa_p(side, 0.5);
    sa_p.solve(ans);

    // Stat st;
    // st.side_occupy(side, ans);
  }
};

void print_ans(vector<Work> ans) {
  cout << (ll)ans.size() << endl;
  for (auto cu : ans) {
    cout << cu.k + 1 << " " << cu.i << " ";
    cout << cu.j << " " << cu.s << endl;
  }

  for (auto cu : ans) {
    stat_ans += D[cu.k] - S[cu.k] + 1;
  }
  stat_ans *= 1e6;
  stat_ans /= T * H * W;
}

int main() {
  srand((unsigned int)time(NULL));
  auto start_time = system_clock::now();

  read_input();

  vector<Work> ans;
  Solver solver;
  solver.solve(ans);

  print_ans(ans);

  cerr << "stat_loop : " << stat_loop << endl;
  cerr << "stat_time : " << stat_time << endl;
  cerr << "stat_ans  : " << stat_ans << endl;
}

// local1.cpp

// 考察
/*
SA_plant.init() : num = min(num, sum)を追加のみ
*/
