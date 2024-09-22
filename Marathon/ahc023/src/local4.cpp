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

typedef vector<Work> vw;
typedef vector<vw> vvw;

bool operator<(const Work &w1, const Work &w2) {
  return (w1.s < w2.s || (w1.s == w2.s && D[w1.k] < D[w2.k]));
}

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

  void rest_num(vw &ans) {
    vvll sum(T, vll(T, 0));
    rep(i, 0, K) sum[S[i] - 1][D[i] - 1]++;
    for (auto cu : ans)
      sum[S[cu.k] - 1][D[cu.k] - 1]--;

    rep(i, 0, T) {
      rep(j, 0, T) cout << sum[i][j] << " ";
      cout << endl;
    }

    ll _sum = 0;
    rep(i, 0, T) rep(j, 0, T) _sum += sum[i][j];
    cerr << "stat_rest : " << _sum << endl;
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

  void side_notuse(vvll &side, vvw &plants) {
    ll N = side.size();
    rep(i, 0, N) {
      ll _N = (ll)plants[i].size();
      ll l = 1e5, r = -1;
      ll num = 0;

      rep(j, 0, _N) {
        l = min(l, plants[i][j].s);
        r = max(r, D[plants[i][j].k]);
        num++;
        int f = 0;
        if (j == _N - 1) f = 1;
        else if (plants[i][j].s != plants[i][j + 1].s) f = 1;

        if (f == 1) {
          if (num != side[i].size()) {
            cout << side[i].size() - num;
            cout << " " << l << " " << r << endl;
          }
          l = 1e5, r = -1, num = 0;
        }
      }
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

  void dfs(ll sta, ll num, ll idx, vvvll &works, vvw &plans, vvll &period) {
    // cout << sta << " " << num << " " << plans[idx].size() << endl;
    if (T <= sta) return;
    // sta日スタート, num個をplans[idx]に格納
    vp range;
    rep(_num, 2, 5) rep(l, 1, 3) rep(r, 1, 3) if (l + r == _num)
        range.push_back({l, r});

    ll best = -1;
    ll sta_j;
    pll best_range;
    for (auto [l, r] : range) {
      vll shu;
      rep(j, sta + 1, T) shu.push_back(j);
      random_device seed_gen;
      mt19937 engine(seed_gen());
      shuffle(shu.begin(), shu.end(), engine);
      // rep(j, sta + 1, T) {
      for (ll j : shu) {
        ll sum = 0;
        rep(_l, 0, l) {
          rep(_r, 0, r) {
            if (T < sta + _l || T < j + _r) break;
            sum += works[sta + _l][j + _r].size();
          }
        }
        if (sum <= num / 2) continue;
        sum = min(sum, num);
        if (chmax(best, sum)) {
          sta_j = j;
          best_range = {l, r};
        }
      }
    }

    if (best == -1) {
      dfs(sta + 2, num, idx, works, plans, period);
      return;
    }

    ll rest = num;
    rep(_r, 0, best_range.second) {
      rep(_l, 0, best_range.first) {
        if (T < sta + _l || T < sta_j + _r) break;
        auto ite = works[sta + _l][sta_j + _r].begin();
        while (ite != works[sta + _l][sta_j + _r].end()) {
          if (rest == 0) break;
          plans[idx].push_back(Work(*ite, side[idx][num - rest] / W,
                                    side[idx][num - rest] % W, sta));
          ite = works[sta + _l][sta_j + _r].erase(ite);
          rest--;
        }
      }
    }
    period[idx].push_back(sta_j + best_range.second);
    // if (rest) cout << idx << " " << num << " " << rest << endl;
    dfs(sta_j + best_range.second, num, idx, works, plans, period);
  }

  void init(vvw &plans, vvll &period) {
    vp at(N);
    rep(i, 0, N) at[i] = {side[i].size(), i};
    gSort(at);

    vvvll works(T + 1, vvll(T + 1));
    rep(i, 0, K) works[S[i]][D[i]].push_back(i);

    rep(_i, 0, N) {
      auto [num, i] = at[_i];
      period[i].push_back(1);
      dfs(1, num, i, works, plans, period);
      period[i].push_back(T + 1);
    }
  }

  void fill(vvw &plans, vvll &period) {
    vll used(K, 0);
    rep(i, 0, N) for (auto cu : plans[i]) used[cu.k] = 1;
    rep(i, 0, N) {
      rep(j, 0, period[i].size() - 1) {
        ll l = period[i][j];
        ll r = period[i][j + 1] - 1;
        rep(k, 0, plans[i].size()) {
          auto cu = plans[i][k];
          if (l <= S[cu.k] && D[cu.k] <= r) plans[i][k].s = l;
        }
      }
    }

    rep(i, 0, N) {
      rep(j, 0, period[i].size() - 1) {
        ll l = period[i][j];
        ll r = period[i][j + 1] - 1;
        // [l,r]日のものをfill_num個補充する必要がある.
        ll fill_num = side[i].size();
        for (auto cu : plans[i]) {
          if (l <= S[cu.k] && D[cu.k] <= r) fill_num--;
        }
        if (fill_num == 0) continue;

        vp score;
        rep(k, 0, K) {
          if (used[k] == 1) continue;
          if (S[k] < l || r < D[k]) continue;
          ll _score = abs(S[k] - l) + abs(D[k] - r);
          score.push_back({_score, k});
        }
        if (score.size() == 0) continue;
        Sort(score);
        // cerr << l << " " << r << " ";
        // cerr << fill_num << " ";
        ll idx = side[i].size() - fill_num;
        rep(k, 0, fill_num) {
          if (score.size() <= k) break;
          plans[i].push_back(
              Work(score[k].second, side[i][idx] / W, side[i][idx] % W, l));
          used[score[k].second] = 1;
          idx++;
        }
        // cerr << idx - (side[i].size() - fill_num) << " " << i << endl;
      }
    }

    rep(i, 0, N) { Sort(plans[i]); }

    rep(i, 0, N) {
      ll num = 0;
      rep(j, 0, (ll)plans[i].size()) {
        plans[i][j].i = side[i][num] / W;
        plans[i][j].j = side[i][num] % W;
        num++;
        if (j != (ll)plans[i].size() - 1) {
          if (plans[i][j].s != plans[i][j + 1].s) num = 0;
        }
      }
    }
  }

  void solve(vector<Work> &ans) {
    auto start_time = system_clock::now();

    vvll period(N);
    vvw plans(N);
    init(plans, period);
    fill(plans, period);

    // Stat st;
    // st.side_notuse(side, plans);

    rep(i, 0, N) {
      for (auto cu : plans[i])
        ans.push_back(cu);
    }
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

  ll calc_score(vw &ans) {
    ll ret = 0;
    for (auto cu : ans) {
      ret += D[cu.k] - S[cu.k] + 1;
    }
    ret *= 1e6;
    ret /= T * H * W;
    return ret;
  }

  void solve(vector<Work> &ans) {
    auto start_time = system_clock::now();

    ll best_score = 0;
    while (1) {
      if (1.5 < ela_times(start_time)) break;
      stat_loop++;

      Graph g;
      vll road;
      construct_road(g, road);

      vvll side;
      SA_side sa_s(0.05);
      sa_s.solve(g, road, side);

      vw _ans;
      SA_plant sa_p(side, 0.5);
      sa_p.solve(_ans);

      ll score = calc_score(_ans);
      if (chmax(best_score, score)) {
        ans = _ans;
      }
      // Stat st;
      // st.rest_num(ans);
      // st.side_occupy(side, ans);
    }
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

// local4.cpp (← local3.cpp)

// 考察
/*
Solve.solve():
  1.5s(loop : 約10回)の最高スコアを採用するものを導入
  ↑ road, side, plantを1回で行う, road, sideの変更でスコアが変動するため.
*/
