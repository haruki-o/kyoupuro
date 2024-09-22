#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef tuple<ll, ll, ll> T;
typedef vector<T> vT;
typedef vector<vT> vvT;
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

double stat_time = 0.0;
double stat_loop = 0;
ll stat_ans = 0;

ll N;
vvll h, v, d;

vll dy = {-1, 0, 1, 0}, dx = {0, -1, 0, 1};
vc move_c = {'U', 'L', 'D', 'R'};

double ela_times(system_clock::time_point &clock) {
  return duration_cast<microseconds>(system_clock::now() - clock).count() *
         1e-6;
}

void read_input() {
  cin >> N;
  h.resize(N - 1);
  rep(i, 0, N - 1) {
    h[i].resize(N);
    string s;
    cin >> s;
    rep(j, 0, N) h[i][j] = s[j] - '0';
  }

  v.resize(N);
  rep(i, 0, N) {
    v[i].resize(N - 1);
    string s;
    cin >> s;
    rep(j, 0, N - 1) v[i][j] = s[j] - '0';
  }

  d.resize(N);
  rep(i, 0, N) {
    d[i].resize(N);
    rep(j, 0, N) cin >> d[i][j];
  }
}

pll to_dir(ll y, ll x, ll dir) { return {y + dy[dir], x + dx[dir]}; }

bool is_out(ll y, ll x) {
  if (y < 0 || y >= N || x < 0 || x >= N) return true;
  return false;
}

bool is_wall(pll cu, ll dir) {
  ll cu_y = cu.first;
  ll cu_x = cu.second;
  ll to_y = cu.first + dy[dir];
  ll to_x = cu.second + dx[dir];
  if (move_c[dir] == 'U' || move_c[dir] == 'D') {
    if (h[min(cu_y, to_y)][cu_x] == 1) return true;
  }
  if (move_c[dir] == 'L' || move_c[dir] == 'R') {
    if (v[cu_y][min(cu_x, to_x)] == 1) return true;
  }
  return false;
}

void check_op_s(string &op_s) {
  if (1e5 < op_s.size()) cerr << "over op_s length" << endl;

  pll cu = {0, 0};
  vvll used(N, vll(N, 0));
  rep(i, 0, op_s.size()) {
    used[cu.first][cu.second] = 1;
    if (op_s[i] == 'U') cu.first -= 1;
    if (op_s[i] == 'D') cu.first += 1;
    if (op_s[i] == 'L') cu.second -= 1;
    if (op_s[i] == 'R') cu.second += 1;
  }
  if (cu.first != 0 || cu.second != 0) cerr << "final not {0,0}" << endl;

  int flag = 1;
  rep(i, 0, N) rep(j, 0, N) if (used[i][j] == 0) flag = 0;
  if (flag == 0) cerr << "exist never visited mass" << endl;
}

ll ans_calc_score(string &op_s) {
  if (1e5 < op_s.size()) return INF;
  ll ret = 0;
  vvvll g(N, vvll(N));
  pll cu = {0, 0};
  rep(i, 0, op_s.size()) {
    g[cu.first][cu.second].push_back(i);
    if (op_s[i] == 'U') cu.first -= 1;
    if (op_s[i] == 'D') cu.first += 1;
    if (op_s[i] == 'L') cu.second -= 1;
    if (op_s[i] == 'R') cu.second += 1;
  }
  g[0][0].push_back(op_s.size());
  rep(i, 0, N) rep(j, 0, N) if (g[i][j].size() == 0) return INF;
  rep(i, 0, N) {
    rep(j, 0, N) {
      if (g[i][j].size() == 0) cerr << "out" << endl;
      ll last_day = g[i][j][g[i][j].size() - 1];
      ll d_L = (op_s.size() - last_day) * d[i][j];
      g[i][j].insert(g[i][j].begin(), 0);
      g[i][j].push_back(op_s.size());
      rep(k, 1, g[i][j].size()) {
        if (k == 1) {
          ll last_d = d_L + (d[i][j]) * (g[i][j][k] - 1);
          ll add = (d_L + last_d) * (g[i][j][k]) / 2;
          ret += add;
        } else {
          // 初項0, 項数 num, 最終項 d[i][j] * (num  - 1)
          ll num = g[i][j][k] - g[i][j][k - 1];
          ret += d[i][j] * (num - 1) * num / 2;
        }
      }
    }
  }
  ret /= (ll)op_s.size();

  return ret;
}

struct Solver {
  Solver() {}

  ll edit_a(pll cu, vvll &a, string &op_s) {
    if (op_s.size() == 0) return -1;

    vvvll g(N, vvll(N));
    rep(i, 0, op_s.size()) {
      if (op_s[i] == 'U') cu.first -= 1;
      if (op_s[i] == 'D') cu.first += 1;
      if (op_s[i] == 'L') cu.second -= 1;
      if (op_s[i] == 'R') cu.second += 1;
      g[cu.first][cu.second].push_back(i);
    }
    ll ret = 0;
    rep(i, 0, N) {
      rep(j, 0, N) {
        if (g[i][j].size() == 0) {
          ret += (a[i][j] + (a[i][j] + d[i][j] * op_s.size())) * op_s.size() /
                 2 / op_s.size();
          a[i][j] += d[i][j] * op_s.size();
          continue;
        }

        g[i][j].push_back(op_s.size());
        ll sum = 0;
        rep(k, 0, g[i][j].size() - 1) {
          if (k == 0) {
            ll l = a[i][j];
            ll r = l + d[i][j] * g[i][j][0];
            ll num = g[i][j][0];
            sum += (l + r) * num / 2;
          } else {
            ll r = (g[i][j][k] - g[i][j][k - 1]) * d[i][j];
            ll num = g[i][j][k] - g[i][j][k - 1];
            sum += r * num / 2;
          }
        }
        ret += sum / op_s.size();

        ll last_day = g[i][j][g[i][j].size() - 2];
        ll d_L = (op_s.size() - last_day) * d[i][j];
        a[i][j] = d_L;
      }
    }
    return ret / (N * N);
  }

  string move_origin(pll &cu_cor) {
    string ret = "";

    vvll dp(N, vll(N, -2));
    dp[cu_cor.first][cu_cor.second] = -1;
    queue<pll> qu;
    qu.push(cu_cor);
    pll to_cor;
    while (!qu.empty()) {
      pll cu = qu.front();
      qu.pop();
      rep(dir, 0, 4) {
        auto [to_y, to_x] = to_dir(cu.first, cu.second, dir);
        if (is_out(to_y, to_x)) continue;
        if (dp[to_y][to_x] != -2) continue;
        if (is_wall(cu, dir)) continue;
        dp[to_y][to_x] = dir;
        qu.push({to_y, to_x});
      }
      if (cu.first == 0 && cu.second == 0) {
        to_cor = cu;
        break;
      }
    }

    cu_cor = to_cor;
    while (dp[cu_cor.first][cu_cor.second] != -1) {
      ret += move_c[dp[cu_cor.first][cu_cor.second]];
      cu_cor = to_dir(cu_cor.first, cu_cor.second,
                      (dp[cu_cor.first][cu_cor.second] + 2) % 4);
    }

    reverse(ret.begin(), ret.end());
    cu_cor = to_cor;

    return ret;
  }

  void construct_dfs(ll c, ll p, vvll &g, vll &visited, vpp &order) {
    visited[c] = 1;
    ll cu_y = c / N, cu_x = c % N;
    ll dir = rand() % 4;

    rep(_dir, 0, 4) {
      dir = (dir + 1) % 4;
      auto [to_y, to_x] = to_dir(cu_y, cu_x, dir);
      if (is_out(to_y, to_x)) continue;
      if (is_wall({cu_y, cu_x}, dir)) continue;
      if (visited[to_y * N + to_x]) continue;
      if (to_y * N + to_x == p) continue;
      ll to = to_y * N + to_x;
      g[to].push_back(c);
      g[c].push_back(to);
      construct_dfs(to, c, g, visited, order);
    }

    if (c != 0 && g[c].size() == 1) order.push_back({{cu_y, cu_x}, 1});
  }

  void construct_g(vvll &g, vpp &order) {
    vll visited(N * N, 0);
    construct_dfs(0, -1, g, visited, order);
  }

  string route_free(pll &cu, vvll &a, vvll &g) {
    // {cost, dir}
    vvp dp(N, vp(N, {0, -2}));
    // {cost, dir, cor}
    priority_queue<tuple<ll, ll, pll>, vector<tuple<ll, ll, pll>>> pq;
    pq.push({0, -1, cu});
    while (!pq.empty()) {
      auto [cu_cost, cu_dir, cu_cor] = pq.top();
      ll cu_y = cu_cor.first, cu_x = cu_cor.second;
      pq.pop();
      if (dp[cu_y][cu_x].second != -2) continue;
      dp[cu_y][cu_x] = {cu_cost, cu_dir};
      rep(dir, 0, 4) {
        auto [to_y, to_x] = to_dir(cu_y, cu_x, dir);
        if (is_out(to_y, to_x)) continue;
        if (is_wall({cu_y, cu_x}, dir)) continue;
        if (dp[to_y][to_x].second != -2) continue;
        pq.push({cu_cost + a[to_y][to_x], dir, {to_y, to_x}});
      }
    }

    ll ma = 0;
    pll to;
    rep(i, 0, N) rep(j, 0, N) if (chmax(ma, dp[i][j].first)) to = {i, j};
    cu = to;

    string ret = "";
    while (dp[to.first][to.second].second != -1) {
      ret += move_c[dp[to.first][to.second].second];
      to =
          to_dir(to.first, to.second, (dp[to.first][to.second].second + 2) % 4);
    }

    reverse(ret.begin(), ret.end());
    return ret;
  }

  string route_add(pll &cu, vvll &a) {
    // {ava, dir}
    vvp dp(N, vp(N, {0, -2}));
    // {ava, sum, num, dir, cor}
    priority_queue<tuple<ll, ll, ll, ll, pll>,
                   vector<tuple<ll, ll, ll, ll, pll>>>
        pq;
    pq.push({a[cu.first][cu.second], a[cu.first][cu.second], 1, -1, cu});
    while (!pq.empty()) {
      auto [cu_ava, cu_sum, cu_num, cu_dir, cu_cor] = pq.top();
      ll cu_y = cu_cor.first, cu_x = cu_cor.second;
      pq.pop();
      if (dp[cu_y][cu_x].second != -2) continue;
      dp[cu_y][cu_x] = {cu_ava, cu_dir};
      rep(dir, 0, 4) {
        auto [to_y, to_x] = to_dir(cu_y, cu_x, dir);
        if (is_out(to_y, to_x)) continue;
        if (is_wall({cu_y, cu_x}, dir)) continue;
        if (dp[to_y][to_x].second != -2) continue;
        ll add_ava = (cu_sum + a[to_y][to_x]) / (cu_num + 1);
        pq.push(
            {add_ava, cu_sum + a[to_y][to_x], cu_num + 1, dir, {to_y, to_x}});
      }
    }
    pll to = {0, 0};
    rep(i, 0, N) {
      rep(j, 0, N) {
        if (dp[to.first][to.second] < dp[i][j]) to = {i, j};
      }
    }
    cu = to;

    string ret = "";
    while (dp[to.first][to.second].second != -1) {
      ret += move_c[dp[to.first][to.second].second];
      to =
          to_dir(to.first, to.second, (dp[to.first][to.second].second + 2) % 4);
    }
    reverse(ret.begin(), ret.end());
    return ret;
  }

  void solve(string &ans) {
    auto start_time = system_clock::now();

    vvll g(N * N);
    vpp order;
    string op_s = "";
    construct_g(g, order);

    pll cu_cor = {0, 0};
    vvll a(N, vll(N, 0));
    rep(i, 0, N) rep(j, 0, N) a[i][j] = d[i][j] * 20;

    while (1) {
      if (1.9 < ela_times(start_time)) break;
      stat_loop++;

      pll best_cu_cor;
      string best_op_s;
      vvll best_a(N, vll(N));
      ll best_score = INF;

      rep(i, 0, 2) {
        vvll _a(N, vll(N, 0));
        _a = a;
        pll _cu_cor = cu_cor;
        string _op_s;
        if (i == 0) _op_s = route_add(_cu_cor, _a);
        if (i == 1) _op_s = route_free(_cu_cor, _a, g);
        ll _score = edit_a(cu_cor, _a, _op_s);
        if (chmin(best_score, _score)) {
          best_cu_cor = _cu_cor;
          best_op_s = _op_s;
          best_a = _a;
        }
      }

      if (1e5 - 500 < best_op_s.size() + op_s.size()) break;
      op_s += best_op_s;
      cu_cor = best_cu_cor;
      a = best_a;
    }

    op_s += move_origin(cu_cor);
    ans = op_s;
  }
};

int main() {
  srand((unsigned int)time(NULL));
  auto start_time = system_clock::now();

  read_input();

  string ans;
  Solver solver;
  solver.solve(ans);

  check_op_s(ans);
  cout << ans << endl;

  // stat_time = ela_times(start_time);
  stat_ans = ans_calc_score(ans);
  cerr << "stat_loop : " << stat_loop << endl;
  cerr << "stat_time : " << stat_time << endl;
  cerr << "all_time  : " << ela_times(start_time) << endl;
  cerr << "stat_ans  : " << stat_ans << endl;
}

// local18.cpp <- local17.cpp
// route_addとroute_freeの選択性にしました. 