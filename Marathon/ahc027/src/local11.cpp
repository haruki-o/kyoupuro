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

  void edit_a(pll cu, vvll &a, string op_s) {
    vvvll g(N, vvll(N));
    rep(i, 0, op_s.size()) {
      if (op_s[i] == 'U') cu.first -= 1;
      if (op_s[i] == 'D') cu.first += 1;
      if (op_s[i] == 'L') cu.second -= 1;
      if (op_s[i] == 'R') cu.second += 1;
      g[cu.first][cu.second].push_back(i);
    }

    rep(i, 0, N) {
      rep(j, 0, N) {
        if (g[i][j].size() == 0) {
          a[i][j] = d[i][j] * op_s.size();
          continue;
        }
        ll last_day = g[i][j][g[i][j].size() - 1];
        ll d_L = (op_s.size() - last_day) * d[i][j];
        a[i][j] = d_L;
      }
    }
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

  void construct_uf(vvll &g) {
    vp ava;
    rep(i, 0, N) {
      rep(j, 0, N) {
        rep(dir, 2, 4) {
          auto [to_y, to_x] = to_dir(i, j, dir);
          if (is_out(to_y, to_x)) continue;
          if (is_wall({i, j}, dir)) continue;
          ava.push_back({i * N + j, to_y * N + to_x});
        }
      }
    }
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    std::shuffle(ava.begin(), ava.end(), engine);

    UnionFind uf(N * N);
    for (auto [cu, to] : ava) {
      if (uf.same(cu, to)) continue;
      uf.unite(cu, to);
      g[cu].push_back(to);
      g[to].push_back(cu);
    }
  }

  void decide_order(ll c, ll p, vvll &g, vll &visited, vpp &order) {
    visited[c] = 1;
    ll cu_y = c / N, cu_x = c % N;
    if (c != 0 && g[c].size() == 1) order.push_back({{cu_y, cu_x}, 1});
    for (auto to : g[c]) {
      ll to_y = to / N, to_x = to % N;
      if (visited[to_y * N + to_x]) continue;
      if (to_y * N + to_x == p) continue;
      decide_order(to, c, g, visited, order);
    }
  }

  void construct_g(vvll &g, vpp &order) {
    vll visited(N * N, 0);
    // construct_uf(g);
    // decide_order(0, -1, g, visited, order);
    construct_dfs(0, -1, g, visited, order);
  }

  string route(pll cu, pll to, vvll &g) {
    vvll dp(N, vll(N, -2));
    queue<pll> qu;
    qu.push(cu);
    dp[cu.first][cu.second] = -1;
    while (!qu.empty()) {
      auto [cu_y, cu_x] = qu.front();
      qu.pop();
      for (ll _to : g[cu_y * N + cu_x]) {
        ll to_y = _to / N, to_x = _to % N;
        if (dp[to_y][to_x] != -2) continue;
        rep(dir, 0, 4) {
          auto [_y, _x] = to_dir(cu_y, cu_x, dir);
          if (_y == to_y && _x == to_x) {
            dp[to_y][to_x] = dir;
          }
        }
        qu.push({to_y, to_x});
      }
    }

    string op_s = "";
    while (dp[to.first][to.second] != -1) {
      op_s += move_c[dp[to.first][to.second]];
      to = to_dir(to.first, to.second, (dp[to.first][to.second] + 2) % 4);
    }
    reverse(op_s.begin(), op_s.end());
    return op_s;
  }

  string route_free(pll cu, pll to, vvll &g, vvll &a) {
    auto start_time = system_clock::now();

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
    string ret = "";
    while (dp[to.first][to.second].second != -1) {
      ret += move_c[dp[to.first][to.second].second];
      to =
          to_dir(to.first, to.second, (dp[to.first][to.second].second + 2) % 4);
    }
    reverse(ret.begin(), ret.end());
    stat_time += ela_times(start_time);
    return ret;
  }

  void decide_op_s(string &op_s, vvll &g, vpp &order, map<pll, string> &store0,
                   map<pll, string> &store1) {
    string l_op_s = "", r_op_s = "";
    int f = 0;
    rep(i, 1, order.size()) {
      ll cu = order[i - 1].first.first * N + order[i - 1].first.second;
      ll to = order[i].first.first * N + order[i].first.second;
      string at_op_s = "";
      if (order[i].second == 0) {
        if (store0.count({cu, to})) at_op_s = store0[{cu, to}];
        else f = 1;
      } else {
        if (store1.count({cu, to})) at_op_s = store1[{cu, to}];
        else {
          at_op_s = route(order[i - 1].first, order[i].first, g);
          store1[{cu, to}] = at_op_s;
        }
      }
      if (f) r_op_s += at_op_s;
      else l_op_s += at_op_s;
    }

    if (!f) {
      op_s = l_op_s + r_op_s;
      return;
    }
    rep(i, 1, order.size() - 1) {
      ll cu = order[i - 1].first.first * N + order[i - 1].first.second;
      ll to = order[i].first.first * N + order[i].first.second;
      if (order[i].second == 0 && !store0.count({cu, to})) {
        ll toto = order[i + 1].first.first * N + order[i + 1].first.second;
        // {cu, to},{to, toto}がどっちも新しい
        if (!store0.count({to, toto}) && order[i + 1].second == 0) {
          vvll a(N, vll(N));
          edit_a(order[i + 1].first, a, r_op_s + l_op_s);
          string _l_op_s = "", _r_op_s = "";
          _l_op_s = route_free(order[i - 1].first, order[i].first, g, a);
          store0[{cu, to}] = _l_op_s;
          _r_op_s = route_free(order[i].first, order[i + 1].first, g, a);
          store0[{to, toto}] = _r_op_s;
          string _op_s = _l_op_s + _r_op_s;
          op_s = l_op_s + _op_s + r_op_s;
        } else {
          vvll a(N, vll(N));
          edit_a(order[i].first, a, r_op_s + l_op_s);
          string _op_s = route_free(order[i - 1].first, order[i].first, g, a);
          store0[{cu, to}] = _op_s;
          op_s = l_op_s + _op_s + r_op_s;
        }

        break;
      }
    }
  }

  // 2,5,4,3 -> 2,2,5,4,3 -> 2,3,2,5,4,3
  void modify(vpp &order) {
    ll si = order.size();
    ll add_idx = rand() % (si - 2);
    add_idx++;
    ll ins_idx = rand() % (si - 1);
    if (order[ins_idx + 1].second == 0 && rand() % 2 == 0) {
      order.insert(order.begin() + ins_idx + 1, {order[add_idx].first, (ll)0});
    } else {
      order.insert(order.begin() + ins_idx + 1, {order[ins_idx].first, (ll)0});
      order.insert(order.begin() + ins_idx + 1, {order[add_idx].first, (ll)0});
    }
  }

  void solve(string &ans) {
    auto start_time = system_clock::now();

    map<pll, string> store0, store1;

    vvll g(N * N);
    vpp order;
    order.push_back({{0, 0}, 1});
    construct_g(g, order);
    order.push_back({{0, 0}, 1});
    string op_s = "";
    decide_op_s(op_s, g, order, store0, store1);

    ll best_score = ans_calc_score(op_s);

    while (1) {
      if (1.5 < ela_times(start_time)) break;
      stat_loop++;

      vpp _order;
      _order = order;
      modify(_order);
      string _op_s = "";
      decide_op_s(_op_s, g, _order, store0, store1);
      ll score = ans_calc_score(_op_s);
      if (chmin(best_score, score)) {
        op_s = _op_s;
        order = _order;
        cerr << best_score << endl;
      }
    }

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

// local11.cpp <- local9.cpp
// edit_a()を変更し, aを変えました. 没