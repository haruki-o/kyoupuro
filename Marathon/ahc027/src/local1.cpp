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
typedef vector<char> vc;
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

pll to(ll y, ll x, ll dir) { return {y + dy[dir], x + dx[dir]}; }

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

  void edit_a(pll &cu, vvll &a, string &op_s) {
    for (auto op : op_s) {
      if (op == 'U') cu.first -= 1;
      if (op == 'D') cu.first += 1;
      if (op == 'L') cu.second -= 1;
      if (op == 'R') cu.second += 1;
      rep(i, 0, N) rep(j, 0, N) a[i][j] += d[i][j];
      a[cu.first][cu.second] = 0;
    }
  }

  ll eval_range(pll &cu, pll &to_cor, vvll &range, vvll &visited, vvll &a) {
    // 塗る範囲, cuから範囲までの距離, 汚れ度, 塗って無い点の数
    ll ret = 0;
    ret -=
        (abs(cu.first - to_cor.first) + abs(cu.second - to_cor.second)) * 100;

    ll sum = 0;
    ll never_visit_sum = 0;
    ll min_paint_sum = N * 3;
    // ll min_paint_sum = 30 + rand() % 30;

    priority_queue<pair<ll, pll>, vector<pair<ll, pll>>> pq;
    pq.push({0, to_cor});
    while (!pq.empty()) {
      auto [add_score, cu_cor] = pq.top();
      pq.pop();
      ll cu_y = cu_cor.first;
      ll cu_x = cu_cor.second;
      if (min_paint_sum <= sum) break;
      sum++;

      if (visited[cu_y][cu_x] == 0) never_visit_sum++;
      ret += add_score;
      range[cu_y][cu_x] = 1;
      rep(dir, 0, 4) {
        auto [to_y, to_x] = to(cu_y, cu_x, dir);
        if (is_out(to_y, to_x)) continue;
        if (is_wall({cu_y, cu_x}, dir)) continue;
        if (range[to_y][to_x]) continue;

        // 汚れ度, 一度も塗ってないか
        ll score = 0;
        if (visited[to_y][to_x] == 0) score += 500;
        score += a[to_y][to_x] / 10;
        pq.push({score, {to_y, to_x}});
      }
    }
    cerr << "{";
    cerr << ret << " ";
    cerr << min_paint_sum << " ";
    cerr << never_visit_sum << " ";
    cerr << (abs(cu.first - to_cor.first) + abs(cu.second - to_cor.second))
         << " ";
    cerr << "}";
    return ret;
  }

  void decide_paint_range(pll &cu, vvll &paint, string &op_s, vvll &a) {
    vvll visited(N, vll(N, 0));
    pll cu_cor = {0, 0};
    visited[0][0] = 1;
    for (auto op : op_s) {
      if (op == 'U') cu_cor.first -= 1;
      if (op == 'D') cu_cor.first += 1;
      if (op == 'L') cu_cor.second -= 1;
      if (op == 'R') cu_cor.second += 1;
      visited[cu_cor.first][cu_cor.second] = 1;
    }

    vp ava;
    rep(i, 0, N) rep(j, 0, N) if (!visited[i][j]) ava.push_back({i, j});
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    std::shuffle(ava.begin(), ava.end(), engine);

    ll width = min((ll)ava.size(), (ll)5);
    ll best_score = -INF;
    vvll best_range(N, vll(N));
    rep(i, 0, width) {
      vvll range(N, vll(N, 0));
      ll score = eval_range(cu, ava[i], range, visited, a);
      if (chmax(best_score, score)) {
        best_range = range;
      }
    }
    cerr << endl;
    paint = best_range;
  }

  string next_range(pll &cu_cor, vvll &paint) {
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
        auto [to_y, to_x] = to(cu.first, cu.second, dir);
        if (is_out(to_y, to_x)) continue;
        if (dp[to_y][to_x] != -2) continue;
        if (is_wall(cu, dir)) continue;
        dp[to_y][to_x] = dir;
        qu.push({to_y, to_x});
      }
      if (paint[cu.first][cu.second]) {
        to_cor = cu;
        break;
      }
    }

    cu_cor = to_cor;
    while (dp[cu_cor.first][cu_cor.second] != -1) {
      ret += move_c[dp[cu_cor.first][cu_cor.second]];
      cu_cor = to(cu_cor.first, cu_cor.second,
                  (dp[cu_cor.first][cu_cor.second] + 2) % 4);
    }

    reverse(ret.begin(), ret.end());
    cu_cor = to_cor;

    return ret;
  }

  void dfs(string &ret, pll cu_cor, vvll &visited, vvll &paint) {
    ll dir = rand() % 4;
    rep(_dir, 0, 4) {
      dir = (dir + 1) % 4;
      auto [y, x] = to(cu_cor.first, cu_cor.second, dir);
      if (is_out(y, x)) continue;
      if (is_wall(cu_cor, dir)) continue;
      if (visited[y][x]) continue;
      if (paint[y][x] == 0) continue;
      visited[y][x] = 1;

      ret += move_c[dir];
      dfs(ret, {y, x}, visited, paint);
      ret += move_c[(dir + 2) % 4];
    }
  }

  void bfs(string &ret, pll cu_cor, vvll &visited, vvll &paint) {
    vvll dp1(N, vll(N, -1));
    queue<pll> qu;
    qu.push(cu_cor);
    dp1[cu_cor.first][cu_cor.second] = 0;
    while (!qu.empty()) {
      pll cu = qu.front();
      qu.pop();
      rep(dir, 0, 4) {
        auto [to_y, to_x] = to(cu.first, cu.second, dir);
        if (is_out(to_y, to_x)) continue;
        if (is_wall(cu, dir)) continue;
        if (dp1[to_y][to_x] != -1) continue;
        if (paint[to_y][to_x] != 0) continue;
        dp1[to_y][to_x] = dp1[cu.first][cu.second] + 1;
        qu.push({to_y, to_x});
      }
    }

    pll end = cu_cor;
    rep(i, 0, N) {
      rep(j, 0, N) if (dp1[end.first][end.second] < dp1[i][j]) end = {i, j};
    }

    vvll dp2(N, vll(N, -1));
    qu.push(end);
    dp2[end.first][end.second] = 0;
    while (!qu.empty()) {
      pll cu = qu.front();
      qu.pop();
      rep(dir, 0, 4) {
        auto [to_y, to_x] = to(cu.first, cu.second, dir);
        if (is_out(to_y, to_x)) continue;
        if (is_wall(cu, dir)) continue;
        if (dp2[to_y][to_x] != -1) continue;
        if (paint[to_y][to_x] != 0) continue;
        dp2[to_y][to_x] = dp1[cu.first][cu.second] + 1;
        qu.push({to_y, to_x});
      }
    }

    priority_queue<pair<pll, pll>, vector<pair<pll, pll>>,
                   greater<pair<pll, pll>>>
        pq;
    pq.push({{0, 0}, cu_cor});
    rep(i, 0, N) {
      rep(j, 0, N) {
        if (paint[i][j] == 0) continue;
        pq.push({{dp1[i][j], INF - dp2[i][j]}, {i, j}});
      }
    }

    vp order;
    while (!pq.empty()) {
      auto [_, cu] = pq.top();
      pq.pop();
      order.push_back(cu);
    }
    cerr << "d";

    rep(i, 0, order.size()) {
      cerr << "cu : " << cu_cor.first << " " << cu_cor.second << endl;
      string _ret = "";
      cerr << "to : " << order[i].first << " " << order[i].second << endl;

      vvll dp(N, vll(N, -2));
      dp[cu_cor.first][cu_cor.second] = -1;
      qu.push(cu_cor);
      pll to_cor;
      while (!qu.empty()) {
        pll cu = qu.front();
        qu.pop();
        int f = 0;
        rep(dir, 0, 4) {
          auto [to_y, to_x] = to(cu.first, cu.second, dir);
          if (is_out(to_y, to_x)) continue;
          if (dp[to_y][to_x] != -2) continue;
          if (is_wall(cu, dir)) continue;
          if (paint[to_y][to_x] == 0) continue;
          dp[to_y][to_x] = dir;
          qu.push({to_y, to_x});
          if (to_y == order[i].first && to_x == order[i].second) {
            to_cor = {to_y, to_x};
            f = 1;
            break;
          }
        }
        if (f) {
          break;
        }
      }
      cerr << "to :  " << to_cor.first << " " << to_cor.second << endl;

      cu_cor = to_cor;
      ll sum1 = 0;
      while (dp[cu_cor.first][cu_cor.second] != -1) {
        sum1++;
        if(sum1 == 10)break;
        cerr << cu_cor.first << " " << cu_cor.second << endl;
        _ret += move_c[dp[cu_cor.first][cu_cor.second]];
        cu_cor = to(cu_cor.first, cu_cor.second,
                    (dp[cu_cor.first][cu_cor.second] + 2) % 4);
      }

      reverse(_ret.begin(), _ret.end());
      cu_cor = to_cor;
      ret += _ret;
    }
  }

  string paint_range(pll &cu_cor, vvll &paint) {
    string ret = "";
    pll dfs_cu = cu_cor;
    string dfs_ret = "";
    vvll visited(N, vll(N, 0));
    dfs(dfs_ret, dfs_cu, visited, paint);

    // pll bfs_cu = cu_cor;
    // string bfs_ret = "";
    // bfs(bfs_ret, bfs_cu, visited, paint);
    // if (dfs_ret.size() < bfs_ret.size()) ret = dfs_ret;
    // else ret = bfs_ret;
    ret = dfs_ret;

    return ret;
  }

  bool check_all_visited(string &op_s) {
    vvll visited(N, vll(N, 0));
    pll cu_cor = {0, 0};
    visited[0][0] = 1;
    for (auto op : op_s) {
      if (op == 'U') cu_cor.first -= 1;
      if (op == 'D') cu_cor.first += 1;
      if (op == 'L') cu_cor.second -= 1;
      if (op == 'R') cu_cor.second += 1;
      visited[cu_cor.first][cu_cor.second] = 1;
    }
    int f = 1;
    rep(i, 0, N) rep(j, 0, N) if (!visited[i][j]) f = 0;
    if (f) return true;
    return false;
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
        auto [to_y, to_x] = to(cu.first, cu.second, dir);
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
      cu_cor = to(cu_cor.first, cu_cor.second,
                  (dp[cu_cor.first][cu_cor.second] + 2) % 4);
    }

    reverse(ret.begin(), ret.end());
    cu_cor = to_cor;

    return ret;
  }

  void solve(string &ans) {
    auto start_time = system_clock::now();

    string op_s = "";
    pll cu = {0, 0};
    vvll a(N, vll(N));
    a = d;
    while (1) {
      stat_loop++;

      vvll paint(N, vll(N, 0));
      decide_paint_range(cu, paint, op_s, a);

      // if (stat_loop == 1) {
      //   rep(i, 0, N) {
      //     rep(j, 0, N) cout << paint[i][j] << " ";
      //     cout << endl;
      //   }
      // }

      pll store_cu = cu;
      string _op_s = next_range(cu, paint);
      _op_s += paint_range(cu, paint);

      edit_a(store_cu, a, _op_s);

      op_s += _op_s;
      if (check_all_visited(op_s)) break;
    }

    op_s += move_origin(cu);

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

  stat_time = ela_times(start_time);
  stat_ans = ans_calc_score(ans);
  cerr << "stat_loop : " << stat_loop << endl;
  cerr << "stat_time : " << stat_time << endl;
  cerr << "stat_ans  : " << stat_ans << endl;
}

