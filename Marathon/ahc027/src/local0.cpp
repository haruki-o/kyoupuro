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

ll calc_score(string &op_s) {
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

  void decide_ban_color(vvll &ban_color) {
    ll c_num = 10 + N / 5;
    vp ava;
    ll sum = c_num;
    rep(i, 0, c_num) {
      while (1) {
        ll y = rand() % N;
        ll x = rand() % N;
        if (ban_color[y][x] == -1) {
          ban_color[y][x] = i;
          ava.push_back({y, x});
          break;
        }
      }
    }

    while (1) {
      vp new_ava;
      for (pll cu : ava) {
        ll dir = rand() % 4;
        rep(_dir, 0, 4) {
          dir = (dir + 1) % 4;
          auto [y, x] = to(cu.first, cu.second, dir);
          if (is_out(y, x)) continue;
          if (is_wall(cu, dir)) continue;
          if (ban_color[y][x] != -1) continue;
          new_ava.push_back(cu);
          new_ava.push_back({y, x});
          ban_color[y][x] = ban_color[cu.first][cu.second];
          sum++;
          break;
        }
      }

      std::random_device seed_gen;
      std::mt19937 engine(seed_gen());
      std::shuffle(new_ava.begin(), new_ava.end(), engine);
      ava.clear();
      ava = new_ava;
      if (N * N <= sum) break;
    }
  }

  void decide_color_order(vll &color_order, vvll &ban_color) {
    ll c_num = 0;
    rep(i, 0, N) rep(j, 0, N) chmax(c_num, ban_color[i][j] + 1);
    vll d_sum(c_num, 0);
    vll n_sum(c_num, 0);
    rep(i, 0, N) {
      rep(j, 0, N) {
        d_sum[ban_color[i][j]] += d[i][j];
        n_sum[ban_color[i][j]]++;
      }
    }

    rep(i, 0, c_num) {
      ll visit_num = max((ll)1, d_sum[i] / n_sum[i] / 50);
      if (ban_color[0][0] == i) visit_num = max((ll)1, visit_num - 2);
      // cerr <<d_sum[i] / n_sum[i] << " " << visit_num << endl;
      rep(j, 0, visit_num) color_order.push_back(i);
    }

    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    std::shuffle(color_order.begin(), color_order.end(), engine);
    color_order.insert(color_order.begin(), ban_color[0][0]);
    color_order.push_back(ban_color[0][0]);
  }

  string next_color(pll &cu_cor, ll color, vvll &ban_color,
                    int finish_flag = 0) {
    string ret = "";
    if (ban_color[cu_cor.first][cu_cor.second] == color && finish_flag == 0)
      return ret;

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
      if (ban_color[cu.first][cu.second] == color && finish_flag == 0) {
        to_cor = cu;
        break;
      }
    }

    if (finish_flag == 1) to_cor = {0, 0};

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

  void dfs(string &ret, pll cu_cor, ll color, vvll &ban_color, vvll &visited) {
    rep(dir, 0, 4) {
      auto [y, x] = to(cu_cor.first, cu_cor.second, dir);
      if (is_out(y, x)) continue;
      if (is_wall(cu_cor, dir)) continue;
      if (visited[y][x]) continue;
      if (ban_color[y][x] != color) continue;
      visited[y][x] = 1;
      ret += move_c[dir];
      dfs(ret, {y, x}, color, ban_color, visited);
      ret += move_c[(dir + 2) % 4];
    }
  }

  string range_color(pll &cu_cor, ll color, vvll &ban_color) {
    string ret = "";
    vvll visited(N, vll(N, 0));
    dfs(ret, cu_cor, color, ban_color, visited);
    return ret;
  }

  void decide_op_s(string &op_s, vll &color_order, vvll &ban_color) {
    pll cu_cor = {0, 0};
    rep(i, 0, (ll)color_order.size()) {
      if (i == (ll)color_order.size() - 1)
        op_s += next_color(cu_cor, color_order[i], ban_color, 1);
      else op_s += next_color(cu_cor, color_order[i], ban_color);

      op_s += range_color(cu_cor, color_order[i], ban_color);
    }
  }

  void modify(vll &color_order) {
    ll ran = rand() % 3;
    ll si = color_order.size();
    if (ran == 0) {
      si -= 2;
      ll idx1 = rand() % si;
      ll idx2 = (idx1 + rand() % si) % si;
      idx1++, idx2++;
      swap(color_order[idx1], color_order[idx2]);
    } else if (ran == 1) {
      ll col = color_order[rand() % si];
      ll ins_idx = rand() % si;
      color_order.insert(color_order.begin() + ins_idx, col);
    } else if (ran == 2) {
      map<ll, ll> c_num;
      rep(i, 0, si) c_num[color_order[i]]++;
      vll ava;
      for (auto cu : c_num) {
        if (2 <= cu.second) ava.push_back(cu.first);
      }
      if (ava.size() == 0) return;
      while (1) {
        ll del_idx = rand() % si;
        if (del_idx == 0 || del_idx == si - 1) continue;
        if (c_num[color_order[del_idx]] < 2) continue;
        color_order.erase(color_order.begin() + del_idx);
        break;
      }
    }
  }

  void solve(string &ans) {
    auto start_time = system_clock::now();

    vvll ban_color(N, vll(N, -1));
    decide_ban_color(ban_color);

    // rep(i, 0, N) {
    //   rep(j, 0, N) cout << ban_color[i][j] << " ";
    //   cout << endl;
    // }

    vll color_order;
    decide_color_order(color_order, ban_color);
    string op_s = "";
    decide_op_s(op_s, color_order, ban_color);

    ll best_score = calc_score(op_s);

    while (1) {
      if (1.0 < ela_times(start_time)) break;
      stat_loop++;

      vll _color_order;
      _color_order = color_order;
      modify(_color_order);
      string _op_s = "";
      decide_op_s(_op_s, _color_order, ban_color);
      ll score = calc_score(_op_s);

      if (chmin(best_score, score)) {
        op_s = _op_s;
        cerr << best_score << " " << stat_loop <<  endl;
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

  stat_time = ela_times(start_time);
  stat_ans = calc_score(ans);

  cerr << "stat_loop : " << stat_loop << endl;
  cerr << "stat_time : " << stat_time << endl;
  cerr << "stat_ans  : " << stat_ans << endl;
}


// local0.cpp
// 盤面にN色塗って, color_order[i]=j := i番目にj色を訪れる。
// を用意し, color_order内を焼きなます
// サンプルよりゴミ