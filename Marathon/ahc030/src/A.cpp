#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef long double ld;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
typedef pair<ll, ll> pll;
typedef vector<pll> vp;
typedef vector<vp> vvp;
typedef vector<ld> vld;
typedef vector<vld> vvld;
typedef vector<string> vs;
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

double stat_time = 0.0;
ll stat_loop = 0;
ll stat_ans = 0;

struct Pos {
  ll y, x;

  Pos() {}
  Pos(ll y, ll x) : y(y), x(x) {}
};
typedef vector<Pos> vpos;

struct Oil {
  ll d;
  vpos pos;
  Pos lower_right;
  Pos ans;

  Oil() {
    lower_right.y = 0, lower_right.x = 0;
  }
};
typedef vector<Oil> voil;

ll N, M;
ld ep;
voil oils;
vvll v;
vld e;

vll dy = {0, 1, 0, -1}, dx = {1, 0, -1, 0};

vs color_code = {"red", "orange", "yellow", "chartreuse", "green", "spring green", "cyan", "azure", "blue", "violet", "magenta", "rose"};
vvll color_rgb = {{255, 0, 0}, {255, 165, 0}, {255, 255, 0}, {127, 255, 0}, {0, 255, 0}, {0, 255, 127}, {0, 255, 255}, {0, 127, 255}, {0, 0, 255}, {127, 0, 255}, {255, 0, 255}, {255, 0, 127}};

bool is_out(ll y, ll x) {
  if (y < 0 || y >= N || x < 0 || x >= N) return true;
  return false;
}

double ela_times(system_clock::time_point &clock) {
  return duration_cast<microseconds>(system_clock::now() - clock).count() * 1e-6;
}

void read_input() {
  cin >> N >> M >> ep;
  oils.resize(M);
  rep(i, 0, M) {
    cin >> oils[i].d;

    oils[i].pos.resize(oils[i].d);
    rep(j, 0, oils[i].d) cin >> oils[i].pos[j].y >> oils[i].pos[j].x;
    rep(j, 0, oils[i].d) {
      chmax(oils[i].lower_right.y, oils[i].pos[j].y);
      chmax(oils[i].lower_right.x, oils[i].pos[j].x);
    }
  }

#ifdef DEBUG
  rep(i, 0, M) cin >> oils[i].ans.y >> oils[i].ans.x;
  v.resize(N);
  rep(i, 0, N) {
    v[i].resize(N);
    rep(j, 0, N) cin >> v[i][j];
  }
  e.resize(2 * N * N);
  rep(i, 0, 2 * N * N) cin >> e[i];
#endif
}

struct Solver {
  const ld Multi_Bound = 0.1;
  const ld Check_Bound = 0.9;
  const ld Det_Bound = 0.99;

  Solver() {}

  void init_prob(vvld &prob) {
    rep(i, 0, M) {
      rep(j, 0, N * N) {
        ll move_y = j / N, move_x = j % N;
        if (is_out(oils[i].lower_right.y + move_y, oils[i].lower_right.x + move_x)) continue;
        prob[i][j] = 1 / (ld)((N - oils[i].lower_right.y) * (N - oils[i].lower_right.x));
      }
    }
  }

  void select_dig_maths(vvld &prob, vvll &selected, vpos &dig_maths, ll &sub_oil, ll select_type) {
    Pos ret(0, 0);

    ld ma = 0;
    rep(i, 0, M) rep(j, 0, N * N) chmax(ma, prob[i][j]);

    // if (Multi_Bound < ma) {
    if (0) {
      int flag = 0;
      rep(i, 0, M) {
        rep(j, 0, N * N) {
          if (ma == prob[i][j] && !flag) {
            flag = 1;
            sub_oil = i;
            for (auto cu : oils[i].pos) {
              dig_maths.push_back(Pos({j / N + cu.y, j % N + cu.x}));
            }
          }
        }
      }
    } else {
      vvld sum(N, vld(N, 0));
      rep(i, 0, M) {
        ld ma = 0;
        rep(j, 0, N * N) chmax(ma, prob[i][j]);
        if (select_type == 0 && Check_Bound < ma) continue;
        if (select_type == 1 && Det_Bound < ma) continue;
        rep(j, 0, N * N) {
          if (prob[i][j] == 0) continue;
          ll move_y = j / N, move_x = j % N;
          for (auto cu : oils[i].pos) {
            ll to_y = cu.y + move_y, to_x = cu.x + move_x;
            if (is_out(to_y, to_x)) continue;
            if (selected[to_y][to_x] != -1) continue;
            sum[to_y][to_x] += prob[i][j];
          }
        }
      }
      ld best = 0;
      rep(i, 0, N * N) {
        if (chmax(best, sum[i / N][i % N])) ret = {i / N, i % N};
      }
      if (best == 0) dig_maths.push_back(Pos({-1, -1}));
      else dig_maths.push_back(Pos(ret));
    }
  }

  void single_update_prob(Pos dig_math, ll &res, vvld &prob) {
    vvld at_prob(M, vld(2, 0));
    vld all_prob(2, 0);
    rep(i, 0, M) {
      rep(j, 0, N * N) {
        if (prob[i][j] == 0) continue;
        ll move_y = j / N, move_x = j % N;
        int flag = 0;
        for (auto cu : oils[i].pos) {
          ll to_y = cu.y + move_y, to_x = cu.x + move_x;
          if (is_out(to_y, to_x)) continue;
          if (to_y == dig_math.y && to_x == dig_math.x) flag = 1;
        }
        // dig_mathがoils[i]に含まれる
        if (flag) at_prob[i][1] += prob[i][j];
        else at_prob[i][0] += prob[i][j];
      }
    }
    rep(i, 0, M) rep(j, 0, 2) all_prob[j] += at_prob[i][j];

    if (res == 0) {
      rep(i, 0, M) {
        rep(j, 0, N * N) {
          if (prob[i][j] == 0) continue;
          ll move_y = j / N, move_x = j % N;
          int flag = 0;
          for (auto cu : oils[i].pos) {
            ll to_y = cu.y + move_y, to_x = cu.x + move_x;
            if (is_out(to_y, to_x)) continue;
            if (to_y == dig_math.y && to_x == dig_math.x) flag = 1;
          }
          // dig_mathがoils[i]に含まれる
          if (flag) prob[i][j] = 0;
        }
        ld all = 0;
        rep(j, 0, N * N) all += prob[i][j];
        rep(j, 0, N * N) prob[i][j] /= all;
      }
    }
    if (1 <= res) {
      rep(i, 0, M) {
        // 0.2でM個の中からiが選ばれると
        ld _prob = at_prob[i][1] / all_prob[1];
        // iが選ばれない(厳密には違う)
        ld sum_prob = 1;
        rep(j, 0, res) sum_prob *= (1 - _prob);
        // res回した後, iが1回でも選ばれる
        _prob = 1 - sum_prob;
        rep(j, 0, N * N) {
          if (prob[i][j] == 0) continue;
          ll move_y = j / N, move_x = j % N;
          int flag = 0;
          for (auto cu : oils[i].pos) {
            ll to_y = cu.y + move_y, to_x = cu.x + move_x;
            if (is_out(to_y, to_x)) continue;
            if (to_y == dig_math.y && to_x == dig_math.x) flag = 1;
          }
          // dig_mathがoils[i]に含まれる
          if (flag) prob[i][j] *= _prob;
          else prob[i][j] *= (1 - _prob);
        }
        ld all = 0;
        rep(j, 0, N * N) all += prob[i][j];
        rep(j, 0, N * N) prob[i][j] /= all;
      }
    }
  }

  void multiple_update_prob(vpos dig_maths, ll &res, vvld &prob, ll sub_oil) {
    // Multi_Bound以下の時0, 以上の時すべて試す.
    rep(j, 0, N * N) {
      if (prob[sub_oil][j] <= Multi_Bound) {
        // prob[sub_oil][j]
      }
    }
  }

  void highprob_expect(vvld &prob, map<pll, ll> &pre_coloring) {
    rep(i, 0, M) {
      if (i == color_code.size()) break;
      rep(j, 0, N * N) {
        if (prob[i][j] < Multi_Bound) continue;
        cout << "#c " << j / N << " " << j % N << " ";
        pre_coloring[{j / N, j % N}] = 2;
        vll rgb(3, 0);
        rep(k, 0, 3) rgb[k] = 255 + prob[i][j] * (color_rgb[i][k] - 255);
        cout << "#";
        rep(k, 0, 3) cout << hex << setw(2) << setfill('0') << rgb[k];
        cout << endl;
        cout << dec;
      }
    }
    rep(i, 0, M) {
      ld ma = 0;
      rep(j, 0, N * N) chmax(ma, prob[i][j]);
      if (ma < Det_Bound) continue;
      rep(j, 0, N * N) {
        if (prob[i][j] == 0) continue;
        if (prob[i][j] != ma) continue;
        ll move_y = j / N, move_x = j % N;
        int flag = 0;
        for (auto cu : oils[i].pos) {
          ll to_y = cu.y + move_y, to_x = cu.x + move_x;
          if (is_out(to_y, to_x)) continue;
          cout << "#c " << to_y << " " << to_x << " " << color_code[i] << endl;
          pre_coloring[{to_y, to_x}] = 2;
        }
      }
    }
    for (auto &cu : pre_coloring) {
      if (cu.second == 2) cu.second = 1;
      else if (cu.second == 1) {
        cout << "#c " << cu.first.first << " " << cu.first.second << " ";
        cout << "white" << endl;
        cu.second = 0;
      }
    }
  }

  ld gen(double mean, double stddev) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(mean, stddev);

    return distribution(gen);
  }

  ll gen_sample(ll k, ll v_S) {
    ld mean = (k - v_S) * ep + v_S * (1 - ep);
    ld sigma = k * ep * (1 - ep);
    sigma = sqrt(sigma);
    ld at = gen(mean, sigma);
    at = max((ll)0, (ll)round(at));
    return at;
  }

  bool check_ans(vvll &ans) {
#ifdef DEBUG
    rep(i, 0, N) rep(j, 0, N) {
      if (ans[i][j] != v[i][j]) return false;
    }
#endif
    return true;
  }

  void solve() {
    auto start_time = system_clock::now();

    vvld prob(M, vld(N * N, 0));
    init_prob(prob);
    map<pll, ll> pre_coloring;
    // highprob_expect(prob, pre_coloring);
    vvll selected(N, vll(N, -1));
    int select_type = 0;
    rep(turn, 0, 2 * N * N) {
      vpos dig_maths;
      ll sub_oil;
      select_dig_maths(prob, selected, dig_maths, sub_oil, select_type);

      vvll ans(N, vll(N, 0));
      if (dig_maths[0].y == -1) {
        rep(i, 0, M) {
          ld ma = 0;
          rep(j, 0, N * N) chmax(ma, prob[i][j]);
          rep(j, 0, N * N) {
            if (ma != prob[i][j]) continue;
            for (auto cu : oils[i].pos) {
              ll to_y = cu.y + j / N, to_x = cu.x + j % N;
              if (is_out(to_y, to_x)) continue;
              ans[to_y][to_x]++;
            }
          }
        }
        ll d = 0;
        rep(i, 0, N) rep(j, 0, N) if (ans[i][j]) d++;
        cout << "a " << d << " ";
        rep(i, 0, N) rep(j, 0, N) if (ans[i][j]) cout << i << " " << j << " ";
        cout << endl;
      } else if (dig_maths.size() == 1) {
        cout << "q 1 " << dig_maths[0].y << " " << dig_maths[0].x << endl;
      } else {
        cout << "q " << dig_maths.size() << " ";
        for (auto cu : dig_maths) cout << cu.y << " " << cu.x << " ";
        cout << endl;
      }
      cout.flush();
      ll res = 0;
#ifdef DEBUG
      if (dig_maths[0].y == -1) {
        res = 1;
        if (!check_ans(ans)) res = 0;
      } else {
        for (auto cu : dig_maths) res += v[cu.y][cu.x];
        if (2 <= dig_maths.size()) {
          res = gen_sample(dig_maths.size(), res);
        }
      }
#else
      cin >> res;
#endif

      if (dig_maths.size() == 1 && dig_maths[0].y != -1) single_update_prob(dig_maths[0], res, prob);
      else multiple_update_prob(dig_maths, res, prob, sub_oil);

#ifdef DEBUG
      rep(i, 0, M) {
        cout << "#" << color_code[i];
        cout << "{" << oils[i].ans.y << " " << oils[i].ans.x << "}";
        cout << prob[i][oils[i].ans.y * N + oils[i].ans.x] << endl;
      }
#endif
      // highprob_expect(prob, pre_coloring);
      if (dig_maths.size() == 1 && dig_maths[0].y != -1) {
        selected[dig_maths[0].y][dig_maths[0].x] = res;
      }

      if (dig_maths[0].y == -1) {
        if (res == 1) {
          cerr << "turn : " << turn << endl;
          break;
        }
        if (res == 0) select_type = 1;
      }
    }
  }
};

int main() {
  srand((unsigned int)time(NULL));
  auto start_time = system_clock::now();

  read_input();
  cerr << "M : " << M << endl;
  Solver solver;
  solver.solve();

  stat_time = ela_times(start_time);

  cerr << "stat_loop : " << stat_loop << endl;
  cerr << "stat_time : " << stat_time << endl;
  cerr << "stat_ans  : " << stat_ans << endl;
}

// local0.cpp

// 考察
/*
TLE
*/
