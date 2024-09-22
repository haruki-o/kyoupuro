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

struct SA_plant {
  ll N;
  vvll side;
  double time;

  SA_plant(vvll &side, double time) : side(side), time(time) {
    N = (ll)side.size();
  }

  void init(vvll &plans) {
    ll sum = 0;
    rep(i, 0, N) sum += (ll)side[i].size();

    vp at;
    rep(i, 0, K) at.push_back({D[i], i});
    Sort(at);
    ll idx = 0;
    rep(i, 0, N) {
      ll rat = (double)side[i].size() / sum * K + (double)0.5;
      rep(j, 0, rat) {
        if (K <= idx) break;
        plans[i].push_back(at[idx].second);
        idx++;
      }
    }

    while (idx < K) {
      plans[0].push_back(at[idx].second);
      idx++;
    }
  }

  void modify(vvll &plans) {
    ll idx1 = rand() % N;
    ll idx2 = (idx1 + rand() % N) % N;
    swap(plans[idx1][rand() % plans[idx1].size()],
         plans[idx2][rand() % plans[idx2].size()]);
  }

  ll calc_score(vvll &plans) {
    ll ret = 0;
    rep(i, 0, N) {
      vll plan((ll)plans[i].size(), 0);
      adopt_plan(plans[i], plan, (ll)side[i].size());
      rep(j, 0, (ll)plans[i].size()) {
        if (plan[j] == 1) {
          ret += D[plans[i][j]] - S[plans[i][j]] + 1;
        }
      }
    }
    return ret;
  }

  void solve(vector<Work> &ans) {
    auto start_time = system_clock::now();

    vvll plans(N);
    init(plans);
    ll best_score = calc_score(plans);
    cerr << best_score << endl;
    double s_temp = 30000, e_temp = 10000;

    while (1) {
      if (time < ela_times(start_time)) break;

      vvll _plans(N);
      _plans = plans;
      modify(_plans);

      ll score = calc_score(_plans);
      // 温度関数
      double temp = s_temp + (e_temp - s_temp) * ela_times(start_time) / time;
      // 遷移確率関数(最大化の場合)
      double prob = exp((score - best_score) / temp);
      // if (prob > (rand() % INF) / (double)INF) { // 確率probで遷移する
      if (best_score < score) {
        plans = _plans;
        best_score = score;
      }
    }

    rep(i, 0, N) {
      vp period((ll)plans[i].size());
      rep(j, 0, (ll)plans[i].size()) period[j] = {D[plans[i][j]], j};
      Sort(period);

      ll sum = 0;
      ll pre = 0;
      rep(j, 0, (ll)plans[i].size()) {
        ll idx = period[j].second;
        if (S[plans[i][idx]] <= pre) continue;
        ans.push_back(
            Work(plans[i][idx], side[i][sum] / W, side[i][sum] % W, pre + 1));

        sum++;
        if (sum == (ll)side[i].size()) {
          sum = 0;
          pre = period[j].first;
        }
      }
    }
  }
};

struct Solver {
  Solver() {}

  void construct_road(Graph &G, vll &road) {
    vll dist(H * W, -1);
    vll pre(H * W, -1);
    queue<ll> qu;
    qu.push(G.sg);
    dist[qu.front()] = 0;
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
    ll sum = 0;
    while (eg != -1) {
      sum++;
      if (5 < sum) road.push_back(eg);
      eg = pre[eg];
    }
  }

  void construct_side(Graph &G, vll &road, vvll &side) {
    vll seen(H * W, 0);
    for (ll cu : road)
      seen[cu] = 1;

    queue<pll> qu;
    rep(i, 0, (ll)road.size()) qu.push({i, road[i]});
    vvll _side((ll)road.size());

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

    rep(i, 0, (ll)road.size()) {
      if (!_side[i].empty()) side.push_back(_side[i]);
    }
  }

  void solve(vector<Work> &ans) {
    Graph g;
    vll road;
    construct_road(g, road);

    vvll side;
    construct_side(g, road, side);

    SA_plant sa_p(side, 1.5);
    sa_p.solve(ans);
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

// local0.cpp

// 考察
/*
・K = HW なら
すべてs[i] = 1とし, 入口から近い順にD[i]が小さい順に植えるとできる.
・K > HW のとき,
道路を作る. 道路(road)には植えない.
道路の各マスから小道(side)を作り, 各マスを入口とみなし上をする.
*/
