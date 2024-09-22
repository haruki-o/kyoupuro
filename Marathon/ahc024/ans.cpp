#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef pair<ll, ll> P;
typedef tuple<ll, ll, ll> T;
typedef tuple<ll, ll, ll, ll> F;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vvi> vvvi;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
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
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
#define INF INT_MAX
#define INFF (9223372036854775800)
#define TIME_LIMIT (5.96)
#define def (101010)
// #define MOD (1000000007)
#define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision() << << endl;

ll N, M;
vvll c;
vll dy = {1, 0, -1, 0}, dx = {0, 1, 0, -1};

bool is_out(ll y, ll x) {
  if (y < 0 || y >= N || x < 0 || x >= N) return true;
  return false;
}

struct SA {
  double time;
  vvll g;
  SA(double time) : time(time), g(g) {}

  void init(vll &num) { rep(i, 0, M) num[i] = max((ll)1, num[i] / 5); }

  void modify(vll &num) {
    ll idx = M % rand();
    if (rand() % 2) {
      num[idx]++;
    } else {
      num[idx] = max((ll)1, num[idx] - 1);
    }
  }

  ll calc_score(vll &num){
    vvll ban(N, vll)
  }

  void solve() {
    auto start_time = system_clock::now();

    vll num(M, 0);
    init(num);
    ll best_score = calc_score(num);

    while (1) {
      if (time < ela_times(start_time)) break;

      vll _num(M, 0);
      _num = num;
      modify(_num);
      ll score = calc_score(_num);

      // 温度関数
      double temp = s_temp + (e_temp - s_temp) * ela_times(start_time) / time;
      // 遷移確率関数(最大化の場合)
      double prob = exp((score - best_score) / temp);
      // if (prob > (rand() % INF) / (double)INF) { // 確率probで遷移する
      if (best_score < score) {
        best_score = score;
        num = _num;
      }
    }
  }
};

int main() {
  srand((unsigned int)time(NULL));

  cin >> N >> M;
  c.resize(N);
  rep(i, 0, N) {
    c[i].resize(N);
    rep(j, 0, N) cin >> c[i][j], c[i][j]--;
  }

  vvll g(N);
  vvi used(N, vll(N, 0));
  rep(i, 0, N) {
    rep(j, 0, N) {
      ll at_c = c[i][j];
      rep(dir, 0, 4) {
        ll toi = i + dy[dir];
        ll toj = j + dx[dir];
        ll to_c = c[toi][toj];
        if (is_out(toi, toj)) continue;
        if (at_c == to_c) continue;
        if (used[at_c][to_c]) continue;
        g[at_c].push_back(to_c);
        g[to_c].push_back(at_c);
        used[at_c][to_c] = 1;
        used[to_c][at_c] = 1;
      }
    }
  }

  SA sa(1.9, g);
  sa.solve();
}