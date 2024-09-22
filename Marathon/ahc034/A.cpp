#include <bits/stdc++.h>
using namespace std;
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

ll N;
vvll h;

vll dx = {0, -1, 0, 1}, dy = {1, 0, -1, 0};
vc dc = {'D', 'L', 'U', 'R'};

void decide_visit_v(ll cu, vll &seen, vll &visit_v) {
  ll cu_cost = 0;
  while (1) {
    vll dp(N * N, -INF);
    priority_queue<tuple<ll, ll, ll, ll>, vector<tuple<ll, ll, ll, ll>>> pq;
    pq.push({0, cu, 0, -1});
    dp[cu] = 0;
    seen[cu] = 1;
    while (!pq.empty()) {
      ll cost, cu2, have, dir;
      auto _ = pq.top();
      pq.pop();
      cost = get<0>(_);
      cu2 = get<1>(_);
      have = get<2>(_);
      dir = get<3>(_);
      if (dp[cu2] != -INF) continue;
      dp[cu2] = cost;

      rep(dir, 0, 4) {
        ll to_i = cu2 / N + dy[dir];
        ll to_j = cu2 % N + dx[dir];
        if (to_i < 0 || to_i >= N || to_j < 0 || to_j >= N) continue;
        ll to2 = to_i * N + to_j;
        if (dp[to2] != -INF) continue;
        pq.push({dp[cu2] + have + 100 + h[to2 / N][to2 % N], to2, dp[to2] + have, dir});
      }
    }
    ll best = INF;
    ll best_next;
    rep(i, 0, N) {
      rep(j, 0, N) {
        if (seen[i * N + j] == 0 && h[i][j] > 0) {
          if (chmin(best, dp[i * N + j])) best = i * N + j;
        }
      }
    }
    if (best == INF) break;
    visit_v.push_back(best);
    cu = best;
  }
}

void decide_all_v(vll visit_v) {
  ll have = 0;
  rep(i, 0, visit_v.size() - 1) {
    ll cu = visit_v[i];
    vll dp(N * N, -INF);
    vll dir_v(N * N, -1);
    priority_queue<tuple<ll, ll, ll, ll>, vector<tuple<ll, ll, ll, ll>>> pq;
    pq.push({0, cu, 0, -1});
    dp[cu] = 0;
    while (!pq.empty()) {
      ll cost, cu2, have, dir;
      auto _ = pq.top();
      pq.pop();
      cost = get<0>(_);
      cu2 = get<1>(_);
      have = get<2>(_);
      dir = get<3>(_);
      if (dp[cu2] != -INF) continue;
      dp[cu2] = cost;
      dir_v[cu2] = dir;
      rep(dir, 0, 4) {
        ll to_i = cu2 / N + dy[dir];
        ll to_j = cu2 % N + dx[dir];
        if (to_i < 0 || to_i >= N || to_j < 0 || to_j >= N) continue;
        ll to2 = to_i * N + to_j;
        if (dp[to2] != -INF) continue;
        pq.push({dp[cu2] + have + 100 + h[to2 / N][to2 % N], to2, dp[to2] + have, dir});
      }
    }

    cu = visit_v[i + 1];
    vll idx;
    cout << visit_v[i] / N << " " << visit_v[i] % N << ',';
    cout << visit_v[i + 1] / N << " " << visit_v[i + 1] % N << endl;
    idx.push_back(cu);
    while (cu != visit_v[i]) {
      cu = (cu / N + dy[(dir_v[cu] + 2) % 4]) * N + (cu % N + dx[(dir_v[cu] + 2) % 4]);
      idx.push_back(cu);
      cout << cu / N << " " << cu % N << endl;
    }
    rep(j, 0, idx.size() - 1) {
      ll _cu = idx[i], _to = idx[i + 1];
      rep(dir, 0, 4) {
        if ((_cu / N + dy[dir]) * N + _cu % N + dx[dir] == _to) {
          if (h[_to / N][_to % N] < 0) {
            cout << dc[dir] << endl;
            cout << '-' << -1 * min(-h[_to / N][_to % N], have) << endl;
            have -= min(-h[_to / N][_to % N], have);
          } else {
            cout << dc[dir] << endl;
            cout << '+' << h[_to / N][_to % N] << endl;
            have += -h[_to / N][_to % N];
          }
        }
      }
    }
  }
}

int main() {
  cin >> N;
  h.resize(N);
  rep(i, 0, N) h[i].resize(N);
  rep(i, 0, N) rep(j, 0, N) cin >> h[i][j];

  vll visit_v;
  vll seen(N * N, 0);
  decide_visit_v(0, seen, visit_v);

  decide_all_v(visit_v);
}
