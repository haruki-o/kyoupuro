#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
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
#define PI (3.14159265359)

double stat_time = 0.0;
double stat_loop = 0;
double stat_ans = 0;

ll L, N, S;
vll Y, X;
vll A;
ll Num = 5;
vll dy = {0, 1, 0, -1}, dx = {1, 0, -1, 0};

double ela_times(system_clock::time_point &clock) {
  return duration_cast<microseconds>(system_clock::now() - clock).count() *
         1e-6;
}

struct Pos {
  ll y, x;

  Pos(ll y, ll x) : y(y), x(x) {}

  Pos to(ll dir, ll num) {
    Pos to(y, x);
    to.y = (to.y + dy[dir] * num + L) % L;
    to.x = (to.x + dx[dir] * num + L) % L;
    return to;
  }
};

struct Solver {
  ll randing_cost(vvll &P) {
    ll ret = 0;
    rep(i, 0, L) rep(j, 0, L) rep(dir, 0, 2) {
      ret += pow(P[i][j] - P[(i + dy[dir]) % L][(j + dx[dir]) % L], 2);
    }
    return ret;
  }

  void optimize_P(vvll &P, vvll &fix) {
    rep(num, 0, 1000) {
      int flag = 1;
      rep(i, 0, L) {
        rep(j, 0, L) {
          if (fix[i][j] == 1)
            continue;
          ll sum = 0;
          rep(dir, 0, 4) sum += P[(i + dy[dir] + L) % L][(j + dx[dir] + L) % L];
          if (P[i][j] != sum / 4)
            flag = 0;
          P[i][j] = sum / 4;
        }
      }
      // cerr << randing_cost(P) << endl;
      if (flag) {
        // cerr << "optimize : " << num << endl;
        break;
      }
    }
  }

  void create_P(vvll &P, double limit) {
    auto start_time = system_clock::now();
    vvll fix(L, vll(L, 0));
    rep(i, 0, N) {
      fix[Y[i]][X[i]] = 1;
      P[Y[i]][X[i]] = max((ll)0, min((ll)1000, 1000 / N * i));
    }
    optimize_P(P, fix);

    ll best_score = randing_cost(P);
    while (1) {
      if (limit < ela_times(start_time))
        break;
      stat_loop++;

      vvll _P(L, vll(L));
      _P = P;
      ll idx1 = rand() % N, idx2 = rand() % N;
      if (idx1 == idx2)
        idx2 = (idx2 + rand() % N) % N;
      swap(_P[Y[idx1]][X[idx1]], _P[Y[idx2]][X[idx2]]);
      optimize_P(_P, fix);
      ll score = randing_cost(_P);
      if (chmin(best_score, score)) {
        P = _P;
      }
    }
  }

  void measure_P(vvll &P, vll &measured) {
    rep(i, 0, N) {
      rep(num, 0, Num) {
        cout << i << " " << 0 << " " << 0 << endl;
        cout.flush();
        ll m;
#ifdef _DEBUG
        random_device seed_gen;
        default_random_engine engine(seed_gen());
        normal_distribution<> dist(0.0, S);
        ll _y = Y[A[i]];
        ll _x = X[A[i]];
        m = max((ll)0, min((ll)1000, P[_y][_x] + (ll)dist(engine)));
#else
        cin >> m;
#endif
        measured[i] += m;
      }
      measured[i] /= Num;
    }
  }

  void predict(vvll &P, vll &measured, vll &E) {
    mcf_graph<int, ll> g(2 * N + 2);
    int s = 2 * N, t = 2 * N + 1;

    rep(i, 0, N) {
      g.add_edge(s, i, 1, 0);
      g.add_edge(N + i, t, 1, 0);
    }

    rep(i, 0, N) {
      rep(j, 0, N) {
        ll sum = abs(measured[i] - P[Y[j]][X[j]]);
        g.add_edge(i, N + j, 1, sum);
      }
    }

    auto result = g.flow(s, t, N);
    auto edges = g.edges();
    for (auto e : edges) {
      if (e.from == s || e.to == t || e.flow == 0)
        continue;
      E[e.from] = e.to - N;

      // ll i = e.from, j = e.to - N;
      // ll sum = 0;
      // Pos cu(Y[j], X[j]);
      // rep(dir, 0, 4) {
      //   rep(num, 0, length) {
      //     Pos to = cu.to(dir, num + 1);
      //     sum += abs(measured[i][dir][num] - P[to.y][to.x]);
      //     cerr << to.y << " " << to.x << " " << measured[i][dir][num] << " "
      //     << P[to.y][to.x] << endl;
      //   }
      // }
      // cerr << (A[i] == j) << " " << i << " " << j << " " << sum << endl;
    }
  }

  void solve() {
    vvll P(L, vll(L, 0));
    create_P(P, 2.0);
    rep(i, 0, L) {
      rep(j, 0, L) cout << P[i][j] << " ";
      cout << endl;
    }
    cout.flush();
    vll measured(N, 0);
    measure_P(P, measured);

    vll E(N);
    predict(P, measured, E);

    cout << "-1 -1 -1" << endl;
    rep(i, 0, N) cout << E[i] << endl;
#ifdef _DEBUG
    rep(i, 0, N) stat_ans += (E[i] == A[i]);
#endif
  }
};

int main() {
  srand((unsigned int)time(NULL));

  cin >> L >> N >> S;
  Num = 1;
  Y.resize(N);
  X.resize(N);
  rep(i, 0, N) cin >> Y[i] >> X[i];
#ifdef _DEBUG
  A.resize(N);
  rep(i, 0, N) cin >> A[i];
#endif

  Solver solver;
  solver.solve();

  cerr << "L N S : " << L << " " << N << " " << S << endl;
  // cerr << "stat_loop : " << stat_loop << endl;
  cerr << "stat_ans : " << stat_ans << " " << N - stat_ans << endl;
}

// lcoa1
// 同じところを何回も掘って平均
// ダメそう