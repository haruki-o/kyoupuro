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
double stat_wa = 0;
ll stat_beam_ac = 0;
ll measure_cost = 0;
ll placement_cost = 0;
ll beam_width = 5;
const double create_P_time = 3.5;
ll L, N, S;
vll Y, X;
vll A;
ll length = 0;

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

  ll eval_score(vvll &P) {
    ll ret = 0;
    ll ma = 0;

    rep(i, 0, N) {
      rep(j, i + 1, N) {
        ll sum = 0;
        Pos cui(Y[i], X[i]);
        Pos cuj(Y[j], X[j]);
        rep(dir, 0, 4) {
          rep(num, 1, length + 1) {
            Pos toi = cui.to(dir, num);
            Pos toj = cuj.to(dir, num);
            sum += pow(P[toi.y][toi.x] - P[toj.y][toj.x], 2);
          }
        }
        chmax(ma, sum);
      }
    }

    return ma;
  }

  void optimize_P(vvll &P, vvll &fix) {
    rep(num, 0, 1000) {
      int flag = 1;
      rep(i, 0, L) {
        rep(j, 0, L) {
          if (fix[i][j] == 1) continue;
          ll sum = 0;
          rep(dir, 0, 4) sum += P[(i + dy[dir] + L) % L][(j + dx[dir] + L) % L];
          if (P[i][j] != sum / 4) flag = 0;
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

  void create_P(vvll &P) {
    auto start_time = system_clock::now();
    vvll fix(L, vll(L, 0));
    rep(i, 0, N) {
      fix[Y[i]][X[i]] = 1;
      P[Y[i]][X[i]] = rand() % 1000;
    }
    optimize_P(P, fix);

    ll best_score = eval_score(P);
    while (1) {
      if (create_P_time < ela_times(start_time)) break;
      stat_loop++;

      vvll _P(L, vll(L));
      _P = P;
      ll idx = rand() % N;
      _P[Y[idx]][X[idx]] += rand() % 2 == 0 ? 50 : 950;
      _P[Y[idx]][X[idx]] %= 1000;
      optimize_P(_P, fix);
      ll score = eval_score(_P);
      // cerr << score << " " << best_score << endl;
      if (chmax(best_score, score)) {
        P = _P;
        // cerr << best_score << endl;
      }
    }
  }

  void measure_P(vvll &P, vvvll &measured) {
    rep(i, 0, N) {
      rep(dir, 0, 4) {
        rep(num, 1, length + 1) {
          cout << i << " " << dy[dir] * num << " " << dx[dir] * num << endl;
          measure_cost +=
              (ll)100 * ((ll)10 + abs(dy[dir] * num) + abs(dx[dir] * num));
          cout.flush();
          ll m;
#ifdef _DEBUG
          random_device seed_gen;
          default_random_engine engine(seed_gen());
          normal_distribution<> dist(0.0, S);
          ll _y = (Y[A[i]] + dy[dir] * num + L) % L;
          ll _x = (X[A[i]] + dx[dir] * num + L) % L;
          m = max((ll)0, min((ll)1000, P[_y][_x] + (ll)dist(engine)));
#else
          cin >> m;
#endif
          measured[i][dir][num - 1] = m;
        }
      }
    }
  }

  void predict(vvll &P, vvvll &measured, vll &E) {
    mcf_graph<int, ll> g(2 * N + 2);
    int s = 2 * N, t = 2 * N + 1;

    rep(i, 0, N) {
      g.add_edge(s, i, 1, 0);
      g.add_edge(N + i, t, 1, 0);
    }

    rep(i, 0, N) {
      rep(j, 0, N) {
        ll sum = 0;
        Pos cu(Y[j], X[j]);
        rep(dir, 0, 4) {
          rep(num, 0, length) {
            Pos to = cu.to(dir, num + 1);
            sum += abs(measured[i][dir][num] - P[to.y][to.x]);
          }
        }
        g.add_edge(i, N + j, 1, sum);
      }
    }

    auto result = g.flow(s, t, N);
    auto edges = g.edges();
    for (auto e : edges) {
      if (e.from == s || e.to == t || e.flow == 0) continue;
      E[e.from] = e.to - N;
    }

#ifdef _DEBUG
    vp ac_diff(N);
    for (auto e : edges) {
      if (e.from == s || e.to == t || e.flow == 0) continue;
      ll i = e.from;
      vp diff(N);
      rep(j, 0, N) {
        ll sum = 0;
        Pos cu(Y[j], X[j]);
        rep(dir, 0, 4) {
          rep(num, 0, length) {
            Pos to = cu.to(dir, num + 1);
            sum += abs(measured[i][dir][num] - P[to.y][to.x]);
          }
        }
        diff[j] = {sum, j};
      }
      Sort(diff);
      // cerr << "wormhole : " << i << " ";
      // cerr << endl;
      rep(j, 0, beam_width) {
        // cerr << "cell : " << diff[j].second << " " << diff[j].first;
        // if (A[i] == diff[j].second) cerr << " *";
        // cerr << endl;
        stat_beam_ac += (A[i] == diff[j].second);
      }
      ac_diff[i] = {diff[1].first - diff[0].first, (A[i] == diff[0].second)};
    }
    gSort(ac_diff);
    // rep(j, 0, N) cerr << ac_diff[j].first << " " << ac_diff[j].second <<
    // endl;
#endif
  }
  void solve() {
    vvll P(L, vll(L, 0));
    create_P(P);
    placement_cost = randing_cost(P);
    rep(i, 0, L) {
      rep(j, 0, L) cout << P[i][j] << " ";
      cout << endl;
    }
    cout.flush();
    vvvll measured(N, vvll(4, vll(length)));
    measure_P(P, measured);

    vll E(N);
    predict(P, measured, E);

    cout << "-1 -1 -1" << endl;
    rep(i, 0, N) cout << E[i] << endl;
#ifdef _DEBUG
    rep(i, 0, N) stat_wa += (E[i] != A[i]);
#endif
  }
};

int main() {
  srand((unsigned int)time(NULL));

  cin >> L >> N >> S;
  length = min(L - 1, (ll)pow(S, 0.5));
  Y.resize(N);
  X.resize(N);
  rep(i, 0, N) cin >> Y[i] >> X[i];
#ifdef _DEBUG
  A.resize(N);
  rep(i, 0, N) cin >> A[i];
#endif

  Solver solver;
  solver.solve();

  ll ans = 1e14 / (measure_cost + placement_cost + 1e5);
  rep(i, 0, stat_wa) ans *= 0.8;

  cerr << "L N S : " << L << " " << N << " " << S << endl;
  cerr << "stat_wa : " << stat_wa << endl;
  cerr << "stat_beam_wa : " << N - stat_beam_ac << endl;
  cerr << "score : " << ans << endl;
  cerr << "m_cost : " << measure_cost << endl;
  cerr << "p_cost : " << placement_cost << endl;
}

// local0 -> local2
// eval_score()からranding_cost()を除外
// i = 16 (S = 196)まではノーミスそう
// i = 20 (S = 400, beam_width = 5)まででミスが5個以内くらい

// 少し変更
// create_P()の初期値でoptimize_P()を追加
// create_P_timem beam_width(個までに正解があるか) 導入
// score, measure_cost, placement_cost, stat_beam_ac 導入
// stat_ans -> stat_wa
