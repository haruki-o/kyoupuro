#include <atcoder/all>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
typedef pair<ll, ll> pll;
typedef tuple<ll, ll, ll, ll> F;
typedef vector<pll> vp;
typedef vector<vp> vvp;
typedef vector<F> vF;
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
double stat_wa = 0;
ll measure_cost = 0;
ll placement_cost = 0;

ll L, N, S;
vll Y, X;
vll A;

vll dy = {0, 1, 0, -1}, dx = {1, 0, -1, 0};

double ela_times(system_clock::time_point &clock) {
  return duration_cast<microseconds>(system_clock::now() - clock).count() *
         1e-6;
}

struct Pos {
  ll y, x;
  ll m_dist;

  Pos() {}

  Pos(ll y, ll x) : y(y), x(x) { m_dist = abs(y) + abs(x); }

  Pos to(ll dir, ll num) {
    Pos to(y, x);
    to.y = (to.y + dy[dir] * num + L) % L;
    to.x = (to.x + dx[dir] * num + L) % L;
    return to;
  }
};

ll randing_cost(vvll &P) {
  ll ret = 0;
  rep(i, 0, L) rep(j, 0, L) rep(dir, 0, 2) {
    ret += pow(P[i][j] - P[(i + dy[dir]) % L][(j + dx[dir]) % L], 2);
  }
  return ret;
}

Pos move(Pos &from, Pos &to) {
  Pos ret;
  ll max_y = max(from.y, to.y);
  ll min_y = min(from.y, to.y);
  // トーラスなし
  if (max_y - min_y <= L / 2) ret.y = to.y - from.y;
  // トーラスあり
  else {
    // 上行って下にワープ
    if (max_y == to.y) ret.y = (max_y - min_y) - L;
    else ret.y = L - (max_y - min_y);
  }

  ll max_x = max(from.x, to.x);
  ll min_x = min(from.x, to.x);
  // トーラスなし
  if (max_x - min_x <= L / 2) ret.x = to.x - from.x;
  // トーラスあり
  else {
    // 上行って下にワープ
    if (max_x == to.x) ret.x = (max_x - min_x) - L;
    else ret.x = L - (max_x - min_x);
  }

  ret.m_dist = abs(ret.y) + abs(ret.x);
  return ret;
}

struct Solver {
  ll length;
  ll stat_beam_ac = 0;
  ll beam_width = 5;
  const double create_P_time = 3.5;

  Solver() { length = min((ll)10000 / 4 / N, min(L - 1, (ll)pow(S, 0.5))); }

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

    ll best_eval_score = eval_score(P);
    ll best_randing_score = randing_cost(P);
    while (1) {
      if (create_P_time < ela_times(start_time)) break;
      stat_loop++;

      vvll _P(L, vll(L));
      _P = P;
      ll idx = rand() % N;
      _P[Y[idx]][X[idx]] += rand() % 2 == 0 ? 50 : 950;
      _P[Y[idx]][X[idx]] %= 1000;
      optimize_P(_P, fix);

      if (S <= 200) {
        ll _eval_score = eval_score(_P);
        ll _randing_score = randing_cost(_P);
        if (best_eval_score <= _eval_score &&
            _randing_score <= best_randing_score) {
          best_eval_score = _eval_score;
          best_randing_score = _randing_score;
          P = _P;
        }
      } else {
        ll _eval_score = eval_score(_P);
        if (chmax(best_eval_score, _eval_score)) P = _P;
      }
    }
    rep(i, 0, L) {
      rep(j, 0, L) cout << P[i][j] << " ";
      cout << endl;
    }
    cout.flush();
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
          ll r = (Y[A[i]] + dy[dir] * num + L) % L;
          ll c = (X[A[i]] + dx[dir] * num + L) % L;
          m = max((ll)0, min((ll)1000, P[r][c] + (ll)dist(engine)));
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

struct Solver_200_P {
  ll Max_N;
  ll Min_P, Max_P;
  vll max_P_Y, max_P_X;
  const double create_P_time = 3.5;
  ll Success_num = 1;
  ll Bio_N;
  ll Bio_K;
  ll Th;

  Solver_200_P() {
    Max_N = max(Success_num, L / 4);
    // Min_P = max((ll)0, (ll)500 - S);
    // Max_P = min((ll)1000, (ll)500 + S);
    Min_P = 0, Max_P = 1000;
    Bio_N = (ll)10000 / (N * N);
    Bio_K = Bio_N;
    Th = 500;
  }

  void init(vll &P_Y, vll &P_X) {
    ll sum = 0;
    vvll seen(L, vll(L, 0));
    while (sum != Max_N) {
      ll i = rand() % L;
      ll j = rand() % L;
      if (seen[i][j] == 0) {
        P_Y.push_back(i);
        P_X.push_back(j);
        seen[i][j] = 1;
        sum++;
      }
    }
  }

  void modify(vll &P_Y, vll &P_X) {
    ll _N = (ll)P_Y.size();
    if (rand() % 2 == 0) {
      ll idx = rand() % _N;
      P_Y[idx] += rand() % 2 == 0 ? 1 : L + 1;
      P_Y[idx] %= L;
    } else {
      ll idx = rand() % _N;
      P_X[idx] += rand() % 2 == 0 ? 1 : L + 1;
      P_X[idx] %= L;
    }
  }

  ll calc_score(vll &P_Y, vll &P_X) {
    ll s_N = (ll)P_Y.size();
    ll ret = 0;
    ret += pow((Max_P - Min_P), 2) * s_N;

    vvp all(N, vp(s_N));
    rep(i, 0, N) {
      rep(j, 0, s_N) {
        Pos pi(Y[i], X[i]);
        Pos pj(P_Y[j], P_X[j]);
        Pos mov = move(pi, pj);
        all[i][j] = {mov.m_dist, j};
      }
      Sort(all[i]);
    }

    // 上位 Success_num個の中で被ると減点
    rep(i, 0, N) {
      rep(j, i + 1, N) {
        ll same_num = 0;
        rep(sn, 0, Success_num) {
          // i -> _i
          Pos pi(Y[i], X[i]);
          Pos _pi(P_Y[all[i][sn].second], P_X[all[i][sn].second]);
          // j -> _j
          Pos pj(Y[j], X[j]);
          Pos _pj(P_Y[all[j][sn].second], P_X[all[j][sn].second]);

          Pos move_i = move(pi, _pi);
          Pos move_j = move(pj, _pj);
          if (move_i.y == move_j.y && move_i.x == move_j.x) {
            same_num++;
            ret += 100 * (10 + move_i.m_dist) * 3;
          }
          if (sn == 0) ret += 100 * (10 + move_i.m_dist);
          if (sn == 0) ret += 100 * (10 + move_j.m_dist);
        }
        if (same_num == Success_num) ret = 1e10;
      }
    }

    return ret;
  }

  void create_P(vvll &P) {
    auto start_time = system_clock::now();

    vll best_P_Y, best_P_X;
    init(best_P_Y, best_P_X);
    ll best_score = calc_score(best_P_Y, best_P_X);
    while (1) {
      if (create_P_time < ela_times(start_time)) break;

      vll P_Y, P_X;
      P_Y = best_P_Y;
      P_X = best_P_X;
      modify(P_Y, P_X);
      ll score = calc_score(P_Y, P_X);
      if (chmin(best_score, score)) {
        best_P_Y = P_Y;
        best_P_X = P_X;
      }
    }

    max_P_Y = best_P_Y;
    max_P_X = best_P_X;
    rep(i, 0, L) P[i].assign(L, Min_P);
    rep(i, 0, Max_N) P[max_P_Y[i]][max_P_X[i]] = Max_P;

    rep(i, 0, L) {
      rep(j, 0, L) cout << P[i][j] << " ";
      cout << endl;
    }
    cout.flush();
  }

  void measure_P(vvll &P, vll &E) {
    ll s_N = (ll)max_P_Y.size();
    vvp all(N, vp(s_N));
    rep(i, 0, N) {
      rep(j, 0, s_N) {
        Pos pi(Y[i], X[i]);
        Pos pj(max_P_Y[j], max_P_X[j]);
        Pos mov = move(pi, pj);
        all[i][j] = {mov.m_dist, j};
      }
      Sort(all[i]);
    }

    /*
    rep(i, 0, N) {
      queue<ll> qu;
      rep(j, 0, N) qu.push(j);
      rep(sn, 0, Success_num) {
        if (qu.empty()) break;

        queue<ll> _qu;
        while (!qu.empty()) {
          ll j = qu.front();
          qu.pop();
          Pos pj = Pos(Y[j], X[j]);
          Pos _pj(max_P_Y[all[j][sn].second], max_P_X[all[j][sn].second]);

          Pos mov = move(pj, _pj);
          ll sum = 0;
          rep(b, 0, Bio_N) {
            cout << i << " " << mov.y << " " << mov.x << endl;
            ll r = (max_P_Y[all[A[i]][sn].second] + mov.y + L) % L;
            ll c = (max_P_X[all[A[i]][sn].second] + mov.x + L) % L;
            ll m = input_m(P[r][c]);
            if (Min_P <= m) sum++;
            cerr << P[max_P_Y[all[j][sn].second]][max_P_X[all[j][sn].second]];
            cerr << " " << Min_P << " " << m << endl;
          }
          if (Bio_K <= sum) _qu.push(j);
        }
        while (!_qu.empty()) {
          ll cu = _qu.front();
          _qu.pop();
          qu.push(cu);
        }
      }
      qu.push(0);
      E[i] = qu.front();
    }
    */
    rep(i, 0, N) {
      ll ma = -1e5;
      rep(j, 0, N) {
        Pos pj = Pos(Y[j], X[j]);
        Pos _pj(max_P_Y[all[j][0].second], max_P_X[all[j][0].second]);

        Pos mov = move(pj, _pj);
        ll sum = 0;
        rep(b, 0, Bio_N) {
          cout << i << " " << mov.y << " " << mov.x << endl;
          measure_cost += 100 * (10 + mov.m_dist);

          ll m;
#ifdef _DEBUG
          random_device seed_gen;
          default_random_engine engine(seed_gen());
          normal_distribution<> dist(0.0, S);
          ll r = (Y[A[i]] + mov.y + L) % L;
          ll c = (X[A[i]] + mov.x + L) % L;
          m = max((ll)0, min((ll)1000, P[r][c] + (ll)dist(engine)));
#else
          cin >> m;
#endif

          sum += m;
          if (Max_P <= m) sum += S;
          if (m <= Min_P) sum -= S;
          cerr << P[r][c] << " " << Min_P << " " << m << endl;
        }
        cerr << i << " " << j << " " << A[i] << " " << sum << endl;
        if (chmax(ma, sum)) E[i] = j;
      }
      cerr << E[i] << endl;
    }
  }

  void solve() {
    vvll P(L, vll(L, 0));
    create_P(P);
    placement_cost = randing_cost(P);

    vll E(N, -1);
    measure_P(P, E);

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
  Y.resize(N);
  X.resize(N);
  rep(i, 0, N) cin >> Y[i] >> X[i];
#ifdef _DEBUG
  A.resize(N);
  rep(i, 0, N) cin >> A[i];
#endif

  if (S < 200) {
    Solver solver;
    solver.solve();
  } else {
    Solver_200_P solver;
    solver.solve();
  }

  ll ans = 1e14 / (measure_cost + placement_cost + 1e5);
  rep(i, 0, stat_wa) ans *= 0.8;

  cerr << "L N S : " << L << " " << N << " " << S << endl;
  cerr << "stat_loop : " << stat_loop << endl;
  cerr << "stat_wa : " << stat_wa << endl;
  cerr << "score : " << ans << endl;
  cerr << "m_cost : " << measure_cost << endl;
  cerr << "p_cost : " << placement_cost << endl;
}

// local3 -> local4
/*
p >= 200 のとき, 数か所 Max_P, その他 Min_Pにする.
どちらか当てる確率を0.95にしたい.
(下記はMin_Pの上側確率 (閾値 Th =Min_Pとする))
Max_P - Min_P = S のとき, 50%, 84% (2 * S, 50%, 98%) 二項分布のN を決めた時,
p = 0.5の下側確率と p = 0.84の上側確率が一緒の時の K を求める
prob.cpp は N, p1, p2を与えると(K, p)を出力します.
N = 15で, K = 11, p = 0.922 (表11回以上の時 Max_P と答えると p で正解) 
N = 20で, K = 14, p = 0.942
*/

// 質問回数の上限10^4を忘れていて, 上をやろうとすると  N * N * 3 * 20 = 4 *
// 10^5ぐらいで破綻 10^4
// (N^2)だとするとどちらか当てるのは1~4回しかmeassureできない 0-   : 0.500 0-σ :
// 0.341 σ-2σ: 0.136 2σ-  : 0.022 ただ平均値を使うようにしました この解法は, N *
// Nがダメそう

// chatgpt
/*
表裏がでる2種類のコインA,Bがあり、それぞれの表の出る確率は50%, 84%です.
種類のわからないコインが渡されます.
N回投げたのちどちらかのコインかを95%以上で当ててください.
以上を満たす最小のNを教えてください。
*/

// chatgpt
/*
平均800, 標準偏差400の正規分布に従う生成器Aと,
平均0, 標準偏差400の正規分布に従う生成器Bがあります.
あなたは, どちらの種類かわからない生成器が渡されます.
5回まで与えられた生成期から生成することができます.
その後どちらが与えられたか出力しなさい.
以上の問題を高精度で判断するアルゴリズムはありますか？
-> 平均値が高い方がいいらしい
間違ってそう
*/