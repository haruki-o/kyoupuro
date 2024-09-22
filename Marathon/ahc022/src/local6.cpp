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

double Combi(int n, int y) {
  double v = 1.0;
  if ((1 <= y) && (y < n))
    for (double i = y; i > 0.99; i -= 1.0)
      v *= (n + 1.0 - i) / i;
  return v;
}

double bio_p(int n, int k, double p) {
  return Combi(n, k) * pow(p, k) * pow(1.0 - p, n - k);
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

bool operator<(const Pos &p1, const Pos &p2) {
  return p1.y < p2.y || (p1.y == p2.y && p1.x <= p2.x);
}

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

void optimize_P(vvll &P, vvll &fix) {
  ll sum = 0;
  rep(i, 0, L) rep(j, 0, L) sum += fix[i][j];
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

ll optimize_K(ll N, double p1, double p2) {
  vector<double> c_p1(N + 1, 0), c_p2(N + 1, 0);
  double sum = 0.0;
  rep(i, 0, N + 1) {
    sum += bio_p(N, i, p1);
    c_p1[i] = sum;
  }
  sum = 0;
  rep(i, 0, N + 1) {
    sum += bio_p(N, i, p2);
    c_p2[i] = sum;
  }

  double ma = 0.0;
  int K = 0;
  rep(i, 0, N + 1) if (chmax(ma, min(c_p1[i - 1], 1 - c_p2[i - 1]))) K = i;
  return K;
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
  ll Min_P, Max_P;
  const double create_P_time = 3.5;
  ll Success_num = 1;
  ll Bio_N;
  ll Bio_K;
  ll Th;

  ll m_N;
  vector<Pos> moves;

  vector<Pos> L_Min, L_Max;

  Solver_200_P() {
    Min_P = max((ll)0, (ll)500 - S / 2 - (ll)(0.25 * S));
    Max_P = min((ll)1000, (ll)500 + S / 2 + (ll)(0.25 * S));
    Th = Min_P;
    rep(i, 2, 12) {
      if (N <= pow(2, i - 1)) {
        m_N = i;
        break;
      }
    }
    Bio_N = 10000 / (m_N * N);
    Bio_K = optimize_K(Bio_N, 0.5, 0.945);
    cerr << Bio_N << " " << Bio_K << endl;

    // L_Min.push_back(Pos(L / 4, L / 4));
    // L_Min.push_back(Pos(L / 4, 3 * L / 4));
    // L_Min.push_back(Pos(3 * L / 4, L / 4));
    // L_Min.push_back(Pos(3 * L / 4, 3 * L / 4));
    // L_Max.push_back(Pos(L / 4, L / 2));
    // L_Max.push_back(Pos(L / 2, L / 4));
    // L_Max.push_back(Pos(L / 2, L / 2));
    // L_Max.push_back(Pos(L / 2, 3 * L / 4));
    // L_Max.push_back(Pos(3 * L / 4, L / 2));
    rep(i, 0, L) {
      rep(j, 0, L) {
        if (L / 6 <= i && i <= L / 3 || 2 * L / 3 <= i && i <= 5 * L / 6)
          L_Max.push_back(Pos(i, j));
        else if (L / 6 <= j && j <= L / 3 || 2 * L / 3 <= j && j <= 5 * L / 6)
          L_Max.push_back(Pos(i, j));
        else L_Min.push_back(Pos(i, j));
      }
    }
  }

  void init(vector<Pos> &_moves) {
    map<pll, ll> ma;
    rep(i, 0, m_N) {
      while (1) {
        Pos p1(rand() % L, rand() % L);
        Pos p2(rand() % L, rand() % L);
        Pos _move = move(p1, p2);
        if (ma.count({_move.y, _move.x}) == 0) {
          ma[{_move.y, _move.x}] = 1;
          _moves.push_back(_move);
          break;
        }
      }
    }
  }

  void modify(vector<Pos> &_moves) {
    ll idx = rand() % m_N;
    map<pll, ll> ma;
    rep(i, 0, m_N) if (i != idx) ma[{_moves[i].y, _moves[i].x}];
    while (1) {
      Pos p1(rand() % L, rand() % L);
      Pos p2(rand() % L, rand() % L);
      Pos _move = move(p1, p2);
      if (ma.count({_move.y, _move.x}) == 0) {
        ma[{_move.y, _move.x}] = 1;
        _moves[idx] = _move;
        break;
      }
    }
    // ll diff = rand() % 2 == 0 ? 1 : -1;
    // if (rand() % 2 == 0) {
    //   _moves[idx].y += diff;
    //   if (L / 2 < _moves[idx].y) _moves[idx].y = -L / 2;
    //   if (_moves[idx].y < -L / 2) _moves[idx].y = L / 2;
    // } else {
    //   _moves[idx].x += diff;
    //   if (L / 2 < _moves[idx].x) _moves[idx].x = -L / 2;
    //   if (_moves[idx].x < -L / 2) _moves[idx].x = L / 2;
    // }
  }

  ll calc_score(vector<Pos> &_moves, vvll &P) {
    ll ret = 0;

    vvvll over(L, vvll(L));
    vvll fix(L, vll(L, 0));
    vector<pair<ll, Pos>> all;

    rep(i, 0, N) {
      rep(j, 0, m_N) {
        Pos to =
            Pos((Y[i] + _moves[j].y + L) % L, (X[i] + _moves[j].x + L) % L);
        over[to.y][to.x].push_back(i);
      }
    }

    rep(i, 0, L) {
      rep(j, 0, L) {
        if ((ll)over[i][j].size() != 0) {
          all.push_back({(ll)over[i][j].size(), Pos(i, j)});
          fix[i][j] = 1;
        }
      }
    }
    gSort(all);

    // fixのPの値を決める
    map<ll, ll> ma;
    ll pre_P = 1;
    vvvll seen(L, vvll(L, vll(m_N, -1)));

    for (auto _to : all) {
      Pos to = _to.second;
      vll ok(2, 1);
      rep(_ok, 0, 2) {
        map<ll, ll> ma2;
        rep(i, 0, (ll)over[to.y][to.x].size()) {
          ll idx = over[to.y][to.x][i];
          Pos cu = Pos(Y[idx], X[idx]);
          ll sum = 0;
          ll bit = 0;
          rep(j, 0, m_N) {
            if (seen[cu.y][cu.x][j] == -1) sum++;
            else bit += pow(2, j) * seen[cu.y][cu.x][j];
          }
          if (sum == 1) {
            rep(j, 0, m_N) {
              if (seen[cu.y][cu.x][j] == -1) {
                bit += pow(2, j) * _ok;
                if (ma.count(bit) || ma2.count(bit)) ok[_ok] = 0;
                ma2[bit] = 1;
              }
            }
          }
        }
      }

      if (ok[0] == 0 && ok[1] == 0) return INFF;

      // ll _ok = ok[1 - pre_P] == 1 ? 1 - pre_P : pre_P;
      ll _ok = ok[0] == 1 ? 0 : 1;
      if (ok[0] == 1 && ok[1] == 1) {
        _ok = ok[1 - pre_P] == 1 ? 1 - pre_P : pre_P;
        ll mi = 5;
        // for (auto _cu : L_Min)
        //   if (chmin(mi, move(_cu, to).m_dist)) _ok = 0;
        // for (auto _cu : L_Max)
        //   if (chmin(mi, move(_cu, to).m_dist)) _ok = 1;
        rep(i, -2, 3) {
          rep(j, -2, 3) {
            Pos _to((to.y + i + L) % L, (to.x + j + L) % L);
            if (chmin(mi, move(to, _to).m_dist)) {
              if (P[_to.y][_to.x] == Min_P) _ok = 0;
              if (P[_to.y][_to.x] == Max_P) _ok = 1;
              // if(P[_to.y][_to.x] == Min_P || P[_to.y][_to.x] == Max_P)cerr <<
              // "a";
            }
          }
        }
      }

      rep(i, 0, (ll)over[to.y][to.x].size()) {
        ll idx = over[to.y][to.x][i];
        Pos cu = Pos(Y[idx], X[idx]);
        ll sum = 0;
        ll bit = 0;
        rep(j, 0, m_N) {
          if (seen[cu.y][cu.x][j] == -1) {
            sum++;
            bit += pow(2, j) * _ok;
          } else bit += pow(2, j) * seen[cu.y][cu.x][j];

          Pos to_((cu.y + _moves[j].y + L) % L, (cu.x + _moves[j].x + L) % L);
          if (to.y == to_.y && to.x == to_.x) {
            seen[cu.y][cu.x][j] = _ok;
          }
        }
        if (sum == 1) ma[bit] = 1;
      }
      pre_P = _ok;
      if (_ok == 1) P[to.y][to.x] = Max_P;
      else P[to.y][to.x] = Min_P;
    }

#ifdef _DEBUG
    map<ll, ll> d_ma;
    rep(i, 0, N) {
      ll bit = 0;
      rep(j, 0, m_N) {
        Pos to =
            Pos((Y[i] + _moves[j].y + L) % L, (X[i] + _moves[j].x + L) % L);
        if (P[to.y][to.x] == Max_P) bit += pow(2, j);
      }
      if (d_ma.count(bit)) cerr << "err 重複があります" << endl;
      d_ma[bit] = 1;
    }
#endif

    // fix以外のPの値を決める
    optimize_P(P, fix);
    // 実際にスコア計算
    ret = randing_cost(P);
    return ret;
  }

  void create_P(vvll &P) {
    auto start_time = system_clock::now();

    vector<Pos> best_moves;
    init(best_moves);
    ll best_score = calc_score(best_moves, P);
    while (1) {
      if (create_P_time < ela_times(start_time)) break;
      // stat_loop++;

      vector<Pos> _moves;
      _moves = best_moves;
      // init(_moves);
      vvll _P(L, vll(L, 0));
      modify(_moves);
      ll score = calc_score(_moves, _P);
      if (score != INFF) stat_loop++;
      if (chmin(best_score, score)) {
        best_moves = _moves;
        P = _P;
        // cerr << best_score << endl;
      }
    }

    moves = best_moves;

    rep(i, 0, L) {
      rep(j, 0, L) cout << P[i][j] << " ";
      cout << endl;
    }
    cout.flush();
  }

  void measure_P(vvll &P, vll &E) {
    rep(i, 0, N) {
      vll E_P(m_N, 0);
      rep(j, 0, m_N) {
        ll sum = 0;
        rep(num, 0, Bio_N) {
          cout << i << " " << moves[j].y << " " << moves[j].x << endl;
          measure_cost += 100 * (10 + moves[j].m_dist);
          cout.flush();

          ll m;
#ifdef _DEBUG
          random_device seed_gen;
          default_random_engine engine(seed_gen());
          normal_distribution<> dist(0.0, S);
          ll r = (Y[A[i]] + moves[j].y + L) % L;
          ll c = (X[A[i]] + moves[j].x + L) % L;
          m = max((ll)0, min((ll)1000, P[r][c] + (ll)dist(engine)));
#else
          cin >> m;
#endif

          if (Th <= m) sum++;
        }
        if (Bio_K <= sum) E_P[j] = 1;
      }

      ll ma = -1;
      rep(j, 0, N) {
        ll sum = 0;
        rep(k, 0, m_N) {
          Pos to((Y[j] + moves[k].y + L) % L, (X[j] + moves[k].x + L) % L);
          if (E_P[k] == 0 && P[to.y][to.x] == Min_P) sum++;
          if (E_P[k] == 1 && P[to.y][to.x] == Max_P) sum++;
        }
        if (chmax(ma, sum)) E[i] = j;
      }

#ifdef _DEBUG
      // if (A[i] == E[i]) cerr << "correct !" << endl;
      // else cerr << "out" << endl;
      // rep(k, 0, m_N) {
      //   Pos to((Y[A[i]] + moves[k].y + L) % L, (X[A[i]] + moves[k].x + L) %
      //   L); if (P[to.y][to.x] == Max_P) cerr << 1; if (P[to.y][to.x] ==
      //   Min_P) cerr << 0;
      // }
      // cerr << endl;

      // rep(j, 0, N) {
      //   ll sum = 0;
      //   ll diff = 0;
      //   rep(k, 0, m_N) {
      //     Pos to((Y[j] + moves[k].y + L) % L, (X[j] + moves[k].x + L) % L);
      //     Pos c_to((Y[A[i]] + moves[k].y + L) % L,
      //              (X[A[i]] + moves[k].x + L) % L);
      //     if (P[to.y][to.x] != P[c_to.y][c_to.x]) diff++;
      //     if (E_P[k] == 0 && P[to.y][to.x] == Min_P) sum++;
      //     if (E_P[k] == 1 && P[to.y][to.x] == Max_P) sum++;
      //     if (P[to.y][to.x] == Max_P) cerr << 1;
      //     if (P[to.y][to.x] == Min_P) cerr << 0;
      //   }
      //   cerr << " " << diff << " " << sum << endl;
      // }
#endif
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
    // return 0;
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

// local5 -> local6

// 概要
/*
どっちでもいいとき周囲を探して同じものにしました
*/