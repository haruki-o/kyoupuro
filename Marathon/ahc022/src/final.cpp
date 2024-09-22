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
ll stat_wa = 0;
ll stat_wa2 = 0;
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

ll optimize_K(ll N, double p1, double p2, double &p) {
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
  p = ma;
  return K;
}

/* void optimize(ll &Bio_K, ll Bio_N, ll m_N, ll &Min_P, ll &Max_P, ll &Th,
              double &_p) {
  double ma = 0.0;
  rep(diff, 2, 1001) {
    // rep(_Th, 500 - diff / 2, 500 - diff / 2 + 1) {
    rep(_Th, 1, 1000) {
      ll _Min_P = 500 - diff / 2;
      ll _Max_P = 500 + diff / 2;
      double z1 = (_Th - _Min_P) / (double)S;
      double prob1 = 0.5 * (1.0 + std::erf(-z1 / std::sqrt(2.0)));
      double z2 = (_Th - _Max_P) / (double)S;
      double prob2 = 0.5 * (1.0 + std::erf(-z2 / std::sqrt(2.0)));
      double p;
      Bio_K = optimize_K(Bio_N, prob1, prob2, p);
      // cerr << prob1 << " " << prob2 << " " << p << " " << z1 << " " << z2
      //  << endl;
      if (0.95 <= p) {
        Min_P = _Min_P;
        Max_P = _Max_P;
        _p = p;
        Th = _Th;
        // Th = Min_P;
        cerr << prob1 << " " << prob2 << " " << p << " " << z1 << " ";
        cerr << z2 << " " << Th << endl;
        return;
      }
    }
  }
  Min_P = 0;
  Max_P = 1000;
  double p;
  Bio_K = optimize_K(Bio_N, 0.5,
                     0.5 * (1.0 + std::erf((ll)1000 / S / std::sqrt(2.0))), p);
  _p = p;
}
*/

void optimize(ll &Bio_K, ll Bio_N, ll m_N, ll &Min_P, ll &Max_P, ll &Th,
              double &_p) {
  double ma = 0.0;
  rep(diff, 2, 1001) {
    // rep(_Th, 500 - diff / 2, 500 - diff / 2 + 1) {
    ll _Min_P = 500 - diff / 2;
    ll _Max_P = 500 + diff / 2;
    double _S = S / sqrt(Bio_N);
    double z = (500 - _Max_P) / _S;
    double prob = 0.5 * (1.0 + std::erf(-z / std::sqrt(2.0)));
    double p;

    if (0.9999 <= prob) {
      Min_P = _Min_P;
      Max_P = _Max_P;
      _p = prob;
      Th = 500;
      Bio_K = 1;
      return;
    }
  }
  Min_P = 0;
  Max_P = 1000;
  double p;
  Bio_K = optimize_K(Bio_N, 0.5,
                     0.5 * (1.0 + std::erf((ll)1000 / S / std::sqrt(2.0))), p);
  _p = p;
}

struct SA_init {
  vector<Pos> moves;
  vvll P;
  ll m_N, Min_P, Max_P;
  double time;

  SA_init(vector<Pos> moves, vvll P, ll m_N, ll Min_P, ll Max_P, double time)
      : moves(moves), P(P), m_N(m_N), Min_P(Min_P), Max_P(Max_P), time(time) {}

  bool init(vector<Pos> &_moves, vvll &P) {
    map<pll, ll> _ma;
    rep(i, 0, m_N) {
      while (1) {
        Pos p1(rand() % L, rand() % L);
        Pos p2(rand() % L, rand() % L);
        Pos _move = move(p1, p2);
        if ((ll)5 <= _move.m_dist) continue;
        if (_ma.count({_move.y, _move.x}) == 0) {
          _ma[{_move.y, _move.x}] = 1;
          _moves.push_back(_move);
          break;
        }
      }
    }

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

      // if (ok[0] == 0 && ok[1] == 0) return false;

      ll _ok = ok[1 - pre_P] == 1 ? 1 - pre_P : pre_P;
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

    map<ll, ll> d_ma;
    rep(i, 0, N) {
      ll bit = 0;
      rep(j, 0, m_N) {
        Pos to =
            Pos((Y[i] + _moves[j].y + L) % L, (X[i] + _moves[j].x + L) % L);
        if (P[to.y][to.x] == Max_P) bit += pow(2, j);
      }
      d_ma[bit]++;
      if (8 < d_ma[bit]) return false;
    }

    return true;
  }

  ll calc_score(vector<Pos> &_moves) {
    ll ret = 0;

    vvll over(L, vll(L, 0));
    rep(i, 0, N) {
      rep(j, 0, m_N) {
        Pos to((Y[i] + _moves[j].y + L) % L, (X[i] + _moves[j].y + L) % L);
        over[to.y][to.x]++;
      }
    }
    rep(i, 0, L) rep(j, 0, L) if (over[i][j] > 1) ret++;

    ll _sum = 0;
    rep(i, 0, L) rep(j, 0, L) if (over[i][j] == 1) _sum++;
    // cerr << L * L - ret - _sum << " " << _sum << " " << ret << endl;

    return ret;
  }

  void solve(vector<Pos> &ret_moves, vvll &ret_P) {
    auto start_time = system_clock::now();
    ll best_score = INFF;
    while (1) {
      if (time < ela_times(start_time)) break;

      vector<Pos> _moves;
      vvll _P(L, vll(L, -1));

      if (init(_moves, _P)) {
        ll score = calc_score(_moves);
        if (chmin(best_score, score)) {
          moves = _moves;
          P = _P;
        }
        break;
      }
    }
    ll _score = calc_score(moves);

    ret_moves = moves;
    ret_P = P;
  }
};

struct SA_place {
  vector<Pos> moves;
  vvll P;
  ll m_N, Min_P, Max_P;
  double time;
  vvll fix;

  SA_place(vector<Pos> moves, vvll P, ll m_N, ll Min_P, ll Max_P, double time)
      : moves(moves), P(P), m_N(m_N), Min_P(Min_P), Max_P(Max_P), time(time) {

    fix.resize(L);
    rep(i, 0, L) fix[i].assign(L, 0);
    rep(i, 0, N) {
      rep(j, 0, m_N) {
        Pos _to((Y[i] + moves[j].y + L) % L, (X[i] + moves[j].x + L) % L);
        fix[_to.y][_to.x] = 1;
      }
    }
    optimize_P(P, fix);
  }

  bool modify(vvll &_P, ll over_ma, ll _stat_loop) {
    ll idx = rand() % N;
    ll bit_idx = rand() % m_N;
    Pos to((Y[idx] + moves[bit_idx].y + L) % L,
           (X[idx] + moves[bit_idx].x + L) % L);

    if (_P[to.y][to.x] == Min_P) _P[to.y][to.x] = Max_P;
    else if (_P[to.y][to.x] == Max_P) _P[to.y][to.x] = Min_P;

    map<ll, ll> d_ma;
    rep(i, 0, N) {
      ll sum = 0;
      rep(j, 0, m_N) {
        Pos _to((Y[i] + moves[j].y + L) % L, (X[i] + moves[j].x + L) % L);
        if (_P[_to.y][_to.x] == Max_P) sum += pow(2, j);
      }
      d_ma[sum]++;
      if (over_ma < d_ma[sum]) return false;
    }

    if (_stat_loop % 10 == 0) optimize_P(_P, fix);
    return true;
  }

  ll calc_score(vvll &P) {
    ll ret = randing_cost(P);
    return ret;
  }

  void solve(vvll &ret_P) {
    auto start_time = system_clock::now();

    ll best_score = calc_score(P);
    double s_temp = 30000, e_temp = 10000;

    while (1) {
      if (time < ela_times(start_time)) break;

      vvll _P(L, vll(L, -1));

      while (1) {
        _P = P;
        if (time < ela_times(start_time)) break;
        if (modify(_P, 8, stat_loop)) {
          stat_loop++;
          break;
        }
      }

      ll score = calc_score(_P);
      // 温度関数
      double temp = s_temp + (e_temp - s_temp) * ela_times(start_time) / time;
      // 遷移確率関数(最小化の場合)
      double prob = exp((best_score - score) / temp);
      if (prob > (rand() % INF) / (double)INF) { // 確率probで遷移する
        P = _P;
        best_score = score;
        // cerr << best_score << endl;
#ifdef _DEBUG
        map<ll, ll> d_ma;
        rep(i, 0, N) {
          ll bit = 0;
          rep(j, 0, m_N) {
            Pos to =
                Pos((Y[i] + moves[j].y + L) % L, (X[i] + moves[j].x + L) % L);
            if (P[to.y][to.x] == Max_P) bit += pow(2, j);
          }
          // if (d_ma.count(bit)) cerr << "err 重複があります" << endl;
          d_ma[bit]++;
        }
        // for(auto cu : d_ma)cerr << cu.second;
        // cerr << endl;
#endif
      }
    }

    ret_P = P;
  }
};

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
  double sa_init_time = 0.5;
  double sa_place_time = 3.5;
  ll Success_num = 1;
  ll Bio_N;
  ll Bio_K;
  ll Th;

  ll m_N;
  vector<Pos> moves;

  Pos pos_0;

  Solver_200_P() {
    pos_0 = Pos(0, 0);
    m_N = 4;
    Bio_N = 10000 / (m_N * N);
    Bio_N += 20;
    double p;
    optimize(Bio_K, Bio_N, m_N, Min_P, Max_P, Th, p);

    // cerr << Min_P << " " << Max_P << endl;
    // cerr << Bio_N << " " << setprecision(5) << " " << p << endl;
  }

  void decide_0_pos(vector<Pos> &_moves, vvll &P) {
    ll best = -1;
    rep(i, 0, L) {
      rep(j, 0, L) {
        if (P[i][j] == Min_P || P[i][j] == Max_P) continue;
        ll sum = 0;
        rep(dir, 0, 4) {
          rep(num, 1, 5) {
            Pos to((i + dy[dir] * num + L) % L, (j + dx[dir] * num + L) % L);
            sum += pow(P[to.y][to.x], 2);
          }
        }
        if (chmax(best, sum)) pos_0 = Pos(i, j);
      }
    }
    P[pos_0.y][pos_0.x] = 0;
  }

  void create_P(vvll &P) {
    auto start_time = system_clock::now();

    SA_init sa_init(moves, P, m_N, Min_P, Max_P, sa_init_time);
    sa_init.solve(moves, P);

    sa_place_time = sa_place_time - ela_times(start_time);
    SA_place sa_p(moves, P, m_N, Min_P, Max_P, sa_place_time);
    sa_p.solve(P);

    decide_0_pos(moves, P);

    rep(i, 0, L) {
      rep(j, 0, L) cout << P[i][j] << " ";
      cout << endl;
    }
    cout.flush();
  }

  void measure_E_P(vvll &P, vvll &E_P, ll &rest) {
    ll m_sum = 0;
    vvll sum_E_P(N, vll(m_N, 0));
    vvvll all_m(N, vvll(m_N));
    rep(i, 0, N) {
      rep(j, 0, m_N) {
        ll sum = 0;
        rep(num, 0, Bio_N) {
          cout << i << " " << moves[j].y << " " << moves[j].x << endl;
          cout.flush();
          m_sum++;
          measure_cost += 100 * (10 + moves[j].m_dist);

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

          sum += m;
          if (m == 1000) sum += 20;
          if (m == 0) sum += -20;
          all_m[i][j].push_back(m);

          double _S, _p;
          _S = S / sqrt(num + 1);
          // P[][] == Max_Pの時, 平均Min_P, 分散_Sの時上側確率が0.001の時
          _p = 0.5 * (1.0 + std::erf((Min_P - (sum / (num + 1))) / _S /
                                     std::sqrt(2.0)));
          if (_p < 0.001) {
            E_P[i][j] = 1;
            // if (P[r][c] == Min_P) cerr << "a";
            // else cerr << "b";
            break;
          }
          // P[][] == Min_Pの時, 平均Max_P, 分散_Sの時下側確率が0.001の時
          _p = 0.5 * (1.0 + std::erf(((sum / (num + 1)) - Max_P) / _S /
                                     std::sqrt(2.0)));
          if (_p < 0.001) {
            E_P[i][j] = 0;
            // if (P[r][c] == Max_P) cerr << "a";
            // else cerr << "b";
            break;
          }
        }

        // if (Bio_K <= sum) E_P[i][j] = 1;
        ll at = P[(Y[i] + moves[j].y + L) % L][(X[i] + moves[j].x + L) % L];
        if (E_P[i][j] == -1) {
          // if (S / 8 < abs(500 - sum / Bio_K)) continue;
          if (Th <= sum / Bio_N) E_P[i][j] = 1;
          else E_P[i][j] = 0;

          // if (E_P[i][j] == 0 && at == Max_P)
          //   cerr << abs(500 - sum / Bio_N) << " ";
          // if (E_P[i][j] == 1 && at == Min_P)
          //   cerr << abs(500 - sum / Bio_N) << " ";
        }
        sum_E_P[i][j] += sum;
      }
    }

#ifdef _DEBUG
    ll ans = 0;
    rep(i, 0, N) {
      rep(j, 0, m_N) {
        ll at =
            P[(Y[A[i]] + moves[j].y + L) % L][(X[A[i]] + moves[j].x + L) % L];
        if (at == Min_P && E_P[i][j] == 0) ans++;
        if (at == Max_P && E_P[i][j] == 1) ans++;
      }
    }
    // cerr << "correct   : " << ans << endl;
    // cerr << "incorrect : " << N * m_N - ans << endl;
    // cerr << "_Bio_N: " << (10000 - m_sum) / (N * 3) << endl;
#endif
    rest -= m_sum;
  }

  void new_moves(vvll &P, vector<Pos> &ret_moves, vll &wholes, vll &exits,
                 ll &_m_N) {
    rep(i, 0, (ll)exits.size()) {
      Pos cu(Y[exits[i]], X[exits[i]]);
      Pos _move = move(cu, pos_0);
      ret_moves.push_back(_move);
    }
    _m_N = (ll)ret_moves.size();
  }

  void measure_P(vvll &P, vll &E, vvll &E_P, ll &rest) {
    ll _Bio_N = rest / (N * 8);
    rep(i, 0, N) {
      if (E[i] != -1) continue;
      vll wholes, exits;
      wholes.push_back(i);
      ll dir_sum = 0;
      rep(k, 0, m_N) dir_sum += E_P[i][k] * pow(2, k);
      rep(j, i + 1, N) {
        if (E[j] != -1) continue;
        ll _dir_sum = 0;
        rep(k, 0, m_N) _dir_sum += E_P[j][k] * pow(2, k);
        if (dir_sum == _dir_sum) {
          wholes.push_back(j);
          E[j] = 1;
        }
      }
      rep(j, 0, N) {
        ll _dir_sum = 0;
        rep(k, 0, m_N) {
          ll at = P[(Y[j] + moves[k].y + L) % L][(X[j] + moves[k].x + L) % L];
          if (at == Max_P) _dir_sum += pow(2, k);
        }
        if (dir_sum == _dir_sum) exits.push_back(j);
      }

#ifdef _DEBUG
      for (ll cu1 : wholes) {
        bool flag = false;
        for (ll cu2 : exits) {
          if (A[cu1] == cu2) flag = true;
        }
        if (!flag) stat_wa2++;
      }
#endif

      vector<Pos> _moves;
      ll _m_N;
      new_moves(P, _moves, wholes, exits, _m_N);
      vvll result_m((ll)wholes.size(), vll(_m_N, 0));
      rep(j, 0, (ll)wholes.size()) {
        rep(k, 0, _m_N) {
          ll sum = 0;
          rep(num, 0, _Bio_N) {
            cout << wholes[j] << " " << _moves[k].y << " " << _moves[k].x
                 << endl;
            cout.flush();
            rest--;
            measure_cost += 100 * (10 + _moves[k].m_dist);

            ll m;
#ifdef _DEBUG
            random_device seed_gen;
            default_random_engine engine(seed_gen());
            normal_distribution<> dist(0.0, S);
            ll r = (Y[A[wholes[j]]] + _moves[k].y + L) % L;
            ll c = (X[A[wholes[j]]] + _moves[k].x + L) % L;
            m = max((ll)0, min((ll)1000, P[r][c] + (ll)dist(engine)));
#else
            cin >> m;
#endif

            sum += m;
            if (m == 1000) sum += 20;
            if (m == 0) sum += -20;
          }
          result_m[j][k] = sum / _Bio_N;
        }
      }

      mcf_graph<int, ll> g(2 * N + 2);
      int s = 2 * N, t = 2 * N + 1;

      rep(j, 0, (ll)wholes.size()) g.add_edge(s, wholes[j], 1, 0);
      rep(j, 0, (ll)exits.size()) g.add_edge(N + exits[j], t, 1, 0);

      rep(j, 0, (ll)wholes.size()) {
        rep(k, 0, (ll)exits.size()) {
          ll sum = 0;
          rep(l, 0, _m_N) {
            Pos to((Y[exits[k]] + _moves[l].y + L) % L,
                   (X[exits[k]] + _moves[l].x + L) % L);
            ll val1 = result_m[j][l];
            ll val2 = P[to.y][to.x];
            sum += pow(val2 - val1, 2);
          }
          g.add_edge(wholes[j], N + exits[k], 1, sum);
        }
      }

      auto result = g.flow(s, t, N);
      auto edges = g.edges();
      for (auto e : edges) {
        if (e.from == s || e.to == t || e.flow == 0) continue;
        E[e.from] = e.to - N;
      }

#ifdef _DEBUG
      rep(j, 0, (ll)wholes.size()) {
        ll best = INFF;
        vll all_sum;
        rep(k, 0, (ll)exits.size()) {
          ll sum = 0;
          rep(l, 0, _m_N) {
            Pos to((Y[exits[k]] + _moves[l].y + L) % L,
                   (X[exits[k]] + _moves[l].x + L) % L);
            ll val1 = result_m[j][l];
            ll val2 = P[to.y][to.x];
            sum += pow(val2 - val1, 2);
          }
          all_sum.push_back(sum);
        }
        Sort(all_sum);
        // for (auto cu : all_sum)
        //   cerr << cu << " ";
        // cerr << A[wholes[j]] << " " << E[wholes[j]] << " ";
        // if (A[wholes[j]] != E[wholes[j]]) cerr << "*";
        // cerr << endl;
      }
#endif
    }
  }

  void solve() {
    auto start_time = system_clock::now();
    vvll P(L, vll(L, 0));
    create_P(P);
    placement_cost = randing_cost(P);

    vll E(N, -1);
    vvll E_P(N, vll(m_N, -1));
    ll rest = 10000;
    measure_E_P(P, E_P, rest);

    // cerr << "rest : " << rest << endl;
    measure_P(P, E, E_P, rest);
    // cerr << "rest : " << rest << endl;

    cout << "-1 -1 -1" << endl;
    rep(i, 0, N) cout << max((ll)0, E[i]) << endl;
    cout.flush();
#ifdef _DEBUG
    rep(i, 0, N) stat_wa += (E[i] != A[i]);
#endif
  }
};

int main() {
  srand((unsigned int)time(NULL));
  auto start_time = system_clock::now();

  cin >> L >> N >> S;
  Y.resize(N);
  X.resize(N);
  rep(i, 0, N) cin >> Y[i] >> X[i];
#ifdef _DEBUG
  A.resize(N);
  rep(i, 0, N) cin >> A[i];
#endif

  if (S <= 256) {
    Solver solver;
    solver.solve();
  } else {
    Solver_200_P solver;
    solver.solve();
  }

  ll ans = 1e14 / (measure_cost + placement_cost + 1e5);
  rep(i, 0, stat_wa) ans *= 0.8;
  stat_time = ela_times(start_time);

  cerr << "L N S     : " << L << " " << N << " " << S << endl;
  // cerr << "stat_loop : " << stat_loop << endl;
  // cerr << "stat_time : " << stat_time << endl;
  // cerr << "stat_wa   : " << stat_wa << endl;
  // cerr << "stat_wa2  : " << stat_wa2 << endl;
  cerr << "score     : " << ans << endl;
  // cerr << "m_cost    : " << measure_cost << endl;
  // cerr << "p_cost    : " << placement_cost << endl;
}

// local10 -> local11

// 概要
/*
m_N方向全てが被っているのを8個まで許容しました.
その後, 被っているwholeに対し_m_N <= 8個を調べフローをした.
_m_Nのうち一つは0にする.
*/
