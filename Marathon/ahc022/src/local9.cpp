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
double stat_loop2 = 0;
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
    // cerr << prob << " " << z << endl;
    if (0.99 <= prob) {
      Min_P = _Min_P;
      Max_P = _Max_P;
      _p = prob;
      Th = 500;
      Bio_K = 1;
      cerr << prob << " " << z << " " << _S << endl;

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

      if (ok[0] == 0 && ok[1] == 0) return false;

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
    cerr << L * L - ret - _sum << " " << _sum << " " << ret << endl;

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
        stat_loop2++;
        if (chmin(best_score, score)) {
          moves = _moves;
          P = _P;
        }
      }
    }
    ll _score = calc_score(moves);

    ret_moves = moves;
    ret_P = P;
  }
};

struct SA_eliminate {
  vector<Pos> moves;
  vvll P;
  ll m_N, Min_P, Max_P;
  double time;

  SA_eliminate(vector<Pos> moves, vvll P, ll m_N, ll Min_P, ll Max_P,
               double time)
      : moves(moves), P(P), m_N(m_N), Min_P(Min_P), Max_P(Max_P), time(time) {
    vvll fix(L, vll(L, 0));
    rep(i, 0, N) {
      rep(j, 0, m_N) {
        Pos _to((Y[i] + moves[j].y + L) % L, (X[i] + moves[j].x + L) % L);
        fix[_to.y][_to.x] = 1;
      }
    }
  }

  bool modify(vector<Pos> &_moves, vvll &_P) {
    map<pll, ll> ma;
    for (auto cu : _moves)
      ma[{cu.y, cu.x}] = 1;
    ll idx = rand() % m_N;
    ll _dy = dy[rand() % 4] * (rand() % 5);
    ll _dx = _dy == 0 ? dx[rand() % 2 * 2] * (rand() % 5)
                      : dx[rand() % 2 * 2] * (rand() % 5);
    Pos _from(0, 0);
    Pos _to((_moves[idx].y + _dy + L) % L, (_moves[idx].x + _dx + L) % L);
    Pos new_move = move(_from, _to);
    if (ma.count({new_move.y, new_move.x})) return false;
    // cerr << _moves[idx].y << " " << _moves[idx].x << ", " << new_move.y;
    // cerr << " " << new_move.x << endl;

    rep(i, 0, N) {
      Pos p_to((Y[i] + _moves[idx].y + L) % L, (X[i] + _moves[idx].x + L) % L);
      Pos a_to((Y[i] + _moves[idx].y + _dy + L) % L,
               (X[i] + _moves[idx].x + _dx + L) % L);
      // cerr << " " << p_to.y << " " << p_to.x << ", " << a_to.y << " " <<
      // a_to.x; cerr << " " << _P[p_to.y][p_to.x] << " " << _P[a_to.y][a_to.x]
      // << endl;
      if (_P[a_to.y][a_to.x] == Min_P && _P[p_to.y][p_to.x] == Max_P) {
        // cerr << i << endl;
        return false;
      }
      if (_P[a_to.y][a_to.x] == Max_P && _P[p_to.y][p_to.x] == Min_P) {
        // cerr << i << endl;
        return false;
      }
      _P[a_to.y][a_to.x] = _P[p_to.y][p_to.x];
    }
    _moves[idx] = new_move;

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
    cerr << L * L - ret - _sum << " " << _sum << " " << ret << endl;

    return ret;
  }

  void solve(vector<Pos> &ret_moves, vvll &ret_P) {
    auto start_time = system_clock::now();

    ll best_score = calc_score(moves);
    while (1) {
      if (time < ela_times(start_time)) break;
      stat_loop2++;

      vector<Pos> _moves;
      vvll _P(L, vll(L, -1));
      _moves = moves;
      _P = P;
      if (!modify(_moves, _P)) continue;

      ll score = calc_score(_moves);
      cerr << score << " ";
      if (chmin(best_score, score)) {
        moves = _moves;
        P = _P;
      }
    }
    ret_moves = moves;
    ret_P = P;
  }
};

struct SA_place {
  vector<Pos> moves;
  vvll dir_bit;
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

    dir_bit.resize(L);
    rep(i, 0, L) dir_bit[i].assign(L, 0);
    rep(i, 0, N) {
      rep(j, 0, m_N) {
        Pos _to((Y[i] + moves[j].y + L) % L, (X[i] + moves[j].x + L) % L);
        if (P[_to.y][_to.x] == Max_P) dir_bit[Y[i]][X[i]] += pow(2, j);
      }
    }
  }

  bool modify(vvll &_P, ll over_ma) {
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
      if (over_ma < d_ma[sum] ) return false;
    }

    optimize_P(_P, fix);
    return true;
  }

  ll calc_score(vvll &P) {
    ll ret = randing_cost(P);
    return ret;
  }

  void solve(vvll &ret_P) {
    auto start_time = system_clock::now();

    ll best_score = calc_score(P);
    cerr << best_score << endl;
    double s_temp = 30000, e_temp = 10000;

    while (1) {
      if (time < ela_times(start_time)) break;

      vvll _P(L, vll(L, -1));

      while (1) {
        _P = P;
        if (time < ela_times(start_time)) break;
        if (modify(_P, 4)) {
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
        cerr << best_score << endl;
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
  const double sa_init_time = 0.5;
  const double sa_eliminate_time = 0.0;
  const double sa_place_time = 3.0;
  const double create_P_time = 3.5;
  ll Success_num = 1;
  ll Bio_N;
  ll Bio_K;
  ll Th;

  ll m_N;
  vector<Pos> moves;

  Solver_200_P() {
    Min_P = max((ll)0, (ll)500 - S / 2);
    Max_P = min((ll)1000, (ll)500 + S / 2);
    rep(i, 2, 12) {
      if (N <= pow(2, i - 1)) {
        m_N = i;
        break;
      }
    }
    Bio_N = 10000 / (m_N * N);
    double p;
    optimize(Bio_K, Bio_N, m_N, Min_P, Max_P, Th, p);

    cerr << Min_P << " " << Max_P << endl;
    cerr << Bio_N << " " << Bio_K << " " << setprecision(5) << " " << p << endl;
  }

  void create_P(vvll &P) {
    auto start_time = system_clock::now();

    SA_init sa_init(moves, P, m_N, Min_P, Max_P, sa_init_time);
    sa_init.solve(moves, P);

    // #ifdef _DEBUG
    //     map<ll, ll> d_ma;
    //     rep(i, 0, N) {
    //       ll bit = 0;
    //       rep(j, 0, m_N) {
    //         Pos to = Pos((Y[i] + moves[j].y + L) % L, (X[i] + moves[j].x + L)
    //         % L); if (P[to.y][to.x] == Max_P) bit += pow(2, j);
    //       }
    //       if (d_ma.count(bit)) cerr << "err 重複があります" << endl;
    //       d_ma[bit] = 1;
    //     }
    // #endif

    // SA_eliminate sa_e(moves, P, m_N, Min_P, Max_P, sa_eliminate_time);
    // sa_e.solve(moves, P);

    SA_place sa_p(moves, P, m_N, Min_P, Max_P, sa_place_time);
    sa_p.solve(P);

    // #ifdef _DEBUG
    //     map<ll, ll> d_ma;
    //     rep(i, 0, N) {
    //       ll bit = 0;
    //       rep(j, 0, m_N) {
    //         Pos to = Pos((Y[i] + moves[j].y + L) % L, (X[i] + moves[j].x + L)
    //         % L); if (P[to.y][to.x] == Max_P) bit += pow(2, j);
    //       }
    //       if (d_ma.count(bit)) cerr << "err 重複があります" << endl;
    //       d_ma[bit] = 1;
    //     }
    // #endif

    rep(i, 0, L) {
      rep(j, 0, L) cout << P[i][j] << " ";
      cout << endl;
    }
    cout.flush();
  }

  void measure_P(vvll &P, vll &E) {
    ll m_sum = 0;
    vvll E_P(N, vll(m_N, -1));
    vvll sum_E_P(N, vll(m_N, 0));
    vvvll all_m(N, vvll(m_N));
    rep(i, 0, N) {
      rep(j, 0, m_N) {
        ll sum = 0;
        rep(num, 0, Bio_N) {
          cout << i << " " << moves[j].y << " " << moves[j].x << endl;
          m_sum++;
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

          sum += m;
          if (m == 1000) sum += 20;
          if (m == 0) sum += -20;
          all_m[i][j].push_back(m);

          double _S, _p;
          _S = S / sqrt(num + 1);
          // P[][] == Max_Pの時, 平均Min_P, 分散_Sの時上側確率が0.001の時
          _p = 0.5 * (1.0 + std::erf((Min_P - (sum / (num + 1))) / _S /
                                     std::sqrt(2.0)));
          if (_p < 0.0005) {
            E_P[i][j] = 1;
            // if (P[r][c] == Min_P) cerr << "a";
            // else cerr << "b";
            break;
          }
          // P[][] == Min_Pの時, 平均Max_P, 分散_Sの時下側確率が0.001の時
          _p = 0.5 * (1.0 + std::erf(((sum / (num + 1)) - Max_P) / _S /
                                     std::sqrt(2.0)));
          if (_p < 0.0005) {
            E_P[i][j] = 0;
            // if (P[r][c] == Max_P) cerr << "a";
            // else cerr << "b";
            break;
          }
        }

        // if (Bio_K <= sum) E_P[i][j] = 1;
        ll at = P[(Y[i] + moves[j].y + L) % L][(X[i] + moves[j].x + L) % L];
        if (E_P[i][j] == -1) {
          if (S / 4 < abs(500 - sum / Bio_K)) continue;
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

    ll det_num = 0;
    rep(i, 0, N) rep(j, 0, m_N) if (E_P[i][j] != -1) det_num++;
    cerr << "det   : " << det_num << endl;
    cerr << "rest  : " << N * m_N - det_num << endl;
    if (N * m_N != det_num)
      cerr << "_Bio_N: " << (10000 - m_sum) / (N * m_N - det_num) << endl;

    ll ans1 = 0, ans2 = 0;
    rep(i, 0, N) {
      rep(j, 0, m_N) {
        if (E_P[i][j] != -1 || det_num == N * m_N) continue;
        ll sum = 0;
        ll rest_N = (10000 - m_sum) / (N * m_N - det_num);
        rep(num, 0, rest_N) {
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

          all_m[i][j].push_back(m);
          sum += m;
          if (m == 1000) sum += 20;
          if (m == 0) sum += -20;
        }
        sum_E_P[i][j] = sum;

        ll at = P[(Y[i] + moves[j].y + L) % L][(X[i] + moves[j].x + L) % L];
        ll si = (ll)all_m[i][j].size();
        if (500 <= all_m[i][j][si / 2] && at == Max_P) ans1++;
        if (500 > all_m[i][j][si / 2] && at == Min_P) ans1++;

        if (500 <= sum_E_P[i][j] / (Bio_N + rest_N) && at == Max_P) ans2++;
        if (500 > sum_E_P[i][j] / (Bio_N + rest_N) && at == Min_P) ans2++;

        if (500 <= sum_E_P[i][j] / (Bio_N + rest_N)) E_P[i][j] = 1;
        else E_P[i][j] = 0;
        // if (E_P[i][j] == 1 && at == Min_P)
        //   cerr << abs(sum_E_P[i][j] / (Bio_N + rest_N) - 500) << " ";
        // if (E_P[i][j] == 0 && at == Max_P)
        //   cerr << abs(sum_E_P[i][j] / (Bio_N + rest_N) - 500) << " ";
      }
    }
    cerr << ans1 << " " << ans2 << endl;
    rep(i, 0, N) {
      ll mi = INFF;
      rep(j, 0, N) {
        ll sum = 0;
        rep(k, 0, m_N) {
          Pos to((Y[j] + moves[k].y + L) % L, (X[j] + moves[k].x + L) % L);
          if (E_P[i][k] == 0 && P[to.y][to.x] == Max_P) sum++;
          if (E_P[i][k] == 1 && P[to.y][to.x] == Min_P) sum++;
        }
        if (chmin(mi, sum)) E[i] = j;
      }
    }

#ifdef _DEBUG
    vll _stat(4, 0);

    rep(i, 0, N) {
      vp all;
      rep(j, 0, N) {
        ll sum = 0;
        rep(k, 0, m_N) {
          Pos to((Y[j] + moves[k].y + L) % L, (X[j] + moves[k].x + L) % L);
          if (E_P[i][k] == 0 && P[to.y][to.x] == Max_P) sum++;
          if (E_P[i][k] == 1 && P[to.y][to.x] == Min_P) sum++;
        }
        all.push_back({sum, j});
      }
      _stat[min((ll)3, all[A[i]].first)]++;
      Sort(all);
      if (all[0].second != A[i]) stat_wa2++;
      rep(j, 0, N) if (A[i] == all[j].second) cerr << all[j].first << " ";
      rep(j, 0, N) { cerr << all[j].first; }
      cerr << endl;
    }
    rep(j, 0, 4) cerr << _stat[j] << " ";
    cerr << "(";
    rep(j, 0, 4) cerr << _stat[j] / (double)N << " ";
    cerr << ")";
    cerr << endl;
#endif
  }

  void solve() {
    vvll P(L, vll(L, 0));
    create_P(P);
    placement_cost = randing_cost(P);

    vll E(N, -1);
    // measure_P(P, E);

    cout << "-1 -1 -1" << endl;
    rep(i, 0, N) cout << max((ll)0, E[i]) << endl;
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

  // if (S < 200) {
  //   Solver solver;
  //   solver.solve();
  // } else {
  Solver_200_P solver;
  solver.solve();
  // }

  ll ans = 1e14 / (measure_cost + placement_cost + 1e5);
  rep(i, 0, stat_wa) ans *= 0.8;

  cerr << "L N S     : " << L << " " << N << " " << S << endl;
  cerr << "stat_loop : " << stat_loop << endl;
  cerr << "stat_loo2 : " << (ll)stat_loop2 << endl;
  cerr << "stat_time : " << stat_time << endl;
  cerr << "stat_wa   : " << stat_wa << endl;
  cerr << "stat_wa2  : " << stat_wa2 << " (local5 ver.)" << endl;
  cerr << "score     : " << ans << endl;
  cerr << "m_cost    : " << measure_cost << endl;
  cerr << "p_cost    : " << placement_cost << endl;
}

// local8 -> local9

// 概要
/*
create_P()で焼きなましました.
S >= 500時の確率を変更しました.
TH <= mからTH < mにして, Min_P = 0の時に対応しました.
*/
