#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef pair<ll, ll> P;
typedef pair<P, ll> Pll;
typedef pair<ll, P> llP;
typedef tuple<ll, ll, ll> T;
typedef tuple<ll, ll, ll, ll> F;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vvi> vvvi;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
typedef vector<vvvll> vvvvll;
typedef vector<vvvvll> vvvvvll;
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
#define Copy(from, to) copy(from.begin(), from.end(), to.begin())
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
#define Unique(a) sort(a.begin(), a.end())
#define dPQ priority_queue<ll, vll, greater<ll>>
#define PQ priority_queue<ll, vll>
#define INF INT_MAX
#define INFF (9223372036854775800)
#define TIME_LIMIT (1.99)
#define def (20101)
// #define MOD (1000000007)
// #define MAX (2147483647)
#define MAX (1073741823)
#define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision(10) << << endl;
//偏角ソートはlong ddouble!
// auto ite = s.lower_bound("B");
// cout << __builtin_popcount(N) << endl;

void out_V(string H, ll N) {
  vll V_H(N, 0);
  ll idx = 0;
  rep(i, 0, N) {
    rep(j, i + 1, N) {
      if (H[idx] == '1') {
        V_H[i]++;
        V_H[j]++;
      }
      idx++;
    }
  }
  gSort(V_H);
  for (ll cu : V_H)
    cout << cu << " ";
}

string test(string &H, double ep) {
  ll N = (ll)H.size();
  string ret = "";
  rep(i, 0, N) {
    if (rand() < (ll)(MAX * 2 * ep))
      ret += (char)('1' - H[i] + '0');
    else
      ret += H[i];
  }
  return ret;
}

//[0 : size - 1] : 完全グラフ
//[size : N) : 辺なし
string construct_complete(ll N, ll size) {
  string ret = "";
  rep(j, 0, N) {
    rep(k, j + 1, N) {
      if (j < size && k < size)
        ret += '1';
      else
        ret += '0';
    }
  }
  return ret;
}

//[0 : size - 1] : すべてスターグラフ
//[size : N) : スターグラフにつながる
string construct_star(ll N, ll size) {
  string ret = "";
  rep(j, 0, N) {
    rep(k, j + 1, N) {
      if (j < size || k < size)
        ret += '1';
      else
        ret += '0';
    }
  }
  return ret;
}

//すべての次数がsizeのグラフ
string construct_same(ll N, ll size) {
  string ret = "";
  ll step = size % 2 == 0 ? size / 2 : size / 2 + 1;
  vvll g(N, vll(N, 0));
  vll V_H(N, 0);
  rep(i, 0, N) {
    ll ad = 0;
    if (size % 2 == 1 && i % 2 == 0)
      ad = 1;
    rep(j, i + 1 + ad, i + 1 + step) {
      ll mi = min(i, j % N);
      ll ma = max(i, j % N);
      g[mi][ma] = 1;
      V_H[mi]++;
      V_H[ma]++;
    }
  }
  rep(j, 0, N) {
    rep(k, j + 1, N) {
      if (g[j][k] == 1)
        ret += '1';
      else
        ret += '0';
    }
  }
  // cout << size << endl;
  // out_V(ret, N);
  // cout << endl;
  return ret;
}

void construct_G(vs &G, ll N, ll M, double ep) {
  ll same_l = M - min((ll)M / 2, max((ll)5, N / 4));
  if (0.38 <= ep || ep <= 0.03)
    same_l = M - min(2 * M / 3, max((ll)5, N - 1));
  ll complete_l = 0;
  ll star_l = same_l / 2;
  // 完全グラフ
  rep(i, 0, star_l) G[i] = construct_complete(N, N - 1 - i);

  // スターグラフ
  rep(i, star_l, same_l) G[i] = construct_star(N, i - star_l + 1);

  //同じ次数のグラフ
  rep(i, same_l, M) {
    ll step = max(ll(1), N / (M - same_l - 1));
    if (i == M - 1)
      G[i] = construct_same(N, (N - 1));
    else
      G[i] = construct_same(N, step * (i - same_l));
  }

  // rep(i, 0, M) {
  //   cout << "G[i] : " << i << endl;
  //   out_V(G[i], N);
  //   cout << endl;
  // }
}

ll select_N(ll M, double ep) {
  ll best_N = 0;
  long double best_score = -1;
  vs G(M, "");
  rep(N, 1, 101) {
    G.assign(M, "");
    construct_G(G, N, M, ep);
    ll correct = 0;
    vll V_H(N, 0);
    rep(s, 0, 500) {
      string H = test(G[s % M], ep);
      ll ma = -1;
      ll idx = 0;
      V_H.assign(N, 0);
      rep(i, 0, N) {
        rep(j, i + 1, N) {
          if (H[idx] == '1') {
            V_H[i]++;
            V_H[j]++;
          }
          idx++;
        }
      }
      gSort(V_H);
      rep(i, 0, N - 1) ma = max(ma, V_H[i] - V_H[i + 1]);
      ll sum = 0;
      int flag = 0;
      rep(i, 0, N) {
        if (flag)
          sum++;
        if (i != N - 1) {
          if (V_H[i] - V_H[i + 1] == ma)
            flag = 1;
        }
      }
      ll t = rand() % M;
      ll same_l = M - min((ll)M / 2, max((ll)5, N / 4));
      if (0.38 <= ep || ep <= 0.03)
        same_l = M - min(2 * M / 3, max((ll)5, N - 1));
      ll complete_l = 0;
      ll star_l = same_l / 2;
      if (s == 0) {
        // cout << "N, M : " << N << " " << M << endl;
        // cout << "complete_l : " << complete_l << endl;
        // cout << "star_l : " << star_l << endl;
        // cout << "same_l : " << same_l << endl;
      }
      ll ad = 1 + ep / 0.1;
      // cout << V_H[0] << " " << V_H.back() << endl;
      if (N / 3 < V_H[0] - V_H.back()) {
        // sumの共通部分
        ll sum_l = N - (same_l - star_l);
        ll sum_r = star_l;
        if (sum_l <= sum && sum <= sum_r) {
          ll E = (1 - ep) * (N - 1) + ep;
          // cout << "E , V_H[0] : " << E << " " << V_H[0] << endl;
          if (E - 2 < V_H[0])
            t = star_l + (N - 1 - sum);
          else
            t = sum - 1;
        } else {
          if (sum < sum_l)
            t = sum - 1;
          else
            t = star_l + (N - 1 - sum);
        }
      } else {
        ll mi = INFF;
        rep(i, same_l, M) {
          ll size = 0;
          rep(j, 0, N - 1) size += G[i][j] - '0';
          long double E = (1 - ep) * size + ep * (N - size - 1);
          long double _sum = 0.0;
          rep(j, 0, N) _sum += (E - V_H[j]) * (E - V_H[j]);
          if (_sum < mi) {
            mi = _sum;
            t = i;
          }
          // cout << i << " " << size << " " << _sum << endl;
        }
      }
      t = min(M - 1, max((ll)0, t));

      // cout << (s % M == t) << " ma : " << ma << " ";
      // cout << "s,t : " << s % M << " " << t << " ";
      // cout << "sum : " << sum << " ";
      // cout << "diff : " << V_H[0] - V_H.back() << " ";
      // cout << endl;
      if (s % M == t)
        correct++;
    }

    long double _score = 1.0;
    correct /= 5;
    rep(_i, 0, 100 - correct) _score *= 0.9;
    _score /= N;
    if (best_score < _score) {
      best_N = N;
      best_score = _score;
    }
    // cerr << "N : " << N << " correct : " << correct << endl;
  }
  if (best_N == 1)
    best_N = 100;
  return best_N;
}

void solve(ll M, double ep, vs &G, ll N, string H, ll &t) {
  ll ma = -1;
  vll V_H(N, 0);
  ll idx = 0;
  rep(i, 0, N) {
    rep(j, i + 1, N) {
      if (H[idx] == '1') {
        V_H[i]++;
        V_H[j]++;
      }
      idx++;
    }
  }
  gSort(V_H);
  rep(i, 0, N - 1) ma = max(ma, V_H[i] - V_H[i + 1]);
  ll sum = 0;
  int flag = 0;
  rep(i, 0, N) {
    if (flag)
      sum++;
    if (i != N - 1) {
      if (V_H[i] - V_H[i + 1] == ma)
        flag = 1;
    }
  }

  t = rand() % M;
  ll same_l = M - min((ll)M / 2, max((ll)5, N / 4));
  if (0.38 <= ep || ep <= 0.03)
    same_l = M - min(2 * M / 3, max((ll)5, N - 1));
  ll complete_l = 0;
  ll star_l = same_l / 2;
  ll ad = 1 + ep / 0.1;
  if (N / 3 < V_H[0] - V_H.back()) {
    // sumの共通部分
    ll sum_l = N - (same_l - star_l);
    ll sum_r = star_l;
    if (sum_l <= sum && sum <= sum_r) {
      ll E = (1 - ep) * (N - 1) + ep;
      if (E - 2 < V_H[0])
        t = star_l + (N - 1 - sum);
      else
        t = sum - 1;
    } else {
      if (sum < sum_l)
        t = sum - 1;
      else
        t = star_l + (N - 1 - sum);
    }
  } else {
    ll mi = INFF;
    rep(i, same_l, M) {
      ll size = 0;
      rep(j, 0, N - 1) size += G[i][j] - '0';
      long double E = (1 - ep) * size + ep * (N - size - 1);
      long double _sum = 0.0;
      rep(j, 0, N) _sum += (E - V_H[j]) * (E - V_H[j]);
      if (_sum < mi) {
        mi = _sum;
        t = i;
      }
    }
  }
  t = min(M - 1, max((ll)0, t));
}

void solve_large_ep(ll M, double ep, vs &G, ll N, string H, ll &t) {
  ll ma = -1;
  vll V_H(N, 0);
  ll idx = 0;
  rep(i, 0, N) {
    rep(j, i + 1, N) {
      if (H[idx] == '1') {
        V_H[i]++;
        V_H[j]++;
      }
      idx++;
    }
  }
  gSort(V_H);
  rep(i, 0, N - 1) ma = max(ma, V_H[i] - V_H[i + 1]);
  ll sum = 0;
  int flag = 0;
  rep(i, 0, N) {
    if (flag)
      sum++;
    if (i != N - 1) {
      if (V_H[i] - V_H[i + 1] == ma)
        flag = 1;
    }
  }
  // 0: G[0:49], 1: G[50:99]
  // [0 : 49] 完全グラフ1つと無向グラフ
  // [50 : 99] [1 : 40]の頂点数を選んでウニグラフ
  t = rand() % M;
  ll E;
  if (N < M) {
    // sum共通部分
    ll l = M % 2 == 0 ? N - M / 2 : N - M / 2 - 1;
    ll r = M % 2 == 0 ? M / 2 : M / 2 + 1;
    if (l <= sum && sum <= r) {
      ll E = (N - 1) * (1 - ep) + ep;
      // ll ad = M % 2 == 1 ? -1 : 1;
      if (E - 4 < V_H[0])
        t = M - 1 - (sum - l);
      else
        t = sum - 1;
    } else if (M / 2 < sum)
      t = M / 2 + (N - sum - 1);
    else
      t = sum - 1;
  } else {
    if (M / 2 < sum)
      t = M / 2 + (N - sum - 1);
    else
      t = sum - 1;
  }
  t = min(M - 1, max((ll)0, t));
}

struct solve_small_ep {
  ll M, N;
  double ep;
  solve_small_ep(ll M, double ep) : M(M), ep(ep) {}

  void out_V(string H, ll N) {
    vll V_H(N, 0);
    ll idx = 0;
    rep(i, 0, N) {
      rep(j, i + 1, N) {
        if (H[idx] == '1') {
          V_H[i]++;
          V_H[j]++;
        }
        idx++;
      }
    }
    gSort(V_H);
    for (ll cu : V_H)
      cout << cu << " ";
  }

  string test(string &H, double ep) {
    ll N = (ll)H.size();
    string ret = "";
    rep(i, 0, N) {
      if (rand() < (ll)(MAX * 2 * ep))
        ret += (char)('1' - H[i] + '0');
      else
        ret += H[i];
    }
    return ret;
  }

  //頂点数N, 次数0の個数num, N-num個の頂点は次数sta
  string construct_same(ll N, ll num, ll sta) {
    string ret = "";
    ll step = (sta + 1) / 2;
    vvll g(N, vll(N, 0));
    vll V_H(N, 0);
    rep(i, 0, N - num) {
      ll ad = 0;
      if (sta % 2 == 1 && sta != N - num - 1 && i % 2 == 0)
        ad = 1;
      rep(j, i + 1 + ad, i + 1 + step) {
        ll mi = min(i, j % (N - num));
        ll ma = max(i, j % (N - num));
        g[mi][ma] = 1;
      }
    }
    rep(j, 0, N) {
      rep(k, j + 1, N) {
        if (g[j][k] == 1)
          ret += '1';
        else
          ret += '0';
      }
    }
    // out_V(ret, N);
    // cout << endl;
    return ret;
  }

  void construct_G(vs &G, ll N, ll M, double ep, ll step, ll sta, vll &G_idx,
                   bool &flag) {
    //同じ次数のグラフ
    ll idx = 0;
    rep(num, 0, N) {
      if (sta >= N - num) {
        // cout << "sta, N - num :" << sta << " " << N - num << endl;
        continue;
      }
      ll f_sta = num == 0 ? 0 : sta;
      for (ll _sta = f_sta; _sta < N - num; _sta += step) {
        // cout << num << " " << _sta << endl;
        if (idx == M)
          continue;
        if (N - num < _sta + step)
          _sta = N - num - 1;
        G[idx] = construct_same(N, num, _sta);
        if (G_idx[num] == -1)
          G_idx[num] = idx;
        idx++;
      }
    }
    if (idx != M) {
      flag = false;
      return;
    }
    // cout << "N : " << N << endl;
    // rep(i, 0, M) {
    //   cout << "G[i] : " << i << endl;
    //   out_V(G[i], N);
    //   cout << endl;
    // }
    // cout << "G_idx " << endl;
    // for (ll cu : G_idx)
    //   if (cu != -1)
    //     cout << cu << " ";
    // cout << endl;
  }

  ll select_N(ll M, double ep, ll &step, ll &sta) {
    ll best_N = 0;
    ll best_step = 0;
    ll best_sta = 0;
    long double best_score = -1;
    ll best_correct = 0;
    vs G(M, "");
    vll G_idx((ll)100, -1);

    ll f_step = 1 / (1 - 3 * ep);
    ll f_sta = 10 / (1 - 3 * ep);
    // cerr << "f : " << f_step << " " << f_sta << endl;
    rep(N, 1, 101) {
      rep(_step, f_step, f_step + 3) {
        // rep(_sta, f_sta, f_sta + 1) {
        for (ll _sta = f_sta - 3; _sta <= f_sta + 3; _sta += 3) {
          // cout << "step : " << _step << " start : " << _sta << endl;
          G.assign(M, "");
          G_idx.assign(100, -1);
          bool can_construct = true;
          construct_G(G, N, M, ep, _step, _sta, G_idx, can_construct);
          if (!can_construct)
            continue;
          ll correct = 0;
          vll V_H(N, 0);
          ll loop_num = 100;
          rep(s, 0, loop_num) {
            string H = test(G[s % M], ep);
            ll ma = -1;
            ll idx = 0;
            V_H.assign(N, 0);
            rep(i, 0, N) {
              rep(j, i + 1, N) {
                if (H[idx] == '1') {
                  V_H[i]++;
                  V_H[j]++;
                }
                idx++;
              }
            }
            gSort(V_H);
            rep(i, 0, N - 1) ma = max(ma, V_H[i] - V_H[i + 1]);
            ll sum = 0;
            int flag = 0;
            rep(i, 0, N) {
              if (flag)
                sum++;
              if (i != N - 1) {
                if (V_H[i] - V_H[i + 1] == ma)
                  flag = 1;
              }
            }
            if (ma < (ll)4)
              sum = 0;
            ll t = rand() % M;
            ll l = G_idx[sum];
            if (l == -1)
              continue;
            ll r = G_idx[sum + 1] == -1 ? M : G_idx[sum + 1];
            ll mi = INFF;
            // cout << "l,r : " << l << " " << r << endl;
            rep(i, l, r) {
              ll size = 0;
              rep(j, N - 1, N + N - 3) size += G[i][j] - '0';
              size += G[i][0] - '0';
              // rep(j,0,N-1) size += G[i][j] - '0';
              long double E = (1 - ep) * size + (N - size - 1) * ep;
              long double _sum = 0.0;
              rep(j, N / 5, N - sum - N / 5) _sum +=
                  (E - V_H[j]) * (E - V_H[j]);
              if (_sum < mi) {
                mi = _sum;
                t = i;
              }
              // cout << "i : " << i << " size :" << size << " ";
              // cout << endl;
            }
            t = min(M - 1, max((ll)0, t));
            // cout << "V_H " << endl;
            // out_V(H, N);
            // cout << endl;
            // cout << "G" << endl;
            // out_V(G[s % M], N);
            // cout << endl;
            // cout << (s % M == t) << " ma : " << ma << " ";
            // cout << "s,t : " << s % M << " " << t << " ";
            // cout << "sum : " << sum << " ";
            // cout << endl;
            if (s % M == t)
              correct++;
          }
          long double _score = 1.0;
          correct /= (loop_num / 100);
          rep(_i, 0, 100 - correct) _score *= 0.9;
          _score /= N;
          if (best_score < _score) {
            best_N = N;
            best_step = _step;
            best_sta = _sta;
            best_correct = correct;
            best_score = _score;
          }
          // cout << "correct : " << correct << endl;
          // cerr << "N : " << N << " correct : " << correct << endl;
        }
      }
      // cerr << "N : " << N << " correct : " << best_correct << endl;
    }
    if (best_N == 1)
      best_N = 100;
    sta = best_sta;
    step = best_step;
    // cerr << "b : " << step << " " << sta << endl;
    return best_N;
  }

  void _solve(ll M, double ep, vs &G, ll N, string H, ll &t, vll &G_idx) {
    vll V_H(N, 0);
    ll ma = -1;
    ll idx = 0;
    V_H.assign(N, 0);
    rep(i, 0, N) {
      rep(j, i + 1, N) {
        if (H[idx] == '1') {
          V_H[i]++;
          V_H[j]++;
        }
        idx++;
      }
    }
    gSort(V_H);
    rep(i, 0, N - 1) ma = max(ma, V_H[i] - V_H[i + 1]);
    ll sum = 0;
    int flag = 0;
    rep(i, 0, N) {
      if (flag)
        sum++;
      if (i != N - 1) {
        if (V_H[i] - V_H[i + 1] == ma)
          flag = 1;
      }
    }
    if (ma < (ll)4)
      sum = 0;
    t = rand() % M;
    ll l = G_idx[sum];
    if (l == -1)
      return;
    ll r = G_idx[sum + 1] == -1 ? M : G_idx[sum + 1];
    ll mi = INFF;
    // cout << "l,r : " << l << " " << r << endl;
    rep(i, l, r) {
      ll size = 0;
      rep(j, N - 1, N + N - 3) size += G[i][j] - '0';
      size += G[i][0] - '0';
      // rep(j,0,N-1) size += G[i][j] - '0';
      long double E = (1 - ep) * size + (N - size - 1) * ep;
      long double _sum = 0.0;
      rep(j, N / 5, N - sum - N / 5) _sum += (E - V_H[j]) * (E - V_H[j]);
      if (_sum < mi) {
        mi = _sum;
        t = i;
      }
    }
    t = min(M - 1, max((ll)0, t));
  }

  void solve() {
    vs G(M, "");
    // vvll G_idx;
    ll step = 0;
    ll sta = 0;
    ll N = select_N(M, ep, step, sta);
    vll G_idx((ll)100, -1);

    bool construct_flag = true;
    construct_G(G, N, M, ep, step, sta, G_idx, construct_flag);

    cout << N << endl;
    rep(i, 0, M) cout << G[i] << endl;
    cout.flush();

    ll correct = 0;
    rep(k, 0, 100) {
      string H;
      cin >> H;

      // ll s;
      // cin >> s;
      // string H = test(G[s], ep);

      ll t = 0;
      _solve(M, ep, G, N, H, t, G_idx);

      cout << t << endl;
      cout.flush();

      // cout << (s == t) << " ";
      // cout << "s,t : " << s << " " << t << " ";
      // cout << endl;
      // if (s % M == t)
      //   correct++;
    }

    // cerr << N << " " << correct << endl;
  }
};

struct solve_ep0 {
  ll M, N;
  double ep;

  solve_ep0(ll M, double ep) : M(M), ep(ep) {
    if (31 < M)
      N = 6;
    else if (11 < M)
      N = 5;
    else
      N = 4;
  }

  void out_V(string H, ll N) {
    vll V_H(N, 0);
    ll idx = 0;
    rep(i, 0, N) {
      rep(j, i + 1, N) {
        if (H[idx] == '1') {
          V_H[i]++;
          V_H[j]++;
        }
        idx++;
      }
    }
    gSort(V_H);
    for (ll cu : V_H)
      cout << cu << " ";
  }

  string test(string &H, double ep) {
    ll N = (ll)H.size();
    string ret = "";
    rep(i, 0, N) {
      if (rand() < (ll)(MAX * 2 * ep))
        ret += (char)('1' - H[i] + '0');
      else
        ret += H[i];
    }
    return ret;
  }

  void construct_G(ll N, vs &G) {
    ll V = N * (N - 1) / 2;
    map<vll, string> ma;
    rep(i, 0, (1ll << V)) {
      ll u = 0;
      ll v = 1;
      string s_G = "";
      vll V_H(N, 0);
      rep(j, 0, V) {
        if ((1ll << j) & i) {
          V_H[u]++;
          V_H[v]++;
          s_G += '1';
        } else
          s_G += '0';
        v++;
        if (v == N) {
          u++;
          v = u + 1;
        }
      }
      Sort(V_H);
      if (!ma.count(V_H))
        ma[V_H] = s_G;
    }
    ll sum = 0;
    for (auto cu : ma) {
      if (sum == M)
        break;
      G[sum] = cu.second;
      sum++;
    }
    // rep(i,0,M){
    //   cout << "G[i] :" << i << endl;
    //   out_V(G[i],N);
    //   cout << endl;
    // }
  }

  void _solve(ll M, vs &G, ll N, string H, ll &t) {
    vll V_H(N, 0);
    ll idx = 0;
    rep(i, 0, N) {
      rep(j, i + 1, N) {
        if (H[idx] == '1') {
          V_H[i]++;
          V_H[j]++;
        }
        idx++;
      }
    }
    Sort(V_H);
    rep(m, 0, M) {
      vll _V_H(N, 0);
      ll _idx = 0;
      rep(i, 0, N) {
        rep(j, i + 1, N) {
          if (G[m][_idx] == '1') {
            _V_H[i]++;
            _V_H[j]++;
          }
          _idx++;
        }
      }
      Sort(_V_H);
      int flag = 1;
      // cout << "G[s] " << endl;
      // out_V(G[m],N);
      // cout << endl;
      // cout << "H " << endl;
      // out_V(H,N);
      // cout << endl;
      rep(j, 0, N) if (V_H[j] != _V_H[j]) flag = 0;
      if (flag == 1)
        t = m;
    }
  }

  void solve() {
    vs G(M, "");
    construct_G(N, G);
    cout << N << endl;
    rep(i, 0, M) cout << G[i] << endl;
    cout.flush();

    ll correct = 0;
    rep(k, 0, 100) {
      string H;
      cin >> H;

      // ll s;
      // cin >> s;
      // string H = test(G[s], ep);

      ll t = 0;
      _solve(M, G, N, H, t);

      cout << t << endl;
      cout.flush();
      // cout << (s == t) << " ";
      // cout << "s,t : " << s << " " << t << " ";
      // cout << endl;
      // if (s % M == t)
      //   correct++;
    }

    // cerr << N << " " << correct << endl;
  }
};

int main() {
  srand((unsigned int)time(NULL));
  ll M;
  cin >> M;
  double ep;
  cin >> ep;

  if (ep == 0.0) {
    solve_ep0 ret(M, ep);
    ret.solve();
    return 0;
  } else if (ep <= 0.1) {
    solve_small_ep ret(M, ep);
    ret.solve();
    return 0;
  }
  vs G(M, "");

  ll N = select_N(M, ep);
  if (ep >= 0.30)
    N = 100;
  construct_G(G, N, M, ep);

  cout << N << endl;
  rep(i, 0, M) cout << G[i] << endl;
  cout.flush();

  ll correct = 0;
  rep(k, 0, 100) {
    string H;
    cin >> H;

    // ll s;
    // cin >> s;
    // string H = test(G[s], ep);

    ll t = 0;
    solve(M, ep, G, N, H, t);

    cout << t << endl;
    cout.flush();

    // cout << (s == t) << " ";
    // cout << "s,t : " << s << " " << t << " ";
    // cout << endl;
    // if (s % M == t)
    //   correct++;
  }

  cerr << N << " " << correct << endl;
}
