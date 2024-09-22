#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef long double ld;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
typedef pair<ll, ll> pll;
typedef vector<pll> vp;
typedef vector<vp> vvp;
typedef vector<vvp> vvvp;
typedef vector<ld> vld;
typedef vector<vld> vvld;
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
#define TIME_LIMIT (1.9)
#define def (20101)
// #define MOD (1000000007)
// #define MAX (2147483647)
#define MAX (1073741823)
#define MOD (998244353)

ld ela_times(system_clock::time_point &clock) {
  return duration_cast<microseconds>(system_clock::now() - clock).count() *
         1e-6;
}

ld stat_time = 0.0;

ll n, m;
vvll b;
ll K = 5000;

struct SA {
  double time;
  SA(ld time) : time(time) {}

  vp init(vvll b) {
    vp state;
    ll min_v = 1;
    rep(i, 0, K) {
      ll pre_box = -1;
      ll operate = 1;
      ll idx = -1;
      rep(j, 0, m) {
        if (b[j].size() == 0) continue;
        rep(k, 0, (ll)b[j].size()) {
          if (b[j][k] == min_v) {
            pre_box = j;
            idx = k;
            if (k == (ll)b[j].size() - 1) {
              operate = 2;
            }
          }
        }
      }

      ll aft_box = -1;
      ll best_score = -1;
      rep(j, 0, m) {
        if (pre_box == j) continue;
        // 山jの中で最小値
        ll score = INF;
        rep(k, 0, (ll)b[j].size()) {
          if (b[j][k] < score) score = b[j][k];
        }
        if (best_score < score) {
          best_score = score;
          aft_box = j;
        }
      }

      rep(j, idx + 2, (ll)b[pre_box].size() - 2) {
        if (b[pre_box][j] < best_score) {
          idx = j;
          best_score = b[pre_box][j];
          // cout << min_v << " x  " << idx << endl;
        }
      }

      if (operate == 1) {
        rep(j, idx + 1, (ll)b[pre_box].size()) {
          b[aft_box].push_back(b[pre_box][j]);
          if (j == idx + 1) {
            state.push_back({b[pre_box][j], aft_box + 1});
          }
        }
        ll num = (ll)b[pre_box].size() - idx - 1;
        rep(j, 0, num) b[pre_box].pop_back();
      }

      if (operate == 2) {
        state.push_back({b[pre_box][idx], 0});
        b[pre_box].pop_back();
        min_v++;
        if (min_v == 201) break;
      }
    }

    return state;
  }

  void modify(vp &state, vvll b) {
    ll idx = rand() % (state.size() - 5);
    rep(i, 0, idx) {
      auto cu = state[i];
      if (cu.second == 0) {
        rep(i, 0, m) {
          if (b[i][b[i].size() - 1] == cu.first) b[i].pop_back();
        }
      } else {
        ll pre_box = 0, aft_box = cu.second - 1;
        int f = 0;
        rep(i, 0, m) {
          ll num = 0;
          rep(j, 0, b[i].size()) {
            if (b[i][j] == cu.first) {
              pre_box = i;
              f = 1;
            }
            if (f == 1) {
              num++;
              b[aft_box].push_back(b[i][j]);
            }
          }
          rep(j, 0, num) b[pre_box].pop_back();

          if (f == 1) {
            break;
          }
        }
      }
    }
    ll num = state.size() - idx - 1;
    rep(i, 0, num) state.pop_back();
    cerr << "b";
    // ここから
    ll min_v = INF;
    int flag = 1;
    rep(i, 0, m) for (ll cu : b[i]) min_v = min(min_v, cu);
    cerr << min_v << " ";
    rep(i, 0, K) {
      if (195 < min_v) cerr << min_v << " ";
      ll pre_box = -1;
      ll operate = 1;
      ll idx = -1;
      rep(j, 0, m) {
        if (b[j].size() == 0) continue;
        rep(k, 0, (ll)b[j].size()) {
          if (b[j][k] == min_v) {
            pre_box = j;
            idx = k;
            if (k == (ll)b[j].size() - 1) {
              operate = 2;
            }
          }
        }
      }

      ll aft_box = -1;
      ll best_score = -1;
      rep(j, 0, m) {
        if (pre_box == j) continue;
        // 山jの中で最小値
        ll score = INF;
        rep(k, 0, (ll)b[j].size()) {
          if (b[j][k] < score) score = b[j][k];
        }
        if (best_score < score) {
          best_score = score;
          aft_box = j;
        }
      }

      rep(j, idx + 2, (ll)b[pre_box].size() - 2) {
        if (b[pre_box][j] < best_score) {
          idx = j;
          best_score = b[pre_box][j];
          // cout << min_v << " x  " << idx << endl;
        }
      }

      if (flag && operate == 1) {
        aft_box = (pre_box + rand() % (m - 1) + 1) % m;
        flag = 0;
      }

      if (196 <= min_v) {
        cerr << "u : " << operate << endl;
        cerr << pre_box << " " << aft_box << endl;
        rep(i, 0, m) {
          for (ll cu : b[i])
            cerr << cu << " ";
        }
        cerr << endl;
      }
      if (pre_box == -1) break;

      if (operate == 1) {
        rep(j, idx + 1, (ll)b[pre_box].size()) {
          b[aft_box].push_back(b[pre_box][j]);
          if (j == idx + 1) {
            state.push_back({b[pre_box][j], aft_box + 1});
          }
        }
        ll num = (ll)b[pre_box].size() - idx - 1;
        rep(j, 0, num) b[pre_box].pop_back();
      }

      if (operate == 2) {
        state.push_back({b[pre_box][idx], 0});
        b[pre_box].pop_back();
        min_v++;
        if (min_v == 201) break;
      }
    }

    cerr << "c";
  }

  ll calc_score(vp &state, vvll b) {
    ll ret = 0;
    for (auto cu : state) {
      if (cu.second == 0) {
        rep(i, 0, m) {
          if (b[i][b[i].size() - 1] == cu.first) b[i].pop_back();
        }
      } else {
        ll pre_box = 0, aft_box = cu.second - 1;
        ll idx = -1;
        rep(i, 0, m) {
          rep(j, 0, b[i].size()) {
            if (b[i][j] == cu.first) {
              pre_box = i;
              idx = j;
            }
          }
        }

        ll num = 0;
        rep(i, idx, b[pre_box].size()) {
          b[aft_box].push_back(b[pre_box][i]);
          num++;
        }

        rep(i, 0, num) b[pre_box].pop_back();
        ret += num + 1;
      }
    }
    return max((ll)1, 10000 - ret);
  }

  void solve(vp &ans, vvll b) {
    auto start_time = system_clock::now();

    vp state = init(b);

    ll best_score = calc_score(state, b);
    cerr << best_score << endl;

    while (1) {
      if (time < ela_times(start_time)) break;
      cerr << "u";
      vp _state = state;
      cerr << "m";
      modify(_state, b);
      cerr << endl;
      cerr << "y";
      state.clear();
      state = _state;
      cerr << "a";

      // ll score = calc_score(_state, b);
      // if (score != best_score) cerr << score << endl;

      // if (best_score < score) { // 確率probで遷移する
      //   state = _state;
      //   best_score = score;
      // }
    }
    cerr << "z";
    ans = state;
  }
};

int main() {
  cin >> n >> m;
  b.resize(m);
  rep(i, 0, m) b[i];
  rep(i, 0, m) {
    rep(j, 0, n / m) {
      ll _b;
      cin >> _b;
      b[i].push_back(_b);
    }
  }

  vp ans;
  SA sa(1.9);
  sa.solve(ans, b);

  for (auto cu : ans)
    cout << cu.first << " " << cu.second << endl;
}