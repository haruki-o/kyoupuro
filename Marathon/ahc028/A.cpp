#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef pair<ll, ll> pll;
typedef pair<pll, ll> ppll;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
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
#define INFF (9223372036854775800)
#define INF (10000000000000000)
#define TIME_LIMIT (1.9)
#define def (20101)
// #define MOD (1000000007)
// #define MAX (2147483647)
#define MAX (1073741823)
#define MOD (998244353)

double stat_time = 0.0;
double stat_loop = 0;
ll stat_ans = 0;

ll N, M;
ll si, sj;
vs A, t;

double ela_times(system_clock::time_point &clock) {
  return duration_cast<microseconds>(system_clock::now() - clock).count() * 1e-6;
}

ll dist(ll pre, ll to) {
  return abs(pre / N - to / N) + abs(pre % N - to % N) + 1;
}

// store_path[0] = {{13, 31, 51, 130, 51}}, {}, {}}
map<ll, vector<pair<ll, vll>>> store_path;

struct SA {
  void init(vp &que) {
    vvll all(26);
    rep(i, 0, M) all[t[i][0] - 'A'].push_back(i);

    ll cu = all[A[si][sj] - 'A'][0];
    que.push_back({cu, 0});
    all[A[si][sj] - 'A'].erase(all[A[si][sj] - 'A'].begin());

    while (1) {
      ll score = 0;
      ll best = -1;
      char cu_end_c = t[cu][4];
      if (all[cu_end_c - 'A'].size() == 0) {
        int f = 0;
        rep(i, 0, 26) {
          if (all[i].size() != 0) {
            que.push_back({all[i][0], 0});
            cu = all[i][0];
            all[i].erase(all[i].begin());
            f = 1;
            break;
          }
        }
        if (f == 0) break;
      } else {
        que.push_back({all[cu_end_c - 'A'][0], 0});
        cu = all[cu_end_c - 'A'][0];
        all[cu_end_c - 'A'].erase(all[cu_end_c - 'A'].begin());
      }
    }
  }

  ll modify(vp &que) {
    ll idx1 = rand() % M;
    ll idx2 = (idx1 + rand() % M) % M;

    ll idx1_pre = -1, idx1_aft = -1;
    ll idx2_pre = -1, idx2_aft = -1;
    ll at = 0;
    if (idx1 != 0) {
      for (auto cu1 : store_path[que[idx1 - 1].first]) {
        if (at != que[idx1 - 1].second) {
          at++;
          continue;
        }
        idx1_pre = cu1.second[4];
        break;
      }
    }
    at = 0;
    if (idx1 != M - 1) {
      for (auto cu1 : store_path[que[idx1 + 1].first]) {
        if (at != que[idx1 + 1].second) {
          at++;
          continue;
        }
        idx1_aft = cu1.second[0];
        break;
      }
    }
    at = 0;
    if (idx2 != 0) {
      for (auto cu1 : store_path[que[idx2 - 1].first]) {
        if (at != que[idx2 - 1].second) {
          at++;
          continue;
        }
        idx2_pre = cu1.second[4];
        break;
      }
    }
    at = 0;
    if (idx2 != M - 1) {
      for (auto cu1 : store_path[que[idx2 + 1].first]) {
        if (at != que[idx2 + 1].second) {
          at++;
          continue;
        }
        idx2_aft = cu1.second[0];
        break;
      }
    }

    ll idx1_fro, idx1_bac;
    ll idx2_fro, idx2_bac;
    at = 0;
    for (auto cu1 : store_path[que[idx1].first]) {
      if (at != que[idx1].second) {
        at++;
        continue;
      }
      idx1_fro = cu1.second[0];
      idx1_bac = cu1.second[4];
      break;
    }
    at = 0;
    for (auto cu2 : store_path[que[idx2].first]) {
      if (at != que[idx2].second) {
        at++;
        continue;
      }
      idx2_fro = cu2.second[0];
      idx2_bac = cu2.second[4];
      break;
    }

    ll pre_score = 0;
    if (idx1_pre != -1) pre_score += dist(idx1_pre, idx1_fro);
    if (idx1_aft != -1) pre_score += dist(idx1_bac, idx1_aft);
    if (idx2_pre != -1) pre_score += dist(idx2_pre, idx2_fro);
    if (idx2_aft != -1) pre_score += dist(idx2_bac, idx2_aft);
    swap(que[idx1], que[idx2]);

    ll aft_score = 0;
    if (idx1_pre != -1) aft_score += dist(idx1_pre, idx2_fro);
    if (idx1_aft != -1) aft_score += dist(idx2_bac, idx1_aft);
    if (idx2_pre != -1) aft_score += dist(idx2_pre, idx1_fro);
    if (idx2_aft != -1) aft_score += dist(idx1_bac, idx2_aft);

    return aft_score - pre_score;
  }

  void modify2(vp &que) {
    ll idx = rand() % M;
    ll idx_pre = -1, idx_aft = -1;
    ll at = 0;
    if (idx != 0) {
      for (auto cu1 : store_path[que[idx - 1].first]) {
        if (at != que[idx - 1].second) {
          at++;
          continue;
        }
        idx_pre = cu1.second[4];
        break;
      }
    }
    at = 0;
    if (idx != M - 1) {
      for (auto cu1 : store_path[que[idx + 1].first]) {
        if (at != que[idx + 1].second) {
          at++;
          continue;
        }
        idx_aft = cu1.second[0];
        break;
      }
    }

    ll best_score = INF;
    ll best = 0;
    at = 0;
    for (auto cu1 : store_path[que[idx].first]) {
      ll sum = 0;
      sum += cu1.first;
      if (idx_pre != -1) sum += dist(idx_pre, cu1.second[0]);
      if (idx_aft != -1) sum += dist(cu1.second[4], idx_aft);
      if (chmin(sum, best_score)) best = at;
      at++;
    }
    que[idx].second = best;
  }

  ll calc_score(vp &que) {
    ll ret = 0;
    ll cui = si, cuj = sj;
    rep(i, 0, M) {
      string sub_t = t[que[i].first];
      ll sub_t_idx = que[i].first;
      vll best_move(5);
      ll mi = INF;
      ll at = 0;
      for (auto cu1 : store_path[sub_t_idx]) {
        if (at != que[i].second) {
          at++;
          continue;
        }
        ll sum = 0;
        rep(i, 0, 5) {
          if (i == 0) sum += dist(cui * N + cuj, cu1.second[i]);
          else sum += dist(cu1.second[i - 1], cu1.second[i]);
        }
        if (chmin(mi, sum)) best_move = cu1.second;
        break;
      }
      cui = best_move[4] / N;
      cuj = best_move[4] % N;
      ret += mi;
    }
    return ret;
  }

  void print(vp &que) {
    ll cui = si, cuj = sj;
    rep(i, 0, M) {
      string sub_t = t[que[i].first];
      ll sub_t_idx = que[i].first;
      vll best_move(5);
      ll mi = INF;
      ll at = 0;
      for (auto cu1 : store_path[sub_t_idx]) {
        if (at != que[i].second) {
          at++;
          continue;
        }
        ll sum = 0;
        rep(i, 0, 5) {
          if (i == 0) sum += dist(cui * N + cuj, cu1.second[i]);
          else sum += dist(cu1.second[i - 1], cu1.second[i]);
        }
        if (chmin(mi, sum)) best_move = cu1.second;
        break;
      }
      cui = best_move[4] / N;
      cuj = best_move[4] % N;
      for (auto cu1 : best_move) cout << cu1 / N << " " << cu1 % N << endl;
    }
  }

  void solve() {
    auto start_time = system_clock::now();
    vp que;
    init(que);

    while (1) {
      if (1.9 < ela_times(start_time)) break;

      vp new_que = que;

      ll score;
      if (rand() % 3 == 0) {
        // modify2(new_que);
        score = -1;
      } else score = modify(new_que);

      if (score <= 0) {
        que = new_que;
        cerr << setprecision(3) << ela_times(start_time) << " ";
        cerr << calc_score(new_que) << endl;
      }
    }
    print(que);
  }
};

int main() {
  srand((unsigned int)time(NULL));

  cin >> N >> M;
  cin >> si >> sj;
  A.resize(N);
  rep(i, 0, N) cin >> A[i];
  t.resize(M);
  rep(i, 0, M) cin >> t[i];

  vvll char_all(26);
  rep(i, 0, N) rep(j, 0, N) char_all[A[i][j] - 'A'].push_back(i * N + j);

  rep(i, 0, M) {
    rep(init, 0, char_all[t[i][0] - 'A'].size()) {
      ll ma = 0;
      rep(j, 0, 5) chmax(ma, (ll)char_all[t[i][j] - 'A'].size());
      vvp dp(ma, vp(5, {INF, -1}));
      dp[0][0] = {0, -1};
      rep(_j, 1, 5) {
        rep(_i, 0, char_all[t[i][_j] - 'A'].size()) {
          rep(_k, 0, char_all[t[i][_j - 1] - 'A'].size()) {
            ll pre_cor = char_all[t[i][_j - 1] - 'A'][_k];
            if (_j == 1) pre_cor = char_all[t[i][0] - 'A'][init];
            ll cu_cor = char_all[t[i][_j] - 'A'][_i];
            ll add = dist(pre_cor, cu_cor);
            if (dp[_k][_j - 1].first + add < dp[_i][_j].first) {
              if (_j == 1) {
                dp[_i][_j] = pair(dp[_k][_j - 1].first + add, init);
              } else {
                dp[_i][_j] = pair(dp[_k][_j - 1].first + add, _k);
              }
            }
          }
        }
      }
      ll mi_idx = 0;
      rep(_i, 0, char_all[t[i][4] - 'A'].size()) {
        if (dp[_i][4].first < dp[mi_idx][4].first) {
          mi_idx = _i;
        }
      }
      vll path;
      path.push_back(char_all[t[i][4] - 'A'][mi_idx]);
      ll cu_i = mi_idx;
      ll cu_j = 4;
      while (1) {
        cu_i = dp[cu_i][cu_j].second;
        cu_j--;
        path.push_back(char_all[t[i][cu_j] - 'A'][cu_i]);
        if (cu_j == 0) break;
      }
      reverse(path.begin(), path.end());
      ll sum = 0;
      rep(j, 0, 4) sum += dist(path[j], path[j + 1]);

      store_path[i].push_back({sum, path});
    }
  }

  for (auto &cu : store_path) {
    Sort(cu.second);
  }
  // for (auto cu : store_path) {
  //   cout << cu.first.first << " " << cu.first.second << endl;
  //   for (auto cu2 : cu.second) {
  //     cout << "dist : " << cu2.first << endl;
  //     for (auto cu3 : cu2.second) cout << cu3 << " ";
  //     for (auto cu3 : cu2.second) cout << A[cu3 / N][cu3 % N];
  //     cout << endl;
  //   }
  // }
  SA sa;
  sa.solve();
}
