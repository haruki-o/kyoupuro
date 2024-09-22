#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
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

ll beam_size = 6;
ll N, M, K;
vvll a;
vvvll s;

void dfs(ll cu_i, ll cu_j, vvll &ban, vll state, ll score, map<ll, vll, greater<ll>> &ma) {
  vll cu_ban(N, 0);
  cu_ban = ban[cu_i];
  rep(i, 0, cu_j) {
    rep(j, 0, 3) {
      cu_ban[i + j] += s[state[i]][0][j];
      cu_ban[i + j] %= MOD;
    }
  }
  vp candi(M);
  rep(i, 0, M) {
    ll sum = cu_ban[cu_j] + s[i][0][0];
    sum %= MOD;
    candi[i] = {sum, i};
  }
  gSort(candi);

  if (cu_j != 5) {
    rep(i, 0, beam_size) {
      state[cu_j] = candi[i].second;
      dfs(cu_i, cu_j + 1, ban, state, score + candi[i].first, ma);
    }
  } else {
    rep(i, 0, beam_size) {
      state[cu_j] = candi[i].second;
      ma[score + candi[i].first] = state;
    }
  }
}

struct Solver1 {
  void solve(vT &ans) {
    vvll ban(N, vll(N));
    rep(i, 0, N) rep(j, 0, N) ban[i][j] = a[i][j];

    rep(i, 0, N - 3) {
      map<ll, vll, greater<ll>> ma;
      vll state(N, -1);
      dfs(i, 0, ban, state, 0, ma);

      // cerr << ma.size() << endl;
      // ll idx = 0;
      // for(auto cu : ma){
      //   if(idx < 5)cerr << cu.first << endl;
      //   idx++;
      // }
      // cerr << endl;

      for (auto cu : ma) {
        rep(j, 0, N) {
          if (cu.second[j] == -1) break;
          rep(_i, 0, 3) {
            rep(_j, 0, 3) {
              ban[i + _i][j + _j] += s[cu.second[j]][_i][_j];
              ban[i + _i][j + _j] %= MOD;
            }
          }
        }
        rep(j, 0, N) {
          if (cu.second[j] == -1) break;
          ans.push_back({cu.second[j], i, j});
        }
        break;
      }
    }

    // 下の3行
    rep(j, 0, 6) {
      vll best(3, -1);
      ll best_score = -1;
      rep(x, -1, M) {
        rep(y, -1, M) {
          rep(z, -1, M) {
            vvll sum_s(3, vll(3, 0));
            rep(_i, 0, 3) {
              rep(_j, 0, 3) {
                if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
                if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
                if (z != -1) sum_s[_i][_j] += s[z][_i][_j];
                sum_s[_i][_j] %= MOD;
              }
            }
            ll score = 0;
            rep(i, 6, 9) score += ((ban[i][j] + sum_s[i - 6][0]) % MOD);
            if (chmax(best_score, score)) {
              best[0] = x, best[1] = y, best[2] = z;
            }
          }
        }
      }
      rep(i, 0, 3) if (best[i] != -1) ans.push_back({best[i], 6, j});

      vvll sum_s(3, vll(3, 0));
      rep(_i, 0, 3) {
        rep(_j, 0, 3) {
          ll x = best[0], y = best[1], z = best[2];
          if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
          if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
          if (z != -1) sum_s[_i][_j] += s[z][_i][_j];
          ban[6 + _i][j + _j] += sum_s[_i][_j];
          ban[6 + _i][j + _j] %= MOD;
        }
      }
    }

    // 右の3列
    rep(i, 0, 6) {
      vll best(3, -1);
      ll best_score = -1;
      rep(x, -1, M) {
        rep(y, -1, M) {
          rep(z, -1, M) {
            vvll sum_s(3, vll(3, 0));
            rep(_i, 0, 3) {
              rep(_j, 0, 3) {
                if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
                if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
                if (z != -1) sum_s[_i][_j] += s[z][_i][_j];
                sum_s[_i][_j] %= MOD;
              }
            }
            ll score = 0;
            rep(j, 6, 9) score += ((ban[i][j] + sum_s[0][j - 6]) % MOD);
            if (chmax(best_score, score)) {
              best[0] = x, best[1] = y, best[2] = z;
a            }
          }
        }
      }

      rep(j, 0, 3) if (best[j] != -1) ans.push_back({best[j], i, 6});

      vvll sum_s(3, vll(3, 0));
      rep(_i, 0, 3) {
        rep(_j, 0, 3) {
          ll x = best[0], y = best[1], z = best[2];
          if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
          if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
          if (z != -1) sum_s[_i][_j] += s[z][_i][_j];
          ban[i + _i][6 + _j] += sum_s[_i][_j];
          ban[i + _i][6 + _j] %= MOD;
        }
      }
    }

    // 最後の1マス
    vll best(7, -1);
    ll best_score = -1;
    rep(m, -1, M) {
      rep(n, m, M) {
        rep(o, n, M) {
          rep(p, o, M) {
            rep(x, p, M) {
              rep(y, x, M) {
                rep(z, y, M) {
                  vvll sum_s(3, vll(3, 0));
                  rep(_i, 0, 3) {
                    rep(_j, 0, 3) {
                      if (m != -1) sum_s[_i][_j] += s[m][_i][_j];
                      if (n != -1) sum_s[_i][_j] += s[n][_i][_j];
                      if (o != -1) sum_s[_i][_j] += s[o][_i][_j];
                      if (p != -1) sum_s[_i][_j] += s[p][_i][_j];
                      if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
                      if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
                      if (z != -1) sum_s[_i][_j] += s[z][_i][_j];
                      sum_s[_i][_j] %= MOD;
                    }
                  }
                  ll score = 0;
                  rep(i, 6, 9) rep(j, 6, 9) score += ((ban[i][j] + sum_s[i - 6][j - 6]) % MOD);
                  if (chmax(best_score, score)) {
                    best[0] = x, best[1] = y, best[2] = z;
                    best[3] = m, best[4] = n, best[5] = o, best[6] = p;
                  }
                }
              }
            }
          }
        }
      }
    }
    rep(j, 0, 7) if (best[j] != -1) ans.push_back({best[j], 6, 6});
  }
};

struct Solver2 {
  void solve(vT &ans, ll lower) {
    vvll ban(N, vll(N));
    rep(i, 0, N) rep(j, 0, N) ban[i][j] = a[i][j];

    // 左上 6 * 6
    rep(i, 0, 6) {
      rep(j, 0, 6) {
        vll best(2, -1);
        ll best_score = -1;
        rep(x, -1, M) {
          rep(y, -1, M) {
            vvll sum_s(3, vll(3, 0));
            rep(_i, 0, 3) {
              rep(_j, 0, 3) {
                if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
                if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
                sum_s[_i][_j] %= MOD;
              }
            }
            ll score = 0;
            score += ((ban[i][j] + sum_s[0][0]) % MOD);
            if (chmax(best_score, score)) {
              best[0] = x, best[1] = y;
            }
          }
          if (x == -1 && lower < best_score) break;
        }
        if (best[0] != -1) ans.push_back({best[0], i, j});
        if (best[1] != -1) ans.push_back({best[1], i, j});

        vvll sum_s(3, vll(3, 0));
        rep(_i, 0, 3) {
          rep(_j, 0, 3) {
            ll x = best[0], y = best[1];
            if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
            if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
            ban[i + _i][j + _j] += sum_s[_i][_j];
            ban[i + _i][j + _j] %= MOD;
          }
        }
      }
    }

    // 下の3行
    rep(j, 0, 6) {
      vll best(3, -1);
      ll best_score = -1;
      rep(x, -1, M) {
        rep(y, -1, M) {
          rep(z, -1, M) {
            vvll sum_s(3, vll(3, 0));
            rep(_i, 0, 3) {
              rep(_j, 0, 3) {
                if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
                if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
                if (z != -1) sum_s[_i][_j] += s[z][_i][_j];
                sum_s[_i][_j] %= MOD;
              }
            }
            ll score = 0;
            rep(i, 6, 9) score += ((ban[i][j] + sum_s[i - 6][0]) % MOD);
            if (chmax(best_score, score)) {
              best[0] = x, best[1] = y, best[2] = z;
            }
          }
        }
      }
      rep(i, 0, 3) if (best[i] != -1) ans.push_back({best[i], 6, j});

      vvll sum_s(3, vll(3, 0));
      rep(_i, 0, 3) {
        rep(_j, 0, 3) {
          ll x = best[0], y = best[1], z = best[2];
          if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
          if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
          if (z != -1) sum_s[_i][_j] += s[z][_i][_j];
          ban[6 + _i][j + _j] += sum_s[_i][_j];
          ban[6 + _i][j + _j] %= MOD;
        }
      }
    }

    // 右の3列
    rep(i, 0, 6) {
      vll best(3, -1);
      ll best_score = -1;
      rep(x, -1, M) {
        rep(y, -1, M) {
          rep(z, -1, M) {
            vvll sum_s(3, vll(3, 0));
            rep(_i, 0, 3) {
              rep(_j, 0, 3) {
                if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
                if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
                if (z != -1) sum_s[_i][_j] += s[z][_i][_j];
                sum_s[_i][_j] %= MOD;
              }
            }
            ll score = 0;
            rep(j, 6, 9) score += ((ban[i][j] + sum_s[0][j - 6]) % MOD);
            if (chmax(best_score, score)) {
              best[0] = x, best[1] = y, best[2] = z;
            }
          }
        }
      }

      rep(j, 0, 3) if (best[j] != -1) ans.push_back({best[j], i, 6});

      vvll sum_s(3, vll(3, 0));
      rep(_i, 0, 3) {
        rep(_j, 0, 3) {
          ll x = best[0], y = best[1], z = best[2];
          if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
          if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
          if (z != -1) sum_s[_i][_j] += s[z][_i][_j];
          ban[i + _i][6 + _j] += sum_s[_i][_j];
          ban[i + _i][6 + _j] %= MOD;
        }
      }
    }

    // 最後の1マス
    vll best(7, -1);
    ll best_score = -1;
    rep(m, -1, M) {
      rep(n, m, M) {
        rep(o, n, M) {
          rep(p, o, M) {
            rep(x, p, M) {
              rep(y, x, M) {
                rep(z, y, M) {
                  vvll sum_s(3, vll(3, 0));
                  rep(_i, 0, 3) {
                    rep(_j, 0, 3) {
                      if (m != -1) sum_s[_i][_j] += s[m][_i][_j];
                      if (n != -1) sum_s[_i][_j] += s[n][_i][_j];
                      if (o != -1) sum_s[_i][_j] += s[o][_i][_j];
                      if (p != -1) sum_s[_i][_j] += s[p][_i][_j];
                      if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
                      if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
                      if (z != -1) sum_s[_i][_j] += s[z][_i][_j];
                      sum_s[_i][_j] %= MOD;
                    }
                  }
                  ll score = 0;
                  rep(i, 6, 9) rep(j, 6, 9) score += ((ban[i][j] + sum_s[i - 6][j - 6]) % MOD);
                  if (chmax(best_score, score)) {
                    best[0] = x, best[1] = y, best[2] = z;
                    best[3] = m, best[4] = n, best[5] = o, best[6] = p;
                  }
                }
              }
            }
          }
        }
      }
    }
    rep(j, 0, 7) if (best[j] != -1) ans.push_back({best[j], 6, 6});
  }
};

ll calc_score(vT &ans) {
  vvll _ban(N, vll(N));
  rep(i, 0, N) _ban[i] = a[i];
  for (auto cu : ans) {
    ll i = get<1>(cu), j = get<2>(cu);
    rep(_i, 0, 3) {
      rep(_j, 0, 3) {
        _ban[i + _i][j + _j] += s[get<0>(cu)][_i][_j];
        _ban[i + _i][j + _j] %= MOD;
      }
    }
  }
  ll ret = 0;
  rep(i, 0, N) rep(j, 0, N) ret += _ban[i][j];
  return ret;
}

void last_one_math(vT &ans) {
  vvll ban(N, vll(N));
  rep(i, 0, N) ban[i] = a[i];
  for (auto cu : ans) {
    ll i = get<1>(cu), j = get<2>(cu);
    rep(_i, 0, 3) {
      rep(_j, 0, 3) {
        ban[i + _i][j + _j] += s[get<0>(cu)][_i][_j];
        ban[i + _i][j + _j] %= MOD;
      }
    }
  }
  if (K - ans.size() == 1) {
    vll best(1, -1);
    ll best_score = -1;
    rep(x, -1, M) {
      vvll sum_s(3, vll(3, 0));
      rep(_i, 0, 3) {
        rep(_j, 0, 3) {
          if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
          sum_s[_i][_j] %= MOD;
        }
      }
      ll score = 0;
      rep(i, 6, 9) rep(j, 6, 9) score += ((ban[i][j] + sum_s[i - 6][j - 6]) % MOD);
      if (chmax(best_score, score)) {
        best[0] = x;
      }
    }
    rep(j, 0, 1) if (best[j] != -1) ans.push_back({best[j], 6, 6});
  }
  if (K - ans.size() == 2) {
    vll best(2, -1);
    ll best_score = -1;
    rep(x, -1, M) {
      rep(y, x, M) {
        vvll sum_s(3, vll(3, 0));
        rep(_i, 0, 3) {
          rep(_j, 0, 3) {
            if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
            if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
            sum_s[_i][_j] %= MOD;
          }
        }
        ll score = 0;
        rep(i, 6, 9) rep(j, 6, 9) score += ((ban[i][j] + sum_s[i - 6][j - 6]) % MOD);
        if (chmax(best_score, score)) {
          best[0] = x, best[1] = y;
        }
      }
    }

    rep(j, 0, 2) if (best[j] != -1) ans.push_back({best[j], 6, 6});
  }
  if (K - ans.size() == 3) {
    vll best(3, -1);
    ll best_score = -1;
    rep(x, -1, M) {
      rep(y, x, M) {
        rep(z, y, M) {
          vvll sum_s(3, vll(3, 0));
          rep(_i, 0, 3) {
            rep(_j, 0, 3) {
              if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
              if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
              if (z != -1) sum_s[_i][_j] += s[z][_i][_j];
              sum_s[_i][_j] %= MOD;
            }
          }
          ll score = 0;
          rep(i, 6, 9) rep(j, 6, 9) score += ((ban[i][j] + sum_s[i - 6][j - 6]) % MOD);
          if (chmax(best_score, score)) {
            best[0] = x, best[1] = y, best[2] = z;
          }
        }
      }
    }

    rep(j, 0, 3) if (best[j] != -1) ans.push_back({best[j], 6, 6});
  }
  if (K - ans.size() == 4) {
    vll best(4, -1);
    ll best_score = -1;
    rep(o, -1, M) {
      rep(x, o, M) {
        rep(y, x, M) {
          rep(z, y, M) {
            vvll sum_s(3, vll(3, 0));
            rep(_i, 0, 3) {
              rep(_j, 0, 3) {
                if (o != -1) sum_s[_i][_j] += s[o][_i][_j];
                if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
                if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
                if (z != -1) sum_s[_i][_j] += s[z][_i][_j];
                sum_s[_i][_j] %= MOD;
              }
            }
            ll score = 0;
            rep(i, 6, 9) rep(j, 6, 9) score += ((ban[i][j] + sum_s[i - 6][j - 6]) % MOD);
            if (chmax(best_score, score)) {
              best[0] = x, best[1] = y, best[2] = z;
              best[3] = o;
            }
          }
        }
      }
    }
    rep(j, 0, 4) if (best[j] != -1) ans.push_back({best[j], 6, 6});
  }

  if (K - ans.size() == 5) {
    vll best(5, -1);
    ll best_score = -1;
    rep(n, -1, M) {
      rep(o, n, M) {
        rep(x, o, M) {
          rep(y, x, M) {
            rep(z, y, M) {
              vvll sum_s(3, vll(3, 0));
              rep(_i, 0, 3) {
                rep(_j, 0, 3) {
                  if (n != -1) sum_s[_i][_j] += s[n][_i][_j];
                  if (o != -1) sum_s[_i][_j] += s[o][_i][_j];
                  if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
                  if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
                  if (z != -1) sum_s[_i][_j] += s[z][_i][_j];
                  sum_s[_i][_j] %= MOD;
                }
              }
              ll score = 0;
              rep(i, 6, 9) rep(j, 6, 9) score += ((ban[i][j] + sum_s[i - 6][j - 6]) % MOD);
              if (chmax(best_score, score)) {
                best[0] = x, best[1] = y, best[2] = z;
                best[3] = n, best[4] = o;
              }
            }
          }
        }
      }
    }
    rep(j, 0, 5) if (best[j] != -1) ans.push_back({best[j], 6, 6});
  }
  if (6 <= K - ans.size()) {
    vll best(6, -1);
    ll best_score = -1;
    rep(m, -1, M) {
      rep(n, m, M) {
        rep(o, n, M) {
          rep(x, o, M) {
            rep(y, x, M) {
              rep(z, y, M) {
                vvll sum_s(3, vll(3, 0));
                rep(_i, 0, 3) {
                  rep(_j, 0, 3) {
                    if (m != -1) sum_s[_i][_j] += s[m][_i][_j];
                    if (n != -1) sum_s[_i][_j] += s[n][_i][_j];
                    if (o != -1) sum_s[_i][_j] += s[o][_i][_j];
                    if (x != -1) sum_s[_i][_j] += s[x][_i][_j];
                    if (y != -1) sum_s[_i][_j] += s[y][_i][_j];
                    if (z != -1) sum_s[_i][_j] += s[z][_i][_j];
                    sum_s[_i][_j] %= MOD;
                  }
                }
                ll score = 0;
                rep(i, 6, 9) rep(j, 6, 9) score += ((ban[i][j] + sum_s[i - 6][j - 6]) % MOD);
                if (chmax(best_score, score)) {
                  best[0] = x, best[1] = y, best[2] = z;
                  best[3] = m, best[4] = n, best[5] = o;
                }
              }
            }
          }
        }
      }
    }
    rep(j, 0, 6) if (best[j] != -1) ans.push_back({best[j], 6, 6});
  }
}

int main() {
  cin >> N >> M >> K;

  a.resize(N);
  rep(i, 0, N) a[i].resize(N);
  rep(i, 0, N) rep(j, 0, N) cin >> a[i][j];

  s.resize(M);
  rep(i, 0, M) {
    s[i].resize(3);
    rep(j, 0, 3) s[i][j].resize(3);
  }
  rep(i, 0, M) rep(j, 0, 3) rep(k, 0, 3) cin >> s[i][j][k];

  vT ans1;
  Solver1 solver1;
  solver1.solve(ans1);

  vT ans2;
  Solver2 solver2;
  ll best_score = -1;
  rep(i, 0, 10) {
    vT _ans;
    // 998244353
    // 920000000 から順に
    //   3000000 ずつ下限を下げていく
    ll lower = 920000000 - i * 3000000;
    solver2.solve(_ans, lower);
    cerr << _ans.size() << " " << lower <<  " " << calc_score(_ans) << endl;
    if (_ans.size() <= K) {
      if (chmax(best_score, calc_score(_ans))) {
        ans2 = _ans;
      }
    }
  }

  vT ans;
  if (calc_score(ans1) < calc_score(ans2)) ans = ans2;
  else ans = ans1;

  cerr << calc_score(ans1) << endl;
  cerr << calc_score(ans2) << endl;

  // 最後右下に
  last_one_math(ans);

  cerr << calc_score(ans) << endl;

  cout << ans.size() << endl;
  for (auto cu : ans) cout << get<0>(cu) << " " << get<1>(cu) << " " << get<2>(cu) << endl;
}
