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

void construct_G(vs &G, ll N, ll M) {
  // 完全グラフ1つと無向グラフ
  rep(i, 0, M / 2) {
    rep(j, 0, N) {
      rep(k, j + 1, N) {
        if (j < N - 1 - i && k < N - 1 - i)
          G[i] += '1';
        else
          G[i] += '0';
      }
    }
  }

  // [50 : 99] [0 : 49]の頂点数を選んでウニグラフ
  rep(i, M / 2, M) {
    rep(j, 0, N) {
      rep(k, j + 1, N) {
        if (j <= (i - M / 2) || k <= (i - M / 2))
          G[i] += '1';
        else
          G[i] += '0';
      }
    }
  }

  // rep(i, 0, M) {
  //   ll idx = 0;
  //   vll V_H(N, 0);
  //   rep(j, 0, N) {
  //     rep(k, j + 1, N) {
  //       if (G[i][idx] == '1') {
  //         V_H[j]++;
  //         V_H[k]++;
  //       }
  //       idx++;
  //     }
  //   }
  //   cout << "G[i] : " << i << "  V_H" << endl;
  //   for (ll cu : V_H)
  //     cout << cu << " ";
  //   cout << endl;
  // }
}

ll select_N(ll M, double ep) {
  ll best_N = 0;
  long double  best_score = -1;
  vs G(M, "");
  rep(N, 1, (ll)101) {
    // cout << "select_N N : " << N << endl;
    G.assign(M, "");
    construct_G(G, N, M);

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

      // cout << (s % M == t) << " ma : " << ma << " ";
      // cout << "s,t : " << s % M << " " << t << " ";
      // cout << "sum : " << sum << " ";
      // cout << endl;
      if (s % M == t)
        correct++;
    }

    long double _score = 1.0;
    correct /= 5;
    rep(_i, 0, 100 - correct)_score *= 0.9;
    _score /= N;
    if (best_score < _score) {
      best_N = N;
      best_score = _score;
    }
    // cerr << "N : " << N << " correct : " << correct << endl;
  }
  if(best_N == 1)best_N = 100;

  return best_N;
}

int main() {
  srand((unsigned int)time(NULL));
  ll M;
  cin >> M;
  double ep;
  cin >> ep;

  ll N = select_N(M, ep);

  // [0 : 49] 完全グラフ1つと無向グラフ
  // [50 : 99] [0 : 49]の頂点数を選んでウニグラフ
  vs G(M, "");
  construct_G(G, N, M);

  // cout << N << endl;
  // rep(i, 0, M) cout << (int)G[i].size() << " ";
  // rep(i, 0, M) cout << G[i] << endl;
  cout.flush();

  ll correct = 0;
  ll f_correct = 0;
  ll s_correct = 0;
  rep(k, 0, 100) {
    //  string H;
    //  cin >> H;

    ll s;
    cin >> s;
    string H = test(G[s], ep);

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
    ll t = rand() % M;
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

    // cout << t << endl;
    // cout.flush();

    // cout << (s == t) << " ma : " << ma << " ";
    // cout << "s,t : " << s << " " << t << " ";
    // cout << endl;
    // cout << "V_H" << endl;
    // for (ll cu : V_H)
    //   cout << cu << " ";
    // cout << endl;
    if (s == t)
     correct++;
  }

  // cerr << "correct : " << correct << endl;
  cerr << N << " " << correct << endl;
}
