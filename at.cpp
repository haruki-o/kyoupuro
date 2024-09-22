#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef pair<ll, ll> P;
typedef tuple<ll, ll, ll> T;
typedef tuple<ll, ll, ll, ll> F;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vvi> vvvi;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
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
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
#define INF INT_MAX
#define INFF (9223372036854775800)
#define TIME_LIMIT (5.96)
#define def (2010101)
// #define MOD (1000000007)
#define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision() << << endl;

int main() {
  ll N, M;
  cin >> N >> M;
  vll A(M), B(M), C(M);
  rep(i, 0, M) cin >> A[i] >> B[i] >> C[i];

  vvP g(N);
  rep(i, 0, M) A[i]--, B[i]--;
  rep(i, 0, M) {
    g[A[i]].push_back({B[i], C[i]});
    g[B[i]].push_back({A[i], C[i]});
  }

  ll ans = INFF;

  rep(i, 1, (1 << N)) {
    vvll dp(N, vll(N, 1e15));
    rep(j, 0, M) {
      if (((1 << A[j]) & i) && ((1 << B[j]) & i)) {
        dp[A[j]][B[j]] = C[j];
        dp[B[j]][A[j]] = C[j];
      }
    }
    // if (i == ((1 << N) - 1)) {
    //   rep(j, 0, N) {
    //     rep(k, 0, N) cout << dp[j][k] << " ";
    //     cout << endl;
    //   }
    // }
    rep(j, 0, N) {
      rep(k, 0, N) {
        rep(l, 0, N) { dp[j][k] = min(dp[j][k], dp[j][l] + dp[l][k]); }
      }
    }
    rep(j, 0, N) {
      rep(k, 0, N) {
        if (dp[j][k] < ans) ans = dp[j][k];
      }
    }

    if (i == ((1 << N) - 1)) {
      rep(j, 0, N) {
        rep(k, 0, N) cout << dp[j][k] << " ";
        cout << endl;
      }
    }
  }

  cout << ans << endl;
}