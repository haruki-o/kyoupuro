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
typedef vector<pll> vp;
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
#define PI (3.14159265359)

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

int main() {
  int N;
  double p1, p2;
  cin >> N >> p1 >> p2;

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

  rep(i, 0, N + 1) {
    cout << setprecision(3) << fixed;
    cout << i << " " << bio_p(N, i, p1) << " " << bio_p(N, i, p2) << " ";
    cout << c_p1[i] << " " << c_p2[i] << endl;
  }

  double ma = 0.0;
  int X = 0;
  rep(i, 0, N + 1) if (chmax(ma, min(c_p1[i - 1], 1 - c_p2[i - 1]))) X = i;
  cout << X << " " << ma << endl;
}