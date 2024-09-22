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
typedef vector<vp> vvp;
typedef vector<vvp> vvvp;
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
ll stat_loop = 0;
ll stat_ans = 0;

ll N, D, Q;
vll w;

bool input_c(char &c, vll &l, vll &r) {
  ll suml = 0, sumr = 0;
  for (ll cu : l)
    suml += w[cu];
  for (ll cu : r)
    sumr += w[cu];

  map<ll, ll> ma;
  for (ll cu : l)
    ma[cu] = 1;
  for (ll cu : r)
    if (ma.count(cu)) return false;

  if (suml < sumr) c = '<';
  else if (suml > sumr) c = '>';
  else c = '=';

  return true;
}

int main() {
  cin >> N >> D >> Q;
#ifdef _DEBUG
  w.resize(N);
  rep(i, 0, N) cin >> w[i];
#endif

  vll A;
  A.push_back(0);
  ll _l = -1, _r = 1;

  ll sort_num = 0;
  while ((ll)A.size() != N) {
    sort_num++;
    ll _mid = (_l + _r) / 2;

    ll nl = 1, nr = 1;
    vll l(nl), r(nr);
    if ((ll)A.size() == 1) {
      l[0] = A[0];
      r[0] = A.size();
    } else {
      l[0] = A[_mid];
      r[0] = A.size();
    }

    cout << nl << " " << nr << " ";
    for (ll cu : l)
      cout << cu << " ";
    for (ll cu : r)
      cout << cu << " ";
    cout << endl;
    cout.flush();

    char c;
#ifdef _DEBUG
    if (!input_c(c, l, r)) {
      cerr << "exist same element" << endl;
      return 0;
    }
#else
    cin >> c;
#endif

    int flag = 0;
    if (c == '<') _l = _mid;
    else if (c = '>') _r = _mid;
    else {
      _l = _mid;
      flag = 1;
    }

    if (_r - _l <= 1) flag = 1;

    if (flag) {
      A.insert(A.begin() + _l + 1, (ll)A.size());
      _l = -1, _r = (ll)A.size();
    }

    if (flag) {
      cout << "sort_num: " << sort_num << endl;
      for (ll cu : A)
        cout << w[cu] << " ";
      cout << endl;
    }
  }

  cerr << "N, Q     : " << N << " " << Q << endl;
  cerr << "sort_num: " << sort_num  << " " << 2 *  N<< endl;
}