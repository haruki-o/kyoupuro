#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;

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
typedef vector<vector<ll>> vvll;
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
typedef unordered_set<ll> usi;
typedef unordered_set<ll> usll;
#define rep(i, l, n) for (ll i = (ll)(l); i < (ll)(n); i++)
#define repd(i, n, l) for (ll i = (ll)(n); i > (ll)(l); i--)
#define Copy(from, to) copy(from.begin(), from.end(), to.begin())
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
#define Unique(a) sort(a.begin(), a.end())
#define dPQ priority_queue<ll, vll, greater<ll>>
#define PQ priority_queue<ll, vll>
#define INF INT_MAX
#define INFF (9223372036854775805)
#define TIME_LIMIT (1.99)
#define def (10000)
// #define MOD (1000000007)
#define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision(10) << << endl;
//偏角ソートはlong double!
// auto ite = s.lower_bound("B");

// 状態の初期化
void init(vll &X) {
  map<ll, int> ma;
  rep(i, 0, 100) {
    ll at = rand();
    if (!ma.count(at % (def * 2) - def)) {
      X[i] = at % (def * 2) - def;
      ma[rand() % (def * 2) - def] = 1;
    } else
      i--;
  }
  gSort(X);
}

// 状態遷移
// void modify(STATE &state) {}

// 状態のスコア計算
long double calc_score(vll &X, vector<tuple<long double, ll, ll>> &ta, vll &a) {
  ll ind = 0;
  vll iti(10, 0);
  rep(i, 0, 100) {
    ll sum = 0;
    ll x1 = X[i], y1 = def - abs(X[i]);
    if (x1 < 0) {
      y1 *= -1;
      x1 *= -1;
    }
    // cout << i << " " << x1 << " " << y1 << endl;
    while (1) {
      ll x2 = get<1>(ta[ind]), y2 = get<2>(ta[ind]);
      if (x2 < 0) {
        y2 *= -1;
        x2 *= -1;
      }
      // cout << x2 << " " << y2 << endl;
      if (y2 * x1 >= y1 * x2)
        break;
      ind++;
      sum++;
      if (ind == (int)ta.size())
        break;
    }
    // cout << sum << endl;
    if (ind == (int)ta.size())
      break;
    if (0 < sum && sum <= 10)
      iti[sum - 1]++;
  }
  // cout << "iti[] : ";
  // rep(i, 0, 10) cout << iti[i] << " ";
  // cout << endl;

  ll score = 0;
  rep(i, 0, 10) score += min(a[i], iti[i]);
  ll sum = 0;
  rep(i, 0, 10) sum += a[i];

  return (double)score / sum * 1000000;
}

int main() {
  srand((unsigned int)time(NULL));

  int N, K;
  cin >> N >> K;
  vll a(10);
  rep(i, 0, 10) cin >> a[i];
  //(偏角,x,y)
  vector<tuple<long double, ll, ll>> ta(N);
  rep(i, 0, N) {
    ll x, y;
    cin >> x >> y;
    ta[i] = {atan2(y, x), x, y};
    if (y < 0)
      get<0>(ta[i]) += PI;
  }
  // rep(i, 0, N) cout << ta[i] << endl;
  Sort(ta);

  vll X(100);
  init(X);
  // cout << calc_score(X, ta, a) << endl;
  cout << 100 << endl;
  rep(i, 0, 100){
    cout << rand() % def << " " << rand() % def << " " << rand() % def << " " << rand() % def << endl;
  }
}