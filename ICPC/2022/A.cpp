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
#define INFF (9223372036854775800)
#define TIME_LIMIT (1.99)
#define def (4010101)
// #define MOD (1000000007)
#define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision(10) << << endl;
//偏角ソートはlong double!
// auto ite = s.lower_bound("B");

int main() {
  while (1) {
    ll N;
    cin >> N;
    if (N == 0)
      break;
    vvi g(N);
    vvi _g(N);
    rep(i, 0, N) {
      ll c1, c2;
      cin >> c1 >> c2;
      c1--;
      c2--;
      g[c1].push_back(i);
      g[c2].push_back(i);
      _g[i].push_back(c1);
      _g[i].push_back(c2);
    }
    rep(i, 0, N) Sort(g[i]);
    rep(i, 0, N) Sort(_g[i]);
    vvi seen(N, vi(2, 0));
    rep(i, 0, N) {
      if (_g[i][0] == _g[i][1]) {
        seen[i][0] = 1;
        seen[i][1] = 1;
      }
    }
    ll sta = 0;
    ll card = g[0][0];
    rep(i, 0, N) {
      if (_g[i][0] != _g[i][1]) {
        sta = i;
        card = g[i][0];
        break;
      }
    }
    ll ans = 0;
    multiset<ll> se;
    rep(i, 0, N) {
      if (_g[i][0] != _g[i][1]) {
        rep(j, 0, 2) { se.insert(i); }
      }
    }
    ll ind = 0;

    while (!se.empty()) {
      ll at = (sta + 1) % N;
      while (1) {
        cout << "at " << at << endl;
        if (seen[at][0] == 0 || seen[at][1] == 0)
          ans++;
        if (card == _g[at][0] || card == _g[at][1])
          break;
        at = (at + 1) % N;
      }
      rep(j, 0, 2) {
        if (card == _g[g[card][0]][j])
          seen[g[card][0]][j] = 1;
        if (card == _g[g[card][1]][j])
          seen[g[card][1]][j] = 1;
      }
      se.erase(se.find(g[card][0]));
      se.erase(se.find(g[card][1]));

      cout << "ans :" << ans << endl;
      for (ll cu : se)
        cout << cu << " ";
      cout << endl;

      // sta更新
      sta = (at + 1)% N;
      if (se.lower_bound(sta) == se.end()) {
        sta = *se.begin();
      } else {
        sta = *se.lower_bound(sta);
      }
      if (se.lower_bound(sta) == se.end()) {
        sta = *se.begin();
      } else {
        sta = *se.lower_bound(sta);
      }
      if (seen[sta][0] == 0)
        card = _g[sta][0];
      else
        card = _g[sta][1];
      cout << "next " << sta << " " << card << endl;
    }
    cout << ans << endl;
  }
}
