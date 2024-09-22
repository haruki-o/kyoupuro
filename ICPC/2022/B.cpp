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
    vector<set<ll>> _g(N);
    rep(i, 0, N) {
      ll c1, c2;
      cin >> c1 >> c2;
      c1--;
      c2--;
      g[c1].push_back(i);
      g[c2].push_back(i);
      _g[i].insert(c1);
      _g[i].insert(c2);
    }
    rep(i, 0, N) Sort(g[i]);
    vvi seen(N, vi(2, 0));
    ll sta = -1;
    ll card = g[0][0];
    multiset<ll> se;

    rep(i, 0, N) {
      auto ite = _g[i].begin();
      ite++;
      if (*_g[i].begin() == *ite) {
        seen[i][0] = 1;
        seen[i][1] = 1;
      } else {
        se.insert(i);
        se.insert(i);

        if (sta == -1) {
          sta = i;
          card = g[i][0];
        }
      }
    }
    if(se.empty()){
      cout << 0 << endl;
      continue;
    }
    ll ans = 0;
    ll ind = 0;
    auto ite = se.begin();
    auto sta_ite = ite;
    // for (ll cu : se)
    //   cout << cu << " ";
    // cout << endl;
    ll sto = 0;
    while (!se.empty()) {
      ite = se.lower_bound(sto);
      // cout << "ite " << *ite << "  sto : " << sto << endl;
      sta_ite = ite;
      auto ite1 = ite;
      while (1) {
        while (*ite1 == *ite) {
          ite++;
          if (ite == se.end())
            ite = se.begin();
        }
        ll fro = *_g[*ite1].begin();
        ans++;
        // cout << "*ite " << *ite << " ite1 " << *ite1 << " fro "  << fro << endl;
        // cout << "_g[]" << endl;
        // rep(i,0,N){
        //   for(ll cu : _g[i])cout << cu <<  " ";
        //   cout << endl;
        // }
        if (_g[*ite].find(fro) != _g[*ite].end()) {
          // cout << "b";
          _g[*ite].erase(fro);
          _g[*ite1].erase(fro);
          sto = *ite;
          se.erase(sta_ite);
          se.erase(ite);
          // cout << "se[]";
          // for (ll cu : se)
          //   cout << cu << " ";
          // cout << endl;
          // cout << ans << endl;
          break;
        } else {
          // cout << "a";
          _g[*ite].insert(fro);
          _g[*ite1].erase(fro);
        }
        ite1 = ite;
      }
    }
      cout << ans << endl;
  }
}
