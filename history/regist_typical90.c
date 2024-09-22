#include <bits/stdc++.h>
// #include <atcoder/all>
#include <iostream>

using namespace std;
// using namespace atcoder;

// using mint = modint998244353;
// using mint = modint1000000007;

typedef long long ll;
typedef pair<ll, ll> P;
typedef tuple<ll, ll, ll> T;
typedef tuple<ll, ll, ll, ll> F;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vvi> vvvi;
typedef vector<ll> vll;
typedef vector<vector<ll>> vvll;
typedef vector<vvll> vvvll;
typedef vector<vvvll> vvvvll;
typedef vector<vector<double>> vvd;
typedef vector<double> vd;
typedef vector<P> vP;
typedef vector<vP> vvP;
typedef vector<T> vT;
typedef unordered_set<ll> usi;
typedef unordered_set<ll> usll;
#define rep(i, l, n) for (ll i = (ll)(l); i < (ll)(n); i++)
#define repd(i, n, l) for (ll i = (ll)(n); i > (ll)(l); i--)
#define Copy(from, to) copy(from.begin(), from.end(), to.begin())
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater<ll>())
#define Unique(a) sort(a.begin(), a.end())
#define dPQ prioritySo_queue<ll, vll, greater<ll>>
#define PQ priority_queue<ll, vll>
#define INF INT_MAX
#define INFF LLONG_MAX
#define TIME_LIMIT (2.0)
#define def (1000000)
// #define MOD (1000000007)
#define MOD (998244353)
#define PI (3.14159265359)

int main() {
  map<int,ll> ma;
  int N;cin>>N;
  vi a(N);
  rep(i,0,N)cin>>a[i];
  rep(i,0,N){
    if(ma.count(a[i]) != 0)ma[a[i]]++;
    else ma[a[i]] = 1;
  }

  ll ans = 0;
  rep(i,0,N){
    if(a[i] == 1){
      ans += ma[1] * (ma[1]);
      continue;
    }
    rep(j,1,a[i] + 1){
      if(j * j > a[i])break;
      if(a[i] % j == 0){
        if(j == 1){
          if(ma.count(1) != 0 && ma.count(a[i]) != 0){
            ans += ma[1] * (ma[a[i]]) * 2;
          }
        }
        else if(j == a[i] / j){
          if(ma[j] >= 2)ans += ma[j] * (ma[j]);
        }
        else if(ma.count(j) != 0 && ma.count(a[i] / j) != 0){
          ans += ma[j] * ma[a[i] / j] * 2;
        }
      }
    }
    // cout << i << " " << a[i] << " " << ans << endl;
  }
  cout << ans << endl;

}