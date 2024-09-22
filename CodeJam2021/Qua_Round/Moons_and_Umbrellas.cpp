#include <bits/stdc++.h>
//#include <atcoder/all>
#include <iostream>

using namespace std;
//using namespace atcoder;

//using mint = modint998244353;

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
#define Copy(from, to) copy(from.begin(), from.end(),to.begin())
#define Sort(a) sort(a.begin(),a.end())
#define gSort(a) sort(a.begin(),a.end(),greater<ll>())
#define Unique(a) sort(a.begin(),a.end())
#define dPQ prioritySo_queue<ll,vll,greater<ll>>
#define PQ priority_queue<ll,vll> 
#define INF INT_MAX
#define INFF LLONG_MAX
#define TIME_LIMIT (2.0)
#define def (201010)
//1000000007 998244353
#define MOD (998244353)
#define PI (3.14159265359)

int main(){
  int T;cin>>T;
  rep(i,0,T){
    int X,Y;cin>>X>>Y;
    string s;cin>>s;
    int N = (int)s.size();
    vvi dp(2,vi(N+2,0));
    if(s[0] == 'C')dp[1][0] = def;
    if(s[0] == 'J')dp[0][0] = def;
    rep(j,1,N){
      if(s[j] == '?'){
        dp[0][j] = min(dp[0][j-1],dp[1][j-1] + Y);
        dp[1][j] = min(dp[1][j-1],dp[0][j-1] + X);
      }
      if(s[j] == 'C'){
        dp[0][j] = min(dp[0][j-1],dp[1][j-1] + Y);
        dp[1][j] = def;
      }
      if(s[j] == 'J'){
        dp[0][j] = def;
        dp[1][j] = min(dp[1][j-1],dp[0][j-1] + X);
      }
    }
    cout << "Case #" << i+1 << ": " << min(dp[0][N-1], dp[1][N-1]) << endl;
  }
}