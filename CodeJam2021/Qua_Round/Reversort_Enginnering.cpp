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
    cout << "Case #" << i+1 << ": ";
    int N,C;cin>>N>>C;
    if(C < N-1 || (N+1)*N/2 -1 < C){
      cout << "IMPOSSIBLE" << endl;
      continue;
    }
    vi ans(N);
    //0で前、1で後ろ
    int pre = 0;
    int l = 0, r = N-1;
    C = C - N + 1;
    rep(j,0,N){
      if(j == 0){
        if(C >= N-j-1){
          ans[N-1] = j+1;
          r = N-2;
          pre = 1;
          C -= N-j-1;
        }
        else{
          ans[0] = j+1;
          l = 1;
          pre = 0;
        }
      }
      else{
        if(C >= N-j-1){
          if(pre == 0){
            ans[r] = j+1;
            r--;
          }
          else{
            ans[l] = j+1;
            l++;
          }
          C -= N-j-1;
          pre = 1 - pre;
        }
        else{
          if(pre == 1){
            ans[r] = j+1;
            r--;
          }
          else{
            ans[l] = j+1;
            l++;
          }
        }
      }
    }
    rep(j,0,N)cout << ans[j] << " ";
    cout << endl;


  }
}