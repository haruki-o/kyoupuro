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
  int T,N,Q;cin>>T>>N>>Q;
  rep(i,0,T){
    vi at;
    //j = 3の時を除いてjを挿入する.
    rep(j,3,N+1){
      if(j == 3){
        rep(k,0,1){
          cout << k+1 << " " << k+2 << " " << k+3 << endl;
          cout.flush();
        }
        int ci;cin >> ci;

        rep(k,0,3)at.push_back((ci + k + 1)%3 + 1);
        continue;
      }
      int l = 0, r = (int)at.size() -1;
      while(r- l > 2){
        int mi = (l + r) /2;
        cout << at[mi] << " " << at[mi+1] << " " << j << endl;
        cout.flush();

        int ci; cin >> ci;
        if(ci == at[mi]){
          r = mi-1;
        }
        else if(ci == at[mi+1]){
          l = mi+2;
        }
        else{
          at.insert(at.begin() + mi+1, j);
          break;
        }
      }

      if((int)at.size() != j){
        if(r - l == 2){
          cout << at[r-1] << " " << at[r] << " " << j << endl;
          cout.flush();
          
          int ci; cin >> ci;
          if(ci == at[r])at.insert(at.begin() + r+1, j);
          if(ci == j)at.insert(at.begin() + r, j);
          if(ci == at[r-1]){
            cout << at[l] << " " << at[r-1] << " " << j << endl;
            cout.flush();

            cin >> ci;
            if(ci == j)at.insert(at.begin() + l+1, j);
            if(ci == at[l])at.insert(at.begin() + l, j);
          }
        }
        else if(r - l == 1){
          cout << at[l] << " " << at[r] << " " << j << endl;
          cout.flush();

          int ci; cin >> ci;
          if(ci == at[r])at.insert(at.begin() + r+1, j);
          else if(ci == at[l]) at.insert(at.begin() + l, j);
          else if(ci == j)at.insert(at.begin() + l+1, j);
        }
        else if(r == l){
          if(l == 0){
            cout << at[0] << " " << at[1] << " " << j << endl;
            cout.flush();
            int ci; cin >> ci;
            if(ci == at[0])at.insert(at.begin(), j);
            else if(ci == at[1])at.insert(at.begin() + 2, j);
            else if(ci == j)at.insert(at.begin() + 1, j);
          }
          else if(l == (int)at.size()-1){
            int si = (int)at.size() -1;
            cout << at[si] << " " << at[si-1] << " " << j << endl;
            cout.flush();
            int ci; cin >> ci;
            if(ci == at[si])at.insert(at.begin() + si + 1, j);
            else if(ci == at[si-1])at.insert(at.begin() + si -1, j);
            else if(ci == j)at.insert(at.begin() + si, j);
          }
          else{
            cout << at[l+1] << " " << at[l] << " " << j << endl;
            cout.flush();
            int ci; cin >> ci;
            if(ci == at[l])at.insert(at.begin() + l, j);
            else if(ci == j)at.insert(at.begin() + l+1, j); 
          }
        }
      }
      // cout << "at : ";
      // rep(k,0,(int)at.size())cout << at[k] << " ";
      // cout <<endl;
    }
    rep(i,0,N)cout << at[i] << " ";
    cout << endl;
    cout.flush();
    int re;cin >> re;
  }
}