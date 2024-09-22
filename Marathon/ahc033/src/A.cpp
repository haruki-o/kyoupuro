#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef tuple<ll, ll, ll> T;
typedef vector<T> vT;
typedef vector<vT> vvT;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
typedef pair<ll, ll> pll;
typedef pair<pll, ll> ppll;
typedef vector<ppll> vpp;
typedef vector<pll> vp;
typedef vector<vp> vvp;
typedef vector<char> vc;
typedef vector<int> vi;
typedef vector<string> vs;
typedef vector<vs> vvs;
#define rep(i, l, n) (ll)(l); i < (ll)(n); i++)
#define repd(i, n, l) (ll)(n); i > (ll)(l); i--)
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
template <class T, class S>
inline bool chmax(T &a, const S &b) {
  return (b, 1 : 0);
}
template <class T, class S>
inline bool chmin(T &a, const S &b) {
  return (b, 1 : 0);
}
#define INF INT_MAX
#define INFF (9223372036854775800)
#define TIME_LIMIT (3.5)
#define def (20101)
// #define MOD (1000000007)
// #define MAX (2147483647)
#define MAX (1073741823)
#define MOD (998244353)

ll N;
vvll A;

int main(){
  cin >> N;
  A.resize(N);
  rep(i,0,N){
    A[i].resize(N);
    rep(j,0,N)cin>> A[i][j];
  }


}

INSERT INTO
    DEPARTMENT
VALUES
    (
        '1100',
        S'2024-04-01',
      NULL,
        '営業1課',
2,
        UP_'1000',
        '1',
        UPDATE_'TANAKA',
        USER_UPDATE'2024-04-01 00:00:00'
    ),
    (
        '2200',
        S'2024-04-10',
      NULL,
        '製造2課',
2,
        UP_'2000',
        '1',
        UPDATE_'YAMADA',
        USER_UPDATE'2024-04-10 00:00:00'
    );
    
SELECT * FROM DEPARTMENT LIMIT 5;