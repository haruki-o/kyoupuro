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
typedef pair<pll, ll> ppll;
typedef tuple<ll, ll, ll> tll;
typedef vector<tll> vt;
typedef vector<ppll> vpp;
typedef vector<pll> vp;
typedef vector<vp> vvp;
typedef vector<char> vc;
typedef vector<int> vi;
typedef vector<string> vs;
typedef vector<vs> vvs;
#define rep(i, l, n) for (ll i = (ll)(l); i < (ll)(n); i++)
#define repd(i, n, l) for (ll i = (ll)(n); i > (ll)(l); i--)
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
template <class T, class S>
inline bool chmax(T &a, const S &b) {
  return (a < b ? a = b, 1 : 0);
}
template <class T, class S>
inline bool chmin(T &a, const S &b) {
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

struct UnionFind {
  //根のとき-1,iの親はpar[i]
  vi par;
  vll sum;
  UnionFind(ll N) {
    par.resize(N);
    sum.assign(N, 1);
    rep(i, 0, N) par[i] = i;
  }
  // cの根っこを返す
  ll root(ll c) {
    if (par[c] == c)
      return c;
    ll a = root(par[c]);
    par[c] = a;
    return a;
  }
  // aをbに結合(違う時)
  bool unite(ll a, ll b) {
    ll ra = root(a);
    ll rb = root(b);
    if (ra == rb)
      return 0;
    par[ra] = rb;
    sum[rb] += sum[ra];
    return 1;
  }
  bool same(ll a, ll b) {
    if (root(a) == root(b))
      return 1;
    return 0;
  }
};

ll N, _N, allN;
vll A, B;

int main() {
  _N = 0;
  srand((unsigned int)time(NULL));
  auto startTime = system_clock::now();

  cin >> N;
  A.resize(N);
  B.resize(N);
  rep(i, 0, N) cin >> A[i] >> B[i];

  A.push_back(0);
  B.push_back(0);
  _N++;
  map<pll, ll> ma;
  rep(i, 0, N) ma[{A[i], B[i]}] = 1;
  rep(i, 0, 4 * N - 10) {
    while (1) {
      ll _x = rand() % ((ll)pow(10, 9) + 1);
      ll _y = rand() % ((ll)pow(10, 9) + 1);
      if (ma.count({_x, _y}) == 0) {
        ma[{_x, _y}] = 1;
        A.push_back(_x);
        B.push_back(_y);
        _N++;
        break;
      }
    }
  }
  allN = N + _N;
  vvll g(allN);
  rep(i, 0, allN) {
    rep(j, 0, allN) {
      if (A[i] <= A[j] && B[i] <= B[j]) {
        g[i].push_back(j);
      }
    }
  }
  
  vvll ansG(5 *N);
  UnionFind uf(5 * N);
  priority_queue<tll,vt, greater<tll>> pq;
  rep(cu,0,allN){
    for(ll to: g[cu]){
      ll dist = abs(A[cu] - A[to]) + abs(B[cu] - B[to]);
      pq.push({dist, cu, to});
    }
  }
  while(!pq.empty()){
    auto cu = pq.top();
    pq.pop();
    ll p1 = get<1>(cu);
    ll p2 = get<2>(cu);
    if(uf.unite(p1,p2)){
      ansG[p1].push_back(p2);
    }
  }

  vector<pair<pll, pll>> ans;
  vll seen(allN, 0);
  seen[N] = 1;
  queue<ll> qu;
  qu.push(N);
  while (!qu.empty()) {
    ll cu = qu.front();
    qu.pop();
    for (ll to : ansG[cu]) {
      if (seen[to] == 1) continue;
      seen[to] = 1;
      qu.push(to);
      ans.push_back({{A[cu], B[cu]}, {A[to], B[to]}});
    }
  }

  cout << ans.size() << endl;
  rep(i, 0, ans.size()) {
    cout << ans[i].first.first << " " << ans[i].first.second << " ";
    cout << ans[i].second.first << " " << ans[i].second.second << endl;
  }
}