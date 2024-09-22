//#include </usr/include/c++/v1/>
#include <bits/stdc++.h>

using namespace std;

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
#define INFF LLONG_MAX
#define TIME_LIMIT (1.99)
#define def (4010101)
// 1000000007
#define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision(10) << << endl;

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
void out(int N,vvi g) {
  cout << N << " " << 2 * N * N << endl;
  rep(i, 0, N) {
    rep(j, 0, N) {
      if (g[i][j] < 10)
        cout << g[i][j];
      else
        cout << (char)(g[i][j] - 10 + 'a');
    }
    cout << endl;
  }
}

int main() {
  vll sum(16,0);
  cout << "Qを入力" << endl;
  int Q;cin>>Q;
  rep(i,0,Q){
  cout << "Nを入力" << endl;
    int N;
    cin >> N;
    vector<pair<int,int>> a;
    rep(i, 0, N - 1) {
      rep(j, 0, N) if (i != N - 2 || j != N - 1) a.push_back(pair(i*N + j, (i+1)*N + j));
    }
    rep(i, 0, N) {
      rep(j, 0, N - 1) if (i != N - 1 || j != N - 2) a.push_back(pair(i*N+j,i*N+j+1));
    }

    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    shuffle(a.begin(), a.end(), engine);

    UnionFind uf(N*N);
    vector<pair<int,int>> apt;
    rep(i, 0, (int)a.size()) {
      if (uf.unite(a[i].first, a[i].second)) {
        apt.push_back(pair(a[i].first, a[i].second));
      }
    }

    vvi g(N, vi(N, 0));
    rep(i, 0, (int)apt.size()) {
      int fi = apt[i].first;
      int se = apt[i].second;
      // cout << '{' << fi/N << "," << fi%N << "},{" << se/N << " " << se%N << '}' << endl;
      //同じ行
      if (fi / N == se / N) {
        g[fi/N][fi%N] += 4;
        g[se/N][se%N] += 1;
      }
      //同じ列
      if (fi % N == se % N) {
        g[fi/N][fi%N] += 8;
        g[se/N][se%N] += 2;
      }
    }
    rep(i,0,N)rep(j,0,N)sum[g[i][j]]++;
  out(N,g);

  }
  // ll all = 0;
  // rep(i,0,16)all +=sum[i];
  // cout << "all " << all << endl;
  // rep(i,0,16)cout << sum[i] << endl;

}

//入力 N ([6:10])
//出力 char chg[N-1][N-1]