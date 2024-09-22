#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;

typedef long long ll;
// typedef pair<ll, ll> P;
// typedef pair<P, ll> Pll;
// typedef pair<ll, P> llP;
typedef tuple<ll, ll, ll> T;
typedef tuple<ll, ll, ll, ll> F;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vvi> vvvi;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
typedef vector<vvvll> vvvvll;
typedef vector<vvvvll> vvvvvll;
typedef vector<vector<double>> vvd;
typedef vector<double> vd;
// typedef vector<P> vP;
// typedef vector<vP> vvP;
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
#define def (301010)
// #define MOD (1000000007)
#define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision(10) << << endl;
//偏角ソートはlong double!
// auto ite = s.lower_bound("B");
// cout << __builtin_popcount(N) << endl;

struct Calendar {
  int X, M;
  vector<string> P;

  Calendar(int X, int M, vector<string> &P) : M(M), X(X), P(P) {}
  Calendar() {}

  void init(int X, int M) {
    X = X;
    M = M;
    rep(i, 0, M) rep(j, 0, 2 * X) P[i] += '9';
  }

  void output() { rep(i, 0, (int)P.size()) cout << P[i] << endl; }
};

struct Order {
  int X, M, C, E;
  ll d;
  int J;
  vll r,t;
  
  Order(int X,int M,int C,int E) : X(X),M(M),C(C),E(E){
    
  }
};
struct Solver {
  int X, M, C, E;
  ll costA[20][9], costB[20][9];
  ll bestScore = 0;

  void input_init() {
    cin >> X >> M >> C >> E;
    rep(i, 0, X) rep(j, 0, 9) cin >> costA[i][j] >> costB[i][j];
  }
  void input_reactive(int iter) {
    ll score;
    int V,D;
    cin >> score >> V >> D;
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < X; j++) {
        double cap;
        int letCnt;
        cin >> cap >> letCnt;
      }
    }
    bestScore = max(bestScore, score);
    cerr << "turn " << iter << " Score = " << score << ' '
         << ",  chVioN = " << V << ' ' << ",  letOpN = " << D << endl;
  }

  void output_() {

    for (int i = 0; i < M; i++) {
      string calendar;
      for (int j = 0; j < X; j++) {
        calendar.push_back('9');
        calendar.push_back('9');
      }
      cout << calendar << endl;
    }
  }
};

int main(int argc, char *argv[]) {
  Solver S;
  S.input_init();

  for (int i = 0; i < S.E; i++) {
    S.output_();
    S.input_reactive(i);
  }
  cerr << "Best Score = " << S.bestScore << endl;

  return 0;
}
