#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

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
#define TIME_LIMIT (0.01)
#define ALL_TIME_LIMIT (0.25)
#define def (301010)
#define MOD (1000000007)
// #define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision(10) << << endl;
//偏角ソートはlong double!
// auto ite = s.lower_bound("B");
// cout << __builtin_popcount(N) << endl;

int main() {
  int sum = 0;
  rep(K, 2, 6) {
    rep(N, 3 * K + 9, 3 * K + 34) {
      rep(i, 0, 20) {
        string _s = ("0000" + to_string(sum));
        _s = _s.substr((int)_s.size() - 4, 4);

        string s = "";
        rep(j, 1, K + 1) rep(k, 0, 100) s += (char)('0' + j);
        rep(j, 0, N * N - 100 * K) s += '0';

        std::random_device seed_gen;
        std::mt19937 engine(seed_gen());
        std::shuffle(s.begin(), s.end(), engine);

        fstream output_fstream;
        output_fstream.open(_s + ".txt", std::ios_base::out);

        output_fstream << N << " " << K << endl;
        rep(j, 0, N) {
          rep(k, 0, N) output_fstream << s[j * N + k];
          output_fstream << endl;
        }
        sum++;
      }
    }
  }
}