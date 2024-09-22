//#include </usr/include/c++/v1/>
#include <bits/stdc++.h>

using namespace std;
using namespace chrono;

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
#define TIME_LIMIT (2.0)
#define def (201010)
// 1000000007
#define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision(10) << << endl;

auto startClock = system_clock::now();

int N, initT;
vvi initG;
vvi G;
vvi maG;
vi rest(16, 0);
vvi seen;
vi dy = {0, -1, 0, 1}, dx = {-1, 0, 1, 0};
vi pow2 = {1, 2, 4, 8};

void init() {
  cin >> N >> initT;
  initG.resize(N);
  G.resize(N);
  maG.resize(N);
  rep(i, 0, N) {
    string s;
    cin >> s;
    initG[i].resize(N);
    G[i].assign(N, 0);
    maG[i].resize(N);
    rep(j, 0, N) {
      if (s[j] <= '9')
        initG[i][j] = s[j] - '0';
      else
        initG[i][j] = s[j] - 'a' + 10;
    }
  }
  // rep(i, 0, N) {
  //   rep(j, 0, N)printf(" %2d ", initG[i][j]);
  //   printf("\n");
  // }
  rep(i, 0, N) rep(j, 0, N) rest[initG[i][j]]++;

  seen.resize(N);
  rep(i, 0, N) seen[i].assign(N, 0);
}

void out(vvi g) {
  rep(i, 0, N) {
    rep(j, 0, N) {
      if (G[i][j] < 10)
        cout << g[i][j];
      else
        cout << (char)(G[i][j] - 10 + 'a');
    }
    cout << endl;
  }
}

int ma = 1;
// di: {左、上、右、下}

bool dfs(int h, int w, int di, int sum) {
  // cout << h << " " << w << " " << di << " " << sum << endl;
  const double time =
      duration_cast<microseconds>(system_clock::now() - startClock).count() *
      1e-6;
  if (time > TIME_LIMIT) {
    cout << "Time limit " << endl;
    cout << ma << endl;
    out(maG);
    return;
  }

  if (sum > ma) {
    ma = sum;
    maG = G;
  }
  if (sum > 1 && (G[3][3] == 2 && G[2][3] == 13 && G[2][2] == 6))
    out(G);

  // k方向にj種類目のタイルを置く
  rep(k, 0, 4) {
    if (h + dy[k] < 0 || h + dy[k] >= N || w + dx[k] < 0 || w + dx[k] >= N)
      continue;
    if (seen[h + dy[k]][w + dx[k]] == 1)
      continue;
    //置く方向に、今いるタイルが(h,w)方向に伸びているか
    if (!((1 << k) & G[h][w]))
      continue;
    // cout << "k " << k << endl;
    rep(j, 1, 16) {
      if (rest[j] == 0)
        continue;
      // j種類目のタイルに（ｋ＋２）％２方向の道があるか
      if (!((1 << ((k + 2) % 4)) & j))
        continue;
      // cout << h << " " << w << " j " << j <<  " "<<endl;
      int ju = 1;
      rep(i, 0, 4) {
        if ((1 << i) & j) {
          int ath = h + dy[i];
          int atw = w + dx[i];
          if (ath < 0 || ath >= N || atw < 0 || atw >= N) {
            ju = 0;
            continue;
          }
          if (seen[ath][atw] == 1)
            ju = 0;
        }
      }
      //残りもあり&&無駄な辺がない＆＆はみ出してない
      if (ju == 1) {
        rest[j]--;
        seen[h + dy[k]][w + dx[k]] = 1;
        G[h + dy[k]][w + dx[k]] = j;
        cout << j << endl;
        dfs(h + dy[k], w + dx[k], k, sum + 1);
        seen[h + dy[k]][w + dx[k]] = 0;
        G[h + dy[k]][w + dx[k]] = 0;
        rest[j]++;
      }
    }
  }
}

int main() {
  init();
  rep(i, 0, 4) {
    if (rest[(1 << i)] != 0) {
      seen[N / 2][N / 2] = 1;
      G[N / 2][N / 2] = (1 << i);
      dfs(N / 2, N / 2, i, 1);
      seen[N / 2][N / 2] = 0;
      G[N / 2][N / 2] = 0;
    }
  }
  cout << ma << endl;
  out(maG);
}

// 6
// 2 2 0 3
// 000000
// 000000
// 006d00
// 000200
// 000000
// 000000