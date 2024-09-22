//#include </usr/include/c++/v1/>
#include <bits/stdc++.h>

using namespace std;
using namespace chrono;

typedef long long ll;
typedef pair<int, int> P;
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
#define TIME_LIMIT (0.002)
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

int N, initT;
vvi initG;
vvi G;
vvi G2;
vi parts(16, 0);
vi parts2(16, 0);
vi dy = {0, -1, 0, 1}, dx = {-1, 0, 1, 0};
vi pow2 = {1, 2, 4, 8};
int pre_delete_h = -1, pre_delete_w = -1;
int all_sco = 0;
vvi add_E_sum;
vvi delete_E_sum;



bool range_out(int h, int w) {
  if (h < 0 || h >= N || w < 0 || w >= N)
    return true;
  return false;
}

void out_1vec(vi g) {
  rep(i, 0, (int)g.size()) {
    cout << g[i] << " ";
    if (i % N == N - 1)
      cout << endl;
  }
}

void out_parts(vi g) {
  // rep(i, 0, 16) printf("%2d ", (int)i);
  // cout << endl;
  rep(i, 0, 16) printf("%2d ", g[i]);
  cout << endl;
}

bool same_g() {
  bool ju = true;
  rep(i, 0, 16) if (parts[i] != parts2[i]) ju = false;
  return ju;
}

void init() {
  cin >> N >> initT;
  initG.resize(N);
  G.resize(N);
  G2.resize(N);
  rep(i, 0, N) {
    string s;
    cin >> s;
    initG[i].resize(N);
    G[i].resize(N);
    G2[i].resize(N);
    rep(j, 0, N) {
      if (s[j] <= '9')
        initG[i][j] = s[j] - '0';
      else
        initG[i][j] = s[j] - 'a' + 10;
    }
  }
  G = initG;
  rep(i, 0, N) rep(j, 0, N) parts[initG[i][j]]++;
}

vvi generate() {
  vector<pair<int, int>> a;
  rep(i, 0, N - 1) {
    rep(j, 0, N) if (i != N - 2 || j != N - 1)
        a.push_back({i * N + j, (i + 1) * N + j});
  }
  rep(i, 0, N) {
    rep(j, 0, N - 1) if (i != N - 1 || j != N - 2)
        a.push_back({i * N + j, i * N + j + 1});
  }

  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());
  shuffle(a.begin(), a.end(), engine);

  UnionFind uf(N * N);
  vector<pair<int, int>> apt;
  rep(i, 0, (int)a.size()) {
    if (uf.unite(a[i].first, a[i].second)) {
      apt.push_back({a[i].first, a[i].second});
    }
  }

  vvi g(N, vi(N, 0));
  rep(i, 0, (int)apt.size()) {
    int fi = apt[i].first;
    int se = apt[i].second;
    // cout << '{' << fi/N << "," << fi%N << "},{" << se/N << " " << se%N << '}'
    // << endl;
    //同じ行
    if (fi / N == se / N) {
      g[fi / N][fi % N] += 4;
      g[se / N][se % N] += 1;
    }
    //同じ列
    if (fi % N == se % N) {
      g[fi / N][fi % N] += 8;
      g[se / N][se % N] += 2;
    }
  }

  rep(i, 0, N) rep(j, 0, N) parts2[g[i][j]]++;

  rep(i, 0, N) {
    rep(j, 0, N) {
      rep(k, 0, 4) {
        if (!((1 << k) & g[i][j]))
          continue;
        int ath = g[i][j] + dy[k];
        int atw = g[i][j] + dx[i];
        if (ath < 0 || ath >= N || atw < 0 || atw >= N)
          continue;
        if (!(1 << ((k + 2) % 4) & g[ath][atw]))
          continue;
      }
    }
  }
  return g;
}

P add_edge() {
  int b_h = -1, b_w = -1, b_di = -1;
  int sco = 0; // 0,1,2,3,4
  vvP all(5);
  rep(i, 0, N * N) {
    int cuh = i / N;
    int cuw = i % N;
    //右,下のみ
    rep(j, 2, 4) {
      //辺があったらダメ
      if ((G2[cuh][cuw] & (1 << j)))
        continue;
      int ath = cuh + dy[j], atw = cuw + dx[j];
      if (range_out(ath, atw))
        continue;
      // N-1,N-1はだめ
      if (ath == N - 1 && atw == N - 1)
        continue;
      //辺(cuh,cuw) -> (ath,atw)を加えた時
      int at_sco = 0;
      //パーツは変更前(変更前パーツは1つ減る)
      if (parts2[G2[ath][atw]] > parts[G2[ath][atw]])
        at_sco++;

      if (parts2[G2[cuh][cuw]] > parts[G2[cuh][cuw]])
        at_sco++;

      //パーツ変更後
      if (ath == cuh) {
        if (parts2[G2[cuh][cuw] + 4] < parts[G2[cuh][cuw] + 4])
          at_sco++;
        if (parts2[G2[ath][atw] + 1] < parts[G2[ath][atw] + 1])
          at_sco++;
      } else {
        if (parts2[G2[cuh][cuw] + 8] < parts[G2[cuh][cuw] + 8])
          at_sco++;
        if (parts2[G2[ath][atw] + 2] < parts[G2[ath][atw] + 2])
          at_sco++;
      }
      all[at_sco].push_back({cuh * N + cuw, j});
    }
  }
  repd(i, 4, -1) {
    if ((int)all[i].size() != 0) {
      P at = all[i][rand() % (int)all[i].size()];
      b_h = at.first / N;
      b_w = at.first % N;
      b_di = at.second;
      break;
    }
  }
  parts2[G2[b_h][b_w]]--;
  parts2[G2[b_h + dy[b_di]][b_w + dx[b_di]]]--;
  G2[b_h][b_w] += pow2[b_di];
  G2[b_h + dy[b_di]][b_w + dx[b_di]] += pow2[(b_di + 2) % 4];
  parts2[G2[b_h][b_w]]++;
  parts2[G2[b_h + dy[b_di]][b_w + dx[b_di]]]++;

  add_E_sum[b_h][b_w]++;
  all_sco += sco;
  return {b_h, b_w};
}

vi seen_detect_dfs;
void detect_dfs(int c, int p, int sh, int sw) {
  int cuh = c / N, cuw = c % N;
  rep(i, 0, 4) {
    int ath = cuh + dy[i], atw = cuw + dx[i];
    if (ath == p / N && atw == p % N && p != -1)
      continue;

    if (range_out(ath, atw))
      continue;
    if (!(G2[cuh][cuw] & (1 << i)))
      continue;
    if (seen_detect_dfs[ath * N + atw] != -1)
      continue;
    seen_detect_dfs[ath * N + atw] = i;

    if (ath == sh && atw == sw)
      return;

    detect_dfs(ath * N + atw, c, sh, sw);
  }
}

void delete_edge(int sh, int sw) {
  seen_detect_dfs.assign(N * N, -1);
  detect_dfs(sh * N + sw, -1, sh, sw);

  queue<int> qu;
  qu.push(sh * N + sw);
  int b_h = -1, b_w = -1, b_di = -1;
  int sco = 0;
  qu.push(sh * N + sw);
  vi ran;
  while (!qu.empty()) {
    ran.push_back(qu.front());
    int cuh = qu.front() / N;
    int cuw = qu.front() % N;
    qu.pop();
    int ath = cuh + dy[(seen_detect_dfs[cuh * N + cuw] + 2) % 4];
    int atw = cuw + dx[(seen_detect_dfs[cuh * N + cuw] + 2) % 4];

    //辺(cuh,cuw) -> (ath,atw)を削除した時
    int at_sco = 0;

    //パーツは変更前(変更前パーツは1つ減る)
    if (parts2[G2[ath][atw]] > parts[G2[ath][atw]])
      at_sco++;

    if (parts2[G2[cuh][cuw]] > parts[G2[cuh][cuw]])
      at_sco++;

    //パーツ変更後(変更後パーツは1つ増える)
    if (ath == cuh) {
      if (parts2[G2[cuh][cuw] - 4] < parts[G2[cuh][cuw] - 4])
        at_sco++;
      if (parts2[G2[ath][atw] - 1] < parts[G2[ath][atw] - 1])
        at_sco++;
    } else {
      if (parts2[G2[cuh][cuw] - 8] < parts[G2[cuh][cuw] - 8])
        at_sco++;
      if (parts2[G2[ath][atw] - 2] < parts[G2[ath][atw] - 2])
        at_sco++;
    }

    if (sco < at_sco) {
      sco = at_sco;
      b_h = cuh;
      b_w = cuw;
      b_di = (seen_detect_dfs[b_h * N + b_w] + 2) % 4;
    }

    if (ath == sh && atw == sw)
      break;
    qu.push(ath * N + atw);
  }
  //追加した辺と同じ || 以前消したのと同じ
  if ((b_h == sh && b_w == sw) ||
      (b_h == pre_delete_h && b_w == pre_delete_w) || (sco == 0)) {
    // cout << "a" << endl;
    int at = ran[rand() % (int)ran.size()];
    b_h = at / N;
    b_w = at % N;
    b_di = (seen_detect_dfs[b_h * N + b_w] + 2) % 4;
  }
  pre_delete_h = b_h, pre_delete_w = b_w;
  delete_E_sum[b_h][b_w]++;
  //変更前パーツの個数を減らす
  parts2[G2[b_h][b_w]]--;
  parts2[G2[b_h + dy[b_di]][b_w + dx[b_di]]]--;
  //パーツを変更する
  G2[b_h][b_w] -= pow2[b_di];
  G2[b_h + dy[b_di]][b_w + dx[b_di]] -= pow2[(b_di + 2) % 4];
  //変更後パーツの個数を増やす
  parts2[G2[b_h][b_w]]++;
  parts2[G2[b_h + dy[b_di]][b_w + dx[b_di]]]++;
}

void out(vvi g) {
  cout << N << " " << 2 * N * N * N << endl;
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
  srand((int)time(0));
  init();
  int success_sum = 0;
  int Q;cin>>Q;
  rep(k,0,Q){
    parts2.assign(16,0);
    auto startClock = system_clock::now();
    G2 = generate();
    cout << "generate G2" << endl;
    out(G2);
    int roop_num = 0;
    while (!same_g()) {
      roop_num++;
      const double time =
          duration_cast<microseconds>(system_clock::now() - startClock).count() *
          1e-6;
      if (time > TIME_LIMIT)
        break;
      P add_E = add_edge();
      delete_edge(add_E.first, add_E.second);
    }
    cout << roop_num << " " << same_g()<< endl;
    // cout << endl << "roop_num: " << roop_num << endl;
    // cout << "average_sco: " << all_sco / (double)roop_num << endl;
    // out_parts(parts);
    // out_parts(parts2);
    // out(G2);
    // cout << "add_E_sum" << endl;
    // rep(i, 0, N) {
    //   rep(j, 0, N) cout << add_E_sum[i][j] << " ";
    //   cout << endl;
    // }
    // cout << "delete_E_sum" << endl;
    // rep(i, 0, N) {
    //   rep(j, 0, N) cout << delete_E_sum[i][j] << " ";
    //   cout << endl;
    // }
    if(same_g())success_sum++;
    // cout << same_g() << endl;
  }
  cout << Q << " " << success_sum << endl;
}

//入力 N ([6:10])
//出力 char chg[N-1][N-1]