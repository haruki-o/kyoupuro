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
#define all_TIME_LIMIT (2.98)
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

vvi _G;
vvi _made_G;
vi _route;
int vh, vw;
string ans = "";
string _ans = "";
// vh,vwがdi([0,3])にいく
vi dc = {'L', 'U', 'R', 'D'};

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

void init2() {
  _G.resize(N);
  rep(i, 0, N) _G[i].resize(N);

  _made_G.resize(N);
  rep(i, 0, N) _made_G[i].resize(N);
}
void all_reset() {

  _G = G;
  rep(i, 0, N) _made_G[i].assign(N, 0);

  rep(i, 0, N) {
    rep(j, 0, N) {
      if (G[i][j] == 0) {
        vh = i;
        vw = j;
      }
    }
  }

  _ans = "";
}

P best_parts(int h, int w, int si) {
  int _h, _w;
  int mi_dist = INF;
  int fi_parts = G2[h][w];
  if (si == 1) {
    if (h == 0)
      fi_parts = G2[1][w];
    if (h == N - 1)
      fi_parts = G2[N - 2][w];
  }
  if (si == 2) {
    if (h == 0)
      fi_parts = G2[0][w - 1];
    if (h == N - 1)
      fi_parts = G2[N - 1][w - 1];
  }
  if (3 <= si && si <= 4)
    fi_parts = G2[h - (si - 3)][w + (4 - si)];
  if (si == 6)
    fi_parts = G2[N - 3][N - 1];
  if (si == 7)
    fi_parts = G2[N - 3][N - 2];
  if (si == 8)
    fi_parts = G2[N - 1][N - 3];
  if (si == 9)
    fi_parts = G2[N - 2][N - 3];
  if (si == 10) {
    if (h == 0)
      fi_parts = G2[1][w - 3];
    if (h == N - 1)
      fi_parts = G2[N - 2][w - 3];
  }
  if (si == 11) {
    if (h == 0)
      fi_parts = G2[0][w - 4];
    if (h == N - 1)
      fi_parts = G2[N - 1][w - 4];
  }
  rep(i, 0, N) {
    rep(j, 0, N) {
      if (fi_parts != _G[i][j])
        continue;
      if (_made_G[i][j] == 1)
        continue;
      int dist = abs(h - i) + abs(w - j);
      dist += dist * 2;
      dist += abs(vh - i) + abs(vw - j);
      if (dist < mi_dist) {
        mi_dist = dist;
        _h = i;
        _w = j;
      }
    }
  }

  return {_h, _w};
}

void va_move(int h, int w, int outh, int outw) {
  vvP dp(N, vP(N, {-1, -1}));
  dp[vh][vw].second = 0;
  queue<P> qu;
  qu.push({vh, vw});
  while (!qu.empty()) {
    auto [cuh, cuw] = qu.front();
    qu.pop();
    rep(i, 0, 4) {
      int ath = cuh + dy[i];
      int atw = cuw + dx[i];
      if (range_out(ath, atw))
        continue;
      if (dp[ath][atw].second != -1)
        continue;
      if (_made_G[ath][atw] == 1)
        continue;
      if (ath == outh && atw == outw)
        continue;
      dp[ath][atw] = {(i + 2) % 4, dp[cuh][cuw].second + 1};
      qu.push({ath, atw});
    }
    if (dp[h][w].second != -1)
      break;
  }
  while (!qu.empty())
    qu.pop();

  //(h,w) -> (vh,vw)の復元
  qu.push({h, w});
  string at_ans = "";
  while (!qu.empty()) {
    auto [cuh, cuw] = qu.front();
    qu.pop();
    int _di = dp[cuh][cuw].first;
    if (cuh == vh && cuw == vw)
      break;
    at_ans += dc[_di];
    qu.push({cuh + dy[_di], cuw + dx[_di]});
  }
  //_Gの更新も行う
  repd(i, (int)at_ans.size() - 1, -1) {
    rep(j, 0, 4) if (at_ans[i] == dc[j]) {
      _ans += dc[(j + 2) % 4];
      int nh = vh + dy[(j + 2) % 4];
      int nw = vw + dx[(j + 2) % 4];
      swap(_G[vh][vw], _G[nh][nw]);
      vh = nh;
      vw = nw;
    }
  }
}

void move() {
  rep(i, 1, (int)_route.size()) {
    int ath = _route[i] / N;
    int atw = _route[i] % N;
    int preh = _route[i - 1] / N;
    int prew = _route[i - 1] % N;
    va_move(ath, atw, preh, prew);
    if (preh == ath) {
      if (atw < prew) {
        _ans += dc[2];
      } else
        _ans += dc[0];
    } else {
      if (ath < preh) {
        _ans += dc[3];
      } else
        _ans += dc[1];
    }
    swap(_G[ath][atw], _G[preh][prew]);
    vh = preh;
    vw = prew;
  }
}

void find_best_route(int _h, int _w, int h, int w) {
  vvP dp(N, vP(N, {-1, INF}));
  _route.clear();

  priority_queue<pair<int, T>, vector<pair<int, T>>, greater<pair<int, T>>> pq;
  pq.push({0, {_h, _w, -1}});
  while (!pq.empty()) {
    pair<int, T> cu = pq.top();
    auto [cuh, cuw, di] = cu.second;
    pq.pop();
    if (dp[cuh][cuw].second != INF)
      continue;
    dp[cuh][cuw].second = cu.first;
    dp[cuh][cuw].first = (di + 2) % 4;
    // int si = _h < h ? 3 : 1;
    int si = 0;
    rep(i, si , si + 4) {
      int ati = i %4;
      int ath = cuh + dy[ati];
      int atw = cuw + dx[ati];
      if (range_out(ath, atw))
        continue;
      if (dp[ath][atw].second != INF)
        continue;
      if (_made_G[ath][atw] == 1)
        continue;
      if (ati % 2 == dp[cuh][cuw].first % 2) {
        if (!range_out(ath + dy[ati], atw + dx[ati])) {
          if (_made_G[ath + dy[ati]][atw + dx[ati]] == 1)
            pq.push({dp[cuh][cuw].second + 5, {ath, atw, ati}});
          else
            pq.push({dp[cuh][cuw].second + 4, {ath, atw, ati}});
        } else
          pq.push({dp[cuh][cuw].second + 4, {ath, atw, ati}});
      } else {
        if (!range_out(ath + dy[ati], atw + dx[ati])) {
          if (_made_G[ath + dy[ati]][atw + dx[ati]] == 1)
            pq.push({dp[cuh][cuw].second + 5, {ath, atw, ati}});
          else
            pq.push({dp[cuh][cuw].second + 2, {ath, atw, ati}});
        } else
          pq.push({dp[cuh][cuw].second + 2, {ath, atw, ati}});
      }
    }
    if (dp[h][w].second != INF)
      break;
  }

  queue<P> qu;
  //(h,w) -> (_h,_w)の復元
  qu.push({h, w});
  while (!qu.empty()) {
    auto [cuh, cuw] = qu.front();
    qu.pop();
    _route.push_back(cuh * N + cuw);
    if (cuh == _h && cuw == _w)
      break;
    int _di = dp[cuh][cuw].first;
    qu.push({cuh + dy[_di], cuw + dx[_di]});
  }
  reverse(_route.begin(), _route.end());
}

void set_parts(int h, int w, int si) {
  //(_h, _w) -> (h, w)
  auto [_h, _w] = best_parts(h, w, si);
  if (_h == h && _w == w) {
    _made_G[h][w] = 1;
    return;
  }
  find_best_route(_h, _w, h, w);
  move();
  _made_G[h][w] = 1;
}

bool set_3_3() {
  set_parts(N - 3, N - 3, 0);

  //ここから3 * 3の1行目(2,3)の完成
  set_parts(N - 3, N - 2, 6);
  va_move(N - 2, N - 1, -1, -1);
  // 1 3 2
  // * * vの時
  if (_G[N - 3][N - 1] == G2[N - 3][N - 2]) {
    // 1 v 3
    //* * 2
    //* * *にする
    _made_G[N - 3][N - 1] = 1;
    swap(_G[N - 2][N - 1], _G[N - 3][N - 1]);
    swap(_G[N - 3][N - 1], _G[N - 3][N - 2]);
    _made_G[N - 2][N - 1] = 1;
    _made_G[N - 3][N - 2] = 0;
    vh = N - 3;
    vw = N - 2;
    _ans += "UL";

    // 1 * 3
    //* * v
    //* * 2にする
    va_move(N - 1, N - 1, -1, -1);
    swap(_G[N - 1][N - 1], _G[N - 2][N - 1]);
    _made_G[N - 1][N - 1] = 1;
    _made_G[N - 2][N - 1] = 0;
    vh = N - 2;
    vw = N - 1;
    _ans += "U";

    // 1 3 v
    //* * *
    //* * 2にする
    va_move(N - 3, N - 2, -1, -1);
    swap(_G[N - 3][N - 2], _G[N - 3][N - 1]);
    _made_G[N - 3][N - 2] = 1;
    _made_G[N - 3][N - 1] = 0;
    vh = N - 3;
    vw = N - 1;
    _ans += "R";

    // 1 3 *
    //* * 2
    //* * vにする
    va_move(N - 2, N - 1, -1, -1);
    swap(_G[N - 2][N - 1], _G[N - 1][N - 1]);
    _made_G[N - 2][N - 1] = 1;
    _made_G[N - 1][N - 1] = 0;
    vh = N - 1;
    vw = N - 1;
    _ans += "D";

    // 1 3 *
    //* 2 v
    //* * *にする
    va_move(N - 2, N - 2, -1, -1);
    swap(_G[N - 2][N - 2], _G[N - 2][N - 1]);
    _made_G[N - 2][N - 2] = 1;
    _made_G[N - 2][N - 1] = 0;
    vh = N - 2;
    vw = N - 1;
    _ans += "R";

  } else {
    set_parts(N - 2, N - 2, 7);
  }
  va_move(N - 3, N - 1, -1, -1);
  swap(_G[N - 3][N - 1], _G[N - 3][N - 2]);
  swap(_G[N - 3][N - 2], _G[N - 2][N - 2]);
  _made_G[N - 3][N - 1] = 1;
  _made_G[N - 2][N - 2] = 0;
  vh = N - 2;
  vw = N - 2;
  _ans += "LD";

  //ここから3 * 3の1列目(4,7)の完成
  set_parts(N - 2, N - 3, 8);
  va_move(N - 1, N - 2, -1, -1);
  // 7 * *
  // 4 v * のとき
  if (_G[N - 1][N - 3] == G2[N - 2][N - 3]) {
    // v * *
    // 7 4 *にする
    swap(_G[N - 1][N - 2], _G[N - 1][N - 3]);
    swap(_G[N - 1][N - 3], _G[N - 2][N - 3]);
    _made_G[N - 1][N - 2] = 1;
    _made_G[N - 1][N - 3] = 1;
    _made_G[N - 2][N - 3] = 0;
    vh = N - 2;
    vw = N - 3;
    _ans += "LU";

    //* * *
    // 7 v 4にする
    va_move(N - 1, N - 1, -1, -1);
    swap(_G[N - 1][N - 1], _G[N - 1][N - 2]);
    _made_G[N - 1][N - 1] = 1;
    _made_G[N - 1][N - 2] = 0;
    vh = N - 1;
    vw = N - 2;
    _ans += "L";

    // 7 * *
    // v * 4にする
    va_move(N - 2, N - 3, -1, -1);
    swap(_G[N - 2][N - 3], _G[N - 1][N - 3]);
    _made_G[N - 2][N - 3] = 1;
    _made_G[N - 1][N - 3] = 0;
    _made_G[N - 1][N - 1] = 0;
    vh = N - 1;
    vw = N - 3;
    _ans += "D";

    // 7 * *
    //* 4 vにする
    va_move(N - 1, N - 2, -1, -1);
    swap(_G[N - 1][N - 2], _G[N - 1][N - 1]);
    _made_G[N - 1][N - 2] = 1;
    _made_G[N - 1][N - 1] = 0;
    vh = N - 1;
    vw = N - 1;
    _ans += "R";

    // 7 4 *
    //* v *
    va_move(N - 2, N - 2, -1, -1);
    swap(_G[N - 2][N - 2], _G[N - 1][N - 2]);
    _made_G[N - 2][N - 2] = 1;
    _made_G[N - 1][N - 2] = 0;
    vh = N - 1;
    vw = N - 2;
    _ans += "D";

  } else {
    // 7 4 *
    // v * *を目指す
    set_parts(N - 2, N - 2, 9);
  }
  va_move(N - 1, N - 3, -1, -1);

  // 4 v *
  // 7 * *
  swap(_G[N - 1][N - 3], _G[N - 2][N - 3]);
  swap(_G[N - 2][N - 3], _G[N - 2][N - 2]);
  _made_G[N - 1][N - 3] = 1;
  _made_G[N - 2][N - 2] = 0;
  _ans += "UR";

  vh = N - 2;
  vw = N - 2;
  //右下2 * 2の処理
  // v *
  // 5 *
  if (_G[N - 1][N - 2] == G2[N - 2][N - 2]) {
    swap(_G[N - 2][N - 2], _G[N - 1][N - 2]);
    swap(_G[N - 1][N - 2], _G[N - 1][N - 1]);
    _ans += "DR";
    if (_G[N - 1][N - 2] == G2[N - 1][N - 2] &&
        _G[N - 2][N - 1] == G2[N - 2][N - 1])
      return true;
  }
  // v 5
  //* *
  if (_G[N - 2][N - 1] == G2[N - 2][N - 2]) {
    swap(_G[N - 2][N - 2], _G[N - 2][N - 1]);
    swap(_G[N - 2][N - 1], _G[N - 1][N - 1]);
    _ans += "RD";
    if (_G[N - 1][N - 2] == G2[N - 1][N - 2] &&
        _G[N - 2][N - 1] == G2[N - 2][N - 1])
      return true;
  }

  return false;
}

bool set_2_2() {
  va_move(N - 2, N - 2, -1, -1);
  //右下2 * 2の処理
  // v *
  // 5 *
  if (_G[N - 1][N - 2] == G2[N - 2][N - 2]) {
    swap(_G[N - 2][N - 2], _G[N - 1][N - 2]);
    swap(_G[N - 1][N - 2], _G[N - 1][N - 1]);
    _ans += "DR";
    if (_G[N - 1][N - 2] == G2[N - 1][N - 2] &&
        _G[N - 2][N - 1] == G2[N - 2][N - 1])
      return true;
  }
  // v 5
  //* *
  if (_G[N - 2][N - 1] == G2[N - 2][N - 2]) {
    swap(_G[N - 2][N - 2], _G[N - 2][N - 1]);
    swap(_G[N - 2][N - 1], _G[N - 1][N - 1]);
    _ans += "RD";
    if (_G[N - 1][N - 2] == G2[N - 1][N - 2] &&
        _G[N - 2][N - 1] == G2[N - 2][N - 1])
      return true;
  }

  return false;
}

int main() {
  auto allStartClock = system_clock::now();
  srand((int)time(0));
  init();
  init2();

  int success_sum = 0;
  int complete_sum = 0;
  int same_g_sum = 0;
  int min_index = -1;
  int roop_num = 0;

  int Q = 100000000;
  // cin >> Q;
  rep(k, 0, Q) {
    all_reset();
    const double allTime =
        duration_cast<microseconds>(system_clock::now() - allStartClock)
            .count() *
        1e-6;
    if (allTime > all_TIME_LIMIT) {
      roop_num = k;
      break;
    }
    parts2.assign(16, 0);
    auto startClock = system_clock::now();
    G2 = generate();
    while (!same_g()) {
      const double time =
          duration_cast<microseconds>(system_clock::now() - startClock)
              .count() *
          1e-6;
      if (time > TIME_LIMIT)
        break;
      P add_E = add_edge();
      delete_edge(add_E.first, add_E.second);
    }
    // cout << roop_num << " " << same_g() << endl;

    if (same_g())
      success_sum++;
    else
      continue;

    //左[0, N-3]列
    int break_ju = 0;
    rep(j, 0, N - 2) {
      // if (j == 2 && N <= 7)
      //   if ((int)_ans.size() > 150 + 50 * (N - 6)) {
      //     break_ju = 1;
      //     break;
      //   }
      if (j % 2 == 0) {
        rep(i, 0, N - 2) set_parts(i, j, 0);
        if (vw == j)
          va_move(vh, j + 1, -1, -1);
        // 2ますが同じなら飛ばせる
        if (_G[N - 2][j] == G2[N - 2][j] && _G[N - 1][j] == G2[N - 1][j]) {
          _made_G[N - 2][j] = 1;
          _made_G[N - 1][j] = 1;
          continue;
        }

        set_parts(N - 1, j, 1);
        if (vh == N - 2 && vw == j)
          va_move(N - 2, j + 1, -1, -1);
        if (_G[N - 2][j] == G2[N - 1][j]) {
          break_ju = 1;
          break;
        }
        set_parts(N - 1, j + 1, 2);
        va_move(N - 2, j, -1, -1);
        swap(_G[N - 2][j], _G[N - 1][j]);
        swap(_G[N - 1][j], _G[N - 1][j + 1]);
        vh = N - 1, vw = j + 1;
        _made_G[N - 1][j + 1] = 0;
        _made_G[N - 2][j] = 1;
        _ans += "DR";

      } else {
        repd(i, N - 1, 1) set_parts(i, j, 0);
        if (vw == j)
          va_move(vh, j + 1, -1, -1);
        // 2ますが同じなら飛ばせる
        if (_G[1][j] == G2[1][j] && _G[0][j] == G2[0][j]) {
          _made_G[1][j] = 1;
          _made_G[0][j] = 1;
          continue;
        }

        set_parts(0, j, 1);
        if (vh == 1 && vw == j)
          va_move(1, j + 1, -1, -1);
        if (_G[1][j] == G2[0][j]) {
          break_ju = 1;
          break;
        }
        set_parts(0, j + 1, 2);
        va_move(1, j, -1, -1);
        swap(_G[1][j], _G[0][j]);
        swap(_G[0][j], _G[0][j + 1]);
        vh = 0, vw = j + 1;
        _made_G[0][j + 1] = 0;
        _made_G[1][j] = 1;
        _ans += "UR";
      }
      // cout << k << " " << j << " " << (int)_ans.size() << endl;
    }
    if (break_ju == 1)
      continue;

    //右2列
    rep(i, 0, N - 2) {
      // 2個全部同じなら飛ばせる
      if (vh == i)
        va_move(i + 1, vw, -1, -1);
      int ju = 1;
      rep(j, 0, 2) if (_G[i][N - 2 + j] != G2[i][N - 2 + j]) ju = 0;
      if (ju) {
        rep(j, 0, 2) _made_G[i][N - 2 + j] = 1;
        continue;
      }

      //縦2つ配置
      set_parts(i, N - 2, 3);
      if (vh == i && vw == N - 1)
        va_move(i + 1, N - 1, -1, -1);
      if (_G[i][N - 1] == G2[i][N - 2]) {
        break_ju = 1;
        break;
      }
      set_parts(i + 1, N - 2, 4);
      //ずらす
      va_move(i, N - 1, -1, -1);
      swap(_G[i][N - 1], _G[i][N - 2]);
      swap(_G[i][N - 2], _G[i + 1][N - 2]);
      vh = i + 1;
      vw = N - 2;
      _made_G[i][N - 1] = 1;
      _made_G[i + 1][N - 2] = 0;
      _ans += "LD";
    }
    if (break_ju == 1)
      continue;

    //右下 [2][2]
    int ju = set_2_2();
    if (ju)
      complete_sum++;
    if (ans == "" || (int)_ans.size() < (int)ans.size()) {
      if (ju) {
        ans = _ans;
        min_index = k;
      }
    }
  }
  cerr << fixed << setprecision(3) << 2 - (double)ans.size() / initT << endl;
  // cerr << success_sum << " " << complete_sum << endl;
  // cerr << "min_index " << min_index << endl;
  // cerr << "roop_num " << roop_num << endl;
  // cerr << (int)ans.size() << endl;
  cout << ans << endl;
}

//入力 N ([6:10])
//出力 char chg[N-1][N-1]

// TimeLimit 0.002
//  0.0000460000 out(G2)なし (1000,832)
//  0.0027880000 out(G2)あり (1000,14)

// TimeLimte 0.01
// out()なし (1000,797)
// out()あり (1000,780)

// for f in in/*.txt; do ./a.out < $f > output.txt; done

// find_best_route()の改善(壁に当たるときの + 6) 143.184 -> 143.406
// find_best_route()の改善(壁に当たるときの + 5) 143.444

//上記 + best_parts()の改善(マンハッタン距離 * 4 ->  * 1) 143.997, 143.587, 143.452
//上記 + best_parts()の改善(マンハッタン距離 * 4 ->  * 2) 143.513, 143.98, 143.807
//上記 + best_parts()の改善(マンハッタン距離 * 4 ->  * 3) 143.745, 143.492, 143.496


// find_best_route()の改善(rep(i,0,4)を_h,hの大小で1,3スタート決める) 143.374,143.284



