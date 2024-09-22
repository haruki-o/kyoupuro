#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef pair<ll, ll> pll;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vvi> vvvi;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvll;
typedef vector<pll> vP;
#define rep(i, l, n) for (ll i = (ll)(l); i < (ll)(n); i++)
#define repd(i, n, l) for (ll i = (ll)(n); i > (ll)(l); i--)
#define Copy(from, to) copy(from.begin(), from.end(), to.begin())
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
#define INF INT_MAX
#define INFF (9223372036854775800)
#define TIME_LIMIT (1.9)
#define def (2010101)
// #define MOD (1000000007)
#define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision(10) << << endl;
// 偏角ソートはlong ddouble!
// auto ite = s.lower_bound("B");
// cout << __builtin_popcount(N) << endl;

ll N, M, K;
vll x, y, u, v, w, a, b;
vvp g;
map<pair<ll, ll>, ll> ma;

double stat_time = 0.0;
double stat_loop = 0;

struct UnionFind {
  // 根のとき-1,iの親はpar[i]
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

double ela_times(system_clock::time_point &clock) {
  return duration_cast<microseconds>(system_clock::now() - clock).count() *
         1e-6;
}

ll dist(pair<ll, ll> p1, pair<ll, ll> p2) {
  return (p1.first - p2.first) * (p1.first - p2.first) +
         (p1.second - p2.second) * (p1.second - p2.second);
}

// 状態の初期化
void init(vll &P, vll &B) {
  vector<pair<ll, ll>> house(K);
  rep(i, 0, K)
      house[i] = {(abs(a[i]) + abs(b[i])) * (1 + 0.1 * (rand() % 5)), i};
  gSort(house);
  vi after(K, 0);
  vll at = {50, 100, 300, 600, 1000, 5000};
  for (ll num : at) {
    rep(_i, 0, K) {
      ll i = house[_i].second;
      ll idx = 0;
      ll mi = 1e15;
      int flag = 0;
      rep(j, 1, P.size()) {
        // 電源と家の距離
        ll nec = sqrt(dist({x[j], y[j]}, {a[i], b[i]})) + 1;
        // 既に被覆されている
        if (nec < P[j])
          flag = 1;
        if (num <= nec)
          continue;
        if (nec * nec - P[j] * P[j] < mi) {
          mi = nec * nec - P[j] * P[j];
          idx = j;
        }
      }
      if (!flag && mi != 1e15) {
        ll nec = sqrt(dist({x[idx], y[idx]}, {a[i], b[i]})) + 1;
        P[idx] = nec;
      }
    }
  }

  // {cost, from , to}
  priority_queue<tuple<ll, ll, ll>, vector<tuple<ll, ll, ll>>,
                 greater<tuple<ll, ll, ll>>>
      pq;
  pq.push({0, -1, 0});
  vll dp(N, INFF);
  vll dir(N, 0);
  while (!pq.empty()) {
    auto cu = pq.top();
    pq.pop();
    if (dp[get<2>(cu)] < get<0>(cu))
      continue;
    dp[get<2>(cu)] = get<0>(cu);
    dir[get<2>(cu)] = get<1>(cu);
    for (auto to : g[get<2>(cu)]) {
      if (dp[to.first] != INFF)
        continue;
      pq.push({dp[get<2>(cu)] + to.second, get<2>(cu), to.first});
    }
  }
  queue<ll> qu;
  vi seen(N, 0);
  rep(i, 0, N) {
    if (P[i] != 0) {
      qu.push(i);
      seen[i] = 1;
    }
  }

  while (!qu.empty()) {
    ll cu = qu.front();
    qu.pop();
    ll to = dir[cu];
    if (B[ma[{min(cu, to), max(cu, to)}]] == 1)
      continue;
    B[ma[{min(cu, to), max(cu, to)}]] = 1;
    if (seen[to] == 1 || to == -1)
      continue;
    seen[to] = 1;
    qu.push(to);
  }
}

// 状態遷移
void modify(vll &P, vll &B) {
  vector<pair<ll, ll>> house(K);
  rep(i, 0, K) house[i] = {abs(a[i]) + abs(b[i]), i};
  gSort(house);
  rep(_i, 0, K) {
    ll i = house[_i].second;
    ll idx = 0;
    ll mi = 1e15;
    int flag = 0;
    rep(j, 1, P.size()) {
      // 電源と家の距離
      ll nec = sqrt(dist({x[j], y[j]}, {a[i], b[i]})) + 1;
      // 既に被覆されている
      if (nec < P[j])
        flag = 1;
      // 5000以上の時ダメ
      if (5000 <= nec)
        continue;
      if (nec * nec - P[j] * P[j] < mi) {
        mi = nec * nec - P[j] * P[j];
        idx = j;
      }
    }
    if (!flag) {
      ll nec = sqrt(dist({x[idx], y[idx]}, {a[i], b[i]})) + 1;
      P[idx] = nec;
    }
  }
}

void bfs(vll &seen, vll &B) {
  UnionFind uf(N);
  rep(i, 0, M) if (B[i] == 1) uf.unite(u[i], v[i]);
  ll root = uf.root(0);
  rep(i, 0, N) if (root == uf.root(i)) seen[i] = 1;
}

// 状態のスコア計算
ll calc_score(vll &P, vll &B) {
  ll ret = 0;
  vll seen(N, 0);
  bfs(seen, B);
  vll ava(K, 0);
  rep(i, 0, N) {
    if (seen[i] == 0)
      continue;
    rep(j, 0, K) {
      if (dist({x[i], y[i]}, {a[j], b[j]}) <= P[i] * P[i])
        ava[j] = 1;
    }
  }

  ll sum = 0;
  rep(i, 0, K) sum += ava[i];

  if (sum != K)
    return (ll)1e6 * (sum + 1) / K;
  else {
    rep(i, 0, N) ret += P[i] * P[i];
    rep(i, 0, M) {
      if (B[i] == 1)
        ret += w[i];
    }
  }
  return 1e6 * (1 + (1e8 / (ret + 1e7)));
}

// 山登り
void sa(vll &P, vll &B) {
  auto start_time = system_clock::now();
  while (1) {
    stat_loop++;
    auto now_time = system_clock::now();
    if (TIME_LIMIT < ela_times(start_time))
      break;

    vll new_P(N), new_B(M);
    init(new_P, new_B);
    int new_score = calc_score(new_P, new_B);
    int pre_score = calc_score(P, B);

    if (pre_score < new_score) {
      P = new_P;
      B = new_B;
    }
  }
}

void ans_output(vll &P, vll &B) {
  rep(i, 0, N) cout << P[i] << " ";
  cout << endl;
  rep(i, 0, M) cout << B[i] << " ";
  cout << endl;
}

int main() {
  srand((unsigned int)time(NULL));

  cin >> N >> M >> K;
  x.resize(N), y.resize(N);
  u.resize(M), v.resize(M), w.resize(M);
  a.resize(K), b.resize(K);
  g.resize(N);
  rep(i, 0, N) cin >> x[i] >> y[i];
  rep(i, 0, M) cin >> u[i] >> v[i] >> w[i];
  rep(i, 0, M) u[i]--, v[i]--;
  rep(i, 0, M) {
    g[u[i]].push_back({v[i], w[i]});
    g[v[i]].push_back({u[i], w[i]});
    ma[{min(u[i], v[i]), max(u[i], v[i])}] = i;
  }
  rep(i, 0, K) cin >> a[i] >> b[i];

  vll P(N, 0), B(M, 0);
  init(P, B);
  sa(P, B);
  ans_output(P, B);

  cerr << stat_loop << endl;
  cerr << calc_score(P, B) << endl;
}