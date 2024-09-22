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
typedef vector<vll> vvll;
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
#define rep(i, l, n) for (ll i = (ll)(l); i < (ll)(n); i++)
#define repd(i, n, l) for (ll i = (ll)(n); i > (ll)(l); i--)
#define Copy(from, to) copy(from.begin(), from.end(), to.begin())
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
#define Unique(a) sort(a.begin(), a.end())
#define dPQ priority_queue<ll, vll, greater<ll>>
#define PQ priority_queue<ll, vll>
template <class T, class S> inline bool chmax(T &a, const S &b) {
  return (a < b ? a = b, 1 : 0);
}
template <class T, class S> inline bool chmin(T &a, const S &b) {
  return (a > b ? a = b, 1 : 0);
}
#define INF INT_MAX
#define INFF (9223372036854775800)
#define TIME_LIMIT (5.96)
#define def (10101010)
// #define MOD (1000000007)
#define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision(10) << << endl;
// 偏角ソートはlong ddouble!
// auto ite = s.lower_bound("B");
// cout << __builtin_popcount(N) << endl;

// 下,右,上,左
vll dy = {1, 0, -1, 0};
vll dx = {0, 1, 0, -1};
ll N, W, K, C;
vll a, b, c, d;
vvll S;

void dist_break_rocks(vll houses, ll y, ll x, ll &sum, ll &best_dist) {
  if ((int)houses.size() == 0 || best_dist < sum)
    return;

  vll dir_sum(4, 0);
  vll max_dist_dir(4, (ll)200);
  for (ll house : houses) {
    ll dy = c[house] - y;
    ll dx = d[house] - x;
    if (dy > 0)
      chmin(max_dist_dir[0], abs(c[house] - y));
    if (dy < 0)
      chmin(max_dist_dir[2], abs(c[house] - y));
    if (dx > 0)
      chmin(max_dist_dir[1], abs(d[house] - x));
    if (dx < 0)
      chmin(max_dist_dir[3], abs(d[house] - x));
  }

  for (ll house : houses) {
    ll dy = c[house] - y;
    ll dx = d[house] - x;
    if (dy > 0)
      dir_sum[0] += max_dist_dir[0];
    if (dy < 0)
      dir_sum[2] += max_dist_dir[2];
    if (dx > 0)
      dir_sum[1] += max_dist_dir[1];
    if (dx < 0)
      dir_sum[3] += max_dist_dir[3];
  }

  ll max_dir = 0;
  rep(i, 0, 4) {
    dir_sum[i] -= max_dist_dir[i];
    if (dir_sum[max_dir] < dir_sum[i])
      max_dir = i;
  }
  if (dir_sum[max_dir] == 0) {
    rep(i, 0, 4) {
      if (max_dist_dir[i] != 200) {
        dir_sum[i] = max_dist_dir[i];
        max_dir = i;
      }
    }
  }

  //( {y, x}, {toy, tox} ]
  ll toy = y + dy[max_dir] * max_dist_dir[max_dir];
  ll tox = x + dx[max_dir] * max_dist_dir[max_dir];
  sum += abs(toy - y) + abs(tox - x);

  vll l_houses, r_houses;
  for (ll house : houses) {
    ll dy = c[house] - y;
    ll dx = d[house] - x;
    if (c[house] == toy && d[house] == tox)
      continue;
    if ((dy > 0 && max_dir == 0) || (dy < 0 && max_dir == 2) ||
        (dx > 0 && max_dir == 1) || (dx < 0 && max_dir == 3))
      l_houses.push_back(house);
    else
      r_houses.push_back(house);
  }

  dist_break_rocks(l_houses, toy, tox, sum, best_dist);
  dist_break_rocks(r_houses, y, x, sum, best_dist);
}

void select_break_rocks(vP &breack_rocks, vll houses, ll y, ll x) {
  if ((int)houses.size() == 0)
    return;
  // cout << endl << "from : " << y << " " << x << endl;
  // cout << "分割前 : ";
  // for (ll house : houses)
  //   cout << house << " ";
  // cout << endl;
  vll dir_sum(4, 0);
  vll max_dist_dir(4, (ll)200);
  for (ll house : houses) {
    ll dy = c[house] - y;
    ll dx = d[house] - x;
    if (dy > 0)
      chmin(max_dist_dir[0], abs(c[house] - y));
    if (dy < 0)
      chmin(max_dist_dir[2], abs(c[house] - y));
    if (dx > 0)
      chmin(max_dist_dir[1], abs(d[house] - x));
    if (dx < 0)
      chmin(max_dist_dir[3], abs(d[house] - x));
  }

  for (ll house : houses) {
    ll dy = c[house] - y;
    ll dx = d[house] - x;
    if (dy > 0)
      dir_sum[0] += max_dist_dir[0];
    if (dy < 0)
      dir_sum[2] += max_dist_dir[2];
    if (dx > 0)
      dir_sum[1] += max_dist_dir[1];
    if (dx < 0)
      dir_sum[3] += max_dist_dir[3];
  }

  ll max_dir = 0;
  rep(i, 0, 4) {
    dir_sum[i] -= max_dist_dir[i];
    if (dir_sum[max_dir] < dir_sum[i])
      max_dir = i;
  }
  if (dir_sum[max_dir] == 0) {
    rep(i, 0, 4) {
      if (max_dist_dir[i] != 200) {
        dir_sum[i] = max_dist_dir[i];
        max_dir = i;
      }
    }
  }

  //( {y, x}, {toy, tox} ]
  ll toy = y;
  ll tox = x;
  rep(_i, 0, max_dist_dir[max_dir]) {
    toy += dy[max_dir];
    tox += dx[max_dir];
    breack_rocks.push_back({toy, tox});
  }
  vll l_houses, r_houses;
  for (ll house : houses) {
    ll dy = c[house] - y;
    ll dx = d[house] - x;
    // cout << house << " {" << c[house] << ", " << d[house] << "} dy : " << dy
    //      << " dx : " << dx << endl;
    if (c[house] == toy && d[house] == tox)
      continue;
    if ((dy > 0 && max_dir == 0) || (dy < 0 && max_dir == 2) ||
        (dx > 0 && max_dir == 1) || (dx < 0 && max_dir == 3))
      l_houses.push_back(house);
    else
      r_houses.push_back(house);
  }

  // cout << "dir : " << max_dir << " " << dir_sum[max_dir] << endl;
  // cout << "分割後l : ";
  // for (ll house : l_houses)
  //   cout << house << " ";
  // cout << endl;
  // cout << "分割後r : ";
  // for (ll house : r_houses)
  //   cout << house << " ";
  // cout << endl;

  select_break_rocks(breack_rocks, l_houses, toy, tox);
  select_break_rocks(breack_rocks, r_houses, y, x);
}

void local_solve_mining(vP &break_rocks) {
  ll set_power = 50;
  ll ans = 0;
  ll all_mine_sum = 0;
  for (auto rock : break_rocks) {
    ll mine_sum = 0;
    while (1) {
      cout << rock.first << " " << rock.second << " " << set_power << endl;
      cout.flush();
      ans += set_power + C;
      all_mine_sum++;
      ll r;
      S[rock.first][rock.second] -= set_power;
      if (S[rock.first][rock.second] <= 0) {
        r = 1;
      } else
        r = 0;

      if (r == 0) {
        mine_sum++;
      }
      if (r == 1) {
        if (mine_sum == 0)
          set_power /= 3;
        else
          set_power *= (mine_sum);
        chmax(set_power, (ll)10);
        chmin(set_power, (ll)500);
        break;
      }
    }
  }
  // cerr << ans << endl;
  cerr << all_mine_sum << endl;
}

void solve_mining(vP &break_rocks) {
  ll set_power = 50;
  for (auto rock : break_rocks) {
    ll mine_sum = 0;
    while (1) {
      cout << rock.first << " " << rock.second << " " << set_power << endl;
      cout.flush();

      ll r;
      cin >> r;

      if (r == -1 || r == 2)
        return;

      if (r == 0) {
        mine_sum++;
      }
      if (r == 1) {
        if (mine_sum == 0)
          set_power /= 3;
        else
          set_power *= (mine_sum);
        chmax(set_power, (ll)10);
        chmin(set_power, (ll)500);
        break;
      }
    }
  }
}

int main() {
  cin >> N >> W >> K >> C;

  S.resize(N);
  rep(i, 0, N) {
    S[i].resize(N);
    rep(j, 0, N) cin >> S[i][j];
  }

  a.resize(W);
  b.resize(W);
  c.resize(K);
  d.resize(K);
  rep(i, 0, W) cin >> a[i] >> b[i];
  rep(i, 0, K) cin >> c[i] >> d[i];

  ll best_dist = 1e10;
  ll best_house_W = 0;

  rep(bit, 0, (ll)(pow(W, K))) {
    // house_W[i]=j := 家iは水源jとつながる
    // K_houses[i]=[] := 水源iは家[]とつながる
    vll house_W(K);
    vvll K_houses(W);
    ll at = bit;
    rep(i, 0, K) {
      house_W[i] = at % W;
      K_houses[house_W[i]].push_back(i);
      at /= W;
    }
    ll at_dist = 0;
    rep(i, 0, W) {
      dist_break_rocks(K_houses[i], a[i], b[i], at_dist, best_dist);
    }
    if (chmin(best_dist, at_dist))
      best_house_W = bit;
  }

  vP break_rocks;
  // house_W[i]=j := 家iは水源jとつながる
  // K_houses[i]=[] := 水源iは家[]とつながる
  vll house_W(K);
  vvll K_houses(W);
  ll bit = best_house_W;
  rep(i, 0, K) {
    house_W[i] = bit % W;
    K_houses[house_W[i]].push_back(i);
    bit /= W;
  }
  rep(i, 0, W) {
    break_rocks.push_back({a[i], b[i]});
    // cout << "i : " << i << endl;
    // for (ll house : K_houses[i])
    //   cout << house << " ";
    // cout << endl << endl;
    select_break_rocks(break_rocks, K_houses[i], a[i], b[i]);
  }

  local_solve_mining(break_rocks);
  // solve_mining(break_rocks);
}

// chmin(set_power, (ll)500);を追加して手元1.15倍
