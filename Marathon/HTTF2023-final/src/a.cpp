#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;

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
#define INF INT_MAX
#define INFF (9223372036854775800)
#define TIME_LIMIT (1.99)
#define def (20101)
// #define MOD (1000000007)
#define MOD (998244353)
#define PI (3.14159265359)
#define R ((ll)1000000000)
#define D ((ll)1000000)

bool is_out_circle(ll x, ll y) {
  ll _R = x * x + y * y;
  if (_R > R * R)
    return true;
  return false;
}

void out_point(ll x, ll y) { cout << x << " " << y << endl; }

void print_lattice(vP &candi_center) {
  ll interval = R / 6;
  // 左上から
  for (ll y = R;; y -= interval) {
    if (y < -R)
      break;
    for (ll x = -R;; x += interval) {
      if (R < x)
        break;
      if (is_out_circle(x, y))
        continue;
      candi_center.push_back({x, y});
    }
  }
}

int main() {
  srand((unsigned int)time(NULL));
  long double sd;
  cin >> sd;
  vector<P> candi_center;
  print_lattice(candi_center);

  vll dx = {-1, -1, 1, 1}, dy = {1, -1, -1, 1};
  vector<long double> l_t = {3 * PI / 2, 0, PI / 2, PI};
  vector<long double> r_t = {2 * PI, PI / 2, PI, 3 * PI / 2};
  int ans = 0;
  int cout_sum = 0;
  vi seen(candi_center.size(), 0);
  while (1) {
    ll cx = rand() % (2 * R) - R;
    ll cy = rand() % (2 * R) - R;
    vi distinate;
    rep(i, 0, (int)candi_center.size()) {
      if (seen[i] == 0)
        distinate.push_back(i);
    }
    if ((int)distinate.size() == 0) {
      rep(i, 0, (int)candi_center.size()) {
        distinate.push_back(i);
        seen[i] = 0;
      }
    }
    ll ran = rand() % distinate.size();
    seen[distinate[ran]] = 1;
    cx = candi_center[distinate[ran]].first;
    cy = candi_center[distinate[ran]].second;

    ll len = R / (6 - ans / 10);
    int flag = 1;
    rep(i, 0, 4) {
      ll px = cx + len / 2 * dx[i];
      ll py = cy + len / 2 * dy[i];
      if (is_out_circle(px, py))
        flag = 0;
    }
    if (!flag)
      continue;

    int out_sum = 0;
    int dis_flag = 0;
    rep(i, 0, 4) {
      if (out_sum == 3)
        continue;
      ll px = cx + len / 2 * dx[i];
      ll py = cy + len / 2 * dy[i];
      out_point(px, py);
      cout_sum++;
      cout.flush();
      ll type;
      cin >> type;
      if (type == 0) {
        long double theta;
        cin >> theta;
        if (l_t[i] > theta || theta > r_t[i])
          out_sum++;
      }
      if (type == 1) {
        ll xk, yk;
        cin >> xk >> yk;
        ans++;
        dis_flag = 1;
        break;
      }
      if (type == 2) {
        ans++;
        return 0;
      }
      if (cout_sum == 1000)
        return 0;
    }
    if (cout_sum == 1000)
      return 0;
    if (dis_flag)
      continue;
    if (3 <= out_sum)
      continue;

    // ランダムウォーク
    ll wx = cx, wy = cy;
    ll f_step = len / 2;
    ll step = f_step;
    dis_flag = 0;
    rep(i, 0, 20) {
      out_point(wx, wy);
      cout.flush();
      cout_sum++;
      ll type;
      cin >> type;
      if (type == 0) {
        long double theta;
        cin >> theta;
        wx += step * cos(theta);
        wy += step * sin(theta);
        if (is_out_circle(wx, wy)) {
          wx = max(cx - len / 2, wx);
          wx = min(cx + len / 2, wx);
          wy = max(cy - len / 2, wy);
          wy = min(cy + len / 2, wy);
        }
        step /= 1.3;
      }
      if (type == 1) {
        ll xk, yk;
        cin >> xk >> yk;
        ans++;
        dis_flag = 1;
        break;
      }
      if (type == 2) {
        ans++;
        return 0;
      }
      if (cout_sum == 1000)
        return 0;
    }
    if (cout_sum == 1000)
      return 0;
    if (dis_flag)
      continue;
  }
}