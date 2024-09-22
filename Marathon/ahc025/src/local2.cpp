#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef long double ld;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
typedef pair<ll, ll> pll;
typedef vector<pll> vp;
typedef vector<vp> vvp;
typedef vector<vvp> vvvp;
typedef vector<ld> vld;
typedef vector<vld> vvld;
#define rep(i, l, n) for (ll i = (ll)(l); i < (ll)(n); i++)
#define repd(i, n, l) for (ll i = (ll)(n); i > (ll)(l); i--)
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
template <class T, class S> inline bool chmax(T &a, const S &b) {
  return (a < b ? a = b, 1 : 0);
}
template <class T, class S> inline bool chmin(T &a, const S &b) {
  return (a > b ? a = b, 1 : 0);
}
#define INF INT_MAX
#define INFF (9223372036854775800)
#define TIME_LIMIT (1.9)
#define def (20101)
// #define MOD (1000000007)
// #define MAX (2147483647)
#define MAX (1073741823)
#define MOD (998244353)

ld ela_times(system_clock::time_point &clock) {
  return duration_cast<microseconds>(system_clock::now() - clock).count() *
         1e-6;
}

ld stat_time = 0.0;
ld lam = 0.00001;
ll stat_loop = 0;
ll stat_ans = 0;
ll stat_sort_num = 0;

ll N, D, Q;
vll w;

char input_c(vll &l, vll &r) {
  char c;

#ifdef _DEBUG
  ll suml = 0, sumr = 0;
  for (ll cu : l)
    suml += w[cu];
  for (ll cu : r)
    sumr += w[cu];
  map<ll, ll> ma;
  for (ll cu : l)
    ma[cu] = 1;
  for (ll cu : r)
    if (ma.count(cu)) cerr << "err same element " << cu << endl;
  if (suml < sumr) c = '<';
  else if (suml > sumr) c = '>';
  else c = '=';

#else
  cin >> c;
#endif

  return c;
}

void out(ll nr, ll nl, vll l, vll r) {
  cout << nl << " " << nr << " ";
  for (ll cu : l)
    cout << cu << " ";
  for (ll cu : r)
    cout << cu << " ";
  cout << endl;
  cout.flush();
}

void sort_group(vll &sorted_d, ll idx, vvll &group) {
  if (sorted_d.size() == 0) {
    sorted_d.push_back(idx);
    return;
  }

  ll si = (ll)sorted_d.size();
  ll _l = -1, _r = si;

  while (1) {
    if (Q <= stat_sort_num) break;
    ll _mid = (_l + _r) / 2;

    ll nl = group[sorted_d[_mid]].size(), nr = group[idx].size();
    vll l(nl), r(nr);
    l = group[sorted_d[_mid]], r = group[idx];

    stat_sort_num++;
    out(nl, nr, l, r);
    char c = input_c(l, r);

    int flag = 0;
    if (c == '<') _l = _mid;
    else if (c == '>') _r = _mid;
    else {
      _l = _mid;
      flag = 1;
    }

    if (_r - _l <= 1) flag = 1;

    if (flag) {
      sorted_d.insert(sorted_d.begin() + _l + 1, idx);
      break;
    }
  }
}

void group_sub_elem(ll g_idx, ll idx, vll &sorted_d, vvll &group) {
  for (auto ite = group[g_idx].begin(); ite != group[g_idx].end();) {
    if (*ite == idx) ite = group[g_idx].erase(ite);
    else ite++;
  }
  for (auto ite = sorted_d.begin(); ite != sorted_d.end();) {
    if (*ite == g_idx) ite = sorted_d.erase(ite);
    else ite++;
  }
}

void group_add_elem(ll g_idx, ll idx, vll &sorted_d, vvll &group) {
  group[g_idx].push_back(idx);
  for (auto ite = sorted_d.begin(); ite != sorted_d.end();) {
    if (*ite == g_idx) ite = sorted_d.erase(ite);
    else ite++;
  }
}

void debug_calc_ans(vll &d) {
#ifdef _DEBUG
  vll sum(D, 0);
  rep(i, 0, N) sum[d[i]] += w[i];
  ll all_sum = 0;
  rep(i, 0, N) all_sum += w[i];
  all_sum /= D;
  rep(i, 0, D) stat_ans += pow(sum[i] - all_sum, 2);
  stat_ans /= D;
  stat_ans = 1 + 100 * sqrt(stat_ans);

  vll diff(D, 0);
  vll correct_sum(D, 0);
  rep(i, 0, N) correct_sum[d[i]] += w[i];
  // rep(i, 0, D) cerr << all_sum - correct_sum[i] << endl;
#endif
}

int main() {
  auto start_time = system_clock::now();
  srand((unsigned int)time(NULL));

  cin >> N >> D >> Q;
#ifdef _DEBUG
  w.resize(N);
  rep(i, 0, N) cin >> w[i];
#endif

  vll E(N, 0);
  vvll group(D);
  rep(i, 0, N) group[i % D].push_back(i);

  vll sorted_d;
  rep(i, 0, D) { sort_group(sorted_d, i, group); }

  while (1) {
    if (Q <= stat_sort_num) break;

    ll g_idx = sorted_d[D - 1 - (rand() % 10 == 0)];
    if (group[g_idx].size() == 1) continue;

    vll candi;
    for (ll cu : group[g_idx])
      candi.push_back(cu);
    ll idx = candi[rand() % (ll)group[g_idx].size()];

    group_sub_elem(g_idx, idx, sorted_d, group);
    sort_group(sorted_d, g_idx, group);

    g_idx = sorted_d[(rand() % 10 == 0)];
    group_add_elem(g_idx, idx, sorted_d, group);
    sort_group(sorted_d, g_idx, group);
  }

#ifdef _DEBUG
  rep(i, 0, D) {
    ll sum = 0;
    ll _sum = 0;
    for (ll cu : group[sorted_d[i]]) {
      sum += w[cu];
      _sum += E[cu];
    }
    cerr << "sum, _sum : " << sum << " " << _sum << endl;
  }
#endif

  vll d(N);
  rep(i, 0, D) for (ll cu : group[i]) d[cu] = i;
  rep(i, 0, N) cout << d[i] << " ";
  cout << endl;
  debug_calc_ans(d);
  stat_time = ela_times(start_time);

  cerr << "N, Q         : " << N << " " << Q << endl;
  cerr << "stat_time    : " << stat_time << endl;
  cerr << "stat_sort_num: " << stat_sort_num << endl;
  cerr << "stat_ans     : " << stat_ans << endl;
}

