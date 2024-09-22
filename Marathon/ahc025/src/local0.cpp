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

void sort_f(vll &at) {
  ll si = (ll)at.size();
  vll A;
  A.push_back(at[0]);
  ll _l = -1, _r = 1;

  while ((ll)A.size() != si) {
    stat_sort_num++;
    ll _mid = (_l + _r) / 2;

    ll nl = 1, nr = 1;
    vll l(nl), r(nr);
    if ((ll)A.size() == 1) {
      l[0] = A[0];
      r[0] = A.size();
    } else {
      l[0] = A[_mid];
      r[0] = A.size();
    }

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
      A.insert(A.begin() + _l + 1, at[A.size()]);
      _l = -1, _r = (ll)A.size();
    }
  }

  at.clear();
  for (ll cu : A)
    at.push_back(cu);
}

ld large_k_prb(ll N, ll k, ll _x) {
  ld x = _x * 1e3;
  ld Fx = 1 - exp(-lam * x);
  ld fx = lam * exp(-lam * x);
  ld prb = k * pow(Fx, k - 1) * pow(1 - Fx, N - k) * fx;
  rep(i, 1, k + 1) prb /= i;
  rep(i, 1, N - k + 1) prb /= i;
  rep(i, 1, N + 1) prb *= i;
  prb *= 1000;
  // if (prb > lam) cout << _x << " " << setprecision(3) << prb << endl;
  return prb;
}

void update_prb(char c, vll l, vll r, vvld &prob) {
  vvld c_prob(2, vld(1001, 0));
  rep(i, 0, 1001) c_prob[0][i] = c_prob[0][max((ll)0, i - 1)] + prob[l[0]][i];
  rep(i, 0, 1001) c_prob[1][i] = c_prob[1][max((ll)0, i - 1)] + prob[r[0]][i];

  vvld aft_prob(2, vld(1001, 0));
  rep(i, 0, 1001) {
    if (c == '<') {
      aft_prob[0][i] = prob[l[0]][i] * (1 - c_prob[1][i]);
      aft_prob[1][i] = prob[r[0]][i] * c_prob[0][max((ll)0, i - 1)];
    }
    if (c == '>') {
      aft_prob[0][i] = prob[l[0]][i] * c_prob[1][max((ll)0, i - 1)];
      aft_prob[1][i] = prob[r[0]][i] * (1 - c_prob[0][i]);
    }
    if (c == '=') {
      aft_prob[0][i] = prob[l[0]][i] * c_prob[1][i];
      aft_prob[1][i] = prob[r[0]][i] * c_prob[0][i];
    }
  }

  ld sum_l = 0.0, sum_r = 0.0;
  rep(i, 0, 1001) {
    sum_l += aft_prob[0][i];
    sum_r += aft_prob[1][i];
  }
  rep(i, 0, 1001) {
    prob[l[0]][i] = aft_prob[0][i] / sum_l;
    prob[r[0]][i] = aft_prob[1][i] / sum_r;
  }
}

ll calc_score(vll &E, vll &d) {
  ll t_bar = 0;
  rep(i, 0, N) t_bar += E[i];
  t_bar /= D;

  vll sum(D, 0);
  rep(i, 0, N) sum[d[i]] += E[i];
  ll score = 0;
  rep(i, 0, D) score += pow(sum[i] - t_bar, 2);
  score /= D;
  score = 1 + 100 * sqrt(score);

  return score;
}

void determine_group(vll &d, vvld &prob, ld time) {
  auto start_time = system_clock::now();

  vll E(N, 0);
  rep(i, 0, N) rep(j, 0, 1001) E[i] += prob[i][j] * j;

  rep(i, 0, N) d[i] = i % D;
  ll best_score = calc_score(E, d);

  while (1) {
    if (time < ela_times(start_time)) break;

    vll _d(N);
    _d = d;
    if (rand() % 2) {
      ll idx1 = rand() % N;
      ll idx2 = (idx1 + rand() % N) % N;
      swap(_d[idx1], _d[idx2]);
    } else {
      ll idx = rand() % N;
      _d[idx] = rand() % D;
    }

    ll score = calc_score(E, _d);
    if (chmin(best_score, score)) {
      d = _d;
    }
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

  vvld prob(N, vld(1001, 0));
  vvll used(N, vll(N, 0));
  vll at;
  ll num = 7;
  ll lp = 0;
  rep(i, 0, N) {
    vll at;
    rep(j, 0, num + (lp < N % num)) {
      at.push_back(i);
      i++;
      if (i == N) break;
    }
    i--;
    lp++;


    sort_f(at);
    rep(j, 0, (ll)at.size()) {
      ll cu = at[j];
      rep(_x, 0, 1001) prob[cu][_x] = large_k_prb((ll)at.size(), j + 1, _x);
    }
    rep(j, 0, (ll)at.size()) {
      rep(k, j + 1, (ll)at.size()) {
        used[at[j]][at[k]] = 1;
        used[at[k]][at[j]] = 1;
      }
    }
  }

  vp candi;
  rep(i, 0, N) rep(j, i + 1, N) {
    if (used[i][j] == 0) candi.push_back({i, j});
  }
  random_device seed_gen;
  mt19937 engine(seed_gen());
  shuffle(candi.begin(), candi.end(), engine);

  ll prev = 0;
  rep(i, 0, N) {
    ll sum = 0;
    rep(j, 0, 1001) { sum += prob[i][j] * j * 1e3; }
    prev += abs(w[i] - sum);
  }

  rep(q, 0, Q - stat_sort_num) {
    if (candi.size() < q) break;
    ll nl = 1, nr = 1;
    vll l(1), r(1);
    l[0] = candi[q].first;
    r[0] = candi[q].second;
    out(nl, nr, l, r);

    char c = input_c(l, r);
    update_prb(c, l, r, prob);
  }

  ll aft = 0;
  rep(i, 0, N) {
    ll sum = 0;
    rep(j, 0, 1001) { sum += prob[i][j] * j * 1e3; }
    aft += abs(sum - w[i]);
  }

  vll d(N);
  determine_group(d, prob, TIME_LIMIT - ela_times(start_time));
  rep(i, 0, N) cout << d[i] << " ";
  cout << endl;
  debug_calc_ans(d);
  stat_time = ela_times(start_time);

  cerr << "N, Q         : " << N << " " << Q << endl;
  cerr << "stat_time    : " << stat_time << endl;
  cerr << "stat_sort_num: " << stat_sort_num << endl;
  cerr << "stat_ans     : " << stat_ans << endl;
}