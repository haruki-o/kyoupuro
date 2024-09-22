#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef pair<ll, ll> pll;
typedef tuple<ll, ll, ll> tll;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vvi> vvvi;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
typedef vector<pll> vp;
typedef vector<vp> vvp;
#define rep(i, l, n) for (ll i = (ll)(l); i < (ll)(n); i++)
#define repd(i, n, l) for (ll i = (ll)(n); i > (ll)(l); i--)
#define Copy(from, to) copy(from.begin(), from.end(), to.begin())
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
#define INF INT_MAX
#define INFF (9223372036854775800)
#define TIME_LIMIT (1.95)
#define def (2010101)
// #define MOD (1000000007)
#define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision(10) << << endl;
// 偏角ソートはlong ddouble!
// auto ite = s.lower_bound("B");
// cout << __builtin_popcount(N) << endl;

ll N = 30;
vvll b(30, vll(30, -1));
vll dx = {-1, -1, 0, 0, 1, 1};
vll dy = {-1, 0, -1, 1, 0, 1};
// 上に行くのは偉い
vll dsum = {2, 3, 1, 1, 1, 1};
vll eval_dx = {-1, -1, 0, 0, 1, 1};
vll eval_dy = {-1, 1, -1, 1, -1, 1};

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

// 状態の初期化
void init(vvll &ans_b) {
  ll sum = 0;
  rep(i, 0, 30) rep(j, 0, i + 1) ans_b[i][j] = sum, sum++;
}

// 状態遷移
// void modify(vvll &ans_b) {
//   ll idx1 = rand() % 465;
//   ll idx2 = rand() % 465;
//   pll p1, p2;
//   ll sum = 0;
//   rep(i, 0, 30) {
//     rep(j, 0, i + 1) {
//       if (sum == idx1)
//         p1 = {i, j};
//       if (sum == idx2)
//         p2 = {i, j};
//       sum++;
//     }
//   }
//   // cout << p1.first << " " << p1.second << " -> " << p2.first << " " <<
//   p2.second << endl; swap(ans_b[p1.first][p1.second],
//   ans_b[p2.first][p2.second]);
// }

void modify(vvll &ans_b) {
  pll p1, p2;
  ll idx1 = rand() % 465;
  ll sum = 0;
  rep(i, 0, 30) {
    rep(j, 0, i + 1) {
      if (sum == idx1)
        p1 = {i, j}, p2 = {i, j};
      sum++;
    }
  }

  ll dist = p1.first + 1;
  p2 = {p1.first, (p1.second + rand() % dist) % dist};
  if (rand() % 2 == 0) {
    if (p2.second < p2.first / 2)
      p2.first--;
  }
  swap(ans_b[p1.first][p1.second], ans_b[p2.first][p2.second]);
}

// 状態のスコア計算
ll eval_score(vvll &ans_b) {
  ll ret = 0;

  ll sum = 0;
  rep(i, 0, 29) {
    rep(j, 0, i + 1) {
      if (ans_b[i][j] > ans_b[i + 1][j] || ans_b[i][j] > ans_b[i + 1][j + 1])
        sum++;
    }
  }
  if (sum != 0)
    return 50000 - 50 * sum;

  ll K = 0;
  map<ll, pll> ma;
  rep(i, 0, 30) rep(j, 0, i + 1) ma[b[i][j]] = {i, j};
  rep(i, 0, 30) rep(j, 0, i + 1) {
    ll x = ma[ans_b[i][j]].first;
    ll y = ma[ans_b[i][j]].second;
    K += min(abs(x - i), abs(y - j)) + abs(abs(x - i) - abs(y - j));
  }

  return 100000 - 5 * K;
}

ll raw_score(vvll &ans_b) {
  ll sum = 0;
  rep(i, 0, 29) {
    rep(j, 0, i + 1) {
      if (ans_b[i][j] > ans_b[i + 1][j] || ans_b[i][j] > ans_b[i + 1][j + 1])
        sum++;
    }
  }
  if (sum != 0)
    return 50000 - 50 * sum;

  sum = 0;
  vvll _b(30, vll(30, 0));
  _b = b;
  vvll used(30, vll(30, 0));
  // _bを変更していく(_bをans_bに近づける)
  map<ll, pll> ma;
  rep(i, 0, 30) rep(j, 0, i + 1) ma[ans_b[i][j]] = {i, j};
  rep(num, 0, 20) {
    priority_queue<tll, vector<tll>, greater<tll>> pq;
    rep(i, 0, 30) rep(j, 0, i + 1) {
      ll _num = min(min(j, i - j), min(i, 29 - i));
      if (_num == num) {
        pq.push(
            {abs(i - ma[_b[i][j]].first) + abs(j - ma[_b[i][j]].second), i, j});
      }
    }
    while (!pq.empty()) {
      auto cu = pq.top();
      pq.pop();
      ll i = get<1>(cu), j = get<2>(cu);
      ll from_i, from_j;
      rep(_i, 0, 30) rep(_j, 0, _i + 1) {
        if (ans_b[i][j] == _b[_i][_j])
          from_i = _i, from_j = _j;
      }
      cout << from_i << " " << from_j << " -> " << i << " " << j << endl;
      while (1) {
        if (from_i == i && from_j == j)
          break;
        cout << from_i << " " << from_j << endl;
        ll best_dir = 2;
        rep(dir, 0, 6) {
          if ((i - from_i) * eval_dx[dir] <= 0 &&
              (j - from_j) * eval_dy[dir] <= 0)
            continue;
          ll to_i = from_i + dx[dir];
          ll to_j = from_j + dy[dir];
          if (to_i < 0 || to_i >= 30 || to_j < 0 || to_j >= 30)
            continue;
          if (to_i < to_j)
            continue;
          if (used[to_i][to_j] == 1)
            continue;
          if (dsum[best_dir] <= dsum[dir])
            best_dir = dir;
          // cout << dir << " " << (i - from_i) * eval_dx[dir]  << " " << (j - from_j) * eval_dy[dir] << endl; 
        }
        ll to_i = from_i + dx[best_dir];
        ll to_j = from_j + dy[best_dir];
        swap(_b[from_i][from_j], _b[to_i][to_j]);
        sum++;
        from_i = to_i;
        from_j = to_j;
      }
      used[i][j] = 1;
    }
  }
  return 100000 - 5 * sum;
}

void ans_output(vvll &ans_b) {
  ll K = 0;
  vp pre, aft;

  vvll _b(30, vll(30, 0));
  _b = b;
  vvll used(30, vll(30, 0));
  // _bを変更していく(_bをans_bに近づける)
  rep(i, 0, 30) {
    rep(j, 0, i + 1) {
      ll from_i, from_j;
      rep(_i, 0, 30) rep(_j, 0, _i + 1) {
        if (ans_b[i][j] == _b[_i][_j])
          from_i = _i, from_j = _j;
      }
      // cout << from_i << " " << from_j << " -> " << i << " " << j << endl;
      while (1) {
        if (from_i == i && from_j == j)
          break;
        // cout << from_i << " " << from_j << endl;
        ll best_dir = 2;
        ll best_sum = -1;
        rep(dir, 0, 6) {
          if ((i - from_i) * eval_dx[dir] <= 0 &&
              (j - from_j) * eval_dy[dir] <= 0)
            continue;
          ll to_i = from_i + dx[dir];
          ll to_j = from_j + dy[dir];
          if (to_i < 0 || to_i >= 30 || to_j < 0 || to_j >= 30)
            continue;
          if (to_i < to_j)
            continue;
          if (used[to_i][to_j] == 1)
            continue;
          ll _sum = 0;
          _sum += (i - from_i) * eval_dx[dir] > 0;
          _sum += (j - from_j) * eval_dy[dir] > 0;
        }
        ll to_i = from_i + dx[best_dir];
        ll to_j = from_j + dy[best_dir];
        swap(_b[from_i][from_j], _b[to_i][to_j]);
        pre.push_back({from_i, from_j});
        aft.push_back({to_i, to_j});
        from_i = to_i;
        from_j = to_j;
      }
      used[i][j] = 1;
    }
  }

  K = pre.size();
  cout << K << endl;
  rep(i, 0, K) {
    cout << pre[i].first << " " << pre[i].second << " " << aft[i].first << " "
         << aft[i].second << endl;
  }
}

// 焼きなまし
void sa(vvll &ans_b) {
  auto start_time = system_clock::now();
  double start_temp = 50, end_temp = 1; // 適当な値を入れる（後述）
  ll best_score = raw_score(ans_b);

  while (1) {
    stat_loop++;
    auto now_time = system_clock::now();
    if (TIME_LIMIT < ela_times(start_time))
      break;

    vvll new_ans_b(30, vll(30));
    new_ans_b = ans_b;
    modify(new_ans_b);

    auto stat_start_time = system_clock::now();
    ll new_score = raw_score(new_ans_b);
    stat_time += ela_times(stat_start_time);

    // 温度関数
    double temp = start_temp +
                  (end_temp - start_temp) * ela_times(start_time) / TIME_LIMIT;

    // 遷移確率関数(最大化)
    double prob = exp((new_score - best_score) / temp);

    // cout << new_score << " " << best_score << " " << new_score - best_score
    //  << endl;

    // if (prob > (rand() % INF) / (double)INF) { // 確率probで遷移する
    if (new_score > best_score) {
      ans_b = new_ans_b;
      best_score = new_score;
    }
  }
}

int main() {
  srand((unsigned int)time(NULL));

  rep(i, 0, 30) rep(j, 0, i + 1) cin >> b[i][j];

  vvll ans_b(30, vll(30));
  init(ans_b);
  sa(ans_b);
  ans_output(ans_b);

  cerr << "loop : " << stat_loop << endl;
  cerr << "time : " << stat_time << endl;
  cerr << "score : " << raw_score(ans_b) << endl;
}