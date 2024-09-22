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
// #define TIME_LIMIT (0.1)
#define ALL_TIME_LIMIT (2.5)
#define def (301010)
#define MOD (1000000007)
// #define MOD (998244353)
#define PI (3.14159265359)
// cout << fixed << setprecision(10) << << endl;
//偏角ソートはlong double!
// auto ite = s.lower_bound("B");
// cout << __builtin_popcount(N) << endl;

struct MoveAction {
  int before_row, before_col, after_row, after_col;
  MoveAction(int before_row, int before_col, int after_row, int after_col)
      : before_row(before_row), before_col(before_col), after_row(after_row),
        after_col(after_col) {}
};

struct ConnectAction {
  int c1_row, c1_col, c2_row, c2_col;
  ConnectAction(int c1_row, int c1_col, int c2_row, int c2_col)
      : c1_row(c1_row), c1_col(c1_col), c2_row(c2_row), c2_col(c2_col) {}
};

struct Result {
  vector<MoveAction> move;
  vector<ConnectAction> connect;
  Result(const vector<MoveAction> &move, const vector<ConnectAction> &con)
      : move(move), connect(con) {}
};

struct Computer {
  int x, y, c;
  Computer(int x, int y, int c) : x(x), y(y), c(c) {}
};

struct UnionFind {
  map<P, P> parent;
  map<P, int> sum;
  UnionFind() : parent(), sum() {}

  P find(P x) {
    // xが初めて登場
    if (parent.find(x) == parent.end()) {
      parent[x] = x;
      sum[x] = 1;
      return x;
    } else if (parent[x] == x) {
      return x;
    } else {
      parent[x] = find(parent[x]);
      return parent[x];
    }
  }

  void unite(P x, P y) {
    x = find(x);
    y = find(y);
    if (x != y) {
      parent[x] = y;
      sum[parent[x]] += sum[x];
    }
  }
};

struct Solver {
  static constexpr char USED = 'x';
  //右,下,左,上
  static constexpr int DR[4] = {0, 1, 0, -1};
  static constexpr int DC[4] = {1, 0, -1, 0};

  int N, K;
  int action_count_limit;
  mt19937 engine;
  vs field;
  vs used_field;
  UnionFind uf;

  vvP type_computers;
  vector<vector<string>> select_fields;
  vector<P> select_seeds;
  vector<vector<MoveAction>> select_moves;
  vector<vector<ConnectAction>> select_connects;

  //連結個数,種類,seed,行動回数
  // F best = {0, -1, -1, INF};

  Solver(int N, int K, const vs &field, const vs &used_field,
         const UnionFind &uf, int action_count_limit, int seed = 0)
      : N(N), K(K), used_field(used_field), uf(uf),
        action_count_limit(action_count_limit), field(field) {
    engine.seed(seed);

    type_computers.resize(K + 1);
    select_seeds.resize(K + 1);
    select_fields.resize(K + 1);
    select_moves.resize(K + 1);
    select_connects.resize(K + 1);

    rep(i, 1, K + 1) select_fields[i] = field;
    rep(i, 0, N) {
      rep(j, 0, N) if (field[i][j] != 'x') {
        type_computers[field[i][j] - '0'].push_back({i, j});
      }
    }
  }

  void print_pos(int row, int col) {
    cerr << "{" << row << ", " << col << "}" << endl;
  }

  bool can_move(int row, int col, int dir) const {
    int nrow = row + DR[dir];
    int ncol = col + DC[dir];
    if (0 <= nrow && nrow < N && 0 <= ncol && ncol < N) {
      // return field[nrow][ncol] == '0';
      return true;
    }
    return false;
  }

  P select_seed(int c) {
    // cerr << (int)type_computers[c].size() << endl;
    if ((int)type_computers[c].size() == 0)
      return {-1, -1};
    return type_computers[c][rand() % (int)type_computers[c].size()];
  }

  void fill_ava_vector(vvi &ava, int row, int col, vector<string> &_field,
                       int c, int &connect_sum) {
    // cerr << "call fill_ava_vector() " << row << " " << col << endl;
    ava[row][col] = 1;
    connect_sum++;
    rep(dir, 0, 4) {
      int atr = row, atc = col;
      while (1) {
        // cerr << atr << " " << atc << endl;
        if (!can_move(atr, atc, dir))
          break;
        atr += DR[dir];
        atc += DC[dir];
        if (ava[atr][atc] == 1)
          break;
        if (_field[atr][atc] - '0' != c && _field[atr][atc] != '0')
          break;
        // cerr << "{" << atr << ", " << atc << "}" << endl;
        ava[atr][atc] = 1;
        if (_field[atr][atc] - '0' == c)
          fill_ava_vector(ava, atr, atc, _field, c, connect_sum);
      }
    }
  }

  int bfs(vvi &ava, int row, int col, vector<string> &_field, int c) {
    vvi seen(N, vi(N, -1));
    queue<int> qu;
    qu.push(row * N + col);
    seen[row][col] = 0;
    while (!qu.empty()) {
      int at = qu.front();
      qu.pop();

      rep(dir, 0, 4) {
        int atr = at / N, atc = at % N;
        if (!can_move(atr, atc, dir))
          continue;
        atr += DR[dir];
        atc += DC[dir];
        if (seen[atr][atc] != -1)
          continue;
        if (_field[atr][atc] != '0')
          continue;
        seen[atr][atc] = seen[atr - DR[dir]][atc - DC[dir]] + 1;
        if (seen[atr][atc] < 10)
          qu.push(atr * N + atc);

        if (ava[atr][atc] == 1)
          return seen[atr][atc];
      }
    }
    return N * N;
  }

  void bfs2(vvi &ava, int row, int col, vector<string> &_field, int c,
            vector<MoveAction> &_ret, int &connect_sum) {
    // cerr << "call bfs2() " << row << " " << col << endl;
    vvi seen(N, vi(N, -1));
    vvi dir_g(N, vi(N, -1));
    queue<int> qu;
    qu.push(row * N + col);
    seen[row][col] = 0;
    //逆順時のスタート
    int star = -1, stac = -1;
    while (!qu.empty()) {
      int at = qu.front();
      qu.pop();
      rep(dir, 0, 4) {
        int atr = at / N, atc = at % N;
        if (!can_move(atr, atc, dir))
          continue;
        atr += DR[dir];
        atc += DC[dir];
        if (seen[atr][atc] != -1)
          continue;
        if (_field[atr][atc] != '0')
          continue;
        seen[atr][atc] = seen[atr - DR[dir]][atc - DC[dir]] + 1;
        dir_g[atr][atc] = dir;
        if (seen[atr][atc] < 10)
          qu.push(atr * N + atc);

        if (ava[atr][atc] == 1) {
          star = atr;
          stac = atc;
          break;
        }
      }
      if (star != -1)
        break;
    }
    // cerr << c << " {" << row << ", " << col << "} -> {" << star << ", " <<
    // stac
    //      << "}" << endl;
    while (!qu.empty())
      qu.pop();

    // cerr << "start fill() :";
    // print_pos(star, stac);
    swap(_field[star][stac], _field[row][col]);
    fill_ava_vector(ava, star, stac, _field, c, connect_sum);

    vector<MoveAction> at_ret;
    while (!(star == row && stac == col)) {
      int ner = star + DR[(dir_g[star][stac] + 2) % 4];
      int nec = stac + DC[(dir_g[star][stac] + 2) % 4];

      at_ret.emplace_back(ner, nec, star, stac);

      star = ner;
      stac = nec;
    }

    repd(i, (int)at_ret.size() - 1, -1) { _ret.push_back(at_ret[i]); }
  }

  void print_answer(const Result &res) {
    // cerr << "call print_answer()" << endl;
    cout << res.move.size() << endl;
    for (auto m : res.move) {
      cout << m.before_row << " " << m.before_col << " " << m.after_row << " "
           << m.after_col << endl;
    }
    cout << res.connect.size() << endl;
    for (auto m : res.connect) {
      cout << m.c1_row << " " << m.c1_col << " " << m.c2_row << " " << m.c2_col
           << endl;
    }
  }

  void move(int move_limit = -1) {
    // for (int i = 0; i < move_limit; i++) {
    //   int row = engine() % N;
    //   int col = engine() % N;
    //   int dir = engine() % 4;
    //   if (field[row][col] != '0' && can_move(row, col, dir)) {
    //     swap(field[row][col], field[row + DR[dir]][col + DC[dir]]);
    //     ret.emplace_back(row, col, row + DR[dir], col + DC[dir]);
    //     action_count_limit--;
    //   }
    // }

    rep(c, 1, K + 1) {
      P seed = select_seed(c);
      select_seeds[c] = seed;

      if (seed.first == -1)
        continue;
      vector<MoveAction> _ret;
      // 0: 何も置いていない, 1: c種類のコンピュータが置かれてある or cが置ける
      vvi ava(N, vi(N, 0));

      int move_sum = 0;
      int connect_sum = 0;

      vector<string> _field(N);
      _field = field;

      fill_ava_vector(ava, seed.first, seed.second, _field, c, connect_sum);

      while (1) {
        int min_dist = N * N;
        int min_pos;
        rep(i, 0, N) {
          rep(j, 0, N) {
            if (_field[i][j] - '0' != c || ava[i][j] == 1)
              continue;
            int _dist = bfs(ava, i, j, _field, c);
            // cerr << i << " " << j << " " << _dist << endl;
            if (_dist < min_dist) {
              min_dist = _dist;
              min_pos = i * N + j;
            }
          }
        }
        if (min_dist == N * N)
          break;
        if (move_limit < move_sum + min_dist)
          break;
        // cerr << "{" << min_pos / N << ", " << min_pos % N << "} " << min_dist
        //  << endl;
        move_sum += min_dist;
        bfs2(ava, min_pos / N, min_pos % N, _field, c, _ret, connect_sum);
      }

      // cerr << c << " " << connect_sum << " move_sum : " << move_sum << endl;

      select_fields[c] = _field;
      select_moves[c] = _ret;
    }
  }

  bool can_connect(int row, int col, int dir, int c) const {
    int nrow = row + DR[dir];
    int ncol = col + DC[dir];
    while (0 <= nrow && nrow < N && 0 <= ncol && ncol < N) {
      if (select_fields[c][nrow][ncol] == select_fields[c][row][col]) {
        return true;
      } else if (select_fields[c][nrow][ncol] != '0') {
        return false;
      }
      nrow += DR[dir];
      ncol += DC[dir];
    }
    return false;
  }

  ConnectAction line_fill(int row, int col, int dir, int c) {
    int nrow = row + DR[dir];
    int ncol = col + DC[dir];
    while (0 <= nrow && nrow < N && 0 <= ncol && ncol < N) {
      if (select_fields[c][nrow][ncol] == select_fields[c][row][col]) {
        return ConnectAction(row, col, nrow, ncol);
      }
      assert(select_fields[c][nrow][ncol] == '0');
      select_fields[c][nrow][ncol] = USED;
      nrow += DR[dir];
      ncol += DC[dir];
    }
    assert(false);
  }

  void connect() {
    // int min_c = -1;
    // int max_connect_sum = 0;
    // int min_move_sum = N * N;
    rep(c, 1, K + 1) {
      // print_answer(Result(select_moves[c], select_connects[c]));
      int connect_sum = 0;
      queue<P> qu;
      map<P, int> ma;
      if (select_seeds[c].first == -1)
        continue;
      // cerr << c << " ";
      // print_pos(select_seeds[c].first, select_seeds[c].second);
      ma[select_seeds[c]] = 1;
      qu.push(select_seeds[c]);
      while (!qu.empty()) {
        int row = qu.front().first;
        int col = qu.front().second;
        qu.pop();

        // if (action_count_limit <= connect_sum + (int)select_moves[c].size())
        //   break;
        connect_sum++;
        // print_pos(row, col);
        for (int dir = 0; dir < 4; dir++) {
          if (can_connect(row, col, dir, c)) {
            ConnectAction re = line_fill(row, col, dir, c);
            if (!ma.count({re.c2_row, re.c2_col}) &&
                action_count_limit > (int)select_connects[c].size() +
                                         (int)select_moves[c].size()) {
              qu.push({re.c2_row, re.c2_col});
              ma[{re.c2_row, re.c2_col}] = 1;
              select_connects[c].emplace_back(re);
            }
          }
        }
      }

      // if ((max_connect_sum < connect_sum) ||
      //     (max_connect_sum == connect_sum &&
      //      (int)select_moves[c].size() < min_move_sum)) {
      //   max_connect_sum = connect_sum;
      //   min_move_sum = (int)select_moves[c].size();
      //   min_c = c;
      // }
      // cerr << c << " " << connect_sum << " " << (int)select_moves[c].size()
      //      << endl;
    }
    // return {min_c, min_move_sum, max_connect_sum};
  }

  // int connect_forest_bfs(int row, int col, int c, vvi &to_forest) {
  //   int ret = 0;
  //   queue<P> qu;
  //   qu.push({row, col});
  //   to_forest[row][col] = 1;
  //   while (!qu.empty()) {
  //     P cu = qu.front();
  //     qu.pop();
  //     rep(dir, 0, 4) {
  //       if (!can_move(cu.first, cu.second, dir))
  //         continue;
  //       int tor = cu.first + DR[dir];
  //       int toc = cu.second + DC[dir];
  //       if (select_fields[c][tor][toc] != USED)
  //         continue;
  //       if (to_forest[tor][toc] == 1)
  //         continue;
  //       qu.push({tor, toc});
  //       to_forest[tor][toc] = 1;
  //       if (used_field[tor][toc] - c == '0')
  //         ret++;
  //     }
  //   }
  //   return ret;
  // }

  P forest_connect() {
    int best_c = -1;
    int max_score = -1;
    rep(c, 1, K + 1) {
      if (select_connects[c].size() == 0)
        continue;
      int _score = (int)select_connects[c].size() *
                   ((int)select_connects[c].size() + 1) / 2;
      _score /= ((int)select_connects[c].size() + (int)select_moves[c].size());
      // cerr << "first :" << _score << " -> ";
      if (max_score < _score) {
        max_score = _score;
        best_c = c;
      }
      if (action_count_limit <
          (int)select_connects[c].size() + (int)select_moves[c].size() + 2)
        continue;
      P from_pos = {-1, -1};
      P mid_pos = {-1, -1};
      P to_pos = {-1, -1};

      vvP link_computer(N, vP(N, {-1, -1}));
      map<P, int> from_forest;
      for (auto m : select_connects[c]) {
        from_forest[{m.c1_row, m.c1_col}] = 1;
        from_forest[{m.c2_row, m.c2_col}] = 1;
      }

      //中継点を探す
      queue<P> qu;
      for (auto m : from_forest) {
        rep(dir, 0, 4) {
          int ner = m.first.first;
          int nec = m.first.second;
          while (can_move(ner, nec, dir)) {
            ner += DR[dir];
            nec += DC[dir];
            if (select_fields[c][ner][nec] == USED)
              break;
            // from_forestに当たった
            if (from_forest.count({ner, nec}))
              break;
            if (select_fields[c][ner][nec] == '0')
              continue;
            if (link_computer[ner][nec].first == -1) {
              link_computer[ner][nec] = m.first;
              qu.push({ner, nec});
            }
            break;
          }
        }
      }

      vvi to_forest(N, vi(N, 0));
      //中継点から森へ
      while (!qu.empty()) {
        P cu = qu.front();
        qu.pop();
        rep(dir, 0, 4) {
          int ner = cu.first;
          int nec = cu.second;
          while (can_move(ner, nec, dir)) {
            ner += DR[dir];
            nec += DC[dir];
            if (select_fields[c][ner][nec] == USED) {
              if (!from_forest.count({ner, nec}) &&
                  used_field[ner][nec] - c == '0') {
                int to_connect_size = uf.sum[uf.find({ner, nec})];
                int sum = 0;
                rep(_i, to_connect_size,
                    to_connect_size + (int)select_connects[c].size()) sum += _i;
                sum -= to_connect_size + (int)select_connects[c].size();
                sum /= ((int)select_connects[c].size() +
                        (int)select_moves[c].size() + 2);
                if (_score < sum) {
                  _score = sum;
                  from_pos = link_computer[cu.first][cu.second];
                  mid_pos = cu;
                  to_pos = {ner, nec};
                }
              }
              break;
            }
            if (select_fields[c][ner][nec] != '0')
              break;
          }
        }
      }

      if (max_score < _score) {
        max_score = _score;
        best_c = c;
        if (from_pos.first != -1) {
          select_connects[c].push_back(ConnectAction(
              from_pos.first, from_pos.second, mid_pos.first, mid_pos.second));
          select_connects[c].push_back(ConnectAction(
              mid_pos.first, mid_pos.second, to_pos.first, to_pos.second));
        }
      }
      // cerr << c << " " << _score << endl;
    }
    // cerr << "best_c :" << best_c << " max_score : " << max_score << endl;
    return {best_c, max_score};
  }

  Result solve() {
    // create random moves
    int max_score = 0;
    vector<MoveAction> ret1;
    vector<ConnectAction> ret2;

    auto startClock = system_clock::now();
    double TIME_LIMIT = 0.1;
    if (N * N / 3 > K * 100)
      TIME_LIMIT = ALL_TIME_LIMIT / (K + 3);
    while (1) {
      const double time =
          duration_cast<microseconds>(system_clock::now() - startClock)
              .count() *
          1e-6;
      if (time > TIME_LIMIT)
        break;
      move((int)(action_count_limit / (double)2.2));
      // cerr << "after move()" << endl;
      connect();
      // cerr << "after connect()" << endl;
      P cu = forest_connect();
      if (cu.first == -1)
        continue;
      if (max_score < cu.second) {
        max_score = cu.second;
        ret1 = select_moves[cu.first];
        ret2 = select_connects[cu.first];
      }
      rep(j, 1, K + 1) {
        select_connects[j].clear();
        select_moves[j].clear();
        select_fields[j] = field;
      }
      // print_answer(Result(ret1,ret2));
      // cerr << "after solve " << endl;
      // cerr << "select_seeds ";
      // rep(j, 1, K + 1) print_pos(select_seeds[j].first,
      // select_seeds[j].second); cerr << endl;
    }
    // cerr << max_connect_num << " " << min_move_num << endl;
    return Result(ret1, ret2);
  }
};

int calc_score(int N, vector<string> field, const Result &res) {
  for (auto r : res.move) {
    assert(field[r.before_row][r.before_col] != '0');
    assert(field[r.after_row][r.after_col] == '0');
    swap(field[r.before_row][r.before_col], field[r.after_row][r.after_col]);
  }

  UnionFind uf;
  for (auto r : res.connect) {
    P p1(r.c1_row, r.c1_col), p2(r.c2_row, r.c2_col);
    uf.unite(p1, p2);
  }

  vector<P> computers;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (field[i][j] != '0') {
        computers.emplace_back(i, j);
      }
    }
  }

  int score = 0;
  for (int i = 0; i < (int)computers.size(); i++) {
    for (int j = i + 1; j < (int)computers.size(); j++) {
      auto c1 = computers[i];
      auto c2 = computers[j];
      if (uf.find(c1) == uf.find(c2)) {
        score +=
            (field[c1.first][c1.second] == field[c2.first][c2.second]) ? 1 : -1;
      }
    }
  }

  return max(score, 0);
}

void print_answer(const Result &res) {
  // cerr << "call print_answer()" << endl;
  cout << res.move.size() << endl;
  for (auto m : res.move) {
    cout << m.before_row << " " << m.before_col << " " << m.after_row << " "
         << m.after_col << endl;
  }
  cout << res.connect.size() << endl;
  for (auto m : res.connect) {
    cout << m.c1_row << " " << m.c1_col << " " << m.c2_row << " " << m.c2_col
         << endl;
  }
}

void change_field(vs &field, const Result &res, vs &used_field, UnionFind &uf) {
  static constexpr int DR[4] = {0, 1, 0, -1};
  static constexpr int DC[4] = {1, 0, -1, 0};

  for (auto m : res.move) {
    swap(field[m.before_row][m.before_col], field[m.after_row][m.after_col]);
  }
  char c = field[res.connect[0].c1_row][res.connect[0].c1_col];
  for (auto m : res.connect) {
    used_field[m.c1_row][m.c1_col] = c;
    used_field[m.c2_row][m.c2_col] = c;
    uf.unite({m.c1_row, m.c1_col}, {m.c2_row, m.c2_col});
    rep(dir, 0, 4) {
      if ((m.c2_row - m.c1_row) * DR[dir] <= 0 &&
          (m.c2_col - m.c1_col) * DC[dir] <= 0)
        continue;
      field[m.c1_row][m.c1_col] = 'x';
      while (!(m.c1_row == m.c2_row && m.c1_col == m.c2_col)) {
        m.c1_row += DR[dir];
        m.c1_col += DC[dir];
        field[m.c1_row][m.c1_col] = 'x';
      }
    }
  }
}

int main() {
  int N, K;
  cin >> N >> K;
  vector<string> field(N), inisial_field(N), used_field(N);
  rep(i, 0, N) cin >> field[i];
  rep(i, 0, N) rep(j, 0, N) used_field[i] += '0';

  inisial_field = field;

  rep(i, 0, 1) {
    field = inisial_field;
    srand((unsigned int)time(NULL));
    auto all_startClock = system_clock::now();

    vector<MoveAction> ans_move;
    vector<ConnectAction> ans_connect;
    UnionFind uf;
    while (1) {
      const double all_time =
          duration_cast<microseconds>(system_clock::now() - all_startClock)
              .count() *
          1e-6;
      if (all_time > ALL_TIME_LIMIT)
        break;
      Solver s(N, K, field, used_field, uf,
               K * 100 - (int)ans_connect.size() - (int)ans_move.size());
      auto ret = s.solve();
      // // cerr << "after2" << endl;
      // if ((int)ans_connect.size() + (int)ret.connect.size() +
      //         (int)ans_move.size() + (int)ret.move.size() >
      //     K * 100)
      //   continue;
      if (ret.connect.empty())
        continue;
      change_field(field, ret, used_field, uf);
      ans_move.insert(ans_move.end(), ret.move.begin(), ret.move.end());
      ans_connect.insert(ans_connect.end(), ret.connect.begin(),
                         ret.connect.end());
      // print_answer(Result(ans_move, ans_connect));
      // rep(j, 0, N) {
      //   rep(k, 0, N) {
      //     if (used_field[j][k] != '0') {
      //       if (uf.sum[{j,k}] > 1) {
      //         cerr << uf.sum[{j, k}] << " {" << j << " ," << k  << " } "<<
      //         endl;
      //       }
      //     }
      //   }
      // }
      // rep(j, 0, 5) cerr << endl;
    }

    print_answer(Result(ans_move, ans_connect));
    // cerr << (int)ans_move.size() + (int)ans_connect.size() << endl;
    // cerr << "Score = ";
    cerr << calc_score(N, inisial_field, Result(ans_move, ans_connect)) << endl;
  }
  return 0;
}

// 238332が元
// bfs()をO(100 * N ^2) -> O(100 * 10 * 10)