#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef pair<ll, ll> P;
typedef pair<P, ll> Pll;
typedef pair<ll, P> llP;
typedef tuple<int, int, int> T;
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
template <class T, class S> inline bool chmax(T &a, const S &b) {
  return (a < b ? a = b, 1 : 0);
}
template <class T, class S> inline bool chmin(T &a, const S &b) {
  return (a > b ? a = b, 1 : 0);
}
#define INF INT_MAX
#define INFF (9223372036854775800)
#define TIME_LIMIT (5.5)
#define def (20101)
// #define MOD (1000000007)
// #define MAX (2147483647)
#define MAX (1073741823)
#define MOD (998244353)
#define PI (3.14159265359)

void print_num(ll num) {
  ll sum = 0;
  vc ret;
  while (num) {
    ret.push_back(num % 10 + '0');
    num /= 10;
    sum++;
    if (sum == 3 && num) {
      sum = 0;
      ret.push_back(',');
    }
  }
  repd(i, (ll)ret.size() - 1, -1) cerr << ret[i];
}

double ela_times(system_clock::time_point &clock) {
  return duration_cast<microseconds>(system_clock::now() - clock).count() *
         1e-6;
}

int D;
vi dx = {1, 0, 0, -1, 0, 0};
vi dy = {0, 1, 0, 0, -1, 0};
vi dz = {0, 0, 1, 0, 0, -1};

double stat_time = 0.0;

struct Cube {
  vvvi b;

  Cube() { b.assign(D, vvi(D, vi(D, 0))); }

  Cube(vvi &f, vvi &r) {
    b.assign(D, vvi(D, vi(D, 0)));
    rep(x, 0, D) rep(y, 0, D) rep(z, 0, D) {
      if (f[z][x] == 1 && r[z][y] == 1)
        b[x][y][z] = 1;
    }
  }

  Cube(vvi &f, vvi &r, Cube &nec_c) {
    b.assign(D, vvi(D, vi(D, 0)));
    b = nec_c.b;
    while (1) {
      bool loop_flag = 1;
      rep(x, 0, D) rep(y, 0, D) rep(z, 0, D) {
        bool flag = 0;
        if (b[x][y][z]) {
          bool flag1 = 0, flag2 = 0;
          rep(_x, 0, D) if (_x != x && b[_x][y][z]) flag1 = 1;
          rep(_y, 0, D) if (_y != y && b[x][_y][z]) flag2 = 1;
          if (flag1 && flag2) {
            b[x][y][z] = 0;
            loop_flag = 0;
          }
        }
      }
      if (loop_flag)
        break;
    }
  }

  void update(Cube &nec_c) {
    b = nec_c.b;
    vvvi fix_c(D, vvi(D, vi(D, 0)));
    rep(i, 0, 2 * D) {
      int x = rand() % D;
      int y = rand() % D;
      int z = rand() % D;
      if (nec_c.b[x][y][z])
        fix_c[x][y][z] = 1;
    }
    while (1) {
      bool loop_flag = 1;
      rep(x, 0, D) rep(y, 0, D) rep(z, 0, D) {
        bool flag = 0;
        if (b[x][y][z] && !fix_c[x][y][z]) {
          bool flag1 = 0, flag2 = 0;
          rep(_x, 0, D) if (_x != x && b[_x][y][z]) flag1 = 1;
          rep(_y, 0, D) if (_y != y && b[x][_y][z]) flag2 = 1;
          if (flag1 && flag2) {
            b[x][y][z] = 0;
            loop_flag = 0;
          }
        }
      }
      if (loop_flag)
        break;
    }
  }
};

struct Point {
  int x, y, z;

  Point() {}

  Point(int x, int y, int z) : x(x), y(y), z(z){};

  bool out_range() {
    if (x < 0 || x >= D || y < 0 || y >= D || z < 0 || z >= D)
      return true;
    return false;
  }

  bool to() {
    Point to_cor(x, y, z);
    rep(i, 0, 6) {
      to_cor.x += dx[i];
      to_cor.y += dy[i];
      to_cor.z += dz[i];
      if (to_cor.out_range())
        return false;
    }
    x = to_cor.x;
    y = to_cor.y;
    z = to_cor.z;
    return true;
  }

  bool to(int dir) {
    Point to_cor(x, y, z);
    to_cor.x += dx[dir];
    to_cor.y += dy[dir];
    to_cor.z += dz[dir];
    if (to_cor.out_range())
      return false;
    x = to_cor.x;
    y = to_cor.y;
    z = to_cor.z;
    return true;
  }

  bool next(Point to) {
    if (abs(to.x - x) + abs(to.y - y) + abs(to.z - z) == 1)
      return true;
    return false;
  }

  void print() { cout << x << " " << y << " " << z << endl; }
};

bool operator<(const Point &a, const Point &b) noexcept {
  return (a.x < b.x) || (a.x == b.x && a.y < b.y) ||
         (a.x == b.x && a.y == b.y && a.z <= b.z);
}

struct Part {
  int vol = 1;
  vector<Point> b;
  vector<Point> points;

  void rotate(int shaft) {
    // shaft軸を中心に90°回転
    //(ex z軸) (a, b, c) -> ((a + bi) * (i) = (-b, a, c))
    ll ma = 0;
    for (auto &cu : points) {
      if (shaft == 0) {
        swap(cu.y, cu.z);
        chmax(ma, cu.y);
        cu.y *= -1;
      }
      if (shaft == 1) {
        swap(cu.z, cu.x);
        chmax(ma, cu.z);
        cu.z *= -1;
      }
      if (shaft == 2) {
        swap(cu.x, cu.y);
        chmax(ma, cu.x);
        cu.x *= -1;
      }
    }
    for (auto &cu : points) {
      if (shaft == 0)
        cu.y += ma;
      if (shaft == 1)
        cu.z += ma;
      if (shaft == 2)
        cu.x += ma;
    }
    Sort(points);
  }

  void shrink(vector<Part> &parts) {
    // pointsをとりあえず元に戻す
    points.clear();
    for (Point cu : b)
      points.push_back(cu);

    vi ava_idx;
    rep(i, 0, vol) {
      // そのブロックを消した時シルエットの条件を満たすか
      bool flag1 = 0, flag2 = 0;
      for (Part part : parts) {
        for (Point cor : part.b) {
          if (cor.x == b[i].x && cor.z == b[i].z && cor.y != b[i].y)
            flag1 = 1;
          if (cor.y == b[i].y && cor.z == b[i].z && cor.x != b[i].x)
            flag2 = 1;
        }
      }
      if (!flag1 || !flag2)
        continue;
      // そのブロックを消した時分離しないか
      queue<int> qu;
      qu.push((i + 1) % vol);
      vi seen(vol, 0);
      seen[(i + 1) % vol] = 1;
      while (!qu.empty()) {
        int cu = qu.front();
        qu.pop();
        rep(j, 0, vol) {
          if (seen[j] || i == j)
            continue;
          int dist = 0;
          dist += abs(points[cu].x - points[j].x);
          dist += abs(points[cu].y - points[j].y);
          dist += abs(points[cu].z - points[j].z);
          if (dist == 1) {
            qu.push(j);
            seen[j] = 1;
          }
        }
      }
      int flag = 1;
      rep(j, 0, vol) if (i != j && seen[j] == 0) flag = 0;
      if (flag)
        ava_idx.push_back(i);
    }
    if (ava_idx.size() == 0)
      return;

    int idx = ava_idx[rand() % ava_idx.size()];
    vol--;
    b.erase(b.begin() + idx);
    points.erase(points.begin() + idx);
    Sort(points);
  }

  void extend(Cube &nec_c, Cube &sub_c, vector<Part> &parts) {
    // pointsをとりあえず元に戻す
    points.clear();
    for (Point cu : b)
      points.push_back(cu);

    vector<Point> ava_point;
    Cube used_c;
    for (Part cu : parts) {
      for (Point cu_p : cu.b)
        used_c.b[cu_p.x][cu_p.y][cu_p.z] = 1;
    }
    for (Point cu : b) {
      rep(dir, 0, 6) {
        Point to(cu.x, cu.y, cu.z);
        if (!to.to(dir))
          continue;
        if (!nec_c.b[to.x][to.y][to.z])
          continue;
        if (used_c.b[to.x][to.y][to.z])
          continue;
        ava_point.push_back(to);
      }
    }
    if (ava_point.size() == 0)
      return;
    Point e_p = ava_point[rand() % ava_point.size()];
    vol++;
    points.push_back(e_p);
    b.push_back(e_p);
    Sort(points);
  }
};

void input_silhouette(vvvi &f, vvvi &r) {
  rep(i, 0, D) {
    string s;
    cin >> s;
    rep(j, 0, D) f[0][i][j] = s[j] - '0';
  }
  rep(i, 0, D) {
    string s;
    cin >> s;
    rep(j, 0, D) r[0][i][j] = s[j] - '0';
  }
  rep(i, 0, D) {
    string s;
    cin >> s;
    rep(j, 0, D) f[1][i][j] = s[j] - '0';
  }
  rep(i, 0, D) {
    string s;
    cin >> s;
    rep(j, 0, D) r[1][i][j] = s[j] - '0';
  }
}

void init_all_parts(Cube nec_c, Cube sub_c, vector<Part> &parts) {
  Cube fix_c;
  rep(x, 0, D) rep(y, 0, D) rep(z, 0, D) {
    if (sub_c.b[x][y][z] == 0 || fix_c.b[x][y][z] == 1)
      continue;
    bool flag1 = 0, flag2 = 0;
    rep(_x, 0, D) if (_x != x && fix_c.b[_x][y][z]) flag1 = 1;
    rep(_y, 0, D) if (_y != y && fix_c.b[x][_y][z]) flag2 = 1;
    if (flag1 && flag2)
      continue;

    vector<Point> candi;
    rep(dir, 0, 6) {
      Point to(x, y, z);
      if (!to.to(dir))
        continue;
      if (fix_c.b[to.x][to.y][to.z])
        continue;
      if (!nec_c.b[to.x][to.y][to.z])
        continue;
      candi.push_back(to);
    }

    fix_c.b[x][y][z] = 1;
    Part cu_part;
    cu_part.b.push_back(Point(x, y, z));
    cu_part.points.push_back(Point(x, y, z));
    if (candi.empty()) {
      parts.push_back(cu_part);
      continue;
    }

    int loop_num = 0;
    int upper = rand() % 5;
    while (!candi.empty()) {
      if (loop_num == upper)
        break;
      int idx = rand() % (int)candi.size();
      Point cu;
      cu = candi[idx];
      candi.erase(candi.begin() + idx);
      if (fix_c.b[cu.x][cu.y][cu.z])
        continue;

      fix_c.b[cu.x][cu.y][cu.z] = 1;
      cu_part.b.push_back(cu);
      cu_part.points.push_back(cu);
      cu_part.vol++;
      loop_num++;

      rep(dir, 0, 6) {
        Point to(cu.x, cu.y, cu.z);
        if (!to.to(dir))
          continue;
        if (fix_c.b[to.x][to.y][to.z])
          continue;
        if (!nec_c.b[to.x][to.y][to.z])
          continue;
        candi.push_back(to);
      }
    }

    Sort(cu_part.points);
    parts.push_back(cu_part);
  }
}

bool same_part(Part &part1, Part &part2) {
  if (part1.vol != part2.vol)
    return false;

  rep(x, 0, 4) {
    rep(y, 0, 4) {
      rep(z, 0, 4) {
        bool flag = 1;
        Point p1, p2;
        p1 = part1.points[0];
        p2 = part2.points[0];
        int _dx = p2.x - p1.x;
        int _dy = p2.y - p1.y;
        int _dz = p2.z - p1.z;
        rep(j, 0, part2.vol) {
          Point _p1, _p2;
          _p1 = part1.points[j];
          _p2 = part2.points[j];
          if (_p1.x + _dx != _p2.x || _p1.y + _dy != _p2.y ||
              _p1.z + _dz != _p2.z)
            flag = 0;
        }
        if (flag)
          return true;
        part2.rotate(2);
      }
      part2.rotate(1);
    }
    part2.rotate(0);
  }
  return false;
}

ll calc_score(vector<Part> &parts1, vector<Part> &parts2) {
  ll r = 0;
  long double v = 0;
  rep(i, 0, (ll)parts1.size()) { r += parts1[i].vol; }
  rep(i, 0, (ll)parts2.size()) { r += parts2[i].vol; }

  vi used((int)parts2.size(), 0);
  rep(i, 0, (ll)parts1.size()) {
    rep(j, 0, (ll)parts2.size()) {
      if (used[j])
        continue;
      if (same_part(parts1[i], parts2[j])) {
        r -= parts1[i].vol * 2;
        v += (long double)1 / parts1[i].vol;
        used[j] = 1;
        break;
      }
    }
  }

  return 1e9 * (r + v);
}

ll diff_calc_score(vector<Part> &parts1, vector<Part> &parts2, vi &b_idx1,
                   vi &b_idx2, int type) {
  if (type == 0) {
    rep(i, 0, (ll)parts1.size()) {
      rep(j, 0, (ll)parts2.size()) {
        if (b_idx2[j] != 0)
          continue;
        if (same_part(parts1[i], parts2[j])) {
          b_idx2[j] = i + 1;
          break;
        }
      }
    }
  }
  // parts1が変更された(parts1の1つの要素が0になってる)
  if (type == 1) {
    // parts1[i]と対応するparts2[j]も0にする
    rep(i, 0, (ll)parts1.size()) {
      rep(j, 0, (ll)parts2.size()) {
        if ((int)parts1.size() < b_idx2[j])
          b_idx2[j] = 0;
        if (b_idx1[i] == 0 && b_idx2[j] == i + 1)
          b_idx2[j] = 0;
      }
    }
    rep(i, 0, (ll)parts1.size()) {
      if (b_idx1[i] != 0)
        continue;
      rep(j, 0, (ll)parts2.size()) {
        if (b_idx2[j] != 0)
          continue;
        if (same_part(parts1[i], parts2[j])) {
          b_idx2[j] = i + 1;
          break;
        }
      }
    }
  }
  // parts2が変更された(parts2の1つの要素が0になってる)
  if (type == 2) {
    vi used((int)parts1.size() + 1, 0);
    rep(i, 0, (ll)parts2.size()) {
      if (b_idx2[i] <= (int)parts1.size())
        used[b_idx2[i]] = 1;
    }
    rep(i, 0, (ll)parts1.size()) {
      rep(j, 0, (ll)parts2.size()) {
        if (b_idx2[j] != 0 || used[i + 1] == 1)
          continue;
        if (same_part(parts1[i], parts2[j])) {
          b_idx2[j] = i + 1;
          break;
        }
      }
    }
    rep(j, 0, (ll)parts2.size()) {
      if ((int)parts1.size() < b_idx2[j])
        b_idx2[j] = 0;
    }
  }

  if (0 <= type && type <= 2) {
    rep(i, 0, (ll)parts1.size()) { b_idx1[i] = i + 1; }
    int val = (int)parts1.size() + 1;
    rep(j, 0, (ll)parts2.size()) {
      if (b_idx2[j] == 0) {
        b_idx2[j] = val;
        val++;
      }
    }
  }
  ll r = 0;
  long double v = 0;
  rep(i, 0, (ll)parts1.size()) { r += parts1[i].vol; }
  rep(i, 0, (ll)parts2.size()) { r += parts2[i].vol; }
  rep(i, 0, (ll)parts1.size()) {
    rep(j, 0, (ll)parts2.size()) {
      if (b_idx1[i] == b_idx2[j]) {
        r -= parts1[i].vol * 2;
        v += (long double)1 / parts1[i].vol;
      }
    }
  }
  return 1e9 * (r + v);
}

void modify(Cube &nec_c, Cube &sub_c, vector<Part> &parts, vi &b_idx, int idx) {
  if (rand() % 2) {
    parts[idx].extend(nec_c, sub_c, parts);
  } else {
    if (parts[idx].vol == 1)
      return;
    parts[idx].shrink(parts);
    if (parts[idx].vol == 0)
      parts.erase(parts.begin() + idx);
  }
  b_idx[idx] = 0;
}

bool modify_break(Cube &nec_c, Cube &sub_c, vector<Part> &parts,
                  int break_idx) {
  vector<Part> new_parts;
  new_parts = parts;
  // int break_idx = rand() % (int)new_parts.size();

  new_parts.erase(new_parts.begin() + break_idx);
  Cube fix_c;
  for (Part cu : new_parts) {
    for (Point cu_p : cu.b)
      fix_c.b[cu_p.x][cu_p.y][cu_p.z] = 1;
  }

  while (1) {
    int flag = 1;
    int i = rand() % (int)new_parts.size();
    rep(_i, 0, (int)new_parts.size()) {
      i = (i + 1) % (int)new_parts.size();
      queue<Point> qu;
      int upper = rand() % 3 + 1;
      for (Point _p : new_parts[i].b)
        qu.push(_p);
      while (!qu.empty()) {
        rep(dir, 0, 6) {
          if (upper == 0)
            break;
          Point to;
          to = qu.front();
          if (!to.to(dir))
            continue;
          if (fix_c.b[to.x][to.y][to.z])
            continue;
          if (!nec_c.b[to.x][to.y][to.z])
            continue;
          fix_c.b[to.x][to.y][to.z] = 1;
          new_parts[i].b.push_back(to);
          new_parts[i].points.push_back(to);
          new_parts[i].vol++;
          upper--;
          flag = 0;
          qu.push(to);
        }
        qu.pop();
      }
      Sort(new_parts[i].points);
    }
    if (flag)
      break;
  }
  rep(x, 0, D) rep(z, 0, D) {
    int flag = 0;
    rep(_y, 0, D) if (nec_c.b[x][_y][z]) flag = 1;
    if (flag) {
      flag = 0;
      rep(_y, 0, D) if (fix_c.b[x][_y][z]) flag = 1;
      if (!flag)
        return false;
    }
  }
  rep(y, 0, D) rep(z, 0, D) {
    int flag = 0;
    rep(_x, 0, D) if (nec_c.b[_x][y][z]) flag = 1;
    if (flag) {
      flag = 0;
      rep(_x, 0, D) if (fix_c.b[_x][y][z]) flag = 1;
      if (!flag)
        return false;
    }
  }

  parts = new_parts;
  return true;
}

void print(vector<Part> parts1, vector<Part> parts2) {
  vvi ans(2, vi(D * D * D, 0));
  ll idx = 1;
  rep(i, 0, (ll)parts1.size()) {
    Part cu = parts1[i];
    for (Point cu_p : cu.b) {
      ans[0][cu_p.x * D * D + cu_p.y * D + cu_p.z] = idx;
    }
    idx++;
  }
  vi used_idx((int)parts1.size() + 3, 0);
  rep(i, 0, (ll)parts2.size()) {
    int same_idx = -1;
    rep(j, 0, (ll)parts1.size()) {
      if (used_idx[j + 1])
        continue;
      if (same_part(parts1[j], parts2[i])) {
        same_idx = j + 1;
        used_idx[j + 1] = 1;
        break;
      }
    }
    if (same_idx == -1) {
      same_idx = idx;
      idx++;
    }
    Part cu = parts2[i];
    for (Point cu_p : cu.b) {
      ans[1][cu_p.x * D * D + cu_p.y * D + cu_p.z] = same_idx;
    }
  }

  ll sum = 0;
  rep(i, 0, 2) rep(j, 0, D * D * D) chmax(sum, ans[i][j]);
  cout << sum << endl;
  for (ll cu : ans[0])
    cout << cu << " ";
  cout << endl;
  for (ll cu : ans[1])
    cout << cu << " ";
  cout << endl;
}

void solve(Cube &nec_c1, Cube &sub_c1, Cube &nec_c2, Cube &sub_c2,
           vector<Part> &parts1, vector<Part> &parts2) {
  auto startClock = system_clock::now();
  vector<Part> all_best_parts1;
  vector<Part> all_best_parts2;
  ll all_best_score = 1e15;
  // 初期解
  init_all_parts(nec_c1, sub_c1, parts1);
  init_all_parts(nec_c2, sub_c2, parts2);

  // 壊す・修復するの焼きなまし
  while (1) {
    if (TIME_LIMIT < ela_times(startClock))
      break;
    // 壊す
    vi _b_idx1((int)parts1.size(), 0);
    vi _b_idx2((int)parts2.size(), 0);
    diff_calc_score(parts1, parts2, _b_idx1, _b_idx2, 0);
    if (parts1.size() >= parts2.size()) {
      vi ava_idx;
      vi seen((int)parts1.size(), 0);
      for (int idx : _b_idx2)
        if (idx <= (int)parts1.size())
          seen[idx - 1] = 1;
      int i = rand() % (int)parts1.size();
      rep(_i, 0, (int)parts1.size()) {
        i = (i + 1) % (int)parts1.size();
        if (!seen[i])
          ava_idx.push_back(i);
      }
      i = rand() % (int)parts1.size();
      rep(_i, 0, (int)parts1.size()) {
        i = (i + 1) % (int)parts1.size();
        if (seen[i])
          ava_idx.push_back(i);
      }
      rep(i, 0, (int)parts1.size()) {
        if (modify_break(nec_c1, sub_c1, parts1, ava_idx[i]))
          break;
      }
    }
    if (TIME_LIMIT < ela_times(startClock))
      break;
    if (parts1.size() < parts2.size()) {
      vi ava_idx;
      int i = rand() % (int)parts2.size();
      rep(_i, 0, (int)parts2.size()) {
        i = (i + 1) % (int)parts2.size();
        if (_b_idx2[i] > (int)parts1.size())
          ava_idx.push_back(i);
      }
      i = rand() % (int)parts2.size();
      rep(_i, 0, (int)parts2.size()) {
        i = (i + 1) % (int)parts2.size();
        if (_b_idx2[i] <= (int)parts1.size())
          ava_idx.push_back(i);
      }
      rep(i, 0, (int)parts2.size()) {
        if (modify_break(nec_c2, sub_c2, parts2, ava_idx[i]))
          break;
      }
    }

    // 修復する(同形パーツ固定)
    vi b_idx1((int)parts1.size(), 0);
    vi b_idx2((int)parts2.size(), 0);
    ll best_score = diff_calc_score(parts1, parts2, b_idx1, b_idx2, 0);
    print_num(best_score);
    cerr << " ";
    int loop_num = 0;
    int pre_update = -1;
    auto _startClock = system_clock::now();
    while (1) {
      if (0.15 * (D - 4) < ela_times(_startClock))
        break;
      loop_num++;
      if (rand() % 2) {
        vi seen((int)parts1.size(), 0);
        rep(i, 0, (int)parts2.size()) {
          if (b_idx2[i] <= (int)parts1.size())
            seen[b_idx2[i] - 1] = 1;
        }
        vi ava_idx;
        rep(i, 0, (int)parts1.size()) {
          if (seen[i] == 0)
            ava_idx.push_back(i);
        }
        if (ava_idx.empty())
          continue;
        int idx = ava_idx[rand() % (int)ava_idx.size()];
        modify(nec_c1, sub_c1, parts1, b_idx1, idx);
        best_score = diff_calc_score(parts1, parts2, b_idx1, b_idx2, 1);
      } else {
        vi ava_idx;
        rep(i, 0, (int)parts2.size()) {
          if (b_idx2[i] > (int)parts1.size())
            ava_idx.push_back(i);
        }
        if (ava_idx.empty())
          continue;
        int idx = ava_idx[rand() % (int)ava_idx.size()];
        modify(nec_c2, sub_c2, parts2, b_idx2, idx);
        best_score = diff_calc_score(parts1, parts2, b_idx1, b_idx2, 2);
      }
    }

    print_num(best_score);
    cerr << " ";
    // 共通パーツ数
    int _sum = 0;
    rep(i, 0, (int)parts2.size()) {
      if (b_idx2[i] <= (int)parts1.size())
        _sum++;
    }
    cerr << _sum << " ";
    // 修復する(スコア)
    loop_num = 0;
    pre_update = -1;
    while (1) {
      loop_num++;
      if (TIME_LIMIT < ela_times(startClock))
        break;
      if (loop_num - pre_update > 100)
        break;

      vi _b_idx1((int)parts1.size(), -1);
      vi _b_idx2((int)parts2.size(), -1);
      _b_idx1 = b_idx1;
      _b_idx2 = b_idx2;
      if (rand() % 2) {
        vector<Part> _parts1 = parts1;
        auto _startClock = system_clock::now();
        modify(nec_c1, sub_c1, _parts1, _b_idx1);
        ll score = diff_calc_score(_parts1, parts2, _b_idx1, _b_idx2, 1);
        stat_time += ela_times(_startClock);

        if (chmin(best_score, score)) {
          pre_update = loop_num;
          parts1 = _parts1;
          b_idx1 = _b_idx1;
          b_idx2 = _b_idx2;
        }
      } else {
        vector<Part> _parts2 = parts2;
        auto _startClock = system_clock::now();
        modify(nec_c2, sub_c2, _parts2, _b_idx2);
        ll score = diff_calc_score(parts1, _parts2, _b_idx1, _b_idx2, 2);
        stat_time += ela_times(_startClock);

        if (chmin(best_score, score)) {
          pre_update = loop_num;
          parts2 = _parts2;
          b_idx1 = _b_idx1;
          b_idx2 = _b_idx2;
        }
      }
    }

    // 共通パーツが延ばせないか山登り
    loop_num = 0;
    pre_update = -1;
    while (1) {
      loop_num++;
      if (TIME_LIMIT < ela_times(startClock))
        break;
      if (loop_num - pre_update > 20)
        break;
      int idx1 = rand() % (int)parts1.size();
      int idx2 = -1;
      rep(i, 0, (int)parts2.size()) if (idx1 + 1 == b_idx2[i]) idx2 = i;
      if (idx2 == -1)
        continue;
      vector<Part> _parts1, _parts2;
      _parts1 = parts1;
      _parts2 = parts2;
      _parts1[idx1].extend(nec_c1, sub_c1, _parts1);
      _parts2[idx2].extend(nec_c2, sub_c2, _parts2);
      if (same_part(_parts1[idx1], _parts2[idx2])) {
        pre_update = loop_num;
        parts1 = _parts1;
        parts2 = _parts2;
      }
    }
    best_score = calc_score(parts1, parts2);
    print_num(calc_score(parts1, parts2));
    cerr << " ";
    // 成功
    if (chmin(all_best_score, best_score)) {
      all_best_parts1 = parts1;
      all_best_parts2 = parts2;
    } else {
      parts1 = all_best_parts1;
      parts2 = all_best_parts2;
    }
    print_num(calc_score(parts1, parts2));
    cerr << " " << (int)parts1.size() << " " << (int)parts2.size() << " ";
    // 共通パーツ数
    _sum = 0;
    rep(i, 0, (int)parts2.size()) {
      if (b_idx2[i] <= (int)parts1.size())
        _sum++;
    }
    cerr << _sum << " ";
    cerr << endl;
  }
  parts1 = all_best_parts1;
  parts2 = all_best_parts2;
}

int main() {
  srand((unsigned int)time(NULL));
  cin >> D;
  vvvi f(2, vvi(D, vi(D, 0)));
  vvvi r(2, vvi(D, vi(D, 0)));
  input_silhouette(f, r);

  Cube nec_c1(f[0], r[0]), nec_c2(f[1], r[1]);
  Cube sub_c1(f[1], r[1], nec_c1), sub_c2(f[1], r[1], nec_c2);

  vector<Part> parts1, parts2;
  solve(nec_c1, sub_c1, nec_c2, sub_c2, parts1, parts2);

  print_num(calc_score(parts1, parts2));
  cerr << endl;
  print(parts1, parts2);
  cerr << D << endl;
  cerr << (int)parts1.size() << " " << (int)parts2.size() << endl;
  cerr << stat_time << endl;
}

// local6
// 基本はlocal4
// 修復パート(同形固定) -> 修復パート(スコア)
// だめそう