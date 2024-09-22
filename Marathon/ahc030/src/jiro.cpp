#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/distributions/normal.hpp>
#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef long double ld;
#define FOR(i, a, b) for (int i = a; i < (b); i++)
#define RFOR(i, a, b) for (int i = a; i >= (b); i--)
#define range(a) a.begin(), a.end()
#define Yes() cout << "Yes" << endl
#define No() cout << "No" << endl
#define MP make_pair
int dx[8] = {-1, 0, 1, 0, -1, -1, 1, 1};
int dy[8] = {0, -1, 0, 1, -1, 1, -1, 1};
using P = pair<int, int>;
const long long INF = 1LL << 60;
void chmin(long long& a, long long b) {
  if (a > b) a = b;
}
void chmax(long long& a, long long b) {
  if (a < b) a = b;
}

struct xorshift128 {
  const unsigned int ini_x = 123456789, ini_y = 362436069, ini_z = 521288629, ini_w = 88675123;
  unsigned int x, y, z, w;

  xorshift128() {}

  // シードによってx,y,z,wを初期化 ... initialize x,y,z,w by seed
  void set_seed(unsigned int seed) {
    x = ini_x, y = ini_y, z = ini_z, w = ini_w ^ seed;
  }

  unsigned int randint() {
    unsigned int t = x ^ (x << 11);
    x = y;
    y = z;
    z = w;
    return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
  }

  // [0,r)の範囲の整数で乱数発生 ... generate random integer in [0,r)
  unsigned int randint(unsigned int r) {
    assert(r != 0);
    return randint() % r;
  }

  // [l,r)の範囲の整数で乱数発生 ... generate random integer in [l,r)
  unsigned int randint(unsigned int l, unsigned int r) {
    assert(r > l);
    return l + randint(r - l);
  }

  // 長さnの順列をランダムに生成し、その前k個分を出力する ... generate a random permutation of size n, and return the first k
  vector<int> randperm(int n, int k) {
    assert(k >= 0 && k <= n);
    vector<int> ret(n);
    for (int i = 0; i < n; i++) {
      ret[i] = i;
    }
    for (int i = 0; i < k; i++) {
      swap(ret[i], ret[randint(i, n)]);
    }
    return vector<int>(ret.begin(), ret.begin() + k);
  }

  // [0,1]の範囲の実数で乱数発生 ... generate random real number in [0,1]
  double uniform() {
    return double(randint()) * 2.3283064370807974e-10;
  }

  // [0,r]の範囲の実数で乱数発生 ... generate random real number in [0,r]
  double uniform(double r) {
    assert(r >= 0.0);
    return uniform() * r;
  }

  // [l,r]の範囲の実数で乱数発生 ... generate random real number in [l,r]
  double uniform(double l, double r) {
    assert(r >= l);
    return l + uniform(r - l);
  }
};

xorshift128 Random;

const int64_t CYCLES_PER_SEC = 2800000000;

struct Timer {
  int64_t start;
  Timer() {
    reset();
  }
  void reset() {
    start = getCycle();
  }
  void plus(double a) {
    start -= (a * CYCLES_PER_SEC);
  }
  inline double get() {
    return (double)(getCycle() - start) / CYCLES_PER_SEC;
  }
  inline int64_t getCycle() {
    uint32_t low, high;
    __asm__ volatile("rdtsc" : "=a"(low), "=d"(high));
    return ((int64_t)low) | ((int64_t)high << 32);
  }
};

Timer timer;

// UnionFind  ex.)UnionFind uf(N);  という風に定義する
// coding: https://youtu.be/TdR816rqc3s?t=726
// comment: https://youtu.be/TdR816rqc3s?t=6822
struct UnionFind {
  vector<int> d;
  UnionFind(int n = 0) : d(n, -1) {}  // コンストラクタ(n:頂点数)
  int find(int x) {
    if (d[x] < 0) return x;
    return d[x] = find(d[x]);
  }
  bool unite(int x, int y) {
    x = find(x);
    y = find(y);
    if (x == y) return false;
    if (d[x] > d[y]) swap(x, y);
    d[x] += d[y];
    d[y] = x;
    return true;
  }
  bool same(int x, int y) {
    return find(x) == find(y);
  }
  int size(int x) {
    return -d[find(x)];
  }
};

int N, M, restQuery, oilNum;
double eps;
double cost = 0.0;
vector<int> statesNum;
bitset<512> validRange;
ll hashTable[20][400];

std::string floatToColorCode(float value) {
  assert(value >= 0.0 && value <= M);
  int r, g, b;
  double k = 1.0 - log(value + 1.0) / log(M + 1.0);
  // cerr << "k: " << k << endl;

  b = static_cast<int>(k * 255);
  r = 255;
  g = static_cast<int>(k * 255);

  // RGB値を16進数のカラーコードに変換
  std::stringstream ss;
  ss << "#" << std::setfill('0') << std::setw(2) << std::hex << r << std::setw(2) << std::hex << g << std::setw(2) << std::hex << b;

  return ss.str();
}

// 正規分布の確率密度関数を用いて尤度を計算する関数
double normalLikelihood(double mean, double variance, double x) {
  double eps = 1e-8;                   // 分散が0に近い値を避けるための微小値
  if (variance < eps) variance = eps;  // 分散が0または非常に小さい値であれば、微小値に設定

  double sigma = sqrt(variance);  // 標準偏差は分散の平方根
  double coeff = 1.0 / (sigma * sqrt(2.0 * M_PI));
  double exponent = exp(-0.5 * pow((x - mean) / sigma, 2));
  return coeff * exponent;
}

double logNormalLikelihood(double mean, double variance, double x) {
  double eps = 1e-8;                   // 分散が0に近い値を避けるための微小値
  if (variance < eps) variance = eps;  // 分散が0または非常に小さい値であれば、微小値に設定

  double sigma = sqrt(variance);                             // 標準偏差は分散の平方根
  double log_coeff = -log(sigma) - 0.5 * log(2.0 * M_PI);    // 対数係数部分
  double exponent_part = -0.5 * pow((x - mean) / sigma, 2);  // 指数関数部分の対数
  return log_coeff + exponent_part;
}

// 確率密度関数
double normal_pdf(double mu, double x, double sigma) {
  boost::math::normal_distribution<> dist(mu, sigma);
  return boost::math::pdf(dist, x);
}

// 積分関数
double log_integrate_normal_pdf(double a, double b, double x, double sigma) {
  auto f = [x, sigma](double mu) { return normal_pdf(mu, x, sigma); };
  boost::math::quadrature::gauss<double, 15> integrator;
  double result = integrator.integrate(f, a, b);  // 結果を直接受け取る
  return log(result);
}

// 対数加算を行う関数
double logAdd(double log_x, double log_y) {
  // log_xとlog_yのうち大きい方を保持
  double max_log = std::max(log_x, log_y);
  // 数値が非常に小さい時のアンダーフローを防ぐ
  if (max_log == -std::numeric_limits<double>::infinity()) {
    return max_log;
  }
  return max_log + log(exp(log_x - max_log) + exp(log_y - max_log));
}

struct OilBasic {
  int id;
  int d;
  int xsize, ysize;
  bitset<512> position;
  OilBasic() {}
  OilBasic(int id, int d, bitset<512> position, int xmin, int xmax, int ymin, int ymax) : id(id), d(d), position(position) {
    xsize = xmax - xmin + 1;
    ysize = ymax - ymin + 1;
  }
  vector<pair<int, int>> toPairs() {
    vector<pair<int, int>> ret;
    FOR(i, 0, N) FOR(j, 0, N) if (position.test(i * N + j)) ret.push_back({i, j});
    return ret;
  }
  int cntPosition() {
    return (N - xsize + 1) * (N - ysize + 1);
  }
  bitset<512> getPararellPosition(int delX, int delY) {
    int delta = delX * N + delY;
    return position << delta;
  }
};

int postQuery(int t, bitset<512> coords) {
  assert(restQuery > 0);
  assert(t == 1 || t == 2 || t == 3);
  int res;
  int size = coords.count();
  assert(size > 0 && size <= N * N);
  assert(size == coords.count());
  if (t == 1) {
    assert(size == 1);
    FOR(i, 0, N) FOR(j, 0, N) if (coords.test(i * N + j)) {
      cout << "q " << 1 << " " << i << " " << j << endl;
    }
    cost += 1.0;
    cin >> res;
  } else if (t == 2) {
    cout << "q " << size << " ";
    int cnt = 0;
    FOR(i, 0, N) FOR(j, 0, N) if (coords.test(i * N + j)) {
      cout << i << " " << j << " ";
      ++cnt;
    }
    assert(cnt == coords.count());
    cout << endl;
    cost += 1.0 / sqrt(size);
    cin >> res;
  } else {
    cout << "a " << size << " ";
    FOR(i, 0, N) FOR(j, 0, N) if (coords.test(i * N + j)) {
      cout << i << " " << j << " ";
    }
    cout << endl;
    cost += 1.0 / sqrt(size);
    cin >> res;
  }
  restQuery--;
  return res;
}

struct Oil : OilBasic {
  int idx, dx, dy;
  bitset<512> places;
  Oil() {}
  Oil(const OilBasic& oilBasic, int idx) : OilBasic(oilBasic), idx(idx) {
    dx = idx % (N - xsize + 1), dy = idx / (N - xsize + 1);
    places = getPararellPosition(dx, dy);
  }
  Oil(const OilBasic& oilBasic, int dx, int dy) : OilBasic(oilBasic), dx(dx), dy(dy) {
    assert(dx >= 0 && dx <= N - xsize);
    assert(dy >= 0 && dy <= N - ysize);
    idx = dx + dy * (N - xsize + 1);
    places = getPararellPosition(dx, dy);
  }
  bool canMove(int deltaX, int deltaY) {
    if (dx + deltaX < 0 || dx + deltaX > N - xsize) return false;
    if (dy + deltaY < 0 || dy + deltaY > N - ysize) return false;
    return true;
  }
  pair<int, int> getLeftTop() {
    return {dx, dy};
  }
  pair<int, int> getRightBottom() {
    return {dx + xsize - 1, dy + ysize - 1};
  }
  Oil getMoveOil(int deltaX, int deltaY) {
    assert(canMove(deltaX, deltaY));
    return Oil(extractOilBasic(), dx + deltaX, dy + deltaY);
  }
  OilBasic extractOilBasic() const {
    return *this;
  }
  int overlapNumUsingCoords(bitset<512> coords) {
    return (coords & places).count();
  }
  void debug() {
    FOR(i, 0, N) {
      FOR(j, 0, N) {
        if (places.test(i * N + j)) cerr << "#";
        else cerr << ".";
      }
      cerr << endl;
    }
    cerr << endl;
  }
};

struct OilState {
  int id, num;
  bitset<512> oilMap[25];
  vector<Oil> oils;
  ll hashNum;
  OilState() {}
  OilState(int id, vector<Oil> oils) : id(id), oils(oils) {
    int idx = 0;
    for (auto oil : oils) {
      oilMap[idx] = oil.places;
      idx++;
    }
    num = idx;
    FOR(i, 0, M) hashNum ^= hashTable[i][oils[i].idx];
  }
  bitset<512> oilsPosition() {
    bitset<512> ret;
    ret.reset();
    FOR(i, 0, num) ret |= oilMap[i];
    return ret;
  }
  int overlapNumUsingCoords(bitset<512> coords) {
    int ret = 0;
    FOR(i, 0, num) {
      ret += (coords & oilMap[i]).count();
    }
    return ret;
  }
  void changeOil(int cidx, Oil newOil) {
    hashNum ^= hashTable[cidx][oils[cidx].idx];
    oils[cidx] = newOil;
    oilMap[cidx] = newOil.places;
    hashNum ^= hashTable[cidx][oils[cidx].idx];
    return;
  }
  void swapOil(int a, int b) {
    auto Atopleft = oils[a].getLeftTop();
    auto Btopleft = oils[b].getLeftTop();
    auto Abottomright = oils[a].getRightBottom();
    auto Bbottomright = oils[b].getRightBottom();

    int newAX = Btopleft.first + (Bbottomright.first - Btopleft.first + 1) - (Abottomright.first - Atopleft.first + 1);
    int newAY = Btopleft.second + (Bbottomright.second - Btopleft.second + 1) - (Abottomright.second - Atopleft.second + 1);
    int newBX = Atopleft.first + (Abottomright.first - Atopleft.first + 1) - (Bbottomright.first - Btopleft.first + 1);
    int newBY = Atopleft.second + (Abottomright.second - Atopleft.second + 1) - (Bbottomright.second - Btopleft.second + 1);

    newAX += Random.randint(3) - 1;
    newAY += Random.randint(3) - 1;
    newBX += Random.randint(3) - 1;
    newBY += Random.randint(3) - 1;

    newAX = max(0, newAX);
    newAX = min(N - oils[a].xsize, newAX);
    newAY = max(0, newAY);
    newAY = min(N - oils[a].ysize, newAY);
    newBX = max(0, newBX);
    newBX = min(N - oils[b].xsize, newBX);
    newBY = max(0, newBY);
    newBY = min(N - oils[b].ysize, newBY);

    auto newOilA = Oil(oils[a].extractOilBasic(), newAX, newAY);
    auto newOilB = Oil(oils[b].extractOilBasic(), newBX, newBY);

    changeOil(a, newOilA);
    changeOil(b, newOilB);

    return;
  }
  bool isMovePerOil(int id, int deltaX, int deltaY) {
    return oils[id].canMove(deltaX, deltaY);
  }
  Oil getMoveOil(int m, int dx, int dy) {
    return oils[m].getMoveOil(dx, dy);
  }
  int includeNum(pair<int, int> crd) {
    int ret = 0;
    FOR(i, 0, M) if (oilMap[i].test(crd.first * N + crd.second)) ret++;
    return ret;
  }
  void debug() {
    bitset<512> tmp = oilsPosition();
    FOR(i, 0, N) {
      FOR(j, 0, N) {
        if (tmp.test(i * N + j)) cerr << "#";
        else cerr << ".";
      }
      cerr << endl;
    }
  }
  void debugColor() {
    bitset<512> tmp;
    FOR(i, 0, num) tmp |= oilMap[i];
    FOR(i, 0, N) FOR(j, 0, N) {
      cout << "#c " << i << " " << j << " " << floatToColorCode(tmp.test(i * N + j)) << endl;
    }
    return;
  }
};

struct Game1 {
  ll stateALL;
  priority_queue<pair<double, int>> pos;
  vector<OilState> oilStates;
  double AnswerDecide = log(10);
  double DisposeDecide = log(1e5);
  Game1(vector<OilBasic> oils) {
    stateALL = 1;
    FOR(i, 0, M) stateALL *= statesNum[i];
    builtFromScratch(oils, stateALL);
  }

  void builtFromScratch(vector<OilBasic> oils, ll stateALL) {
    function<void(vector<int>)> dfs = [&](vector<int> idxs) {
      if ((int)idxs.size() == M) {
        vector<Oil> tmp;
        FOR(m, 0, M) tmp.push_back(Oil(oils[m], idxs[m]));
        oilStates.push_back(OilState((int)oilStates.size(), tmp));
        pos.push({-log(stateALL), (int)oilStates.size() - 1});
        return;
      }
      int size = (int)idxs.size();
      FOR(i, 0, statesNum[size]) {
        idxs.push_back(i);
        dfs(idxs);
        idxs.pop_back();
      }
      return;
    };
    dfs(vector<int>());
  }

  // posが大きい2つの確率を返す
  pair<int, int> getIndexMaxPos() {
    assert((int)pos.size() >= 2);
    auto maxpos1 = pos.top();
    pos.pop();
    auto maxpos2 = pos.top();
    pos.push(maxpos1);
    return {maxpos1.first, maxpos2.first};
  }

  // 確率を更新する共通関数
  void updatePossibilityUsingCoords(bitset<512> coords, int x) {
    int k = coords.count();
    double log_evidence = -numeric_limits<double>::infinity();
    double log_maxPos = -numeric_limits<double>::infinity();
    queue<pair<double, int>> tmp;
    while (!pos.empty()) {
      auto [pp, id] = pos.top();
      pos.pop();
      auto num = oilStates[id].overlapNumUsingCoords(coords);
      double avg = (k - num) * eps + num * (1.0 - eps);
      double var = k * eps * (1.0 - eps);
      double log_likelihood = logNormalLikelihood(avg, var, x);
      pp += log_likelihood;
      log_evidence = logAdd(log_evidence, pp);
      log_maxPos = max(log_maxPos, pp);
      tmp.push({pp, id});
    }
    log_maxPos -= log_evidence;
    while (!tmp.empty()) {
      auto [pp, id] = tmp.front();
      tmp.pop();
      pp -= log_evidence;
      if (log_maxPos - pp > DisposeDecide) continue;
      pos.push({pp, id});
    }
    return;
  }

  // 最大の解をpostQueryして確率を更新
  void postQueryMaxPosState() {
    int id;
    double pp;
    assert((int)pos.size() >= 2);
    if (Random.uniform() > 0.1) {
      auto ppstr = pos.top();
      pp = ppstr.first;
      id = ppstr.second;
    } else {
      auto maxppstr = pos.top();
      pos.pop();
      auto ppstr = pos.top();
      pos.push(maxppstr);
      pp = ppstr.first;
      id = ppstr.second;
    }
    auto coords = oilStates[id].oilsPosition();
    auto res = postQuery(2, coords);
    updatePossibilityUsingCoords(coords, res);
    return;
  }

  // 確率が十分に大きいか判定
  bool judgeAnswer() {
    if ((int)pos.size() < 2) return true;
    auto [maxValue1, maxValue2] = getIndexMaxPos();
    return maxValue1 - maxValue2 > AnswerDecide;
  }

  // Answerをpost
  bool postQueryAnswer() {
    auto [maxValue, maxId] = pos.top();
    auto coords = oilStates[maxId].oilsPosition();
    auto res = postQuery(3, coords);
    if (res == 0) pos.pop();
    return res == 1;
  }

  // 一番確率が高い状態をデバッグ
  void debugColorMaxPos() {
    auto [maxValue, maxId] = pos.top();
    auto maxState = oilStates[maxId];
    // cerr << "max: " << maxState.pos << endl;
    maxState.debugColor();
    return;
  }
};

struct Game2 {
  int oils[40][40];
  void solve() {
    FOR(i, 0, N) FOR(j, 0, N) {
      bitset<512> coords;
      coords.reset();
      coords.set(i * N + j);
      auto res = postQuery(2, coords);
      oils[i][j] = res;
    }
  }
  void output() {
    bitset<512> coords;
    coords.reset();
    FOR(i, 0, N) FOR(j, 0, N) {
      if (oils[i][j] > 0) coords.set(i * N + j);
    }
    auto res = postQuery(3, coords);
    assert(res == 1);
  }
};

struct Game3 {
  bool possible[20][400];
  vector<OilBasic> oils;
  vector<int> includes[20][400];
  vector<int> oi;
  Game3(vector<OilBasic> oils) : oils(oils) {
    FOR(i, 0, M) FOR(j, 0, statesNum[i]) possible[i][j] = true;
  }
  vector<pair<int, int>> getMustPostPlaces() {
    set<pair<int, int>> coords;
    FOR(i, 0, M) {
      auto oilBasic = oils[i];
      auto vec = oilBasic.toPairs();
      for (auto [x1, y1] : vec) {
        for (auto [x2, y2] : vec) {
          coords.insert({x2 - x1, y2 - y1});
        }
      }
    }

    vector<vector<bool>> field(N, vector<bool>(N, false));
    auto setField = [&](int x, int y) {
      for (auto crd : coords) {
        int nx = x + crd.first, ny = y + crd.second;
        if (nx < 0 || nx >= N || ny < 0 || ny >= N) continue;
        field[nx][ny] = true;
      }
    };
    auto findNextPost = [&]() {
      int minDist = 1000;
      pair<int, int> ret = {-1, -1};
      FOR(i, 0, N) FOR(j, 0, N) {
        if (field[i][j]) continue;
        int dist = abs(i - N / 2) + abs(j - N / 2);
        if (dist < minDist) {
          minDist = dist;
          ret = {i, j};
        }
      }
      return ret;
    };
    vector<pair<int, int>> ret;
    while (true) {
      auto [x, y] = findNextPost();
      if (x == -1) break;
      ret.push_back({x, y});
      setField(x, y);
    }

    return ret;
  }
  void postMustPlaces() {
    auto vec = getMustPostPlaces();
    int k = 0;
    for (auto [x, y] : vec) {
      auto res = postQuery(1, bitset<512>().set(x * N + y));
      if (res == 0) {
        FOR(m, 0, M) FOR(id, 0, statesNum[m]) if (possible[m][id]) {
          auto oil = Oil(oils[m], id);
          if (oil.places.test(x * N + y)) {
            possible[m][id] = false;
          }
        }
      } else {
        FOR(i, 0, res) oi.push_back(k);
        updateInclude(x, y, res, k++);
      }
      if (calcStatesNum() < 200000) break;
    }
    return;
  }
  void updateInclude(int x, int y, int res, int k) {
    FOR(m, 0, M) {
      FOR(i, 0, statesNum[m]) {
        if (!possible[m][i]) continue;
        auto oil = Oil(oils[m], i);
        if (oil.places.test(x * N + y)) {
          includes[m][k].push_back(i);
        }
      }
    }
    return;
  }

  int calcStatesNum() {
    ll cnt = 0LL;

    auto oii = oi;
    while ((int)oii.size() < M) oii.push_back(-1);
    sort(range(oii));

    vector<int> posNumExclude(M, 0);
    FOR(m, 0, M) FOR(id, 0, statesNum[m]) if (possible[m][id]) posNumExclude[m]++;
    FOR(i, 0, 400) {
      ll c = 0LL;
      FOR(m, 0, M) {
        c += (int)includes[m][i].size();
        posNumExclude[m] -= (int)includes[m][i].size();
      }
      if (c == 0) break;
    }

    do {
      ll tmp = 1LL;
      FOR(m, 0, M) {
        if (oii[m] == -1) tmp *= posNumExclude[m];
        else tmp *= (int)includes[m][oii[m]].size();
      }
      cnt += tmp;
    } while (next_permutation(range(oii)));

    return cnt;
  }
  vector<vector<int>> getPossibleIdxs() {
    auto oii = oi;
    while ((int)oii.size() < M) oii.push_back(-1);
    sort(range(oii));

    vector<set<int>> posNumExcludeSt(M);
    FOR(m, 0, M) FOR(id, 0, statesNum[m]) if (possible[m][id]) posNumExcludeSt[m].insert(id);
    FOR(i, 0, 400) {
      int c = 0;
      FOR(m, 0, M) {
        c += (int)includes[m][i].size();
        for (auto s : includes[m][i]) posNumExcludeSt[m].erase(s);
        if (c == 0) break;
      }
    }
    vector<vector<int>> posNumExclude(M);
    FOR(m, 0, M) posNumExclude[m] = vector<int>(range(posNumExcludeSt[m]));

    vector<vector<int>> ret;
    do {
      function<void(int, vector<int>)> dfs = [&](int idx, vector<int> idxs) {
        if (idx == M) {
          ret.push_back(idxs);
          return;
        }
        if (oii[idx] == -1) {
          for (auto id : posNumExclude[idx]) {
            idxs.push_back(id);
            dfs(idx + 1, idxs);
            idxs.pop_back();
          }
        } else {
          for (auto od : includes[idx][oii[idx]]) {
            idxs.push_back(od);
            dfs(idx + 1, idxs);
            idxs.pop_back();
          }
        }
        return;
      };
      dfs(0, vector<int>());
    } while (next_permutation(range(oii)));

    cerr << "ret.size(): " << ret.size() << endl;

    return ret;
  }
};

struct Game5 {
  vector<OilBasic> oilBasic;
  OilState oilState;
  double statePos;
  vector<pair<bitset<512>, int>> query;
  set<ll> hashSt;
  Game5(vector<OilBasic> oilBasic) : oilBasic(oilBasic) {
    vector<Oil> oils;
    FOR(i, 0, M) oils.push_back(Oil(oilBasic[i], 0));
    oilState = OilState(0, oils);
  }
  void postQueryNow() {
    auto coords = oilState.oilsPosition();
    auto res = postQuery(2, coords);
    query.push_back({coords, res});
    statePos = calcPosState(oilState);
    return;
  }
  double calcPosState(OilState state) {
    double sum_log_likelihood = 0.0;
    for (auto [crd, x] : query) {
      int k = crd.count();
      auto num = state.overlapNumUsingCoords(crd);
      double avg = (k - num) * eps + num * (1.0 - eps);
      double var = k * eps * (1.0 - eps);
      double log_likelihood = logNormalLikelihood(avg, var, x);
      sum_log_likelihood += log_likelihood;
    }
    return sum_log_likelihood;
  }
  void solve1() {
    int m = Random.randint(M);
    int id = Random.randint(statesNum[m]);
    auto oil = Oil(oilBasic[m], id);
    auto newOilState = oilState;
    newOilState.changeOil(m, oil);
    if (hashSt.count(newOilState.hashNum)) return;
    double newPos = calcPosState(newOilState);
    if (newPos - statePos > 0.0) {
      // cerr << "[UPDATE1] newPos: " << newPos << " statePos: " << statePos << endl;
      oilState = newOilState;
      statePos = newPos;
      return;
    }
    return;
  }
  void solve2() {
    int a = Random.randint(M);
    int b = Random.randint(M);
    if (a == b) return;
    auto newOilState = oilState;
    newOilState.swapOil(a, b);
    if (hashSt.count(newOilState.hashNum)) return;
    double newPos = calcPosState(newOilState);
    if (newPos - statePos > 0.0) {
      // cerr << "[UPDATE2] newPos: " << newPos << " statePos: " << statePos << endl;
      oilState = newOilState;
      statePos = newPos;
      return;
    }
    return;
  }
  void solve3() {
    int m = Random.randint(M);
    int dx = Random.randint(3) - 1;
    int dy = Random.randint(3) - 1;
    if (!oilState.isMovePerOil(m, dx, dy)) return;
    auto newOil = oilState.getMoveOil(m, dx, dy);
    auto newOilState = oilState;
    newOilState.changeOil(m, newOil);
    if (hashSt.count(newOilState.hashNum)) return;
    double newPos = calcPosState(newOilState);
    if (newPos - statePos > 0.0) {
      // cerr << "[UPDATE3] newPos: " << newPos << " statePos: " << statePos << endl;
      oilState = newOilState;
      statePos = newPos;
      return;
    }
    return;
  }
  bool judgeOutput() {
    if ((int)query.size() < 5) return false;
    int size = (int)query.size();
    if (oilState.oilsPosition() != query[size - 1].first) return false;
    if (query[size - 1].first != query[size - 2].first) return false;
    if (query[size - 2].first != query[size - 3].first) return false;
    if (query[size - 3].first != query[size - 4].first) return false;
    if (query[size - 4].first != query[size - 5].first) return false;
    return true;
  }
  bool outputQuery() {
    bitset<512> coords;
    auto bs = oilState.oilsPosition();
    if (hashSt.count(oilState.hashNum)) return false;
    auto res = postQuery(3, bs);
    hashSt.insert(oilState.hashNum);
    return res == 1;
  }
};

struct Game6 {
  bitset<512> validRange;
  vector<pair<bitset<512>, int>> query;
  double pos;
  bitset<512> vvbit;
  Game6() {
    validRange.reset();
    FOR(i, 0, oilNum) validRange.set(i);
    vvbit.reset();
    FOR(i, 0, N * N) vvbit.set(i);
  }
  void postQueryNow() {
    auto res = postQuery(2, validRange);
    query.push_back({validRange, res});
    auto flip = validRange;
    flip = flip.flip() & vvbit;
    query.push_back({flip, oilNum - res});
    pos = calcPos(validRange);
    cerr << "pos: " << pos << endl;
    return;
  }
  double calcPos(bitset<512> bs) {
    auto flip = bs;
    flip = flip.flip() & vvbit;
    double sum_log_likelihood = 0.0;
    for (auto [crd, x] : query) {
      int k = crd.count();
      auto num = (crd & bs).count();
      double avg = (k - num) * eps + num * (1.0 - eps);
      double var = k * eps * (1.0 - eps);
      double log_likelihood = logNormalLikelihood(avg, var, x);
      sum_log_likelihood += log_likelihood;
    }
    return sum_log_likelihood;
  }
  void expand() {
    auto newValidRange = validRange;
    int x = Random.randint(N);
    int y = Random.randint(N);
    set<pair<int, int>> st;
    st.insert({x, y});
    int size = Random.randint(1, 9);
    int nx = x, ny = y;
    while ((int)st.size() < size) {
      int dir = Random.randint(4);
      if (nx + dx[dir] < 0 || nx + dx[dir] >= N) continue;
      if (ny + dy[dir] < 0 || ny + dy[dir] >= N) continue;
      nx += dx[dir];
      ny += dy[dir];
      st.insert({nx, ny});
    }

    for (auto [nx, ny] : st) {
      newValidRange.set(nx * N + ny);
    }

    UnionFind uf(N * N);
    FOR(i, 0, N) FOR(j, 0, N) if (newValidRange.test(i * N + j)) FOR(dir, 0, 4) {
      int ni = i + dx[dir], nj = j + dy[dir];
      if (ni < 0 || ni >= N) continue;
      if (nj < 0 || nj >= N) continue;
      if (newValidRange.test(ni * N + nj)) uf.unite(i * N + j, ni * N + nj);
    }

    FOR(i, 0, N * N) {
      if (!newValidRange.test(i)) continue;
      if (uf.size(i) < 4) newValidRange.reset(i);
    }

    auto newPos = calcPos(newValidRange);
    if (newPos - pos > 0.0) {
      cerr << "[UPDATE1] newPos: " << newPos << " pos: " << pos << endl;
      validRange = newValidRange;
      pos = newPos;
    }
    return;
  }
  void vanish() {
    auto newValidRange = validRange;
    int x = Random.randint(N);
    int y = Random.randint(N);
    set<pair<int, int>> st;
    st.insert({x, y});
    int size = Random.randint(1, 9);
    int nx = x, ny = y;
    while ((int)st.size() < size) {
      int dir = Random.randint(4);
      if (nx + dx[dir] < 0 || nx + dx[dir] >= N) continue;
      if (ny + dy[dir] < 0 || ny + dy[dir] >= N) continue;
      nx += dx[dir];
      ny += dy[dir];
      st.insert({nx, ny});
    }

    for (auto [nx, ny] : st) {
      newValidRange.reset(nx * N + ny);
    }

    UnionFind uf(N * N);
    FOR(i, 0, N) FOR(j, 0, N) if (newValidRange.test(i * N + j)) FOR(dir, 0, 4) {
      int ni = i + dx[dir], nj = j + dy[dir];
      if (ni < 0 || ni >= N) continue;
      if (nj < 0 || nj >= N) continue;
      if (newValidRange.test(ni * N + nj)) uf.unite(i * N + j, ni * N + nj);
    }

    FOR(i, 0, N * N) {
      if (!newValidRange.test(i)) continue;
      if (uf.size(i) < 4) newValidRange.reset(i);
    }

    auto newPos = calcPos(newValidRange);
    if (newPos - pos > 0.0) {
      cerr << "[UPDATE2] newPos: " << newPos << " pos: " << pos << endl;
      validRange = newValidRange;
      pos = newPos;
    }

    return;
  }
};

int main(void) {
  ios::sync_with_stdio(0);
  cin.tie(0);
  timer.reset();
  double TIMELIMIT = 2.0;
  random_device rnd;  // 非決定的な乱数生成器
  unsigned long long sd = (unsigned long long)rnd();
  Random.set_seed(sd);
  cin >> N >> M >> eps;
  restQuery = N * N * 2;
  validRange.reset();
  FOR(i, 0, N * N) validRange.set(i);
  vector<OilBasic> oils(M);
  for (int i = 0; i < M; i++) {
    int d;
    cin >> d;
    bitset<512> position;
    position.reset();
    int xmin = 1000, xmax = -1, ymin = 1000, ymax = -1;
    for (int j = 0; j < d; j++) {
      int x, y;
      cin >> x >> y;
      xmin = min(xmin, x);
      xmax = max(xmax, x);
      ymin = min(ymin, y);
      ymax = max(ymax, y);
      position.set(x * N + y);
    }
    oilNum += d;
    oils[i] = OilBasic(i, d, position, xmin, xmax, ymin, ymax);
    statesNum.push_back(oils[i].cntPosition());
  }
  FOR(m, 0, M) FOR(i, 0, statesNum[m]) hashTable[m][i] = Random.randint(1000000000);
  ll stateALLINIT = 1;
  FOR(i, 0, M) stateALLINIT = min(statesNum[i] * stateALLINIT, (ll)1e9);

  if (stateALLINIT < 400000) {
    Game1 game1(oils);

    while (restQuery > 0) {
      if (restQuery == 1 || game1.judgeAnswer()) {
        bool flag = game1.postQueryAnswer();
        if (flag || restQuery == 0) break;
      }
      game1.postQueryMaxPosState();
    }
  } else if (M > 10) {
    Game2 game2;
    game2.solve();
    game2.output();
  } else {
    Game5 game5(oils);
    while (restQuery > 0) {
      if (restQuery == 1 || game5.judgeOutput()) {
        auto flag = game5.outputQuery();
        if (flag || restQuery == 0) break;
      }
      game5.postQueryNow();
      FOR(i, 0, 1000) {
        if (Random.uniform() < 0.33) game5.solve1();
        else if (Random.uniform() < 0.5) game5.solve2();
        else game5.solve3();
      }
    }
  }

  return 0;
}
