#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
typedef pair<ll, ll> pll;
typedef pair<pll, ll> ppll;
typedef tuple<ll, ll, ll> tll;
typedef vector<tll> vt;
typedef vector<ppll> vpp;
typedef vector<pll> vp;
typedef vector<vp> vvp;
typedef vector<char> vc;
typedef vector<int> vi;
typedef vector<string> vs;
typedef vector<vs> vvs;
#define rep(i, l, n) for (ll i = (ll)(l); i < (ll)(n); i++)
#define repd(i, n, l) for (ll i = (ll)(n); i > (ll)(l); i--)
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
template <class T, class S>
inline bool chmax(T &a, const S &b) {
  return (a < b ? a = b, 1 : 0);
}
template <class T, class S>
inline bool chmin(T &a, const S &b) {
  return (a > b ? a = b, 1 : 0);
}
#define INF INT_MAX
#define INFF (9223372036854775800)
#define TIME_LIMIT (3.5)
#define def (20101)
// #define MOD (1000000007)
// #define MAX (2147483647)
#define MAX (1073741823)
#define MOD (998244353)

ll N, M, T, LA, LB;
vll u, v;
vll t;
vll x, y;
vvll g, gMat;

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
    if (par[c] == c) return c;
    ll a = root(par[c]);
    par[c] = a;
    return a;
  }
  // aをbに結合(違う時)
  bool unite(ll a, ll b) {
    ll ra = root(a);
    ll rb = root(b);
    if (ra == rb) return 0;
    par[ra] = rb;
    sum[rb] += sum[ra];
    return 1;
  }
  bool same(ll a, ll b) {
    if (root(a) == root(b)) return 1;
    return 0;
  }
};

double elaTimes(system_clock::time_point &clock) {
  return duration_cast<microseconds>(system_clock::now() - clock).count() * 1e-6;
}

class Event {
 public:
  char type;
  ll l, PA, PB;
  ll v;

  void setSignal(ll l, ll PA, ll PB) {
    this->type = 's';
    this->l = l;
    this->PA = PA;
    this->PB = PB;
  }

  void setMove(ll v) {
    this->type = 'm';
    this->v = v;
  }

  void printSignal() {
    cout << this->type << " " << this->l << " ";
    cout << this->PA << " " << this->PB << endl;
  }

  void printMove() {
    cout << this->type << " " << this->v << endl;
  }
};

struct PressG {
  vvll G;
  map<ll, vll> clusterToVertices;
  map<ll, ll> clusterToAIndex;

  // 圧縮後の頂点c1,c2が接続しているとき、圧縮前の接続する頂点を返す。
  pll edgeComposeVertice(ll c1, ll c2) {
    pll ret = {-1, -1};
    for (ll v1 : clusterToVertices[c1]) {
      for (ll v2 : clusterToVertices[c2]) {
        if (gMat[v1][v2] == 1) {
          ret = {v1, v2};
          return ret;
        }
      }
    }
    return ret;
  }

  // 圧縮前の頂点Vが圧縮後の頂点cを構成しているか
  bool isContainVertice(ll c, ll V) {
    for (auto cu : clusterToVertices[c]) {
      if (cu == V) return true;
    }
    return false;
  }

  // 圧縮後の頂点cを構成する圧縮前の頂点Vを取り除いたとき、連結になるか
  bool canRemoveVertice(ll c, ll V) {
    vll vertices = clusterToVertices[c];
    ll size = (ll)vertices.size();
    vll sub;
    rep(i, 0, size) {
      if (vertices[i] != V) {
        sub.push_back(vertices[i]);
      }
    }
    UnionFind uf(size - 1);
    rep(i, 0, size - 1) {
      rep(j, i, size - 1) {
        if (gMat[sub[i]][sub[j]]) uf.unite(i, j);
      }
    }
    rep(i, 0, size - 1) if (!uf.same(0, i)) return false;
    return true;
  }

  // クラスターc内の圧縮前の頂点V1からV2までの経路
  vector<Event> moveInCluster(ll c, ll V1, ll V2) {
    vll vertices = clusterToVertices[c];
    ll size = (ll)vertices.size();
    // bfs
    queue<ll> qu;
    qu.push(V1);
    map<ll, ll> ma;
    ma[V1] = -1;
    while (!qu.empty()) {
      ll cu = qu.front();
      qu.pop();
      rep(i, 0, size) {
        ll to = vertices[i];
        if (gMat[cu][to] == 0) continue;
        if (ma.count(to) > 0) continue;
        ma[to] = cu;
        qu.push(to);
      }
    }

    // 経路復元
    vll path;
    while (V2 != -1) {
      path.insert(path.begin(), V2);
      V2 = ma[V2];
    }

    // 経路の先頭からEventを詰める
    vector<Event> ret;
    rep(i, 1, path.size()) {
      Event addEvent;
      addEvent.setMove(path[i]);
      ret.push_back(addEvent);
    }
    return ret;
  }

  // クラスターcを構成する圧縮前の頂点Vをクラスタ-c'に移動して調整する
  void adjustClusterNum() {
    auto startTime = system_clock::now();
    const double limitTime = 0.1;
    ll size = clusterToVertices.size();
    ll loop1 = 0, loop2 = 0;
    while (elaTimes(startTime) < limitTime) {
      loop1++;
      ll idx1 = rand() % size;
      while (clusterToVertices[idx1].size() <= LB) {
        idx1 = (idx1 + 1) % size;
        if (elaTimes(startTime) < limitTime) break;
      }
      ll idx2 = rand() % size;
      pll edge;
      rep(i, 0, size) {
        idx2 = (idx2 + 1) % size;
        if (idx1 != idx2) {
          edge = edgeComposeVertice(idx1, idx2);
          if (edge.first != -1) break;
        }
      }
      if (edge.first == -1) continue;

      ll idx1Sisze = (ll)clusterToVertices[idx1].size();
      ll idx2Sisze = (ll)clusterToVertices[idx2].size();
      if (idx1Sisze > idx2Sisze) {
        swap(idx1, idx2);
        swap(edge.first, edge.second);
      }

      if (!canRemoveVertice(idx2, edge.second)) continue;
      // cerr << "--------------" << endl;
      // cerr << "object" << endl;
      // cerr << edge.first << " " << edge.second << endl;
      // cerr << "idx1" << endl;
      // for (ll cu : clusterToVertices[idx1]) cerr << cu << " ";
      // cerr << endl;
      // cerr << "idx2" << endl;
      // for (ll cu : clusterToVertices[idx2]) cerr << cu << " ";
      // cerr << endl;

      clusterToVertices[idx1].push_back(edge.second);
      remove(clusterToVertices[idx2].begin(), clusterToVertices[idx2].end(), edge.second);
      clusterToVertices[idx2].resize(clusterToVertices[idx2].size() - 1);
      // cerr << "idx1" << endl;
      // for (ll cu : clusterToVertices[idx1]) cerr << cu << " ";
      // cerr << endl;
      // cerr << "idx2" << endl;
      // for (ll cu : clusterToVertices[idx2]) cerr << cu << " ";
      // cerr << endl;
      loop2++;
      // break;
    }
    // cerr << "loop1, loop2: " << loop1 << " " << loop2 << endl;
  }

  // クラスターに所属してない頂点からクラスターを生成する
  void addCluster() {
    // LR個より大きいクラスターを削除
    map<ll, vll> newMa;
    ll idx = 0;
    for (auto cu : clusterToVertices) {
      if (cu.second.size() <= LB) {
        newMa[idx] = cu.second;
        idx++;
      }
    }
    clusterToVertices.erase(clusterToVertices.begin(), clusterToVertices.end());
    clusterToVertices = newMa;

    vll used(N, 0);
    while (1) {
      for (auto clus : clusterToVertices) {
        for (ll cu : clus.second) used[cu] = 1;
      }
      int flag = 0;
      rep(i, 0, N) if (!used[i]) flag = 1;
      if (!flag) break;

      rep(i, 0, N) {
        if (used[i]) continue;
        vll addClusterVector;
        priority_queue<pll, vp, greater<pll>> pq;
        pq.push({0, i});
        while (!pq.empty() && addClusterVector.size() < LB) {
          pll cu = pq.top();
          pq.pop();
          if (!used[cu.second]) addClusterVector.push_back(cu.second);
          used[cu.second] = 1;
          for (ll to : g[cu.second]) {
            if (used[to]) continue;
            pq.push({used[to], to});
          }
        }
        ll addIndex = clusterToVertices.size();
        clusterToVertices[addIndex] = addClusterVector;
      }
    }
  }

  // PressGを構築する
  void constructPressG(vll &A) {
    vp ok;

    rep(i, 0, A.size()) {
      if (i != 0 && i != A.size() - 1) {
        if (A[i] == 0 && A[i + 1] == 0) break;
      }
      repd(size, LB, 0) {
        if (A.size() < i + size) continue;

        map<ll, ll> ma;
        rep(j, i, i + size) ma[A[j]] = 1;
        queue<ll> qu;
        qu.push(A[i]);
        vll seen(A.size(), 0);
        seen[A[i]] = 1;
        while (!qu.empty()) {
          ll cu = qu.front();
          qu.pop();
          for (ll to : g[cu]) {
            if (ma.count(to) == 0) continue;
            if (seen[to] == 1) continue;
            qu.push(to);
            seen[to] = 1;
          }
        }
        bool flag = 1;
        rep(j, i, i + size) if (seen[A[j]] == 0) flag = 0;
        if (flag) {
          ok.push_back({i, size});
          break;
        }
      }

      vll pre;
      vll cu;
      int flag = 1;
      if (2 <= ok.size()) {
        pll _ = ok[ok.size() - 2];
        rep(j, _.first, _.first + _.second) {
          pre.push_back(A[j]);
        }
        _ = ok[ok.size() - 1];
        rep(j, _.first, _.first + _.second) {
          cu.push_back(A[j]);
        }

        int flag = 0;
        // cuがpreの一部
        if (pre.size() < cu.size()) continue;
        rep(j, 0, cu.size()) {
          if (pre[pre.size() - 1 - j] != cu[cu.size() - 1 - j]) flag = 1;
        }
        if (!flag) {
          ok.erase(ok.begin() + ok.size() - 1);
        }
      }
    }

    // cerr << "ok.size() : " << ok.size() << endl;
    // G.resize((ll)min((ll)ok.size(), (ll)100));
    G.resize((ll)ok.size());
    rep(i, 0, ok.size()) {
      rep(j, ok[i].first, ok[i].first + ok[i].second) {
        clusterToVertices[i].push_back(A[j]);
      }
      clusterToAIndex[i] = ok[i].first;
    }
    rep(i, 0, G.size()) {
      rep(j, i + 1, G.size()) {
        pll _ = edgeComposeVertice(i, j);
        ll flag = 0;
        for (auto cu : clusterToVertices[i]) {
          for (auto to : clusterToVertices[j]) {
            if (cu == to) flag = 1;
          }
        }
        if (flag) continue;
        if (_.first != -1) {
          G[i].push_back(j);
          G[j].push_back(i);
        }
      }
    }
  }
};

struct Solver {
  const ll DEPTH = min(LB, (ll)15);
  void dfs(ll c, ll p, vll &visited, vvll &tree, vll &A, map<pll, ll> &edgeScore) {
    if (visited[c]) return;
    A.push_back(c);
    visited[c] = 1;
    if (p != -1) {
      tree[p].push_back(c);
      tree[c].push_back(p);
    }

    vp vec;
    for (auto cu : g[c]) {
      vec.push_back({edgeScore[{min(c, cu), max(c, cu)}], cu});
    }
    gSort(vec);
    for (auto cu : vec) {
      dfs(cu.second, c, visited, tree, A, edgeScore);
    }
  }

  void initDfs(ll c, ll p, vll &visited, vvll &tree, vll &A, map<pll, ll> &edgeScore, queue<ll> qu) {
    if (visited[c]) return;
    A.push_back(c);
    qu.push(c);
    // if (DEPTH < qu.size()) qu.pop();

    visited[c] = 1;
    if (p != -1) {
      tree[p].push_back(c);
      tree[c].push_back(p);
    }

    vp vec;
    for (auto cu : g[c]) {
      ll dist = abs(x[cu] - x[qu.front()]) + abs(y[cu] - y[qu.front()]);
      dist += edgeScore[{min(c, cu), max(c, cu)}];
      vec.push_back({dist, cu});
    }
    gSort(vec);
    int flag = 0;
    for (auto cu : vec) {
      queue<ll> newQu;
      if (flag) {
        newQu.push(c);
        initDfs(cu.second, c, visited, tree, A, edgeScore, newQu);
      } else initDfs(cu.second, c, visited, tree, A, edgeScore, qu);
      flag = 1;
    }
  }

  void constructTree(vll &A, vvll &tree) {
    map<pll, ll> edgeScore;
    // ダイクストラで何回辺を通ったか
    rep(i, 0, T - 1) {
      vll candi;
      // ダイクストラ
      vll dp(N, INF);
      vll pre(N, -1);
      queue<pll> qu;
      qu.push({t[i], -1});
      dp[t[i]] = 1;
      while (!qu.empty()) {
        auto cu = qu.front();
        qu.pop();
        for (ll to : g[cu.first]) {
          if (to == t[i + 1] && dp[cu.first] + 1 <= dp[to]) {
            candi.push_back(to);
          }
          if (dp[to] != INF) continue;
          qu.push({to, cu.first});
          dp[to] = dp[cu.first] + 1;
          pre[to] = cu.first;
        }
      }

      // 経路復元
      rep(j, 0, candi.size()) {
        ll at = candi[j];
        edgeScore[{min(at, t[i + 1]), max(at, t[i + 1])}]++;
        while (at != -1) {
          ll to = pre[at];
          if (to != -1) edgeScore[{min(at, to), max(at, to)}]++;
          at = to;
        }
      }
    }

    ll root = 0;
    rep(i, 0, N) {
      if (abs(x[i] - 500) + abs(y[i] - 500) < abs(x[root] - 500) + abs(y[root] - 500)) root = i;
    }

    vll visited(N, 0);
    queue<ll> qu;
    initDfs(root, -1, visited, tree, A, edgeScore,qu);
 
    ll addNum = LA - A.size();
    visited.assign(N,0);
    vll _A;
    dfs(root, -1, visited, tree, _A, edgeScore);
    rep(i,0,addNum)A.push_back(_A[i]);

    
  }

  void initAns(vll &A, vector<Event> &ans, PressG &pressG) {
    ll preCluster;
    ll preVertice = 0;
    double stat = 0;
    ll M = 0;
    rep(i, 0, pressG.G.size()) M += pressG.G[i].size();
    // cerr << "N, M: " << pressG.clusterToVertices.size() << " " << M / 2 << endl;
    rep(i, 0, T) {
      // cerr << endl << "t[] index: " << i << endl;
      // pqの初期値
      priority_queue<tll, vt, greater<tll>> pq;
      if (i == 0) {
        rep(j, 0, pressG.G.size()) {
          int flag = 0;
          for (auto cu : g[0]) {
            if (pressG.isContainVertice(j, cu)) {
              // cerr << "prePressVertice: " << cu << endl;
              // cerr << "cluster: " << j << endl;
              // for (auto cu : pressG.clusterToVertices[j]) cerr << cu << " ";
              // cerr << endl;
              flag = 1;
            }
            if (flag) break;
          }
          if (flag) pq.push({1, j, -1});
        }
      } else {
        rep(j, 0, pressG.G.size()) {
          pll edge = pressG.edgeComposeVertice(preCluster, j);
          if (edge.first != -1) pq.push({1, j, -1});
        }
      }

      // ダイクストラ法
      auto startTime = system_clock::now();

      vll dp(pressG.G.size(), INF);
      vll pre(pressG.G.size(), -1);
      while (!pq.empty()) {
        auto cu = pq.top();
        pq.pop();
        if (dp[get<1>(cu)] == INF) {
          dp[get<1>(cu)] = get<0>(cu);
          pre[get<1>(cu)] = get<2>(cu);
        } else continue;
        for (auto to : pressG.G[get<1>(cu)]) {
          if (dp[to] != INF) continue;
          pq.push({get<0>(cu) + 1, to, get<1>(cu)});
        }
      }
      stat += elaTimes(startTime);
      // cerr << "a";
      // 最良の行き先(圧縮後の頂点bestV)を見つける
      ll bestV = -1;
      ll dist = INF;
      ll score = -1;
      rep(j, 0, pressG.G.size()) {
        if (!pressG.isContainVertice(j, t[i])) continue;
        bool flag = 0;
        flag = (dp[j] < dist);
        ll atScore = i;
        if (dp[j] == dist) {
          while (pressG.isContainVertice(j, t[atScore])) {
            atScore++;
            if (atScore == T) break;
          }
          if (chmax(score, atScore)) flag = 1;
        }

        if (flag) {
          bestV = j;
          dist = dp[j];
        }
      }
      // cerr << bestV << endl;
      // cerr << " b ";

      // 経路復元
      ll atV = bestV;
      vll path;
      while (atV != -1) {
        // cerr << atV << " ";
        path.insert(path.begin(), atV);
        atV = pre[atV];
      }
      // cerr << endl << " c " << endl;

      // 経路先頭からEventを詰める
      ll fromInCluster;
      ll toInCluster;
      rep(j, 0, path.size()) {
        // i-1行けるところまで進めた時、帳尻合わせ
        if (i != 0 && j == 0) {
          pll joinEdge = pressG.edgeComposeVertice(preCluster, path[0]);
          vector<Event> moveVec = pressG.moveInCluster(preCluster, preVertice, joinEdge.first);
          for (Event cu : moveVec) ans.push_back(cu);
        }
        // 次のクラスターを照らす
        Event addEvent;
        ll l = (ll)pressG.clusterToVertices[path[j]].size();
        ll PA = pressG.clusterToAIndex[path[j]];
        ll PB = 0;
        addEvent.setSignal(l, PA, PB);
        ans.push_back(addEvent);
        // rep(k, 0, l) cerr << A[PA + k] << " ";
        // cerr << endl;

        if (j != path.size() - 1) {
          pll joinEdge;
          joinEdge = pressG.edgeComposeVertice(path[j], path[j + 1]);
        }
        // 前のクラスターから今のクラスターに移動
        Event addMoveEvent;
        if (i == 0 && j == 0) {
          for (ll cu1 : g[0]) {
            for (ll cu2 : pressG.clusterToVertices[path[j]])
              if (cu1 == cu2) {
                addMoveEvent.setMove(cu1);
                fromInCluster = cu1;
              }
          }
        } else {
          pll joinEdge;
          if (j == 0) {
            joinEdge = pressG.edgeComposeVertice(preCluster, path[0]);
          } else {
            joinEdge = pressG.edgeComposeVertice(path[j - 1], path[j]);
          }
          addMoveEvent.setMove(joinEdge.second);
          fromInCluster = joinEdge.second;
        }
        ans.push_back(addMoveEvent);

        // クラスター内で移動する
        pll joinEdge = {-1, -1};
        if (j != path.size() - 1) {
          joinEdge = pressG.edgeComposeVertice(path[j], path[j + 1]);
        }
        if (j == path.size() - 1) toInCluster = t[i];
        else toInCluster = joinEdge.first;
        // cerr << fromInCluster << " " << toInCluster << endl;
        vector<Event> moveVec = pressG.moveInCluster(path[j], fromInCluster, toInCluster);
        for (Event cu : moveVec) ans.push_back(cu);
        fromInCluster = joinEdge.second;
      }

      // 更新
      preCluster = bestV;
      preVertice = toInCluster;
      if (i == T - 1) break;
      while (pressG.isContainVertice(bestV, t[i + 1])) {
        i++;
        // クラスター内で移動する
        fromInCluster = toInCluster;
        toInCluster = t[i];
        preVertice = toInCluster;
        vector<Event> moveVec = pressG.moveInCluster(bestV, fromInCluster, toInCluster);
        for (Event cu : moveVec) ans.push_back(cu);
        if (i == T) break;
      }
    }
    // cerr << "initAns() ダイクストラ: " << stat << endl;
  }

  void solve(vector<Event> &ans, vll &A) {
    auto startTime = system_clock::now();

    vvll tree(N);
    constructTree(A, tree);

    PressG pressG;
    pressG.constructPressG(A);

    auto startTime2 = system_clock::now();
    initAns(A, ans, pressG);
  }
};

void checkAns(vector<Event> &ans, vll &A) {
  if (A.size() != LA) {
    cerr << "A.size() != LA" << endl;
    return;
  }

  ll at = 0;
  vll ligth(LB, 0);
  for (auto cu : ans) {
    if (cu.type == 's') {
      rep(i, cu.PA, cu.PA + cu.l) {
        if (LB <= cu.PB + i - cu.PA) {
          cerr << "out of Index" << endl;
          return;
        }
        ligth[LB + i - cu.PA] = A[i];
      }
    }
    if (cu.type == 'm') {
      ll flag = 0;
      for (auto cu1 : g[at]) {
        if (cu1 == cu.v) flag = 1;
      }
      if (!flag) {
        cerr << "can't move" << endl;
      }
      at = cu.v;
    }
  }
}

void printAns(vector<Event> &ans, vll &A) {
  for (auto cu : A) cout << cu << " ";
  cout << endl;
  rep(i, 0, ans.size()) {
    if (ans[i].type == 's') ans[i].printSignal();
    else ans[i].printMove();
  }
}

ll calcAns(vector<Event> &ans) {
  ll ret = 0;
  for (auto cu : ans)
    if (cu.type == 's') ret++;
  return ret;
}

int main() {
  srand((unsigned int)time(NULL));
  auto startTime = system_clock::now();

  cin >> N >> M >> T >> LA >> LB;

  u.resize(M);
  v.resize(M);
  rep(i, 0, M) cin >> u[i] >> v[i];
  t.resize(T);
  rep(i, 0, T) cin >> t[i];
  x.resize(N);
  y.resize(N);
  rep(i, 0, N) cin >> x[i] >> y[i];
  g.resize(N);
  rep(i, 0, M) {
    g[u[i]].push_back(v[i]);
    g[v[i]].push_back(u[i]);
  }
  gMat.resize(N);
  rep(i, 0, N) gMat[i].assign(N, 0);
  rep(i, 0, N) {
    for (ll cu : g[i]) {
      gMat[i][cu] = 1;
      gMat[cu][i] = 1;
    }
  }

  vector<Event> ans;
  vll A;
  Solver solver;
  solver.solve(ans, A);

  printAns(ans, A);

  // cerr << "time:  " << elaTimes(startTime) << endl;
  // cerr << "score: " << calcAns(ans) << endl;
  // checkAns(ans, A);
  cerr << calcAns(ans) << endl;
}

// local9 <- local7
// 木の生成方法を遠くの点を採用するようにしました。(ダメ)