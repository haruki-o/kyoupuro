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
      repd(size, LB, 0) {
        if (N < i + size) continue;

        map<ll, ll> ma;
        rep(j, i, i + size) ma[A[j]] = 1;
        queue<ll> qu;
        qu.push(A[i]);
        vll seen(N, 0);
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
    }

    // 必要最低限のok
    vll sum(N, 0);
    rep(i, 0, ok.size()) {
      rep(j, 0, ok[i].second) sum[A[ok[i].first + j]]++;
    }

    vll removeIdxVec;
    rep(i, 0, ok.size()) {
      int flag = 1;
      rep(j, 0, ok[i].second) {
        if (sum[A[ok[i].first + j]] <= 1) flag = 0;
      }
      if (flag) {
        rep(j, 0, ok[i].second) sum[A[ok[i].first + j]]--;
        removeIdxVec.push_back(i);
      }
    }

    // reverse(removeIdxVec.begin(), removeIdxVec.end());
    // cerr << ok.size() <<endl;
    // for (ll removeIdx : removeIdxVec) {
    //   ok.erase(ok.begin() + removeIdx);
    // }

    // std::random_device seed_gen;
    // std::mt19937 engine(seed_gen());
    // std::shuffle(ok.begin(), ok.end(), engine);

    // cerr << "ok.size() : " << ok.size() << endl;
    // G.resize((ll)min((ll)ok.size(), (ll)100));
    G.resize((ll)ok.size());
    rep(i, 0, ok.size()) {
      rep(j, ok[i].first, ok[i].first + ok[i].second) {
        clusterToVertices[i].push_back(A[j]);
      }
      clusterToAIndex[i] = i;
    }
    rep(i, 0, G.size()) {
      rep(j, i + 1, G.size()) {
        pll _ = edgeComposeVertice(i, j);
        if (_.first != -1) {
          G[i].push_back(j);
          G[j].push_back(i);
        }
      }
    }
  }
};

struct Solver {
  void dfs(ll v, vvll &tree, vll &visited) {
    visited[v] = 1;
    for (ll c : g[v]) {
      if (visited[c] == 1) continue;
      visited[c] = 1;
      tree[v].push_back(c);
      tree[c].push_back(v);
      dfs(c, tree, visited);
    }
  }

  void constructTree(vvll &tree) {
    vll visited(N, 0);
    dfs(0, tree, visited);
  }

  void constructADfs(ll v, vll &A, vvll &tree, vll &visited) {
    // cerr << v << " ";
    visited[v] = 1;
    A.push_back(v);
    for (ll c : g[v]) {
      if (visited[c] == 1) continue;
      visited[c] = 1;
      constructADfs(c, A, tree, visited);
    }
  }

  void initConstructA(vll &A, vvll &tree) {
    vll visited(N, 0);
    constructADfs(0, A, tree, visited);
    ll addNum = LA - A.size();
    rep(i, 0, addNum) A.push_back(0);
  }

  void initAns(vll &A, vector<Event> &ans, PressG &pressG) {
    ll preCluster;
    ll preVertice = 0;
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
      // cerr << "a";
      // 最良の行き先(圧縮後の頂点bestV)を見つける
      ll bestV = -1;
      ll dist = INF;
      ll score = INF;
      ll euScore = INF;
      rep(j, 0, pressG.G.size()) {
        if (!pressG.isContainVertice(j, t[i])) continue;
        bool flag = 0;
        flag = (dp[j] < dist);
        ll atScore = i;
        if (dp[j] == dist) {
          // T[i+1]にも進めるか
          while (pressG.isContainVertice(j, t[atScore])) {
            atScore++;
            if (atScore == T) break;
          }
          if (score < atScore) {
            score = atScore;
            flag = 1;
          }
          // T[i + 1]との距離が近い
          ll atBestEuDist = INF;
          if (i != T - 1) {
            for (auto cu : pressG.clusterToVertices[j]) {
              ll atEuDist = (x[t[i + 1]] - x[cu]) * (x[t[i + 1]] - x[cu]);
              atEuDist += (y[t[i + 1]] - y[cu]) * (y[t[i + 1]] - y[cu]);
              chmin(atBestEuDist, atEuDist);
            }
          }
          if (!flag) {
            if (chmin(atBestEuDist, euScore)) {
              flag = 1;
            }
          }
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
  }

  void solve(vector<Event> &ans, vll &A) {
    auto startTime = system_clock::now();

    vvll tree(N);
    constructTree(tree);
    initConstructA(A, tree);

    PressG pressG;
    pressG.constructPressG(A);

    // vll seen(N, 0);
    // for (auto cu : pressG.clusterToVertices) {
    //   for (auto cu1 : cu.second) seen[cu1] = 1;
    // }
    // rep(i, 0, N) if (seen[i] == 0) cerr << "x";
    auto startTime2 = system_clock::now();
    initAns(A, ans, pressG);
    // cerr << "initAns() time: " << elaTimes(startTime2) << endl;
  }
};

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
  cerr << calcAns(ans) << endl;
}

// local5 <- local3
// 木を生成→Aを構築→クラスターを構築(N個)
// local3と変わらず
