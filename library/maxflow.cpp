#include <bits/stdc++.h>
using namespace std;
typedef pair<int, int> P;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<P> vP;
typedef vector<vP> vvP;
#define rep(i, l, n) for (int i = (int)(l); i < (int)(n); i++)
#define INF INT_MAX

// 残余ネットワークを表し、g[i]は始点がiの有向辺の集合で有向辺は(有向辺の終点,重み)で表す。
vvP g;
// 増加可能経路が発見された時逆向きの辺に直接加えるのではなく、それ以降に増加可能経路を発見する際、その辺にアクセスするときに容量を追加する用のmap.
map<P, int> lazy;
vi seen;

int dfs(int u, int t, int f) {
  if (u == t)
    return f;
  seen[u] = 1;
  for (P &cu : g[u]) {
    if (lazy.count(pair(u, cu.first)) != 0) {
      cu.second += lazy[pair(u, cu.first)];
      lazy.erase(pair(u, cu.first));
    }
    if (cu.second == 0)
      continue;
    if (seen[cu.first] == 1)
      continue;
    int r = dfs(cu.first, t, min(f, cu.second));
    if (r > 0) {
      if (lazy.count(pair(cu.first, u)) == 0)
        lazy[pair(cu.first, u)] = r;
      else
        lazy[pair(cu.first, u)] += r;
      cu.second -= r;
      return r;
    }
  }
  return 0;
}

int main() {
  int V, E;
  cin >> V >> E;
  g.resize(V);
  rep(i, 0, E) {
    int u, v, c;
    cin >> u >> v >> c;
    g[u].push_back(pair(v, c));
    g[v].push_back(pair(u, 0));
  }
  int ans = 0;
  // 増加可能経路がなくなるまで繰り返す。
  while (1) {
    seen.assign(V, 0);
    int add = dfs(0, V - 1, INF);
    // 増加可能経路がないことを表す。
    if (add == 0)
      break;
    ans += add;
  }
  cout << ans << endl;
}
