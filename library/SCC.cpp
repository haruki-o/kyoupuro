struct SCC {
  int N;
  vvll G, rG;
  vi seen;
  //トポ順の逆
  vi order;
  // group[i] == group[j]のとき頂点i,jは強連結,[0,N) -> [0,k)
  vi group;
  // group_size[i] := group iのsize, [0,k) -> [0,N-1)
  vi group_size;
  // rebuildGのsize
  int k = 0;
  //[0,k) -> [0,k)の部分集合
  vvi rebuildG;
  vector<set<int>> _rebuildG;

  SCC(vvll &_G) {
    N = _G.size();
    G.resize(N);
    G = _G;
    rG.resize(N);
    rep(i, 0, N) {
      for (int cu : G[i])
        rG[cu].push_back(i);
    }

    seen.assign(N, 0);
    rep(i, 0, N) if (seen[i] == 0) dfs(i, -1);
    seen.assign(N, 0);
    group.resize(N);
    rep(i, 0, N) {
      if (seen[order[N - 1 - i]] == 0) {
        group_size.push_back(0);
        rdfs(order[N - 1 - i], -1, k);
        k++;
      }
    }

    rebuildG.resize(k);
    _rebuildG.resize(k);
    rep(i, 0, N) {
      for (int cu : G[i])
        _rebuildG[group[i]].insert(group[cu]);
    }
    rep(i, 0, k) {
      for (int cu : _rebuildG[i])
        rebuildG[i].push_back(cu);
    }
  }

  void dfs(int c, int p) {
    if (seen[c] == 1)
      return;
    seen[c] = 1;
    for (int cu : G[c]) {
      if (seen[cu] == 1)
        continue;
      dfs(cu, c);
    }
    order.push_back(c);
  }

  void rdfs(int c, int p, int gro) {
    if (seen[c] == 1)
      return;
    seen[c] = 1;
    group[c] = gro;
    if ((int)group_size.size() != 0)
      group_size[(int)group_size.size() - 1]++;
    for (int cu : rG[c]) {
      if (seen[cu] == 1)
        continue;
      rdfs(cu, c, gro);
    }
  }
};
