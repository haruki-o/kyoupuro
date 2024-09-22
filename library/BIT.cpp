struct binary_indexed_tree {
  ll N;
  vector<ll> BIT;

  binary_indexed_tree(ll _N) : N(_N), BIT(_N + 1, 0) {}

  void init(ll _N) {
    N = _N;
    BIT.assign(_N + 1, 0);
  }
  void add(ll i, ll x) {
    i++;
    while (i <= N) {
      BIT[i] += x;
      i += i & -i;
    }
  }

  ll sum(ll i) {
    i++;
    ll ans = 0;
    while (i > 0) {
      ans += BIT[i];
      i -= i & -i;
    }
    return ans;
  }

  ll sum(ll L, ll R) { return sum(R) - sum(L); }
};

binary_indexed_tree BIT(N);
// https://qiita.com/DaikiSuyama/items/7295f5160a51684554a7
// add -> O(logN)
// sum -> O(logN)

// https://atcoder.jp/contests/abc202/submissions/22808803