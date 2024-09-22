//初期化Nは2^iでなくてもよい
// SegmentTree st(N);
// sectionの初期値はあたえないといけない
//葉っぱの一番左はsection[size-1]で接続
struct SegmentTree {
  ll size = 1;
  vll section;
  SegmentTree(ll N) {
    while (size < N) {
      size *= 2;
    }
    section.assign(2 * size - 1, 0);
  }
  //[l,r}にcを加える,kは今見てる節で[x,y} l,r,x,yは[size-1,2size-2]
  void add(ll l, ll r, ll k, ll x, ll y, ll c) {
    if (l <= x && y <= r) {
      section[k] += c;
      return;
    }
    if (r < x || y < l)
      return;
    add(l, r, (k + 1) * 2 - 1, x, y - (y - x) / 2 - 1, c);
    add(l, r, (k + 1) * 2, x + (y - x) / 2 + 1, y, c);
  }
  // c=[0,size}は葉っぱ,葉っぱから根までの総和
  ll current(ll c) {
    c += size - 1;
    ll sum = 0;
    while (c != 0) {
      sum += section[c];
      c = (c - 1) / 2;
    }
    sum += section[0];
    return sum;
  }
};
// max版
struct SegmentTree {
  ll size = 1;
  vll section;
  SegmentTree(ll N) {
    while (size < N) {
      size *= 2;
    }
    section.assign(2 * size - 1, 0);
  }
  void update(ll l, ll r, ll k, ll x, ll y, ll c) {
    if (l <= x && y <= r) {
      section[k] = c > section[k] ? c : section[k];
      return;
    }
    if (r < x || y < l)
      return;
    update(l, r, (k + 1) * 2 - 1, x, y - (y - x) / 2 - 1, c);
    update(l, r, (k + 1) * 2, x + (y - x) / 2 + 1, y, c);
  }
  // c=[0,size)は葉っぱ,葉っぱから根までのmax
  ll current(ll c) {
    c += size - 1;
    ll at_max = 0;
    while (c != 0) {
      at_max = section[c] > at_max ? section[c] : at_max;
      c = (c - 1) / 2;
    }
    at_max = section[0] > at_max ? section[0] : at_max;
    return at_max;
  }
};
struct SegmentTree {
  ll N;
  vector<ll> ST;
  SegmentTree(vector<ll> A) {
    ll N2 = A.size();
    N = 1;
    while (N < N2) {
      N *= 2;
    }
    ST = vector<ll>(N * 2 - 1, -1);
    for (ll i = 0; i < N2; i++) {
      ST[N - 1 + i] = A[i];
    }
    for (ll i = N - 2; i >= 0; i--) {
      ST[i] = max(ST[i * 2 + 1], ST[i * 2 + 2]);
    }
  }

  void re(vector<ll> &A) {
    for (ll i = 0; i < N2; i++) {
      ST[N - 1 + i] = A[i];
    }
    for (ll i = N - 2; i >= 0; i--) {
      ST[i] = max(ST[i * 2 + 1], ST[i * 2 + 2]);
    }
  }
  void update(ll k, ll x) {
    k += N - 1;
    ST[k] = x;
    while (k > 0) {
      k = (k - 1) / 2;
      ST[k] = max(ST[k * 2 + 1], ST[k * 2 + 2]);
    }
  }
  ll range_max(ll L, ll R, ll i, ll l, ll r) {
    if (r <= L || R <= l) {
      return -1;
    } else if (L <= l && r <= R) {
      return ST[i];
    } else {
      ll m = (l + r) / 2;
      return max(range_max(L, R, i * 2 + 1, l, m),
                 range_max(L, R, i * 2 + 2, m, r));
    }
  }
  ll range_max(ll L, ll R) { return range_max(L, R, 0, 0, N); }
  ll max_right(ll L, ll R, ll x, ll i, ll l, ll r) {
    if (r <= L || R <= l) {
      return -1;
    } else if (ST[i] < x) {
      return -1;
    } else if (r - l == 1) {
      return i - (N - 1);
    } else {
      ll m = (l + r) / 2;
      ll tmp = max_right(L, R, x, i * 2 + 1, l, m);
      if (tmp != -1) {
        return tmp;
      } else {
        return max_right(L, R, x, i * 2 + 2, m, r);
      }
    }
  }
  ll max_right(ll L, ll R, ll x) { return max_right(L, R, x, 0, 0, N); }
};
