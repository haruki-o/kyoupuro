//初期化Nは2^iでなくてもよい
// SegmentTree st(N);
// sectionの初期値はあたえないといけない
//葉っぱの一番左はsection[size-1]で接続
struct LazySegmentTree {
  ll size = 1;
  vll section;
  vll lazy;
  LazySegmentTree(ll N) {
    while (size < N) {
      size *= 2;
    }
    section.assign(2 * size - 1, 0);
    lazy.assign(2 * size - 1, 0);
  }
  //一段降下
  void eval(ll k, ll x, ll y) {
    if (lazy[k] != 0)
      section[k] += lazy[k];
    //最下段かどうか
    if (y - x > 0) {
      lazy[(k + 1) * 2 - 1] = lazy[k] / 2;
      lazy[(k + 1) * 2] = lazy[k] / 2;
    }
    lazy[k] = 0;
  }
  //[l,r}にcを加える,kは今見てる節で[x,y} l,r,x,yは[size-1,2size-2]
  void add(ll l, ll r, ll k, ll x, ll y, ll c) {
    //区間が少しでもかぶっていた時用
    eval(k, x, y);
    if (r < x || y < l)
      return;
    if (l <= x && y <= r) {
      lazy[k] = c * (r - l + 1);
      eval(k, x, y);
    } else {
      add(l, r, (k + 1) * 2 - 1, x, y - (y - x) / 2 - 1, c);
      add(l, r, (k + 1) * 2, x + (y - x) / 2 + 1, y, c);
      section[k] = section[(k + 1) * 2 - 1] + section[(k + 1) * 2];
    }
  }
  ll current(ll l, ll r, ll k, ll x, ll y) {
    //区間が少しでもかぶっていた時用
    eval(k, x, y);
    if (r < x || y < l)
      return 0;
    if (l <= x && y <= r) {
      eval(k, x, y);
      return section[k];
    } else {
      ll sl = current(l, r, (k + 1) * 2 - 1, x, y - (y - x) / 2 - 1);
      ll sr = current(l, r, (k + 1) * 2, x + (y - x) / 2 + 1, y);
      return sl + sr;
    }
  }
};
// max版
//お気持ち
// eval:自sectionのみを更新して子どもに渡す
// update:
//完全被覆->evalに更新値を渡して自sectionと比べてだめなら捨てるこれは子供より確実に親が優れているため
//いいときは自sectionのみを更新して子どものlazyに渡すことで子供たちの更新はあとまわしにする
//部分被覆->子供が更新されたとき、親が劣る可能性があるため再帰的に更新する

struct LazySegmentTree {
  ll size = 1;
  vll section;
  vll lazy;
  LazySegmentTree(ll N) {
    while (size < N) {
      size *= 2;
    }
    section.assign(2 * size - 1, 0);
    lazy.assign(2 * size - 1, 0);
  }
  //一段降下
  void eval(ll k, ll x, ll y) {
    if (lazy[k] != 0)
      section[k] = max(lazy[k], section[k]);
    //最下段かどうか
    if (y - x > 0) {
      lazy[(k + 1) * 2 - 1] = max(lazy[(k + 1) * 2 - 1], section[k]);
      lazy[(k + 1) * 2] = max(lazy[(k + 1) * 2], section[k]);
    }
    lazy[k] = 0;
  }
  //[l,r}にcを加える,kは今見てる節で[x,y} l,r,x,yは[size-1,2size-2]
  void update(ll l, ll r, ll k, ll x, ll y, ll c) {
    //区間が少しでもかぶっていた時用
    eval(k, x, y);
    if (r < x || y < l)
      return;
    if (l <= x && y <= r) {
      lazy[k] = c;
      eval(k, x, y);
    } else {
      update(l, r, (k + 1) * 2 - 1, x, y - (y - x) / 2 - 1, c);
      update(l, r, (k + 1) * 2, x + (y - x) / 2 + 1, y, c);
      section[k] = max(section[(k + 1) * 2 - 1], section[(k + 1) * 2]);
    }
  }
  ll current(ll l, ll r, ll k, ll x, ll y) {
    //区間が少しでもかぶっていた時用
    eval(k, x, y);
    if (r < x || y < l)
      return 0;
    if (l <= x && y <= r) {
      eval(k, x, y);
      return section[k];
    } else {
      ll sl = current(l, r, (k + 1) * 2 - 1, x, y - (y - x) / 2 - 1);
      ll sr = current(l, r, (k + 1) * 2, x + (y - x) / 2 + 1, y);
      return max(sl, sr);
    }
  }
};

//abc223 F editorial 参照
//https://atcoder.jp/contests/abc223/editorial/2774
//https://atcoder.github.io/ac-library/document_ja/lazysegtree.html

// op : 区間取得(S), prod, all_prod
// mapping : 区間更新(F), apply

ll op(ll a, ll b){ return min(a, b); }
ll e(){ return 1001001001; }
ll mapping(ll f, ll x){ return f+x; }
ll composition(ll f, ll g){ return f+g; }
ll id(){ return 0; }

lazy_segtree<ll, op, e, ll, mapping, composition, id> seg(v);