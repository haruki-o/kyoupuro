
struct RollingHash {
  int N;

  // hashed[i] := s[0:i]のhash
  // get(l, r) := [l,r]のhashを求める
  vll hashed;
  ll base = 2007;
  vll power;

  RollingHash(string s) {
    N = (int)s.size();

    power.resize(N);
    power[0] = 1;
    rep(i, 1, N) power[i] = (power[i - 1] * base) % MOD;

    hashed.resize(N);
    hashed[0] = (ll)s[0];
    rep(i, 1, N) hashed[i] = (hashed[i - 1] * base + (ll)s[i]) % MOD;
  }

  ll get(ll l, ll r) {
    if (l == 0)
      return hashed[r];
    return (hashed[r] - (hashed[l - 1] * power[r - l + 1] % MOD) + MOD) % MOD;
  }
};
