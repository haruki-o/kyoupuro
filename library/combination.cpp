// fac:p!%MOD
// inv:(a^-1)%MOD
// fact_inv:invの累積
// nCkでn=10^9,k=10^7のとき
vll fac;
vll inv;
vll fact_inv;
//(a^-1)≡-(p%a)^-1*(p/a) (modp)
void COMinit(ll SIZE) {
  fac.resize(SIZE + 4);
  inv.resize(SIZE + 4);
  fact_inv.resize(SIZE + 4);
  fac[1] = 1;
  inv[1] = 1;
  fact_inv[0] = 1;
  fact_inv[1] = 1;
  rep(i, 2, SIZE + 2) {
    fac[i] = fac[i - 1] * i % MOD;
    inv[i] = ((MOD - (MOD / i * inv[MOD % i]) % MOD) + MOD) % MOD;
    fact_inv[i] = fact_inv[i - 1] * inv[i] % MOD;
  }
}

// nCk=n!/k!/(n-k)!=n!*(k!)^-1*((n-k)!)^-1 O(1)
ll COM(ll n, ll k) {
  if (n == 0)
    return 1;
  ll ans = fac[n] * fact_inv[k] % MOD * fact_inv[n - k] % MOD;
  return ans;
}