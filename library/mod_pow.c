//pow(x,n)のmodがわかる
//計算量はlogN

ll mod_pow(ll x,ll n,ll mod){
  x%=mod;
  ll res=1; 
  while(n>0){
    if(n & 1)res*=x;
    x = x*x;
    n >>= 1;
    res%=mod;
    x%=mod;
  }
  return res;
}