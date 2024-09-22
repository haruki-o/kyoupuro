//10進数以外の時はstringとして扱う.

ll baseR_ll(string s,int r){
  ll re = 0;
  rep(i,0,(int)s.size()){
    re = re * r + (s[i] - '0');
  }
  return re;
}

string ll_baseR(ll N,int r){
  string re = "";
  if(N == 0)return "0";
  while(N > 0){
    re = (char)(N % r + '0') + re; 
    N /= r;
  }
  return re;
}
