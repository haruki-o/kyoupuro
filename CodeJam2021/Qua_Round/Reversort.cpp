
  int T;cin>>T;
  rep(i,0,T){
    int N;cin>>N;
    vi a(N);
    rep(j,0,N)cin>>a[j];
    int sum = 0;
    rep(j,0,N-1){
      rep(k,j,N){
        if(a[k] == j+1){
          reverse(a.begin() + j,a.begin() + k + 1);
          break;
        }
        sum++;
      }
      sum++;
      // rep(k,0,N)cout << a[k] << " ";
      // cout << " sum : " << sum << endl;
    }
    cout << "Case #" << i+1 << ": " << sum << endl;
  }