//n!通りの順列
//n==5,5!を例にする
//factoは5!通りの数字保管庫ex)facto[43]=34210;
//


vecit facto;
veci facto_index(10);
//n=0000000111 0.1.2は使ってません
void factorial(int n){
  if(n==0){
    intt u=0;
    rep(i,0,10){
      u+=pow(10,9-i)*(facto_index[i]-1);
    }
    facto.push_back(u);
  }
  rep(i,0,10){
    if((1<<i) & n){
      facto_index[i]=__builtin_popcount(n);
      factorial(n & ~(1<<i));
      facto_index[i]=0;
    }
  }
}
//ここまででfactoには保管終了

//factoをint->stringに変更subject=054213が文字列で扱える
rep(i,0,facto.size()){
    string subject="";
    if(facto[i]<pow(10,9)){
      rep(j,1,10){
      	subject+=facto[i]/(int)pow(10,10-j)%10+'0';
      }
      subject+=facto[i]%10+'0';
    }else{
      rep(j,1,10){
        subject+=facto[i]/(int)pow(10,10-j)%10+'0';
      }
      subject+=facto[i]%10+'0';
    }
}
