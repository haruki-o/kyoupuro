#include <bits/stdc++.h>
using namespace std;
using namespace chrono;

typedef long long ll;
typedef vector<ll> vll;
typedef vector<vll> vvll;
#define rep(i, l, n) for (ll i = (ll)(l); i < (ll)(n); i++)
#define Sort(a) sort(a.begin(), a.end())

// 下,右,上,左
vll dy = {1, 0, -1, 0};
vll dx = {0, 1, 0, -1};

int main() {
  ll num = 100;
  vvll diff(200, vll(200, 0));
  vvll ban(200, vll(200));
  rep(i, 77, 78) {
    ifstream input_fstream;
    string S = ("0000" + to_string(i));
    S = S.substr((int)S.size() - 4, 4);

    input_fstream.open("../in/in0100/" + S + ".txt", std::ios_base::in);

    ll dammy;
    rep(i, 0, 4) input_fstream >> dammy;

    rep(j, 0, 200) {
      rep(k, 0, 200) {
        ll s;
        input_fstream >> s;
        ban[j][k] = s;
      }
    }

    rep(j, 0, 200) {
      rep(k, 0, 200) {
        // 辺の長さ
        ll sur = 10;
        ll sum = 0;

        rep(dir,0,4){
          ll toj = j + dy[dir] * sur;
          ll tok = k + dx[dir] * sur;

          if(toj < 0 || 200 <= toj || tok < 0 || 200 <= tok)continue;
          sum++;
          ll pred = (ban[j][k] + ban[toj][tok]) / 2;
          pred *= abs(j - toj) + abs(k - tok) + 1;

          ll correct = 0;
          ll cuj = j;
          ll cuk = k;
          while(1){
            correct += ban[cuj][cuk];

            if(cuj == toj && cuk == tok)break;

            cuj += dy[dir];
            cuk += dx[dir];
          }

          diff[j][k] += abs(pred - correct);
        }

        diff[j][k] /= sum;

        cout << diff[j][k] << " ";
      }
      cout << endl;
    }
    input_fstream.close();
  }
}

// 2点を調べ、一次関数と見た時、平均 * 辺の長さでずれを確認
// マスmのS_ijを予想した時、どのくらいの範囲なら適用してよいか知るため