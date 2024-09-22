#include <bits/stdc++.h>
using namespace std;
using namespace chrono;

typedef long long ll;
typedef vector<ll> vll;
typedef vector<vll> vvll;
#define rep(i, l, n) for (ll i = (ll)(l); i < (ll)(n); i++)
#define Sort(a) sort(a.begin(), a.end())

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
        // 幅(2の時5 * 5になる)
        ll sur = 5;
        ll sum = 0;

        rep(_i, j - sur, j + sur) {
          rep(_j, k - sur, k + sur) {
            if(_i == j && _j == k)continue;
            if (_i < 0 || 200 <= _i)
              continue;
            if (_j < 0 || 200 <= _j)
              continue;
            sum++;
            diff[j][k] += abs(ban[_i][_j] - ban[j][k]);
          }
        }

        diff[j][k] /= sum;

        cout << diff[j][k] << " ";
      }
      cout << endl;
    }
    input_fstream.close();
  }
}

// マスmの周囲nマスとの差異を調べる
// マスmのS_ijを予想した時、どのくらいの範囲なら適用してよいか知るため