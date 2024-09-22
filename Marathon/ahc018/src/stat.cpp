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
  vvll sum(200, vll(200, 0));
  vll bins(50, 0);
  rep(i, 77, 78) {
    ifstream input_fstream;
    string S = ("0000" + to_string(i));
    S = S.substr((int)S.size() - 4, 4);

    input_fstream.open("../in/in0100/" + S + ".txt", std::ios_base::in);

    vll dammy(4);
    rep(j, 0, 4) input_fstream >> dammy[j];
    vll all;
    ll all_sum = 0;
    rep(j, 0, 200) {
      rep(k, 0, 200) {
        ll s;
        input_fstream >> s;
        sum[j][k] += s;
        bins[s / (ll)100]++;
        // bins[s / (ll)100]++;
        all.push_back(s);
        all_sum += s;
        cout << s << ",";
      }
      cout << endl;
    }
    Sort(all);
    // cout << S << " W, K : " << dammy[1] << " " << dammy[2] << endl;
    // cout << all[200 * 200 / 4] << " " << all[200 * 200 / 2] << " "
    //      << all[200 * 200 * 3 / 4] << endl;
    // cout << all_sum / 200 / 200 << endl;
    input_fstream.close();
  }
  rep(i, 0, 200) rep(j, 0, 200) sum[i][j] /= num;

  // rep(i, 0, 200) {
  //   rep(j, 0, 20) cout << sum[i][j] / 10 << " ";
  //   cout << endl;
  // }

  // rep(i, 0, 50) cout << bins[i] << endl;
}