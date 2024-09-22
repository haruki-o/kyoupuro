#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

int main(int argc, char *argv[]) {
  vector<pair<int, string>> sto(4, pair<int, string>(0, ""));
  int ind = 0;
  int m = 0;
  while (!cin.eof()) {
    string s;
    cin >> s;
    int swi = 0;
    string num_string = "";
    for (int i = 0; i < (int)s.size(); i++) {
      if (s[i] == ':') {
        sto[ind].first = stoi(num_string);
        swi = 1;
      } else if (swi == 0) {
        num_string += s[i];
      } else if (swi == 1) {
        sto[ind].second += s[i];
      }
      if (i == (int)s.size() - 1)
        m = stoi(num_string);
    }
    // cout << sto[ind].first << " " << sto[ind].second<< endl;
    ind++;
  }
  // cout << "m " << m << endl;
  sort(sto.begin(), sto.end());

  string ans = "";
  for (int i = 0; i < 4; i++) {
    if (sto[i].first == 0)
      continue;
    if (m % sto[i].first == 0)
      ans += sto[i].second;
  }
  if (ans != "")
    cout << ans << endl;
  else
    cout << m << endl;
  return 0;
}
