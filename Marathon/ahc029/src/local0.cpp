#include <atcoder/all>
#include <bits/stdc++.h>
using namespace std;
using namespace atcoder;
using namespace chrono;

typedef long long ll;
typedef pair<ll, ll> pll;
typedef pair<pll, ll> ppll;
typedef vector<ll> vll;
typedef vector<vll> vvll;
typedef vector<vvll> vvvll;
typedef vector<ppll> vpp;
typedef vector<pll> vp;
typedef vector<vp> vvp;
typedef vector<char> vc;
typedef vector<int> vi;
typedef vector<string> vs;
typedef vector<vs> vvs;
#define rep(i, l, n) for (ll i = (ll)(l); i < (ll)(n); i++)
#define repd(i, n, l) for (ll i = (ll)(n); i > (ll)(l); i--)
#define Sort(a) sort(a.begin(), a.end())
#define gSort(a) sort(a.begin(), a.end(), greater())
template <class T, class S>
inline bool chmax(T &a, const S &b) {
  return (a < b ? a = b, 1 : 0);
}
template <class T, class S>
inline bool chmin(T &a, const S &b) {
  return (a > b ? a = b, 1 : 0);
}
#define INFF (9223372036854775800)
#define TIME_LIMIT (3.5)
#define def (20101)
// #define MOD (1000000007)
// #define MAX (2147483647)
#define MAX (1073741823)
#define MOD (998244353)

double stat_time = 0.0;
double stat_loop = 0;
ll stat_ans = 0;

double ela_times(system_clock::time_point &clock) {
  return duration_cast<microseconds>(system_clock::now() - clock).count() * 1e-6;
}

ll rand(ll l, ll r) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dis(l, r);
  return (ll)dis(gen);
}

double rand_double(double l, double r) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(l, r);
  return dis(gen);
}

double clamp(double x, double l, double r) {
  return max(l, min(r, x));
}

double gauss(double mu, double sigma) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<double> dis(mu, sigma);
  return dis(gen);
}

const ll MAX_INVEST_LEVEL = 20;

struct Project {
  ll h;
  ll v;
};

enum class CardType {
  WORK_SINGLE = 0,
  WORK_ALL = 1,
  CANCEL_SINGLE = 2,
  CANCEL_ALL = 3,
  INVEST = 4,
};

ostream &operator<<(ostream &os, CardType ct) {
  switch (ct) {
    case CardType::WORK_SINGLE:
      os << "WORK_SINGLE";
      break;
    case CardType::WORK_ALL:
      os << "WORK_ALL";
      break;
    case CardType::CANCEL_SINGLE:
      os << "CANCEL_SINGLE";
      break;
    case CardType::CANCEL_ALL:
      os << "CANCEL_ALL";
      break;
    case CardType::INVEST:
      os << "INVEST";
      break;
  }
  return os;
}

struct Card {
  CardType t;
  ll w;
  ll p;
};

struct Judge {
  ll n;
  ll m;
  ll k;
  ll t;
  ll x_0, x_1, x_2, x_3, x_4;
  vector<Project> t_projects;
  ll t_project_idx;
  vector<vector<Card>> t_cards;

  Judge(ll n, ll m, ll k, ll t) : n(n), m(m), k(k), t(t) {
    x_0 = rand(1, 20);
    x_1 = rand(1, 10);
    x_2 = rand(1, 10);
    x_3 = rand(1, 5);
    x_4 = rand(1, 3);
    t_project_idx = 0;
    t_cards.resize(t);
  }

  ll generateCardType() {
    // 確率分布の総和を計算
    int totalProbability = 0;
    vll probabilities = {x_0, x_1, x_2, x_3, x_4};
    for (int probability : probabilities) {
      totalProbability += probability;
    }

    // 0から総和-1までの一様乱数を生成
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, totalProbability - 1);
    int randomValue = dis(gen);

    // 累積確率を計算して生成された乱数に対応するカードの種類を決定
    int cumulativeProbability = 0;
    for (int i = 0; i < probabilities.size(); ++i) {
      cumulativeProbability += probabilities[i];
      if (randomValue < cumulativeProbability) {
        return (ll)i;
      }
    }

    return -1;
  }

  void init_input(vector<Card> &cards, vector<Project> &projects) {
#ifdef DEBUG
    rep(i, 0, m) {
      ll h, v;
      cin >> h >> v;
      projects.push_back(Project{h, v});
    }
    rep(i, 0, t * m) {
      ll h, v;
      cin >> h >> v;
      t_projects.push_back(Project{h, v});
    }
    rep(i, 0, n) {
      ll t, w;
      cin >> t >> w;
      cards.push_back(Card{(CardType)t, w, 0ll});
    }
    rep(i, 0, t) {
      rep(j, 0, k) {
        ll t, w, p;
        cin >> t >> w >> p;
        t_cards[i].push_back(Card{(CardType)t, w, p});
      }
    }
#else
    rep(i, 0, n) {
      ll t, w;
      cin >> t >> w;
      cards.push_back(Card{(CardType)t, w, 0ll});
    }
    rep(i, 0, m) {
      ll h, v;
      cin >> h >> v;
      projects.push_back(Project{h, v});
    }
#endif
  }

  void after_use_action(ll use_card_i, ll use_target, vector<Card> &cards, vector<Project> &projects, ll &money, ll &invest_level) {
    const Card &use_card = cards[use_card_i];
    if (use_card.t == CardType::INVEST) {
      invest_level++;
    }
#ifdef DEBUG
    ll L = invest_level;
    if (use_card.t == CardType::WORK_SINGLE) {
      projects[use_target].h -= use_card.w;
      if (projects[use_target].h <= 0) {
        money += projects[use_target].v;

        // double b = rand_double(2.0, 8.0);
        // ll h = round(pow(2, b)) * pow(2, L);
        // ll v = round(pow(2, clamp(gauss(b, 0.5), 0.0, 10.0))) * pow(2, L);
        // projects[use_target].h = h;
        // projects[use_target].v = v;
        projects[use_target].h = t_projects[t_project_idx].h * pow(2, L);
        projects[use_target].v = t_projects[t_project_idx].v * pow(2, L);
        t_project_idx++;
      }
    }
    if (use_card.t == CardType::WORK_ALL) {
      for (auto &project : projects) {
        project.h -= use_card.w;
        if (project.h <= 0) {
          money += project.v;

          // double b = rand_double(2.0, 8.0);
          // ll h = round(pow(2, b)) * pow(2, L);
          // ll v = round(pow(2, clamp(gauss(b, 0.5), 0.0, 10.0))) * pow(2, L);
          // project.h = h;
          // project.v = v;
          project.h = t_projects[t_project_idx].h * pow(2, L);
          project.v = t_projects[t_project_idx].v * pow(2, L);
          t_project_idx++;
        }
      }
    }
    if (use_card.t == CardType::CANCEL_SINGLE) {
      // double b = rand_double(2.0, 8.0);
      // ll h = round(pow(2, b)) * pow(2, L);
      // ll v = round(pow(2, clamp(gauss(b, 0.5), 0.0, 10.0))) * pow(2, L);
      // projects[use_target].h = h;
      // projects[use_target].v = v;
      projects[use_target].h = t_projects[t_project_idx].h * pow(2, L);
      projects[use_target].v = t_projects[t_project_idx].v * pow(2, L);
      t_project_idx++;
    }
    if (use_card.t == CardType::CANCEL_ALL) {
      for (auto &project : projects) {
        // double b = rand_double(2.0, 8.0);
        // ll h = round(pow(2, b)) * pow(2, L);
        // ll v = round(pow(2, clamp(gauss(b, 0.5), 0.0, 10.0))) * pow(2, L);
        // project.h = h;
        // project.v = v;
        project.h = t_projects[t_project_idx].h * pow(2, L);
        project.v = t_projects[t_project_idx].v * pow(2, L);
        t_project_idx++;
      }
    }
#else
    projects.clear();
    rep(i, 0, m) {
      ll h, v;
      cin >> h >> v;
      projects.push_back(Project{h, v});
    }
    cin >> money;
#endif
  }

  vector<Card> read_next_cards(ll invest_level, ll turn) {
    vector<Card> cards;
#ifdef DEBUG
    ll L = invest_level;
    rep(j, 0, k) cards.push_back(t_cards[turn][j]);
    rep(i, 0, k) {
      cards[i].w *= pow(2, L);
      cards[i].p *= pow(2, L);
    }
    // rep(i, 0, k) {
    //   if (i == 0) {
    //     cards.push_back(Card{(CardType)0, (ll)pow(2, L), 0});
    //   } else {
    //     ll t = generateCardType();
    //     ll w = 0;
    //     if (t == 0 || t == 1) w = rand(1, 50) * pow(2, L);
    //     ll p;
    //     if (t == 0) {
    //       ll wd = w / pow(2, L);
    //       p = clamp(round(gauss(wd, wd / 3)), 1, 10000) * pow(2, L);
    //     }
    //     if (t == 1) {
    //       ll wd = w / pow(2, L);
    //       p = clamp(round(gauss(wd * m, wd * m / 3)), 1, 10000) * pow(2, L);
    //     }
    //     if (t == 2) {
    //       p = rand(0, 10) * pow(2, L);
    //     }
    //     if (t == 3) {
    //       p = rand(0, 10) * pow(2, L);
    //     }
    //     if (t == 4) {
    //       p = rand(200, 1000) * pow(2, L);
    //     }
    //     cards.push_back(Card{(CardType)t, w, p});
    //   }
#else
    rep(i, 0, k) {
      ll t, w, p;
      cin >> t >> w >> p;
      cards.push_back(Card{(CardType)t, w, p});
    }
#endif

    return cards;
  }
};

struct Solver {
  ll n;
  ll m;
  ll k;
  ll t;
  Judge judge;
  ll turn;
  ll invest_level;
  ll money;
  vector<Project> projects;
  vector<Card> cards;

  Solver(ll n, ll m, ll k, ll t) : n(n), m(m), k(k), t(t), judge(n, m, k, t), turn(0), invest_level(0), money(0) {
  }

  pll select_action(vector<Card> &cards, vector<Project> &projects, ll money) {
    ll c = 0, m = 0;
    ll best_score = 0;
    rep(i, 0, cards.size()) {
      auto use_card = cards[i];
      ll L = invest_level;
      ll score = 0;
      ll _m = 0;
      if (use_card.t == CardType::WORK_SINGLE) {
        score = use_card.w;
        rep(j, 0, projects.size()) {
          auto project = projects[j];
          _m = j;
          if (project.h - use_card.w <= 0) {
            score *= 3;
            break;
          }
        }
      }
      if (use_card.t == CardType::WORK_ALL) {
        score = use_card.w * cards.size();
        rep(j, 0, projects.size()) {
          auto project = projects[j];
          _m = j;
          if (project.h - use_card.w <= 0) {
            score += project.h - use_card.w;
            break;
          }
        }
      }
      if (use_card.t == CardType::CANCEL_SINGLE) {
        score = 0;
      }
      if (use_card.t == CardType::CANCEL_ALL) {
        for (auto &project : projects) {
          score += project.h - project.v;
        }
      }
      if (use_card.t == CardType::INVEST) {
        if (L < 20) score = INFF;
      }

      if (chmax(best_score, score)) {
        c = i;
        m = _m;
      }
    }

    if (cards[c].t != CardType::WORK_SINGLE && cards[c].t != CardType::CANCEL_SINGLE) m = 0;
    return {c, m};
  }

  ll select_next_card(const vector<Card> &next_cards, ll money) {
    ll next_card_idx = 0;
    vll ava;
    string msg = "";
    rep(i, 0, next_cards.size()) {
      msg += "{" + to_string(next_cards[i].w) + " " + to_string(next_cards[i].p) + "}";
      if (money - next_cards[i].p < 0) continue;
      ava.push_back(i);
    }
    cout << "# " << msg << endl;
    assert(ava.size() > 0);
    next_card_idx = ava[rand() % ava.size()];
    assert(next_card_idx < k);

    return next_card_idx;
  }

  ll solve() {
    judge.init_input(cards, projects);

    rep(turn, 0, t) {
      cout << "# current cards ";
      for (auto cu : cards) {
        cout << "{" << cu.t << ", " << to_string(cu.p) << "}";
      }
      cout << endl;

      auto [use_card_i, use_target] = select_action(cards, projects, money);
      cout << "# select_action " << cards[use_card_i].t << " " << to_string(use_target) << endl;

      cout << use_card_i << " " << use_target << endl;
      cout.flush();
      assert(invest_level <= MAX_INVEST_LEVEL);

      judge.after_use_action(use_card_i, use_target, cards, projects, money, invest_level);

      vector<Card> next_cards = judge.read_next_cards(invest_level, turn);
      ll select_card_i = select_next_card(next_cards, money);
      cards[use_card_i] = next_cards[select_card_i];
      cout << "# next_card " << next_cards[select_card_i].t << endl;
      cout << select_card_i << endl;
      cout.flush();

      money -= next_cards[select_card_i].p;
      assert(money >= 0);
      cerr << "turn : " << turn << " " << money << endl;
    }

    cerr << invest_level << endl;
    return money;
  }
};

int main() {
  srand((unsigned int)time(NULL));
  ll n, m, k, t;
  cin >> n >> m >> k >> t;
  Solver solver(n, m, k, t);
  ll score = solver.solve();
  cerr << "score:" << score << endl;
}
