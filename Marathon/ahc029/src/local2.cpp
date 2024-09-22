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
          project.h = t_projects[t_project_idx].h * pow(2, L);
          project.v = t_projects[t_project_idx].v * pow(2, L);
          t_project_idx++;
        }
      }
    }
    if (use_card.t == CardType::CANCEL_SINGLE) {
      projects[use_target].h = t_projects[t_project_idx].h * pow(2, L);
      projects[use_target].v = t_projects[t_project_idx].v * pow(2, L);
      t_project_idx++;
    }
    if (use_card.t == CardType::CANCEL_ALL) {
      for (auto &project : projects) {
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
    rep(i, 0, projects.size()) {
      auto project = projects[i];
      if (projects[m].v / projects[m].h < project.v / project.h) {
        m = i;
      }
    }

    ll best_score = -1;
    rep(i, 0, cards.size()) {
      auto use_card = cards[i];
      ll L = invest_level;
      ll score = 0;
      if (use_card.t == CardType::WORK_SINGLE) {
        score = use_card.w;
      }
      if (use_card.t == CardType::INVEST) {
        if (L < 20) score = INFF;
      }

      if (chmax(best_score, score)) {
        c = i;
      }
    }

    if (cards[c].t != CardType::WORK_SINGLE && cards[c].t != CardType::CANCEL_SINGLE) m = 0;
    return {c, m};
  }

  ll efficiency_project(vector<Project> &projects) {
    ll idx = 0;
    rep(i, 0, projects.size()) {
      if (projects[idx].v / projects[idx].h < projects[i].v / projects[i].h) idx = i;
    }
    return idx;
  }

  ll select_next_card(const vector<Card> &next_cards, vector<Project> &projects, ll money, ll turn, ll invest_level) {
    ll next_card_idx = 0;
    ll m = efficiency_project(projects);
    ll L = invest_level;
    double best_score = 0;
    rep(i, 1, next_cards.size()) {
      if (money - next_cards[i].p < 0) continue;

      if (next_cards[i].t == CardType::WORK_SINGLE) {
        ll all_h = projects[m].h;
        // use
        ll use_day = max((ll)0, (ll)((all_h - next_cards[i].w) / pow(2, L))) + 1;
        double use_efficiency = (projects[m].v - next_cards[i].p) / use_day;
        // not_use
        double not_use_efficiency = projects[m].v / (all_h * pow(2, L));

        if (use_efficiency < not_use_efficiency) continue;
        if (chmax(best_score, next_cards[i].w / (double)next_cards[i].p)) {
          next_card_idx = i;
        }
      }
      if (next_cards[i].t == CardType::WORK_ALL) continue;
      if (next_cards[i].t == CardType::CANCEL_SINGLE) continue;
      if (next_cards[i].t == CardType::CANCEL_ALL) continue;
      if (next_cards[i].t == CardType::INVEST) {
        if (500 < turn) continue;
        next_card_idx = i;
        best_score = INFF;
      }
    }

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
      cout << "# current products ";
      for (auto cu : projects) {
        cout << "{" << cu.h << ", " << cu.v << "}";
      }
      cout << endl;

      auto [use_card_i, use_target] = select_action(cards, projects, money);
      cout << "# select_action " << cards[use_card_i].t << " " << to_string(use_target) << endl;

      cout << use_card_i << " " << use_target << endl;
      cout.flush();
      assert(invest_level <= MAX_INVEST_LEVEL);

      judge.after_use_action(use_card_i, use_target, cards, projects, money, invest_level);

      vector<Card> next_cards = judge.read_next_cards(invest_level, turn);
      ll select_card_i = select_next_card(next_cards, projects, money, turn, invest_level);
      cards[use_card_i] = next_cards[select_card_i];
      cout << "# next_card " << next_cards[select_card_i].t << " " << next_cards[select_card_i].p << endl;
      cout << select_card_i << endl;
      cout.flush();

      money -= next_cards[select_card_i].p;
      assert(money >= 0);
      cout << "# money " << money << endl;
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
