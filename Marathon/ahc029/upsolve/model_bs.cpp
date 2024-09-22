#pragma GCC target("avx2")
// #pragma GCC optimize("O3")
#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")
#include <bits/stdc++.h>
using namespace std;
// #include <atcoder/all>
// using namespace atcoder;
using ll = long long;
#define _overload3(_1,_2,_3,name,...) name
#define _rep(i,n) repi(i,0,n)
#define repi(i,a,b) for(int i=int(a);i<int(b);++i)
#define rep(...) _overload3(__VA_ARGS__,repi,_rep,)(__VA_ARGS__)
#define _rrep(i,n) rrepi(i,0,n)
#define rrepi(i,a,b) for(int i=int(b)-1;i>=int(a);--i)
#define rrep(...) _overload3(__VA_ARGS__,rrepi,_rrep,)(__VA_ARGS__)
#define all(x) begin(x),end(x)
#define rall(x) rbegin(x),rend(x)
template<class T> bool chmax(T &a, const T &b) { if (a<b) { a=b; return 1; } return 0; }
template<class T> bool chmin(T &a, const T &b) { if (b<a) { a=b; return 1; } return 0; }
template<class T, class U> ostream& operator<<(ostream &os, const pair<T, U> &x){os << x.first << ":" << x.second; return os; }
template<class T, class U> ostream& operator<<(ostream &os, const map<T, U> &x){ for(auto &i: x) os << i << " "; return os; }
template<class T> ostream& operator<<(ostream &os, const vector<T> &x){ for(auto &i: x) os << i << " "; return os; }
template<class T> void print(const T &x){ for(auto &i: x) cout << i << endl; }
#define dout(X) cerr << #X << " " << X << " " << endl

template<class T>
string to_str(const T &x){
    stringstream ss;
    ss << x;
    return ss.str();
}

mt19937 rnd;

unsigned yy=1145141919&1919114514;
inline unsigned xorshift(){yy=yy^(yy<<13);yy=yy^(yy>>17);return yy=yy^(yy<<5);}
#define MASK 65535
inline int randInt(){return (int) (xorshift()&MASK);}

std::chrono::high_resolution_clock::time_point now_time;
inline void set_timer() {
    using namespace std::chrono;
    now_time = high_resolution_clock::now();
}
inline ll get_timer() {
    using namespace std::chrono;
    auto ed = high_resolution_clock::now();
    auto t = ed - now_time;
    return duration_cast<milliseconds>(t).count();
}
inline bool check_timer(int lim) {
    return(get_timer() < lim);
}

constexpr int MAX_INVEST_LEVEL = 20;

struct ZobristHash{
    array<uint_fast64_t, 1001> turn;
    array<uint_fast64_t, 2> turn_type;
    array<uint_fast64_t, 21> invest_level;
    uint_fast64_t money;
    uint_fast64_t project_h, project_v;
    array<uint_fast64_t, 6> card_t;
    uint_fast64_t card_w;
    constexpr ZobristHash()
    :turn(), turn_type(), invest_level(), money(0), project_h(0), project_v(0), card_t(), card_w(0){
        uint_fast64_t x = 123456789ULL;
        for(auto &i: turn) i = xorshift64(x);
        for(auto &i: turn_type) i = xorshift64(x);
        money = xorshift64(x);
        project_h = xorshift64(x);
        project_v = xorshift64(x);
        for(auto &i: card_t) i = xorshift64(x);
        card_w = xorshift64(x);
    }
    constexpr uint_fast64_t xorshift64(uint_fast64_t &x) {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        return x;
    }
    static uint_fast64_t hash(uint_fast64_t x, uint_fast64_t base) {
        x += base;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        return x;
    }
};
constexpr ZobristHash zhash;

struct Project {
    int64_t h; //残務量
    int64_t v; //価値
};
ostream& operator<<(ostream& os, Project pj) {
    os << "Project(" << pj.h << ", " << pj.v << ") ";
    return os;
}

enum class CardType {
    WORK_SINGLE = 0,
    WORK_ALL = 1,
    CANCEL_SINGLE = 2,
    CANCEL_ALL = 3,
    INVEST = 4,
    NONE = 5,
};

ostream& operator<<(ostream& os, CardType ct) {
    switch(ct) {
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
    CardType t; //0:通常労働 1:全力労働 2:キャンセル 3:業務転換 4:増資
    int64_t w; //労働力
    int64_t p; //取得コスト
};
ostream& operator<<(ostream& os, Card cd) {
    os << "Card(" << cd.t << ", " << cd.w << ", " << cd.p << ") ";
    return os;
}

struct Judge {
    const int n; //方針カードの数
    const int m; //プロジェクトの数
    const int k; //追加方針カードの提示数

    Judge(int n, int m, int k): n(n), m(m), k(k) {}

    vector<Card> read_initial_cards() const {
        vector<Card> cards; cards.reserve(n);
        for (int i = 0; i < n; i++) {
            int64_t t, w;
            cin >> t >> w;
            cards.push_back(Card{(CardType)t, w, 0ll});
        }
        return cards;
    }

    vector<Project> read_projects() const {
        static vector<Project> projects; projects.resize(0);
        for (int i = 0; i < m; i++) {
            int64_t h, v;
            cin >> h >> v;
            projects.push_back(Project{h, v});
        }
        return projects;
    }

    static void use_card(int c, int m){
        cout << c << " " << m << endl;
    }

    int64_t read_money() const {
        int64_t money;
        cin >> money;
        return money;
    }

    vector<Card> read_next_cards() const {
        static vector<Card> cards; cards.resize(0);
        for (int i = 0; i < k; i++) {
            int64_t t, w, p;
            cin >> t >> w >> p;
            cards.push_back(Card{(CardType)t, w, p});
        }
        return cards;
    }

    static void select_card(int r) {
        cout << r << endl;
    }

    static void comment(const string& message) {
        cout << "# " << message << endl;
    }
};


enum class TurnType {
    USE_CARD = 0,
    SELECT_CARD = 1,
};

struct Act{
    TurnType type;
    int c, x;
};
ostream& operator<<(ostream& os, const Act& a){
    return os;
}

struct Global {
    int n; //方針カードの数
    int m; //プロジェクトの数
    int k; //追加方針カードの提示数
    int t; //ターン数
    vector<int> x; //方針カードの配布割合
    int total_x;
    int turn;
    vector<vector<Card>> next_cards;
};
Global G;

struct State {
    int turn;
    TurnType turn_type;
    int invest_level;
    int64_t money;
    vector<Project> projects;
    vector<Card> cards;
    vector<Card> next_cards;
    int used_card_idx;
    uint_fast64_t hash;
    State():
        turn(0), turn_type(TurnType::USE_CARD),
        invest_level(0), money(0),
        projects(G.m), 
        cards(G.n), next_cards(G.k), used_card_idx(-1),
        hash(0)
    {
        initHash();
    }

    void initHash(){
        hash ^= zhash.turn[turn];
        hash ^= zhash.turn_type[(int)turn_type];
        hash ^= zhash.invest_level[invest_level];
        hash ^= zhash.hash(money, zhash.money);
        rep(i, G.m) hash ^= zhash.hash(projects[i].h, zhash.project_h);
        rep(i, G.m) hash ^= zhash.hash(projects[i].v, zhash.project_v);
        rep(i, G.n) hash ^= zhash.card_t[(int)cards[i].t];
        rep(i, G.n) hash ^= zhash.hash(cards[i].w, zhash.card_w);
    }
    void setTurn(int turn){
        hash ^= zhash.turn[this->turn];
        hash ^= zhash.turn[turn];
        this->turn = turn;
    }
    void setTurnType(TurnType turn_type){
        hash ^= zhash.turn_type[(int)this->turn_type];
        hash ^= zhash.turn_type[(int)turn_type];
        this->turn_type = turn_type;
    }
    void setInvestLevel(int invest_level){
        hash ^= zhash.invest_level[this->invest_level];
        hash ^= zhash.invest_level[invest_level];
        this->invest_level = invest_level;
    }
    void setMoney(int64_t money){
        hash ^= zhash.hash(this->money, zhash.money);
        hash ^= zhash.hash(money, zhash.money);
        this->money = money;
    }
    void setProject(int idx, const Project &project){
        setProject_h(idx, project.h);
        setProject_v(idx, project.v);
    }
    void setProject_h(int idx, int64_t project_h){
        hash ^= zhash.hash(this->projects[idx].h, zhash.project_h);
        hash ^= zhash.hash(project_h, zhash.project_h);
        this->projects[idx].h = project_h;
    }
    void setProject_v(int idx, int64_t project_v){
        hash ^= zhash.hash(this->projects[idx].v, zhash.project_v);
        hash ^= zhash.hash(project_v, zhash.project_v);
        this->projects[idx].v = project_v;
    }
    void setCard(int idx, const Card &card){
        hash ^= zhash.card_t[(int)this->cards[idx].t];
        hash ^= zhash.hash(this->cards[idx].w, zhash.card_w);
        hash ^= zhash.card_t[(int)card.t];
        hash ^= zhash.hash(card.w, zhash.card_w);
        this->cards[idx] = card;
    }
    
    void advance(const Act &act){
        if(act.type == TurnType::USE_CARD){
            use_card(act.c, act.x);
            read_next_cards();
        } else if(act.type == TurnType::SELECT_CARD){
            select_card(act.c);
            if(money < 0) cerr << "money < 0: " << money << endl;
        }
    }
    inline bool isDone() const{
        return turn == G.t;
    }
    inline bool isDead() const{
        return false;
    }
    const vector<Act>& actions() const{
        static vector<Act> actions; actions.resize(0);
        if(turn_type == TurnType::USE_CARD){
            rep(i, G.n) if(cards[i].t == CardType::INVEST && MAX_INVEST_LEVEL > invest_level){
                actions.push_back(Act{TurnType::USE_CARD, i, 0});
                return actions;
            }

            const int max_n_proj = min(2, G.m);
            static vector<pair<double, int>> target_ps; target_ps.resize(0);
            rep(i, G.m) target_ps.emplace_back(-(projects[i].v-projects[i].h)*pow(0.95, projects[i].h/(1LL<<invest_level)), i);
            sort(all(target_ps));

            set<pair<CardType, int64_t>> cp;
            rep(c, G.n){
                if(cp.count({cards[c].t, cards[c].w})) continue;
                switch(cards[c].t){
                    case CardType::WORK_SINGLE:
                        rep(i, max_n_proj) actions.emplace_back(Act{TurnType::USE_CARD, c, target_ps[i].second});
                        break;
                    case CardType::CANCEL_SINGLE:
                        actions.push_back(Act{TurnType::USE_CARD, c, target_ps.back().second});
                        break;
                    case CardType::WORK_ALL:
                    case CardType::CANCEL_ALL:
                        actions.emplace_back(Act{TurnType::USE_CARD, c, 0});
                        break;
                    case CardType::INVEST:
                        if(MAX_INVEST_LEVEL > invest_level) actions.emplace_back(Act{TurnType::USE_CARD, c, 0});
                        break;
                }
                cp.insert({cards[c].t, cards[c].w});
            }
            shuffle(actions.begin()+1, actions.end(), rnd);
        } else if(turn_type == TurnType::SELECT_CARD){
            set<tuple<CardType, int64_t, int64_t>> cp;
            rep(i, G.k) if(next_cards[i].p <= money) {
                if(cp.count({next_cards[i].t, next_cards[i].w, next_cards[i].p})) continue;
                if(next_cards[i].t == CardType::INVEST && MAX_INVEST_LEVEL <= invest_level) continue;
                actions.emplace_back(Act{TurnType::SELECT_CARD, i, 0});
                cp.insert({next_cards[i].t, next_cards[i].w, next_cards[i].p});
            }
            shuffle(actions.begin()+1, actions.end(), rnd);
        }
        return actions;
    }
    inline double evaluate() const {
        double ret = money;
        ret += (G.t-turn)*(1LL<<invest_level);
        rep(i, G.m) ret += max(0.0, (projects[i].v - projects[i].h) * pow(0.999, projects[i].h/(1LL<<invest_level)));
        return ret;
    }
    
    Project genProject() const{
        Project p;
        p.v = 1LL<<(invest_level+5);
        p.h = round(pow(2, 5-0.25+invest_level));
        return p;
    }

    void workProject(int x, int64_t w){
        setProject_h(x, projects[x].h - w);
        if(projects[x].h <= 0){
            setMoney(money + projects[x].v);
            setProject(x, genProject());
        }
    }    

    void cancelProject(int x){
        setProject(x, genProject());
    }

    const vector<Project>& read_projects() const {
        return projects;
    }

    void use_card(int c, int x) {
        switch (cards[c].t){
            case CardType::WORK_SINGLE:
                workProject(x, cards[c].w);
                break;
            case CardType::WORK_ALL:
                rep(i, G.m) workProject(i, cards[c].w);
                break;
            case CardType::CANCEL_SINGLE:
                cancelProject(x);
                break;
            case CardType::CANCEL_ALL:
                rep(i, G.m) cancelProject(i);
                break;
            case CardType::INVEST:
                setInvestLevel(invest_level+1);
                break;
        }
        used_card_idx = c;
        setCard(used_card_idx, {CardType::NONE, 0, 0});
        setTurnType(TurnType::SELECT_CARD);
    }

    int64_t read_money() const {
        return money;
    }

    Card genCard() const{
        auto rand = uniform_int_distribution(0, G.total_x - 1);
        int r = rand(rnd);
        int t = 0;
        while (r >= G.x[t]) {
            r -= G.x[t];
            t++;
        }
        Card card;
        card.t = CardType(t);
        
        int64_t wd = 25;
        auto r0 = [wd](auto &x){ return wd; };
        auto r2 = [](auto &x){ return 5; };
        auto r4 = [](auto &x){ return 600; };

        card.w = wd*(1LL<<invest_level);
        switch(card.t){
            case CardType::WORK_SINGLE:
                card.p = wd*(1LL<<invest_level)*0.916;
                break;
            case CardType::WORK_ALL:
                card.p = wd*G.m*(1LL<<invest_level)*0.916;
                break;
            case CardType::CANCEL_SINGLE:
            case CardType::CANCEL_ALL:
                card.p = 5*(1LL<<invest_level);
                card.w = 0;
                break;
            case CardType::INVEST:
                card.p = 600*(1LL<<invest_level);
                card.w = 0;
                break;
        }
        return card;
    }

    vector<Card> read_next_cards() {
        if(G.next_cards[0].empty()){
            rep(i, G.next_cards.size()){
                G.next_cards[i].resize(G.k);
                rep(j, G.k) G.next_cards[i][j] = genCard();
            }
        }
        next_cards[0] = {CardType::WORK_SINGLE, 1LL<<invest_level, 0};
        rep(i, 1, G.k) next_cards[i] = G.next_cards[turn-G.turn][i];
        return next_cards;
    }

    void select_card(int r) {
        setCard(used_card_idx, next_cards[r]);
        setMoney(money - next_cards[r].p);
        setTurnType(TurnType::USE_CARD);
        setTurn(turn+1);
    }

    void comment(const string& message) const {
    }

};


struct Node : public std::enable_shared_from_this<Node>{
    State state;
    vector<Act> actions;
    double score;
    int last_action;
    weak_ptr<const Node> parent;
    vector<shared_ptr<Node>> child;
    Node(const State &s)
    : state(s), last_action(-1), parent()
    {
        actions = state.actions();
        score = state.evaluate();
        score = score * 100 + last_action;
    }
    shared_ptr<Node> advanced(const int act_id) const{
        auto s = state;
        s.advance(actions[act_id]);
        auto node = make_shared<Node>(s);
        node->parent = shared_from_this();
        node->last_action = act_id;
        return node;
    }
    void expandChild(){
        rep(i, actions.size()){
            child.emplace_back(advanced(i));
        }
    }
};
inline bool operator<(const shared_ptr<Node> &a, const shared_ptr<Node> &b)
{
    return a->score > b->score;
}

const vector<Act>& beamSearchAction(shared_ptr<Node> state, const size_t beam_width, const size_t beam_depth, const int time_limit)
{
    using PNode = shared_ptr<Node>;
    using CPNode = shared_ptr<const Node>;
    priority_queue<PNode> now_beam;
    now_beam.emplace(state);
    CPNode best_node;
    CPNode best_done_node;

    [&](){
        rep(t, beam_depth){
            if(!check_timer(time_limit)) return;
            priority_queue<PNode> next_beam;
            unordered_set<uint_fast64_t> hash_set;
            while(!now_beam.empty()){
                if(!check_timer(time_limit)) return;
                const auto now_node = now_beam.top();
                now_beam.pop();
                now_node->expandChild();
                for(auto &c: now_node->child){
                    // if(!check_timer(time_limit)) return;
                    if(c->state.isDead()){c.reset(); continue;}
                    if(next_beam.size()>=beam_width && next_beam.top()->score > c->score){c.reset(); continue;}
                    if(c->state.isDone()){
                        if(!best_done_node || c->state.money > best_done_node->state.money) best_done_node = c;
                        continue;
                    }
                    // if(!best_node || c->score > best_node->score) best_node = c;
                    if(hash_set.count(c->state.hash)){c.reset(); continue;}
                    next_beam.emplace(c);
                    hash_set.insert(c->state.hash);
                    if(next_beam.size()>beam_width){
                        hash_set.erase(next_beam.top()->state.hash);
                        next_beam.pop();
                    }
                }
            }
            if(best_done_node) break;
            now_beam.swap(next_beam);
        }
    }();

    // cerr << "now_beam_size = " << now_beam.size() << endl;
    if(now_beam.size()>1){
        while(!now_beam.empty()){
            best_node = now_beam.top();
        now_beam.pop();
        }
    }
    if(best_done_node) best_node = best_done_node;
    static vector<Act> actions;
    actions.resize(0);
    if(!best_node){
        actions.emplace_back(state->state.actions()[0]);
        return actions;
    }
    while(auto p = best_node->parent.lock()){
        actions.emplace_back(p->actions[best_node->last_action]);
        best_node = p;
    }
    reverse(actions.begin(), actions.end());
    return actions;
    // return best_node->actions[actions.back()];
}

int BEAM_WIDTH = 100;
constexpr int BEAM_DEPTH = 10;
constexpr int TIME_LIMIT = 2000-20;
// constexpr int TIME_LIMIT = 4000-100;
// constexpr int TIME_LIMIT = 6000-100;
// constexpr int TIME_LIMIT = 20000-100;

struct Solver {
    const int n;
    const int m;
    const int k;
    const int t;
    Judge judge;
    int invest_level;
    int64_t money;
    vector<Project> projects;
    vector<Card> cards;

    State state;
    vector<Act> actions;

    Solver(int n, int m, int k, int t) : 
        n(n), m(m), k(k), t(t), 
        judge(n, m, k), 
        invest_level(0), money(0),
        projects(), cards(),
        state()
    {
    }

    int64_t solve() {
        cards = judge.read_initial_cards();
        state.cards = cards;
        projects = judge.read_projects();
        state.projects = projects;
        cerr << "start" << endl;

        for (G.turn = 0; G.turn < t; ++G.turn) {
            auto [use_card_i, use_target] = select_action();
            const Card& use_card = cards[use_card_i];
            if (use_card.t == CardType::INVEST) invest_level++;

            // example for comments
            stringstream msg;
            // msg << get_timer() << " ms "; 
            msg << "used Card(t=" << use_card.t << " w=" << use_card.w << " p=" << use_card.p << ") to target " << use_target;
            // cerr << "used Card(t=" << use_card.t << " w=" << use_card.w << " p=" << use_card.p << ") to target " << use_target << endl;
            judge.comment(msg.str());
            judge.use_card(use_card_i, use_target);
            state.use_card(use_card_i, use_target);
            assert(invest_level <= MAX_INVEST_LEVEL);

            projects = judge.read_projects();
            // cerr << "projects:" << projects << endl;
            // cerr << "projects:" << state.projects << endl;
            rep(i, m) state.setProject(i, projects[i]);
            // state.projects = projects;
            money = judge.read_money();
            // cerr << "money: " << money << " " << state.money << endl;
            assert(state.money == money);
            if(G.turn == t-1){
                cerr << "time = " << get_timer() << " ms" << endl;
                cerr << "score = " << state.money << endl;
            }

            vector<Card> next_cards = judge.read_next_cards();
            // cerr << "next_cards:" << next_cards << endl;
            state.next_cards = next_cards;
            int select_card_i = select_next_card(next_cards);
            cards[use_card_i] = next_cards[select_card_i];
            judge.select_card(select_card_i);
            state.select_card(select_card_i);
            // cerr << "cards:" << cards << endl;
            // cerr << "cards:" << state.cards << endl;
            money -= next_cards[select_card_i].p;
            assert(money >= 0);

            rep(i, 1, k) G.x[(int)next_cards[i].t]++;
            G.total_x += k-1;

            G.next_cards[0].resize(0);
        }
        return money;
    }

    pair<int, int> select_action() {
        auto node = make_shared<Node>(state);
        if(actions.size()>1){
            return {actions[1].c, actions[1].x};
        }else{
            actions = beamSearchAction(node, BEAM_WIDTH, BEAM_DEPTH, TIME_LIMIT);
            return {actions[0].c, actions[0].x};
        }
    }

    int select_next_card(const vector<Card>& next_cards) {
        auto node = make_shared<Node>(state);
        if(state.turn > G.t*0.3) BEAM_WIDTH = clamp(BEAM_WIDTH + (TIME_LIMIT*state.turn/G.t*0.95 < get_timer() ? -1 : 1), 70, 200);
        actions = beamSearchAction(node, BEAM_WIDTH, BEAM_DEPTH, TIME_LIMIT);
        return actions[0].c;
    }
};

int main() {
    int n, m, k, t;
    cin >> n >> m >> k >> t;
    cerr << n << " " << m << " " << k << " " << t << endl;
    G.n = n;
    G.m = m;
    G.k = k;
    G.t = t;
    G.x = {21, 11, 11, 6, 4};
    G.total_x = accumulate(G.x.begin(), G.x.end(), 0);
    G.next_cards.resize(BEAM_DEPTH);
    set_timer();
    Solver solver(n, m, k, t);
    int64_t score = solver.solve();
    // cerr << "time = " << get_timer() << " ms" << endl;
    // cerr << "score = " << score << endl;
}
