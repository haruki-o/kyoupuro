
constexpr int N = 100;			
constexpr int M_Max = 300;		
constexpr int K_Max = 5000;		

int M = 0;
int K = 0;

template <class T>
using NArr = CapArr<T, N>;

template <class T>
using MArr = CapArr<T, M_Max>;

template <class T>
using KArr = CapArr<T, K_Max>;




struct Result {
	NArr<int> powers;	
	MArr<bool> ons;		
};

struct Edge {
	array<int, 2> uv;
	int w;			
};

struct Around {
	int vi;
	int ei;
	int w;
};

struct IOServer {
	NArr<pint> srcs_;
	MArr<Edge> edges_;
	KArr<pint> dsts_;
	NArr<CapArr<Around, 32>> arounds_;		

	void InitInput(ChronoTimer& timer) {
		istream& is = cin;
		int dummy;
		is >> dummy;	
		timer.Init();		

		is >> M >> K;

		srcs_.resize(N);
		REP(i, N) {
			is >> srcs_[i].x >> srcs_[i].y;
		}

		edges_.resize(M);
		arounds_.resize(N);
		REP(i, M) {
			auto& e = edges_[i];
			is >> e.uv[0] >> e.uv[1] >> e.w;
			--e.uv[0];
			--e.uv[1];

			{
				auto& a = arounds_[e.uv[0]].push();
				a.ei = i;
				a.vi = e.uv[1];
				a.w = e.w;
			}
			{
				auto& a = arounds_[e.uv[1]].push();
				a.ei = i;
				a.vi = e.uv[0];
				a.w = e.w;
			}
		}

		dsts_.resize(K);
		REP(i, K) {
			is >> dsts_[i].x >> dsts_[i].y;
		}
	}

	void Output(const Result& result) {
		ostream& os = cout;
		REP(i, N) {
			os << result.powers[i] << " ";
		}
		os << endl;
		REP(i, M) {
			os << (result.ons[i] ? 1 : 0) << " ";
		}
		os << endl;

	}

	void Finalize() {
	}
};
IOServer server;

#define USE_SA_POINT_FILTER 1
#define USE_SA_ROLLBACK 1

#define DOUBLE_PHASE 0
#define GREADY_UPDATE 0
#define USE_LEVEL_CHANGE_COST 0



struct RandomTable {
	vector<int> table_;

	void push(int value, int count) {
		table_.reserve(table_.size() + count);
		REP(i, count) {
			table_.emplace_back(value);
		}
	}

	template <class ENGINE>
	int operator()(ENGINE& engine) {
		return table_[engine() % table_.size()];
	}
};

#define USE_ACCEPT_SCORE 0



struct SAChecker {
	Xor64* rand_ = nullptr;
	double* totalMaxScore_ = nullptr;

	double temp = 0;		
	double divTemp = 0;

	double currentScore = 0;
	double maxScore = 0;

	int noMaxUpdateCount = 0;				
	int nextRollbackCheckCount = INT_MAX;	

	inline bool operator()(double newScore) {
		++noMaxUpdateCount;

		if (newScore > currentScore) {
			currentScore = newScore;
			if (newScore > maxScore) {
				maxScore = newScore;
				noMaxUpdateCount = 0;

				if (newScore > *totalMaxScore_) {
					*totalMaxScore_ = newScore;
				}
			}

			return true;
		}

		else if (newScore == currentScore) {
			return true;
		}

		else {
			if (exp((newScore - currentScore) * divTemp) * UINT32_MAX > (*rand_)(UINT32_MAX)) {
				currentScore = newScore;
				return true;
			}
			else {
				return false;
			}
		}
	}

};

template <class F>
struct SATransition {
	const char* name;
	F func;
	int weight;
};
template <class F>
auto MakeTransition(const char* name, F&& func, int weight) {
	return SATransition<F>{ name, func, weight };
}
#define MAKE_TRANS(func, weight) MakeTransition(#func, [&](SAChecker& sac, State& state) { func(sa, sac, state); }, weight)

struct SimulatedAnnealing {
	vector<SAChecker> checkers;

	double totalMaxScore = 0;				
	double timeRate = 0;				
	double temp = 0;					
	double divTemp = 0;					


	Xor64 rand_;

	double startTemp_ = 200;					
	double endTemp_ = 1;						
	int stepLoopCount = 1000;					

	double rollbackStartRate_ = 999.0;			
	int rollbackCount_ = INT_MAX;				
	double rollbackNextMulti_ = 1.1;			


public:
	template <class STATE, class... TRANSITION>
	void Run2(ChronoTimer& timer, vector<STATE>& states, tuple<SATransition<TRANSITION>...>& transitions) {

		vector<STATE> maxStates = states;
		checkers.resize(states.size());
		totalMaxScore = states[0].EvalScore();
		REP(pi, checkers.size()) {
			auto& checker = checkers[pi];
			checker.rand_ = &rand_;
			checker.totalMaxScore_ = &totalMaxScore;

			checker.temp = 0;
			checker.divTemp = 0;
			checker.currentScore = states[pi].EvalScore();
			checker.maxScore = checker.currentScore;
			checker.noMaxUpdateCount = 0;
			checker.nextRollbackCheckCount = rollbackCount_;

			chmax(totalMaxScore, checker.maxScore);
		}

		RandomTable randTable;
		TupleLoop(transitions, [&](auto&& e, size_t i) {
			randTable.push((int)i, e.weight);
		});

		const auto startTime = timer.Now();
		const auto endTime = timer.EndTime();
		const double subTimeCountDiv = 1.0 / (double)(endTime - startTime).count();

		vector<int> pis(states.size());
		iota(ALL(pis), 0);

		bool forceEnd = false;
		while (!timer.IsOver()) {
			timeRate = (timer.Now() - startTime).count() * subTimeCountDiv;
			temp = startTemp_ * std::pow(endTemp_ / startTemp_, timeRate);		
			divTemp = 1.0 / temp;
			for (auto& checker : checkers) {
				checker.temp = temp;
				checker.divTemp = divTemp;
			}


			for (int rp = 0; rp < stepLoopCount; ++rp) {
				int ti = (int)randTable(rand_);

				auto exec = [&](auto&& e, size_t i) {
					for (int pi : pis) {
						auto& checker = checkers[pi];
						e.func(checker, states[pi]);

						if (states[pi].RawScore() > maxStates[pi].RawScore()) {
							maxStates[pi] = states[pi];
						}
						else {
							if (timeRate >= rollbackStartRate_) {
								if (checker.noMaxUpdateCount >= checker.nextRollbackCheckCount) {
									states[pi] = maxStates[pi];
									checker.noMaxUpdateCount = 0;
									checker.nextRollbackCheckCount = (int)round(checker.nextRollbackCheckCount * rollbackNextMulti_);
								}
							}
						}
					}
				};

				TupleAccess(transitions, ti, exec);
			}
			if (forceEnd) {
				break;
			}

			{
				constexpr double start = 0.2;
				constexpr double end = 1.0;
				int targetPointCount = (int)states.size();
				if (timeRate >= end) {
					targetPointCount = 1;
				}
				else if (timeRate >= start) {
					double r = 1.0 - (timeRate - start) / (end - start);		
					targetPointCount = 1 + (int)floor(states.size() * r);
				}
				if ((int)pis.size() > targetPointCount) {
					sort(ALL(pis), [&](int a, int b) {
						return checkers[a].maxScore > checkers[b].maxScore;
					});
					pis.resize(targetPointCount);
				}
			}
		}

	}

	void ForceUpdate() {
	}

private:
	template <class Tuple, class Func>
	void TupleLoop(Tuple & t, Func && f) {
		TupleLoop2(t, forward<Func>(f), make_index_sequence<tuple_size<Tuple>::value>{});
	}
	template <class Tuple, class Func, size_t... Indics>
	void TupleLoop2(Tuple & t, Func && f, index_sequence<Indics...>) {
		using swallow = int[];
		(void)swallow {
			(TupleLoop3<Tuple, Func, Indics>(t, f), 0)...
		};
	}
	template <class Tuple, class Func, size_t Index>
	void TupleLoop3(Tuple & t, Func & f) {
		f(get<Index>(t), Index);
	}

	template <class Tuple, class Func>
	void TupleAccess(Tuple & t, int i, Func && f) {
		TupleAccess2(t, i, forward<Func>(f), make_index_sequence<tuple_size<Tuple>::value>{});
	}
	template <class Tuple, class Func, size_t... Indics>
	void TupleAccess2(Tuple & t, int i, Func && f, index_sequence<Indics...>) {
		using swallow = int[];
		(void)swallow {
			(TupleAccess3<Tuple, Func, Indics>(t, i, f), 0)...
		};
	}
	template <class Tuple, class Func, size_t Index>
	void TupleAccess3(Tuple & t, int i, Func & f) {
		if (i == Index) {
			f(get<Index>(t), Index);
		}
	}

};

using VPath = NArr<int>;		
using EPath = NArr<int>;		

struct PathFinder {
	NArr<NArr<int>> costs_;
	NArr<NArr<int>> froms_;
	NArr<NArr<int>> edges_;		

	void Init() {
		struct Node {
			int point;
			int totalCost;
			bool operator < (Node const& n) const {
				return totalCost > n.totalCost;
			}
		};

		PriorityQueue<Node> que;

		costs_.resize(N);
		froms_.resize(N);
		edges_.resize(N);

		REP(start, N) {
			auto& costs = costs_[start];
			auto& froms = froms_[start];
			auto& edges = edges_[start];
			costs.assign(N, INT_MAX);
			froms.assign(N, -1);
			edges.assign(N, -1);

			costs[start] = 0;
			que.push(Node{ start, 0 });

			while (!que.empty()) {
				Node node = que.top();
				que.pop();

				if (costs[node.point] != node.totalCost) {
					continue;
				}

				for (cauto& a : server.arounds_[node.point]) {
					int next = a.vi;
					int cost = a.w;

					VASSERT((s64)node.totalCost + (s64)cost <= INT_MAX);
					int nextCost = node.totalCost + cost;
					if (nextCost < costs[next]) {
						costs[next] = nextCost;
						froms[next] = node.point;
						edges[next] = a.ei;
						que.push(Node{ next, nextCost });
					}
				}
			}
		}
	}

	int Cost(int a, int b) const {
		return costs_[a][b];
	}



	void GetPath(int from, int to, VPath& vis, VPath& eis) const {
		vis.clear();
		eis.clear();

		cauto& froms = froms_[to];
		cauto& edges = edges_[to];
		int cur = from;
		while (cur >= 0) {
			vis.push(cur);
			if (edges[cur] >= 0) {
				eis.push(edges[cur]);
			}
			cur = froms[cur];
		}
	}
	void GetPathNoFrom(int from, int to, VPath& vis, VPath& eis) const {
		vis.clear();
		eis.clear();

		cauto& froms = froms_[to];
		cauto& edges = edges_[to];
		int cur = from;
		while (cur >= 0) {
			if (cur != from) {
				vis.push(cur);
			}
			if (edges[cur] >= 0) {
				eis.push(edges[cur]);
			}
			cur = froms[cur];
		}
	}
};
#define USE_COMPETITORS 1



struct Level2 {
	pint range;	
	int power;
	int powerCost;
};
struct Caster {
	KArr<int> kis;				
	KArr<int> ki2li;			
	KArr<Level2> levels;		
};

struct Area {
	KArr<NArr<pint>> nearLevels_;		
	NArr<Caster> casters_;

	NArr<NArr<int>> competitors_;		


	void Init() {
		array<int, 5001> powerCosts;
		REP(power, 5001) {
			powerCosts[power] = square(power);
		}

		KArr<int> dists;
		nearLevels_.resize(K);
		casters_.resize(N);

		REP(start, N) {
			auto& caster = casters_[start];
			caster.kis.clear();
			caster.levels.clear();

			dists.assign(K, -1);
			cauto& src = server.srcs_[start];
			REP(k, K) {
				cauto& dst = server.dsts_[k];
				s64 d2 = square(src.x - dst.x) + square(src.y - dst.y);
				if (d2 > square(5000)) {
					continue;
				}
				s64 lower = (s64)sqrt(d2);
				FOR(d, lower, 5001) {
					if (d2 <= square(d)) {
						caster.kis.push(k);
						dists[k] = d;
						break;
					}
				}
			}

			sort(ALL(caster.kis), [&](int a, int b) {
				return dists[a] < dists[b];
				});



			caster.ki2li.assign(K, -1);

			int prevD = -1;
			{
				auto& level = caster.levels.push();
				level.power = 0;
				level.powerCost = 0;
				level.range.a = 0;
				level.range.b = 0;
			}
			REP(i, caster.kis.size()) {
				int ki = caster.kis[i];
				int d = dists[ki];
				if (d != prevD) {
					if (prevD >= 0) {
						caster.levels.back().range.b = i;
					}
					auto& level = caster.levels.push();
					level.power = d;
					level.powerCost = square(d);
					level.range.a = i;
					level.range.b = -1;
				}
				nearLevels_[ki].push(pint{ caster.levels.size() - 1, start });
				caster.ki2li[ki] = caster.levels.size() - 1;
				prevD = d;
			}
			if (prevD >= 0) {
				caster.levels.back().range.b = caster.kis.size();
			}
		}

		REP(ki, K) {
			sort(ALL(nearLevels_[ki]), [&](const pint& a, const pint& b) {
				return casters_[a.second].levels[a.first].power < casters_[b.second].levels[b.first].power;
				});
		}


		array<bitset<N>, N> competitorsBit;
		REP(ki, K) {
			cauto& n = nearLevels_[ki];
			REP(i, n.size()) {
				int vi = n[i].second;
				FOR(j, i + 1, n.size()) {
					int vj = n[j].second;
					competitorsBit[vi].set(vj);
					competitorsBit[vj].set(vi);
				}
			}
		}
		competitors_.resize(N);
		REP(vi, N) {
			cauto& bit = competitorsBit[vi];
			REP(vj, N) {
				if (bit.test(vj)) {
					competitors_[vi].push(vj);
				}
			}
		}
	}


	int GetPower(int vi, int li) const {
		return casters_[vi].levels[li].power;
	}
	int GetPowerCost(int vi, int li) const {
		return casters_[vi].levels[li].powerCost;
	}

};


struct ChoiceTable {
private:
	vector<int> table_;
	int m_ = 0;

public:
	void Init(int n) {
		table_.resize(n);
		iota(table_.begin(), table_.end(), 0);
		m_ = 0;
	}
	void Choice(int m, Xor64& rand) {
		VASSERT(m <= (int)table_.size());
		REP(i, m) {
			swap(table_[i], table_[i + rand(table_.size() - i)]);
		}
		m_ = m;
	}

	using iterator = vector<int>::iterator;
	using const_iterator = vector<int>::const_iterator;

	inline const_iterator begin() const {
		return table_.begin();
	}
	inline const_iterator end() const {
		return table_.begin() + m_;
	}
	inline iterator begin() {
		return table_.begin();
	}
	inline iterator end() {
		return table_.begin() + m_;
	}
	inline int operator [](int index) const {
		return table_[index];
	}
};

template <int CAP>
struct ChoiceTableS {
private:
	CapArr<int, CAP> table_;
	int m_ = 0;

public:
	void Init(int n) {
		table_.resize(n);
		iota(table_.begin(), table_.end(), 0);
		m_ = 0;
	}
	void Choice(int m, Xor64& rand) {
		VASSERT(m <= (int)table_.size());
		REP(i, m) {
			swap(table_[i], table_[i + (int)rand(table_.size() - i)]);
		}
		m_ = m;
	}

	using iterator = vector<int>::iterator;
	using const_iterator = vector<int>::const_iterator;

	inline const_iterator begin() const {
		return table_.begin();
	}
	inline const_iterator end() const {
		return table_.begin() + m_;
	}
	inline iterator begin() {
		return table_.begin();
	}
	inline iterator end() {
		return table_.begin() + m_;
	}
	inline int operator [](int index) const {
		return table_[index];
	}
};

struct State {
	NArr<int> levels;	
	KArr<s8> refCounts;		
	s64 edgeCost;
	s64 powerCost;
	NArr<s8> useVis;		

	s64 rawScore;
	double score;

	void Init() {
		levels.assign(N, 0);
		refCounts.assign(K, 0);
		edgeCost = 0;
		powerCost = 0;
		useVis.clear();

		rawScore = 0;
		score = 0;
	}
	double EvalScore() const {
		return score;
	}
	double RawScore() const {
		return score;
	}
};
static constexpr int StateSize = sizeof(State);

struct Solver {
	Xor64 rand_;

	PathFinder pathFinder_;
	Area area_;

	NArr<int> easyCost_;		

	State bestState_;

	void CheckMaxScore(const State& state) {
		if (state.rawScore > bestState_.rawScore) {
			bestState_ = state;
		}
	}

	void Run(ChronoTimer& timer) {
		{

			pathFinder_.Init();
			area_.Init();
			MakeAllTree(easyCost_);


		}

		vector<State> states;
		InitState(states);
		bestState_ = states.front();

		timer.Start(TIME_LIMIT);

		SimulatedAnnealing sa;
		{
			sa.startTemp_ = HP.StartTemp;
			sa.endTemp_ = HP.EndTemp;
			sa.stepLoopCount = 100;

			sa.rollbackStartRate_ = 0.0;
			sa.rollbackCount_ = HP.RollbackCount;
			sa.rollbackNextMulti_ = HP.RollbackNextMulti;

			auto transitions = make_tuple(
				MAKE_TRANS(Transition_DownPower, HP.TransitionDown)
				, MAKE_TRANS(Transition_UpPower, HP.TransitionUp)
			);
			sa.Run2(timer, states, transitions);
		}

		s64 bestEdgeCost = INT64_MAX;
		MArr<bool> bestOns;
		MArr<bool> ons;
		REP(vi, N) {
			if (bestState_.levels[vi] > 0) {
				s64 edgeCost = MakeTreeWithStart(bestState_.levels, ons, vi);
				if (edgeCost < bestEdgeCost) {
					bestEdgeCost = edgeCost;
					bestOns = ons;
				}
			}
		}

		Result result;
		result.powers.resize(N);
		REP(i, N) {
			result.powers[i] = area_.GetPower(i, bestState_.levels[i]);
		}
		result.ons = bestOns;

		server.Output(result);
	}

	void InitState(vector<State>& states) {
		int StateCount = HP.InitialGreedyLevelCount + HP.InitialMaxLevelCount + HP.InitialRandLevelCount;
		states.resize(StateCount);

		int si = 0;
		REP(i, HP.InitialGreedyLevelCount) {
			State& state = states[si];
			state.Init();
			InitLevelGreedy(state.levels);
			InitOther(state);
			++si;
		}

		if (HP.InitialMaxLevelCount > 0) {
			State& baseState = states[si];
			{
				baseState.Init();
				InitLevelMax(baseState.levels);
				InitOther(baseState);
				++si;
			}

			REP(i, HP.InitialMaxLevelCount - 1) {
				states[si] = baseState;
				++si;
			}
		}

		REP(i, HP.InitialRandLevelCount) {
			State& state = states[si];
			state.Init();
			InitLevelRand(state.levels);
			InitOther(state);
			++si;
		}

		VASSERT(si == StateCount);
	}

	s64 MakeTree(const NArr<int>& levels) {
		struct Node {
			int to;
			int cost;

			bool operator < (const Node& r) const {
				return cost > r.cost;
			}
		};

		static NArr<s8> terminals;
		static NArr<int> costs;		
		static NArr<int> froms;
		costs.assign(N, INT_MAX);
		froms.assign(N, -1);

		s64 totalCost = 0;
		int nextVi = 0;
		terminals.clear();
		terminals.push(0);
		FOR(ni, 1, N) {
			if (levels[ni] > 0) {
				terminals.push(ni);
			}
		}
		int nextTi = (int)rand_(terminals.size());

		while (nextTi >= 0) {
			static VPath vis;
			const int nextVi = terminals[nextTi];
			int from = froms[nextVi];
			if (from < 0) {
				vis.clear();
				vis.push(nextVi);
				terminals.RemoveSwap(nextTi);
			}
			else {
				vis.clear();
				cauto& froms = pathFinder_.froms_[from];
				cauto& edges = pathFinder_.edges_[from];
				int cur = nextVi;
				while (cur != from) {
					vis.push(cur);

					int ei = edges[cur];
					VASSERT(ei >= 0);
					totalCost += server.edges_[ei].w;

					cur = froms[cur];
				}
				terminals.RemoveSwap(nextTi);
			}

			int bestCost = INT_MAX;
			nextTi = -1;
			REP(ti, terminals.size()) {
				int ni = terminals[ti];
				int minCost = costs[ni];
				int minCi = -1;
				for (int ci : vis) {
					VASSERT(levels[ni] >= 0);

					int cost = pathFinder_.Cost(ci, ni);
					if (cost < minCost) {
						minCost = cost;
						minCi = ci;
					}
				}
				if (minCi >= 0) {
					costs[ni] = minCost;
					froms[ni] = minCi;
				}
				if (minCost < bestCost) {
					bestCost = minCost;
					nextTi = ti;
				}
			}
			VASSERT(terminals.empty() || nextTi >= 0);
		}

		return totalCost;
	}

	s64 MakeTreeWithStart(const NArr<int>& levels, MArr<bool>& ons, int start) {
		ons.assign(M, false);

		struct Node {
			int to;
			int cost;

			bool operator < (const Node& r) const {
				return cost > r.cost;
			}
		};

		static CapSet<s8, N> lefts;
		static NArr<int> costs;		
		static NArr<int> froms;
		costs.assign(N, INT_MAX);
		froms.assign(N, -1);

		s64 totalCost = 0;
		int nextVi = 0;
		lefts.Clear();
		lefts.Add(0);
		FOR(ni, 1, N) {
			if (levels[ni] > 0) {
				lefts.Add(ni);
			}
		}
		nextVi = start;

		while (!lefts.empty()) {
			static VPath vis;
			int from = froms[nextVi];
			if (from < 0) {
				vis.clear();
				vis.push(nextVi);
				lefts.Remove(nextVi);
			}
			else {
				vis.clear();
				cauto& froms = pathFinder_.froms_[from];
				cauto& edges = pathFinder_.edges_[from];
				int cur = nextVi;
				while (cur != from) {
					vis.push(cur);
					lefts.ForceRemove(cur);

					int ei = edges[cur];
					VASSERT(ei >= 0);
					ons[ei] = true;
					totalCost += server.edges_[ei].w;

					cur = froms[cur];
				}
			}

			int bestCost = INT_MAX;
			nextVi = -1;
			for (int ni : lefts) {
				int minCost = costs[ni];
				int minCi = -1;
				for (int ci : vis) {
					VASSERT(levels[ni] >= 0);

					int cost = pathFinder_.Cost(ci, ni);
					if (cost < minCost) {
						minCost = cost;
						minCi = ci;
					}
				}
				if (minCi >= 0) {
					costs[ni] = minCost;
					froms[ni] = minCi;
				}
				if (minCost < bestCost) {
					bestCost = minCost;
					nextVi = ni;
				}
			}
			VASSERT(lefts.empty() || nextVi >= 0);
		}

		return totalCost;
	}


	void MakeAllTree(NArr<int>& fromCosts) {
		struct Node {
			int to;
			int cost;

			bool operator < (const Node& r) const {
				return cost > r.cost;
			}
		};

		static PriorityQueue<Node> que;
		que.Clear();
		que.push({ 0, 0 });

		NArr<int>& costs = fromCosts;		
		static NArr<int> froms;
		static NArr<bool> used;
		costs.assign(N, INT_MAX);
		froms.assign(N, -1);
		used.assign(N, false);

		while (!que.empty()) {
			Node node = que.top();
			que.pop();

			if (used[node.to]) {
				continue;
			}
			used[node.to] = true;

			fromCosts[node.to] = node.cost;

			int ci = node.to;
			for (cauto& ar : server.arounds_[node.to]) {
				int ni = ar.vi;
				if (!used[ni]) {
					int cost = pathFinder_.Cost(ci, ni);
					if (cost < costs[ni]) {
						costs[ni] = cost;
						que.push(Node{ ni, cost });
					}
				}
			}
		}
	}

	s64 EasyEdgeCost(const NArr<int>& levels) {
		s64 cost = 0;
		REP(i, N) {
			if (levels[i] > 0) {
				cost += easyCost_[i];
			}
		}
		return cost;
	}

	void InitLevelGreedy(NArr<int>& levels) {
		static KArr<int> kis;
		if (kis.empty()) {
			kis.resize(K);
			iota(ALL(kis), 0);
		}
		shuffle(ALL(kis), rand_);

		for (int ki : kis) {
			cauto& srcs = area_.nearLevels_[ki];
			int bestVi = -1;
			int bestLevel = 0;
			int bestAppendCost = INT_MAX;
			for (auto [level, vi] : srcs) {
				int nowPowerCost = area_.GetPowerCost(vi, levels[vi]);
				int newPowerCost = area_.GetPowerCost(vi, level);
				int appendCost = newPowerCost - nowPowerCost;
				if (appendCost < bestAppendCost) {
					bestAppendCost = appendCost;
					bestLevel = level;
					bestVi = vi;
				}
			}
			VASSERT(bestVi >= 0);
			chmax(levels[bestVi], bestLevel);
		}
	}

	void InitLevelMax(NArr<int>& levels) {
		REP(vi, N) {
			cauto& caster = area_.casters_[vi];
			levels[vi] = caster.levels.size() - 1;
		}
	}

	void InitLevelRand(NArr<int>& levels) {
		static KArr<bool> used;
		used.assign(K, false);
		REP(vi, N) {
			cauto& caster = area_.casters_[vi];
			levels[vi] = (int)rand_(caster.levels.size());
			REP(i, caster.levels[levels[vi]].range.b) {
				int ki = caster.kis[i];
				used[ki] = true;
			}
		}


		static KArr<int> kis;
		if (kis.empty()) {
			kis.resize(K);
			iota(ALL(kis), 0);
		}
		shuffle(ALL(kis), rand_);

		for (int ki : kis) {
			if (used[ki]) {
				continue;
			}
			cauto& srcs = area_.nearLevels_[ki];
			int bestVi = -1;
			int bestLevel = 0;
			int bestAppendCost = INT_MAX;
			for (auto [level, vi] : srcs) {
				int nowPowerCost = area_.GetPowerCost(vi, levels[vi]);
				int newPowerCost = area_.GetPowerCost(vi, level);
				int appendCost = newPowerCost - nowPowerCost;
				if (appendCost < bestAppendCost) {
					bestAppendCost = appendCost;
					bestLevel = level;
					bestVi = vi;
				}
			}
			VASSERT(bestVi >= 0);
			chmax(levels[bestVi], bestLevel);
		}
	}

	void InitOther(State& state) {
		state.powerCost = 0;
		state.refCounts.assign(K, 0);
		state.useVis.clear();
		REP(vi, N) {
			state.powerCost += area_.GetPowerCost(vi, state.levels[vi]);

			if (state.levels[vi] > 0) {
				state.useVis.push(vi);
			}

			cauto& caster = area_.casters_[vi];
			REP(i, caster.levels[state.levels[vi]].range.b) {
				int ki = caster.kis[i];
				++state.refCounts[ki];
				VASSERT(state.refCounts[ki] >= 0);
			}
		}
		state.edgeCost = MakeTree(state.levels);
		state.score = CalcScore(state);
		state.rawScore = (s64)round(state.score);
	}

	double CalcScore(s64 totalCost) {
		return 1000000 * (1 + 100000000 / (double)(totalCost + 10000000));
	}
	double CalcScore(const State& state) {
		return CalcScore(state.powerCost + state.edgeCost);
	}

	double ScoreToCost(double score) {
		return 1e14 / (score - 1e6) - 1e7;
	}

	bool TryDown(State& state, int vi) {
		cauto& caster = area_.casters_[vi];
		auto& level = state.levels[vi];

		int startI = 0;
		int newLevel = level;
		{
			int i = caster.levels[level].range.b - 1;
			while (i >= 0) {
				int ki = caster.kis[i];
				VASSERT(state.refCounts[ki] > 0);
				if (state.refCounts[ki] == 1) {
					break;
				}
				--i;
			}
			if (i < 0) {
				startI = 0;
				newLevel = 0;
			}
			else {
				int ki = caster.kis[i];
				int li = caster.ki2li[ki];
				startI = caster.levels[li].range.b;
				newLevel = li;
			}
		}
		if (level == newLevel) {
			return false;
		}
		int endI = caster.levels[level].range.b;
		FOR(i, startI, endI) {
			int ki = caster.kis[i];
			--state.refCounts[ki];
		}
		level = newLevel;
		return true;
	}

	void DownPossible(State& state) {
		if (rand_(2) == 0) {
			REP(vi, N) {
				TryDown(state, vi);
			}
		}
		else {
			RREP(vi, N) {
				TryDown(state, vi);
			}
		}
	}

	void Transition_DownPower(SimulatedAnnealing& sa, SAChecker& checker, State& state) {
		int visI = (int)rand_(state.useVis.size());
		int svi = state.useVis[visI];

		auto oldLevels = state.levels;
		auto oldRefCounts = state.refCounts;

		int newLevel = (int)rand_(state.levels[svi]);

		VASSERT(newLevel != state.levels[svi]);
		state.levels[svi] = newLevel;

		bool levelUp = false;
		{
			cauto& caster = area_.casters_[svi];
			int startI = caster.levels[newLevel].range.b;
			int endI = caster.levels[oldLevels[svi]].range.b;
			FOR(i, startI, endI) {
				int ki = caster.kis[i];

				--state.refCounts[ki];

				if (state.refCounts[ki] == 0) {
					cauto& srcs = area_.nearLevels_[ki];

					int bestVi = -1;
					int bestLevel = 0;
					if (srcs.size() == 1) {
						bestLevel = srcs[0].first;
						bestVi = srcs[0].second;
					}
					else {
						int bestAppendCost = INT_MAX;
						for (auto&& [level, vi] : srcs) {
							if (vi == svi) {
								continue;
							}
							int nowPowerCost = area_.GetPowerCost(vi, state.levels[vi]);
							int newPowerCost = area_.GetPowerCost(vi, level);
							int appendCost = newPowerCost - nowPowerCost;
							if (nowPowerCost == 0) {
								appendCost += easyCost_[vi];
							}

							if (appendCost < bestAppendCost) {
								bestAppendCost = appendCost;
								bestLevel = level;
								bestVi = vi;

								if (bestAppendCost <= 0) {
									break;
								}
							}
						}
					}

					VASSERT(bestVi >= 0);
					if (bestLevel > state.levels[bestVi]) {
						cauto& caster = area_.casters_[bestVi];
						int startI = caster.levels[state.levels[bestVi]].range.b;
						int endI = caster.levels[bestLevel].range.b;
						FOR(i, startI, endI) {
							int ki = caster.kis[i];
							++state.refCounts[ki];
						}
						state.levels[bestVi] = bestLevel;
						levelUp = true;

					}
				}
			}
		}

		if (levelUp) {
			DownPossible(state);
		}

		bool requireUpdateEdge = false;
		REP(i, N) {
			if ((oldLevels[i] == 0 && state.levels[i] != 0) ||
				(oldLevels[i] != 0 && state.levels[i] == 0)) {
				requireUpdateEdge = true;
				break;
			}
		}

		s64 newPowerCost = 0;
		REP(i, N) {
			newPowerCost += area_.GetPowerCost(i, state.levels[i]);
		}

		s64 newEdgeCost = 0;
		if (requireUpdateEdge) {
			newEdgeCost = MakeTree(state.levels);
		}
		else {
			newEdgeCost = state.edgeCost;
		}

		double newScore = CalcScore(newPowerCost + newEdgeCost);

		if (checker(newScore)) {
			if (requireUpdateEdge) {
				state.edgeCost = newEdgeCost;
				state.useVis.clear();
				REP(i, N) {
					if (state.levels[i]) {
						state.useVis.push(i);
					}
				}
			}
			state.powerCost = newPowerCost;
			state.score = newScore;
			state.rawScore = (s64)round(state.score);

			CheckMaxScore(state);
		}
		else {
			state.levels = oldLevels;
			state.refCounts = oldRefCounts;
		}
	}

	void Transition_DownPowerMulti(SimulatedAnnealing& sa, SAChecker& checker, State& state) {
		CapArr<int, 8> targetVis;
		if (state.useVis.size() == 1) {
			targetVis.push(state.useVis[0]);
		}
		else {
			static ChoiceTable tbl;
			tbl.Init(state.useVis.size());
			int targetCount = 2 + (int)rand_(min(4, state.useVis.size()));
			tbl.Choice(targetCount, rand_);
			for (int i : tbl) {
				targetVis.push(state.useVis[i]);
			}
		}

		auto oldLevels = state.levels;
		auto oldRefCounts = state.refCounts;
		bool levelUp = false;

		REP(ti, targetVis.size()) {
			int svi = targetVis[ti];

			int newLevel = (int)rand_(state.levels[svi]);
			chmax(newLevel, 0);
			state.levels[svi] = newLevel;

			cauto& caster = area_.casters_[svi];
			int startI = caster.levels[newLevel].range.b;
			int endI = caster.levels[oldLevels[svi]].range.b;
			FOR(i, startI, endI) {
				int ki = caster.kis[i];

				--state.refCounts[ki];

				if (state.refCounts[ki] == 0) {
					cauto& srcs = area_.nearLevels_[ki];

					int bestVi = -1;
					int bestLevel = 0;

					int bestAppendCost = INT_MAX;
					for (auto&& [level, vi] : srcs) {
						if (targetVis.find(vi) >= 0) {
							continue;
						}
						int nowPowerCost = area_.GetPowerCost(vi, state.levels[vi]);
						int newPowerCost = area_.GetPowerCost(vi, level);
						int appendCost = newPowerCost - nowPowerCost;
						if (nowPowerCost == 0 && newPowerCost > 0) {
							appendCost += (int)easyCost_[vi];
						}
						if (appendCost < bestAppendCost) {
							bestAppendCost = appendCost;
							bestLevel = level;
							bestVi = vi;

							if (bestAppendCost <= 0) {
								break;
							}
						}
					}
					if (bestVi < 0) {
						for (auto&& [level, vi] : srcs) {
							int nowPowerCost = area_.GetPowerCost(vi, state.levels[vi]);
							int newPowerCost = area_.GetPowerCost(vi, level);
							int appendCost = newPowerCost - nowPowerCost;
							if (nowPowerCost == 0 && newPowerCost > 0) {
								appendCost += (int)easyCost_[vi];
							}
							if (appendCost < bestAppendCost) {
								bestAppendCost = appendCost;
								bestLevel = level;
								bestVi = vi;

								if (bestAppendCost <= 0) {
									break;
								}
							}
						}
					}

					VASSERT(bestVi >= 0);
					if (bestLevel > state.levels[bestVi]) {
						cauto& caster = area_.casters_[bestVi];
						int startI = caster.levels[state.levels[bestVi]].range.b;
						int endI = caster.levels[bestLevel].range.b;
						FOR(i, startI, endI) {
							int ki = caster.kis[i];
							++state.refCounts[ki];
						}
						state.levels[bestVi] = bestLevel;
						levelUp = true;
					}
				}
			}
		}

		if (levelUp) {
			DownPossible(state);
		}

		bool requireUpdateEdge = false;
		REP(i, N) {
			if ((oldLevels[i] == 0 && state.levels[i] != 0) ||
				(oldLevels[i] != 0 && state.levels[i] == 0)) {
				requireUpdateEdge = true;
				break;
			}
		}

		s64 newEdgeCost = 0;
		if (requireUpdateEdge) {
			newEdgeCost = MakeTree(state.levels);
		}
		else {
			newEdgeCost = state.edgeCost;
		}

		s64 newPowerCost = 0;
		REP(i, N) {
			newPowerCost += area_.GetPowerCost(i, state.levels[i]);
		}
		double newScore = CalcScore(newPowerCost + newEdgeCost);

		if (checker(newScore)) {
			if (requireUpdateEdge) {
				state.edgeCost = newEdgeCost;
				state.useVis.clear();
				REP(i, N) {
					if (state.levels[i]) {
						state.useVis.push(i);
					}
				}
			}
			state.powerCost = newPowerCost;
			state.score = newScore;
			state.rawScore = (s64)round(state.score);

			CheckMaxScore(state);
		}
		else {
			state.levels = oldLevels;
			state.refCounts = oldRefCounts;
		}
	}

	void Transition_UpPower(SimulatedAnnealing& sa, SAChecker& checker, State& state) {
		int visI = (int)rand_(state.useVis.size());
		int svi = state.useVis[visI];

		auto oldLevels = state.levels;
		auto oldRefCounts = state.refCounts;

		cauto& caster = area_.casters_[svi];
		int maxUp = caster.levels.size() - 1 - state.levels[svi];
		if (maxUp == 0) {
			return;
		}
		int newLevel = state.levels[svi] + 1 + (int)rand_(maxUp);
		VASSERT(newLevel != state.levels[svi]);

		int startI = caster.levels[state.levels[svi]].range.b;
		int endI = caster.levels[newLevel].range.b;
		FOR(i, startI, endI) {
			int ki = caster.kis[i];
			++state.refCounts[ki];
		}
		state.levels[svi] = newLevel;

		bool update = false;
		if (rand_(2) == 0) {
			for (int vi : area_.competitors_[svi]) {
				VASSERT(vi != svi);
				if (TryDown(state, vi)) {
					update = true;
				}
			}
		}
		else {
			cauto& competitors = area_.competitors_[svi];
			RREP(i, competitors.size()) {
				int vi = competitors[i];
				VASSERT(vi != svi);
				if (TryDown(state, vi)) {
					update = true;
				}
			}
		}
		if (!update) {
			return;
		}
		{
			TryDown(state, svi);
		}

		bool requireUpdateEdge = false;
		REP(i, N) {
			if ((oldLevels[i] == 0 && state.levels[i] != 0) ||
				(oldLevels[i] != 0 && state.levels[i] == 0)) {
				requireUpdateEdge = true;
				break;
			}
		}

		s64 newEdgeCost = 0;
		if (requireUpdateEdge) {
			newEdgeCost = MakeTree(state.levels);
		}
		else {
			newEdgeCost = state.edgeCost;
		}

		s64 newPowerCost = 0;
		REP(i, N) {
			newPowerCost += area_.GetPowerCost(i, state.levels[i]);
		}
		double newScore = CalcScore(newPowerCost + newEdgeCost);

		if (checker(newScore)) {
			if (requireUpdateEdge) {
				state.edgeCost = newEdgeCost;
				state.useVis.clear();
				REP(i, N) {
					if (state.levels[i]) {
						state.useVis.push(i);
					}
				}
			}
			state.powerCost = newPowerCost;
			state.score = newScore;
			state.rawScore = (s64)round(state.score);

			CheckMaxScore(state);
		}
		else {
			state.levels = oldLevels;
			state.refCounts = oldRefCounts;
		}
	}

	void Transition_UpPowerMulti(SimulatedAnnealing& sa, SAChecker& checker, State& state) {
		CapArr<int, 8> targetVis;
		if (state.useVis.size() == 1) {
			targetVis.push(state.useVis[0]);
		}
		else {
			static ChoiceTable tbl;
			tbl.Init(state.useVis.size());
			int targetCount = 2 + (int)rand_(min(4, state.useVis.size()));
			tbl.Choice(targetCount, rand_);
			for (int i : tbl) {
				targetVis.push(state.useVis[i]);
			}
		}

		auto oldLevels = state.levels;
		auto oldRefCounts = state.refCounts;
		for (int svi : targetVis) {
			cauto& caster = area_.casters_[svi];
			int maxUp = caster.levels.size() - 1 - state.levels[svi];
			if (maxUp == 0) {
				return;
			}
			int newLevel = state.levels[svi] + 1 + (int)rand_(maxUp);
			VASSERT(newLevel != state.levels[svi]);

			int startI = caster.levels[state.levels[svi]].range.b;
			int endI = caster.levels[newLevel].range.b;
			FOR(i, startI, endI) {
				int ki = caster.kis[i];
				++state.refCounts[ki];
			}
			state.levels[svi] = newLevel;
		}

		bool update = false;
		if (rand_(2) == 0) {
			REP(vi, N) {
				if (targetVis.find(vi) >= 0) {
					continue;
				}
				if (TryDown(state, vi)) {
					update = true;
				}
			}
		}
		else {
			RREP(vi, N) {
				if (targetVis.find(vi) >= 0) {
					continue;
				}
				if (TryDown(state, vi)) {
					update = true;
				}
			}
		}
		if (!update) {
			return;
		}
		for (int vi : targetVis) {
			TryDown(state, vi);
		}

		bool requireUpdateEdge = false;
		REP(i, N) {
			if ((oldLevels[i] == 0 && state.levels[i] != 0) ||
				(oldLevels[i] != 0 && state.levels[i] == 0)) {
				requireUpdateEdge = true;
				break;
			}
		}

		s64 newEdgeCost = 0;
		if (requireUpdateEdge) {
			newEdgeCost = MakeTree(state.levels);
		}
		else {
			newEdgeCost = state.edgeCost;
		}

		s64 newPowerCost = 0;
		REP(i, N) {
			newPowerCost += area_.GetPowerCost(i, state.levels[i]);
		}
		double newScore = CalcScore(newPowerCost + newEdgeCost);

		if (checker(newScore)) {
			if (requireUpdateEdge) {
				state.edgeCost = newEdgeCost;
				state.useVis.clear();
				REP(i, N) {
					if (state.levels[i]) {
						state.useVis.push(i);
					}
				}
			}
			state.powerCost = newPowerCost;
			state.score = newScore;
			state.rawScore = (s64)round(state.score);

			CheckMaxScore(state);
		}
		else {
			state.levels = oldLevels;
			state.refCounts = oldRefCounts;
		}
	}

};















struct Main {
    void Run(int argc, const char* argv[]) {
        ChronoTimer timer;
        server.InitInput(timer);

        static Solver solver;
        timer.StartMs(TIME_LIMIT);

        solver.Run(timer);

        server.Finalize();
    }
};
