//状態
struct STATE {
};

// 状態の初期化
void init (STATE& state) {
}

// 状態遷移
void modify (STATE& state) {
}

// 状態のスコア計算
int calc_score (STATE& state) {
}

// 焼きなまし法
void sa() {
    STATE state;
    init(state);
    double start_temp, end_temp; // 適当な値を入れる（後述）
    double start_time; // 開始時刻
    while (true) { // 時間の許す限り回す
        double now_time; // 現在時刻
        if (now_time - start_time > TIME_LIMIT) break;

        STATE new_state = state;
        modify(new_state);
        int new_score = calc_score(new_state);
        int pre_score = calc_score(state);

        // 温度関数
        double temp = start_temp + (end_temp - start_temp) * (now_time-start_time) / TIME_LIMIT;
        // 遷移確率関数(最大化の場合)
        double prob = exp((new_score-pre_score)/temp);

        if (prob > (rand()%INF)/(double)INF) { // 確率probで遷移する
            state = new_state;
        }
    }
}

// https://gasin.hatenadiary.jp/entry/2019/09/03/162613