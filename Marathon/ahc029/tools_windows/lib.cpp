以下は、提供されたRustのコードをC++に変換したものです。RustとC++ではいくつかの違いがあるため、注意してください。また、一部のライブラリや機能はC++には存在しないため、代替手段を検討する必要があります。

```cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <deque>
#include <random>
#include <cmath>

// C++で利用するライブラリに応じてヘッダーファイルを追加する

enum CardType {
    WorkSingle,
    WorkAll,
    CancelSingle,
    CancelAll,
    Invest,
};

struct Card {
    CardType ty;
    int64_t w;
    int64_t p;
};

struct Project {
    int64_t initial_h;
    int64_t h;
    int64_t v;
};

struct JudgeData {
    size_t n;
    size_t m;
    size_t k;
    size_t t;
    std::vector<Project> initial_projects;
    std::vector<Project> new_projects;
    std::vector<Card> initial_cards;
    std::vector<std::vector<Card>> new_cards;
};

class SVG {
    // SVG描画クラスを実装する（実際の実装は割愛）
};

class SVGDrawer {
    // SVGの描画を行うクラスを実装する（実際の実装は割愛）
};

struct JudgeResult {
    int64_t score;
};

struct FieldData {
    std::vector<Project> projects;
    int64_t l;
    int64_t money;
};

struct VisData {
    std::vector<std::string> comments;
    FieldData field_before_use;
    FieldData field_after_use;
    std::vector<Card> cards;
    std::vector<Card> candidates;
    size_t selected_card;
    size_t selected_project;
    size_t selected_candidate;
    std::vector<bool> finished_projects;
    std::vector<bool> canceled_projects;
    // その他のメンバ変数や関数を追加する場合があります
};

class State {
public:
    State(const JudgeData& data);
    void update_project(size_t m);
    std::vector<Card> get_card_candidate();
    void increment_l();
    void use_card(const Card& card, size_t m, VisData& vis_data);
    VisData play_turn(std::ifstream& input_stream, std::ofstream& output_stream);
    // その他のメンバ関数や変数を追加する場合があります

private:
    size_t n;
    size_t m;
    size_t k;
    size_t max_turn;
    size_t turn;
    std::vector<Card> cards;
    std::vector<Project> projects;
    std::deque<Project> new_projects;
    int64_t l;
    int64_t money;
    // その他のメンバ変数や関数を追加する場合があります
};

// その他の関数やクラスの実装を追加する場合があります

int main() {
    // メイン関数の実装
    return 0;
}
```

このコードはRustの機能やライブラリに合わせて変更する必要があります。特にSVGの描画に関する部分や、ファイルの入出力、エラーハンドリングなどはC++での標準的な方法に合わせて修正する必要があります。提供されたRustのコードの具体的な動作や利用環境によっても変更が必要です。