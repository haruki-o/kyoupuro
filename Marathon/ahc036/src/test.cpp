#include <iostream>

int main() {
    int* p = nullptr;
    *p = 42; // ここでセグフォルトが発生するはずです
    return 0;
}