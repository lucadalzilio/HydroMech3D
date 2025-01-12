#include <iostream>
#include <tbb/tbb.h>

int main() {
    try {
        tbb::task_arena arena;
        arena.execute([] {
            std::cout << "TBB is working correctly!" << std::endl;
        });
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}
