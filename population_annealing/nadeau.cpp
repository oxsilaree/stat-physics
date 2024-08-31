
#include "nadeau.h"
#include <vector>
#include <iostream>
#include <cassert>
void print_pages() {
    static size_t pagesize = sysconf(_SC_PAGESIZE);
    int64_t bytes = getCurrentRSS();
    assert((bytes % pagesize) == 0);
    size_t pages = bytes / pagesize;
    std::cout << "page size: " << pagesize << "\t";
    std::cout << "bytes: " << bytes << "\t";
    std::cout << "pages: " << pages << std::endl;
}