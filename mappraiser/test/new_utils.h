//
// Created by sbiquard on 17/11/23.
//

#ifndef MIDAPACK_NEW_UTILS_H
#define MIDAPACK_NEW_UTILS_H

#include <mapmat/mapmat.h>

#include <iostream>
#include <string>
#include <vector>

void mat_print(Mat *A, const std::string &name = "");

template <typename T>
void vec_print(std::vector<T> &v, const std::string &name) {
    if (!name.empty())
        std::cout << name << " = ";
    std::cout << "[";
    ulong n = v.size();
    for (ulong i = 0; i < n; i++) {
        std::cout << v[i];
        if (i + 1 < n)
            std::cout << ", ";
    }
    std::cout << "]" << std::endl;
}

template <typename T>
__attribute__((unused)) void print_array(T *arr, int n,
                                         const std::string &name) {
    if (!name.empty())
        std::cout << name << " =" << std::endl;
    std::cout << "[";
    for (int i = 0; i < n; i++) {
        std::cout << arr[i];
        if (i + 1 < n)
            std::cout << ", ";
    }
    std::cout << "]" << std::endl;
}

#endif // MIDAPACK_NEW_UTILS_H
