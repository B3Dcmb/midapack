//
// Created by sbiquard on 17/11/23.
//

#include "new_utils.h"

#include <iostream>

void mat_print(Mat *A, const std::string &name) {
    if (!name.empty())
        std::cout << name << " =" << std::endl;
    std::cout << "[";

    int nnz = A->nnz;
    // loop over rows
    for (int t = 0; t < A->m; t++) {
        std::cout << ((t == 0) ? "[" : " [");
        int c = 0;
        // loop over columns
        for (int j = 0; j < nnz * A->lcount; j++) {
            if (c < nnz && j == nnz * A->indices[t] + c) {
                std::cout << A->values[nnz * t + c];
                c++;
            } else {
                std::cout << 0;
            }
            if (j + 1 < nnz * A->lcount)
                std::cout << ", ";
        }
        std::cout << "]";
        if (t + 1 < A->m)
            std::cout << "," << std::endl;
    }
    std::cout << "]" << std::endl;
}
