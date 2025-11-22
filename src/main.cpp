
#include <iostream>
#include <gtest/gtest.h>
#include "hw1.h"

int main(int argc, char **argv)
{
    if (false) // make false to run unit-tests
    {
        size_t n, m;
        std::cin >> n >> m;
        
        algebra::Matrix m1=algebra::random(n ,m,0,2);
        algebra::show(m1);
        algebra::Matrix m2 = algebra::random(n, m, 0, 2);
        algebra::show(m2);
        algebra::Matrix result = algebra::concatenate(m1, m2);
        algebra::show(result);
        // debug section
    }
    else
    {
        ::testing::InitGoogleTest(&argc, argv);
        std::cout << "RUNNING TESTS ..." << std::endl;
        int ret{RUN_ALL_TESTS()};
        if (!ret)
            std::cout << "<<<SUCCESS>>>" << std::endl;
        else
            std::cout << "FAILED" << std::endl;
    }
    return 0;
}