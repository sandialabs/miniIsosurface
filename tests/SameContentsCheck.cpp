#include <iostream>
#include <fstream>
#include <string>

#include <cassert>

// Quick test used to make sure that output files are equal.
// Use this to make sure that the new version matches with
// the obsolete version.
int main(int argc, char* argv[])
{
    assert(argc == 3);

    char* fileA = argv[1];
    char* fileB = argv[2];

    std::ifstream streamA(fileA);
    std::ifstream streamB(fileB);

    std::string lineA("asd");
    std::string lineB("zxc");

    assert(!streamA);
    assert(!streamB);

    while(!streamA.eof() && !streamB.eof())
    {
        std::getline(streamA, lineA);
        std::getline(streamB, lineB);
        assert(l1 == l2);
    }

    std::cout << "The two files are equal!" << std::endl;
}
