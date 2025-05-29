#include <fstream>
#include <iostream>
#include <string>

int main() {
    std::string num = "5"; // Number of test case to read
    std::ifstream infile("test" + num + ".txt");
    if (!infile) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    std::string line;
    std::string seq[3];
    int i = 0;

    while (std::getline(infile, line)) {
        if (line.empty()) continue; // Skip empty lines
        if (line[0] == '>') {
            i++;
            continue;
        }
        seq[i - 1] += line;
    }

    infile.close();

    for (const auto& aln : seq) {
        std::cout << aln.substr(0, 10) + "..." + aln.substr(aln.size() - 10) << std::endl;
    }

    std::ofstream outfile("cpp_test" + num + ".txt");
    outfile << seq[0] << '\n' << seq[1] << '\n' << seq[2] << '\n';
    outfile.close();

    return 0;
}

