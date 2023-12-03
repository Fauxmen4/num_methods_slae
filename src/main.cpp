#include "matrix_t.h"

void PrintAnswer(Matrix_t x) {
    std::cout << '[';
    std::cout << x.matrix[0][0];
    for (int i = 1; i < x.lines; i++) {
        std::cout << ", " << x.matrix[i][0];
    }
    std::cout << "]\n";
}

void PrintDifference(Matrix_t accur, Matrix_t close) {
    std::cout << '[';
    std::cout << abs(accur.matrix[0][0]-close.matrix[0][0]);
    for (int i = 1; i < accur.lines; i++) {
    std::cout << ", " << abs(accur.matrix[i][0]-close.matrix[i][0]);
    }
    std::cout << "] ";
}

const std::string testDir = "/home/parz1val/Documents/amcp/num_methods/tests/";

// int main() {
//     Matrix_t A(3, 3);
//     Matrix_t b(3, 1);
//     A.Enter();
//     b.Enter();
//     Solution(A, b).Print();
// }

int main() {
    std::string testPath = testDir+"i.txt";
    for (int eps = -3; eps > -6; eps--) {
        std::cout << "Gauss Seidel method with accuracy: " << pow(10,eps) << "\n";
        for (int i = 0; i < 5; i++) {
            testPath[testDir.length()] = toascii(i+48);
            std::pair<Matrix_t, Matrix_t> test = EnterSLAE(3, testPath);
            std::cout << i << ": ";
            std::pair<Matrix_t, int> answer = SeidelMethod(test.first, test.second, pow(10,eps));
            std::cout << "k = " << answer.second << ' ';
            PrintDifference(Solution(test.first, test.second), answer.first);
            PrintAnswer(answer.first);
        }    
    }    
    for (int eps = -3; eps > -6; eps--) {
        std::cout << "Simple iteration method with accuracy: " << pow(10,eps) << "\n";
        for (int i = 0; i < 5; i++) {
            testPath[testDir.length()] = toascii(i+48);
            std::pair<Matrix_t, Matrix_t> test = EnterSLAE(3, testPath);
            std::cout << i << ": ";
            std::pair<Matrix_t, int> answer = SimpleIterMethod(test.first, test.second, pow(10,eps));
            std::cout << "k = " << answer.second << ' ';
            PrintDifference(Solution(test.first, test.second), answer.first);
            PrintAnswer(answer.first);
        }    
    }
    std::cout << "LU decomposition\n";
    for (int i = 0; i < 5; i++) {
        testPath[testDir.length()] = toascii(i+48);
        std::pair<Matrix_t, Matrix_t> test = EnterSLAE(3, testPath);
        std::cout << i << ": ";
        PrintAnswer(LU(test.first, test.second));
    } 
    std::cout << "QR decomposition\n";
    for (int i = 0; i < 5; i++) {
        testPath[testDir.length()] = toascii(i+48);
        std::pair<Matrix_t, Matrix_t> test = EnterSLAE(3, testPath);
        std::cout << i << ": ";
        PrintAnswer(SimpleIterMethod(test.first, test.second, pow(10,-4)).first);
    } 
    return 0;
}
