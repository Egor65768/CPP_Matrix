#include "gtest/gtest.h"

using namespace std;

int main() {
  testing::InitGoogleTest();
  if (RUN_ALL_TESTS()) {
    cout << "Tests failed" << endl;
  }
}