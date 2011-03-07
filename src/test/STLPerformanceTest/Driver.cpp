#include <iostream>
#include <ctime>

#include "Beaker.hpp"

int main(int argc, char **argv) {
  static_cast<void>(argc);
  static_cast<void>(argv);

  double result = 0.;

  Beaker *beaker = new Beaker();
  beaker->Setup(10000000,10);

  std::cout << "Starting Test 1...." << std::endl;
  std::clock_t start = std::clock();
  beaker->Test1();
  std::clock_t end = std::clock();
  result = beaker->GetSummationOfResults1();
  std::cout << "Finished...." << result << std::endl;
  std::cout << "Time: " << (end - start) / (double)CLOCKS_PER_SEC << std::endl;

  std::cout << "Starting Test 2...." << std::endl;
  start = std::clock();
  beaker->Test2();
  end = std::clock();
  result = beaker->GetSummationOfResults2or3();
  std::cout << "Finished...." << result << std::endl;
  std::cout << "Time: " << (end - start) / (double)CLOCKS_PER_SEC << std::endl;

  std::cout << "Starting Test 3...." << std::endl;
  start = std::clock();
  beaker->Test3();
  end = std::clock();
  result = beaker->GetSummationOfResults2or3();
  std::cout << "Finished...." << result << std::endl;
  std::cout << "Time: " << (end - start) / (double)CLOCKS_PER_SEC << std::endl;

  //----------------
#if 0
  std::cout << "Starting Test 2...." << std::endl;
  start = std::clock();
  beaker->Test2();
  end = std::clock();
  std::cout << "Finished...." << std::endl;
  std::cout << "Time: " << (end - start) / (double)CLOCKS_PER_SEC << std::endl;

  std::cout << "Starting Test 3...." << std::endl;
  start = std::clock();
  beaker->Test3();
  end = std::clock();
  std::cout << "Finished...." << std::endl;
  std::cout << "Time: " << (end - start) / (double)CLOCKS_PER_SEC << std::endl;

  std::cout << "Starting Test 1...." << std::endl;
  start = std::clock();
  beaker->Test1();
  end = std::clock();
  std::cout << "Finished...." << std::endl;
  std::cout << "Time: " << (end - start) / (double)CLOCKS_PER_SEC << std::endl;
#endif
  delete beaker;

}

