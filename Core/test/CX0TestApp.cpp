//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Test-Application to check c++11 support.
/*!
 * @file CX0TestApp.cpp
 * This tiny application test some features of the new c++11 standard. One can
 * compile and run this to check the c++11 support of the used compiler/system.
 * The features tested are: shared_pointer, raw strings, auto type's, lamdas,
 * threading and regular expressions.
 */

#include <algorithm>
#include <iostream>
#include <memory>
#include <regex>
#include <string>
#include <thread>

#include "Core/Logging.hpp"

using namespace std;

static const unsigned int num_threads = 10;

// This function will be called from a thread

void call_from_thread(unsigned int tid) {
  cout << "Launched by thread " << tid << endl;
}

class Base {
public:
  Base(const unsigned int in) : i(in) {}
  virtual ~Base(){};
  virtual const unsigned int GetInt() { return i; }
  virtual void SetInt(const unsigned int in) = 0;

protected:
  unsigned int i;
};

class Derived : public Base {
public:
  Derived(const unsigned int in) : Base(in) {}
  virtual ~Derived(){};
  virtual void SetInt(const unsigned int in) { i = in * 10; }
};

int main(int argc, char **argv) {
  std::cout << "  ComPWA Copyright (C) 2013  Mathias Michel " << std::endl;
  std::cout << "  This program comes with ABSOLUTELY NO WARRANTY; for details "
               "see license.txt"
            << std::endl;
  std::cout << std::endl;

  // unsigned int a=1,b=2;
  // int c[3][3];
  // int d = c[a,b];

  //=== shared_ptr in a vector: ===
  vector<shared_ptr<Base>> baseVec;
  baseVec.push_back(shared_ptr<Derived>(new Derived(5)));
  baseVec.push_back(shared_ptr<Derived>(new Derived(10)));
  /*"old" way:
   * vector<Base*> baseVec;
   * baseVec.push_back(new Derived(5));
   * baseVec.push_back(new Derived(10));
   */

  //========= c++11 test ==========
  string hello = "say hello to c++11 Basics: ";
  int test = 7;
  auto neu = test;
  shared_ptr<int> test_ptr(new int(8));
  auto *auto_ptr(new int(9));
  cout << "String: " << hello << endl;
  cout << "auto type: " << neu << endl;
  cout << "shared pointer:  " << *test_ptr << endl;
  cout << "auto pointer:  " << *auto_ptr << endl;
  cout << endl;

  //========= Raw Strings ==========
  cout << "Raw Strings: " << endl;
  string normal_str = "First line.\nSecond line.\nEnd of message.\n";
  string raw_str = "(First line.\nSecond line.\nEnd of message.\n)";
  cout << normal_str << endl;
  cout << raw_str << endl;
  cout << endl;

  //=========  Lambdas  ==========
  cout << "Lamda-functions: " << endl;
  int n = 10;
  vector<int> v(n);
  // initialize the vector with values from 10 to 1
  for (int i = n - 1, j = 0; i >= 0; i--, j++)
    v[j] = i + 1;
  // print the unsorted vector
  for (int i = 0; i < n; i++)
    cout << v[i] << " " << endl;
  // sort the vector
  sort(v.begin(), v.end(), [](int i, int j) -> bool { return (i < j); });
  // print the sorted vector
  for (int i = 0; i < n; i++)
    cout << v[i] << " " << endl;
  cout << endl;

  //=========  Threads  ==========
  cout << "Threads:\n";
  std::thread t[num_threads];
  // Launch a group of threads
  for (unsigned int i = 0; i < num_threads; ++i) {
    t[i] = std::thread(call_from_thread, i);
  }
  cout << "Launched from the main\n" << endl;
  // Join the threads with the main thread
  for (unsigned int i = 0; i < num_threads; ++i) {
    t[i].join();
  }
  cout << endl << "Launched from the main after join\n" << endl;
  cout << endl;

  //========= Regular Expr. ==========
  cout << "Regular Expressions: " << endl;
  cout << "Not supported by gcc up to 4.8" << endl;
  string input;
  regex integer("(\\+|-)?[[:digit:]]+");
  // As long as the input is correct ask for another number
  //  while(true)
  //  {
  cout << "Give me an integer!" << endl;
  //    cin>>input;
  //  input = "1"; //we do not want user interaction in test cases
  // Exit when the user inputs q
  //    if(input=="q")
  //      break;
  //    if(regex_match(input,integer))
  //      cout<<"integer"<<endl;
  //    else{
  //      cout<<"Invalid input"<<endl;
  //    }
  //  }
  cout << endl;

  return 0;
}
