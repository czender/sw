/* scp ~/c++/bug.cc esmf.ess.uci.edu:c++
   g++ -o bug bug.cc
   xlC -o bug bug.cc
   como -o bug bug.cc
   CC -LANG:std -64 -mips4 -o bug bug.cc */
#include <iostream>
#include <valarray>
#include <cstdlib>
int main()
{
  std::cout << "Testing integer representation using <climits>..." << std::endl;
  std::valarray<int> foo(73);
  std::cout << "valarray " << (foo.size() == 73 ? "worked" : "failed") << std::endl;
#ifdef SGIMP64
  std::abort(); // [fnc] Produce core dump
#else
  std::exit(EXIT_FAILURE); // [fnc] Exit nicely
#endif // !ABORT_ON_ERROR
}
