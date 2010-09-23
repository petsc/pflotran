#ifndef NEWTON_H_
#define NEWTON_H_

#include "Block.h"

#include <iostream>
using namespace std;

class Newton {
  
public:
  Newton(int);
  virtual ~Newton();

  void solve();
  void solve(double *r, Block *J, double *conc, double *update);
  
private:

  int n;
  double *x;
  double *b;
  Block *J;
  int *indices;

};

#endif /*NEWTON_H_*/
