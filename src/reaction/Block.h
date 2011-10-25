#ifndef BLOCK_H_
#define BLOCK_H_

#include <math.h>

#include <iostream>
#include <iomanip>
using namespace std;

class Block {
  
public:
  Block();
  Block(int);
  virtual ~Block();

  int getSize();
  double **getValues();
  double getRowAbsMax(int);

  void setValue(int, int, double);
  void setValues(double **);
  void setValues(int, int, Block *);

  void addValue(int, int, double);
  void addValues(double **);
  void addValues(int, int, Block *);

  void scaleRow(int, double);
  void scaleColumn(int, double);
  void scale(double);

  void zero(void);
  void setDiagonal(double d);

  void print();

  
private:

  int size;
  double **A;

};

#endif /*BLOCK_H_*/

