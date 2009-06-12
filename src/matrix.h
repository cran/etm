#ifndef MATRIX
#define MATRIX

#include <iostream>
#include <vector>

#include <math.h>


using namespace std;


class Vector;
class Matrix;


typedef vector< double, allocator<double> > dVector;          // vector of doubles (double vectorx)
typedef dVector::iterator itVector;

typedef vector<Vector> dMatrix;          // vector of double vectors (double Matrix)
typedef dMatrix::iterator itMatrix;

typedef vector<Matrix> dArray;          // vector of double matrices (Array of double matrices)
typedef dArray::iterator itArray;

class Vector : public dVector
{
public:
  Vector();
   
  Vector(const int n);

  Vector(double* v, const int n);

  friend ostream& operator <<(ostream& s, const Vector& v);

  Vector operator +(const Vector& v);

  Vector operator -(const Vector& v);

  friend Vector operator*(double x, const Vector& v);

  friend Vector operator*(const Vector& v, double x);

  Vector operator*(const Vector& v);

  friend double scalar(const Vector& v1, const Vector& v2);

  double max()const;

  double min()const;

  double mean()const;

  void zero(void);

  double norm(void)const;

  Vector diff() const;

  void as_double(double* a);

};

class Matrix: public dMatrix
{
  	
public:
	
  //constructor	
  Matrix();

  Matrix(const int n);
	
  Matrix(const int r,const int c);

  Matrix(double* m, const int r,const int c);
	
  friend ostream& operator<<(ostream& s, const Matrix& m);

  Matrix operator*(const Matrix& m);

  Vector operator*(const Vector& v);

  friend Matrix  operator*(const double, const Matrix& m);
	
  friend Matrix operator*(const Matrix& m, const double);
  
  Matrix operator+(const Matrix& m);
	
  Matrix operator-(const Matrix& m);

  void zero(void);
	
  void identity(void);

  void as_double(double* a);

};


class Array : public dArray
{

 public:
  Array();
 
  Array(const int l);

  Array(const int r,const int c, const int l);

  Array(double*a, const int r, const int c, const int l);

  friend Array operator*(const Matrix& m, const Array& a);

  friend Array operator*(const Array& a, const Matrix& m);

  friend ostream& operator<<(ostream& s,Array a);
 
  void as_double(double* a);

};

#endif
