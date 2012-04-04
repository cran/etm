#include "matrix.h"


/*
  Vector class definitions
*/
Vector::Vector():dVector(){
}

Vector::Vector(const int n):dVector(n){
}

Vector::Vector(double* v, const int n):dVector(){
    for(int i=0;i<n;i++){
	push_back(v[i]);
    }
}

// ostream& operator<<(ostream& s, const Vector& v) {
//     if( v.size() > 0 ) s << "(";
    
//     for(int i=0;i<v.size();i++){
	
// 	s << v[i] ;
	
// 	if( i < (v.size()-1) ) {
// 	    s << ", ";
// 	}
// 	else {
// 	    s << ")"<<endl;
// 	}
//     }
    
//     return s;
// }

Vector operator*(double x, const Vector& v) {
    int i;
    Vector ans(v.size());
    for(i=0;i<v.size();i++) ans[i]= x*v[i];
    return ans;
}

Vector operator*(const Vector& v, double x) {
    int i;
    Vector ans(v.size());
    for(i=0;i<v.size();i++) ans[i]= x*v[i];
    return ans;
}

Vector Vector::operator +(const Vector& v) {
    int i;
    
    if( this->size() != v.size()) {
	// cout << "VECTOR Error: You're trying to add vectors of different sizes\n";
	// cout << v << endl;;
	// cout << *this << endl;
	
	return Vector();
    }
    
    Vector sum(this->size());
    for(i=0;i<this->size();i++) sum[i] = this->at(i)+v[i];
    return sum;
}

Vector Vector::operator -(const Vector& v) {
    int i;
    
    if( this->size() != v.size()) {
	// cout << "VECTOR Error: You're trying to subtract vectors of different sizes\n";
	// cout << v << endl;;
	// cout << *this << endl;
	
	return Vector();
    }
    
    Vector sum(this->size());
    for(i=0;i<this->size();i++) sum[i] = this->at(i)-v[i];
    return sum;
    
}

Vector Vector::operator*(const Vector& v) {
    int i;
    
    if( this->size() != v.size()) {
	// cout << "VECTOR Error: You're trying to multiply vectors of different sizes\n";
	// cout << v << endl;;
	// cout << *this << endl;
	
	return Vector();
    }
    
    Vector p(this->size());
    for(i=0;i<this->size();i++) p[i] = this->at(i)*v[i];
    return p;
}

double scalar(const Vector& v1, const Vector& v2) {
    int i;
    
    double p = 0.0;

    if( v1.size() != v2.size()) {
	// cout << "VECTOR Error: You're trying to multiply vectors of different sizes\n";
	// cout << v1 << endl;;
	// cout << v2 << endl;
	
	return p;
    }
    
    for(i=0;i<v1.size();i++) p += v1[i]*v2[i];
    return p;
}

double Vector::max(void)const {
    double max=this->at(0);
    for(int i=1;i<this->size();i++) if(this->at(i) > max) max=this->at(i);
    return max;
}

double Vector::min(void)const {
    double min=this->at(0);
    for(int i=1;i<this->size();i++) if(this->at(i) < min) min=this->at(i);
    return min;
}

double Vector::mean()const {
    double sum=0;
    for(int i=0;i<this->size();i++) sum += this->at(i);
    return sum/(this->size());
}

void Vector::zero(void) {
    for(int i=0;i<this->size();i++) this->at(i)=0.0;
    
    return;
}

double Vector::norm(void) const
{
    double s=0.0;
    for(int i=0;i<this->size();i++) s+= this->at(i)*this->at(i);
    return sqrt(s);
}

Vector Vector::diff() const
{
    int len = this->size();

    if( len > 1 ) {
	Vector v(len-1);

	for(int i=0; i < (len - 1); i++){
	    v[i] = this->at(i+1) - this->at(i);    
	}

	return v;
    }
		
    return Vector();
}

void Vector::as_double(double* a)
{
    for( int i=0; i<this->size(); i++ ) {	
	a[i] = this->at(i);    
    }
}

/*
  Matrix class definitions
*/
Matrix::Matrix():dMatrix(){
}

Matrix::Matrix(const int n):dMatrix(n){
}

Matrix::Matrix(const int r,const int c)
{
    Vector v(c);
	
    for(int i=0;i<r;i++) {
	push_back(v);
    }
}

Matrix::Matrix(double* m, const int r, const int c)
{
    for(int i=0;i<r;i++){
	Vector v(c);
            
	for(int j=0;j<c;j++) {        	 
	    v[j]=m[i+(j*r)];	       
	}
    
	push_back(v);	
    }
}

// ostream& operator<<(ostream& s,const Matrix& m)
// {
//     for(int i=0; i<m.size();i++) s << m[i];
//     return s;
// }


Matrix Matrix::operator*(const Matrix& m)
{
    if(this->size() == 0 || m.size() == 0 ) return Matrix();
    int rows = this->size();
    int cols = (this->at(0)).size();
    int m_rows = m.size();
    int m_cols = (m.at(0)).size();
    
    if(cols != m_rows) {	
	// cout << "MATRIX Error: Matrix Matrix::operator*(const Matrix& m):" << endl;		
	// cout << "matrices are the wrong size: " << cols << ", " << m_rows << endl;
    
	return Matrix();
    }
	
    Matrix ans(rows,m_cols);
	
    for(int i=0;i<rows;i++) {		
	for(int j=0;j<m_cols;j++) {
	    ans[i][j]=0.0;
	    for(int k=0;k<cols;k++){
		ans[i][j] += (this->at(i)).at(k)*m[k][j];
	    }
	}	
    }
	
    return ans;
}

Vector Matrix::operator*(const Vector& v)
{
    if(this->size() == 0 || v.size() == 0 ) return Vector();
    int rows = this->size();
    int cols = (this->at(0)).size();

    
    if(cols != v.size()) {			
	// cout << "MATRIX Error: multiplying matrix times Vector with wrong sizes\n";    
	return Vector();
    }
	
    Vector ans(rows);
	
    for(int i=0;i<rows;i++){		
	ans[i]=0.0;		
	for(int j=0;j<cols;j++) ans[i] += (this->at(i)).at(j)*v[j];	
    }
	
    return ans;
}

Matrix  operator*(const double x, const Matrix& m)
{
    if( m.size() == 0 ) return Matrix();

    int m_rows = m.size();
    int m_cols = (m.at(0)).size();

	
    Matrix a(m_rows,m_cols);
	
    for(int i=0; i<a.size();i++){
	a[i]= x*m[i];	  
    }
	
    return a;
}
	
Matrix operator*(const Matrix& m, const double x)
{
    if( m.size() == 0 ) return Matrix();

    int m_rows = m.size();
    int m_cols = (m.at(0)).size();

	
    Matrix a(m_rows,m_cols);
	
    for(int i=0; i<a.size();i++){
	a[i]= x*m[i];	  
    }
	
    return a;
}


Matrix Matrix::operator+(const Matrix& m)
{
    if(this->size() == 0 || m.size() == 0 ) return Matrix();
    int rows = this->size();
    int cols = (this->at(0)).size();
    int m_rows = m.size();
    int m_cols = (m.at(0)).size();

	
    if(rows!= m_rows || cols != m_cols) {
	// cout << "MATRIX Error: you're trying to add matrices of different sizes\n";
	return Matrix();	
    }
	
    Matrix ans(m_rows, m_cols);
	
    for(int i=0;i<m_rows;i++){		
	for(int j=0;j<m_cols;j++){			
	    ans[i][j] = (this->at(i)).at(j) + m[i][j]; 		
	}	
    }
	
    return ans;
}

Matrix Matrix::operator-(const Matrix& m)
{
    if(this->size() == 0 || m.size() == 0 ) return Matrix();
    int rows = this->size();
    int cols = (this->at(0)).size();
    int m_rows = m.size();
    int m_cols = (m.at(0)).size();

	
    if(rows!= m_rows || cols != m_cols) {
	// cout << "MATRIX Error: you're trying to add matrices of different sizes\n";
	return Matrix();	
    }
	
    Matrix ans(m_rows, m_cols);
	
    for(int i=0;i<m_rows;i++){		
	for(int j=0;j<m_cols;j++){			
	    ans[i][j] = (this->at(i)).at(j) - m[i][j]; 		
	}	
    }
	
    return ans;
}

void Matrix::zero(void)
{
    if(this->size() == 0 ) return;

    int rows = this->size();
    int cols = (this->at(0)).size();

    for(int i=0;i<rows;i++) {
		
	for(int j=0;j<cols;j++) {			
	    (this->at(i)).at(j)=0.0;		
	}
	
    }
	
    return;

}

void Matrix::identity(void)
{
    if(this->size() == 0 ) return;

    int rows = this->size();
    int cols = (this->at(0)).size();

	
    if(rows!=cols) {		
	// cout << "MATRIX Error: Matrix::identity(): Matrix not square\n";	
    }

	
    zero();
	
    for(int i=0;i<rows;i++) (this->at(i)).at(i)=1.0;
	
    return;
}

void  Matrix::as_double(double* a)
{
    int rows = this->size();  
    
    for( int i=0; i<rows; i++ ) {

	int cols = (this->at(i)).size();

	for( int j=0; j<cols; j++ ) {
	
	    a[i+(j*rows)] = (this->at(i)).at(j);
      
	}    
    }
}



/*
  Array class definitions
*/
Array::Array():dArray(){ 
}

Array::Array(const int len):dArray(len){
}

Array::Array(const int rows, const int cols, const int len)
{
    Matrix m(rows, cols);
	
    for(int i=0;i<len;i++) {
	push_back(m);
    }
}

Array::Array(double*a, const int rows, const int cols, const int len)
{

    for( int k=0; k<len; k++ ) {
     
	Matrix m( rows, cols );

	for( int i=0; i<rows; i++ ) {
	    for( int j=0; j<cols; j++ ) {
		m[i][j] = a[(i+(j*rows)) + k*(rows*cols)] ;
	    }
	}
   
	push_back(m);     
    }
}

// ostream& operator<<(ostream& s, Array a)
// {
//     for(int i=0; i<a.size();i++) {         
// 	s << a[i] << endl;
//     }
	
//     return s;
// }

void  Array::as_double(double* a)
{
    int len = this->size();

    for( int k=0; k<len; k++ ) {

	int rows = (this->at(k)).size();

	for( int i=0; i<rows; i++ ) {

	    int cols = ((this->at(k)).at(i)).size();

	    for( int j=0; j<cols; j++ ) {
		a[(i+(j*rows)) + k*(rows*cols)] = ((this->at(k)).at(i)).at(j);
	    }
	}
    }
}

Array operator*(const Matrix& m, const Array& a)
{
    int len = a.size();

    Array aj;

    for( int k=0; k<len; k++ ) {
	aj.push_back((Matrix)m * (Matrix)a[k]);
    }

    return aj;
}

Array operator*(const Array& a, const Matrix& m)
{
    int len = a.size();

    Array aj;

    for( int k=0; k<len; k++ ) {
	aj.push_back((Matrix)a[k] * (Matrix)m);
    }

    return aj;
}
