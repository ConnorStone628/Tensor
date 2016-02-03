#ifndef __TENSOR__
#define __TENSOR__

#include <stdarg.h>
#include <vector>
#include <iostream>
#include <stdexcept>
template<class T1, class T2> T1 min(const T1& a, const T2& b){return a>b ? a:b;}//fixme deleteme

template <class DATATYPE = double, class INDEXTYPE = int>
class tensor {
  //Contains the dimensions of the tensor, lengths of each dimension, limit for the number of elements to print, the total number of dimensions
  INDEXTYPE *shape, *offsets, printlimit, ndimensions;
  
  //Pointer to the beginning of the tensor data
  DATATYPE* origin;

  //if tensor is using the data from some other object, it will not delete the data when it is deleted
  bool usingother;
  
  //Given the requested dimensions, it fills the shape and offsets objects, also creates the tensor in memory and assigns it to the origin
  void InitializeSelf(const std::vector<INDEXTYPE>& D);
  void InitializeSelf(const tensor<INDEXTYPE,INDEXTYPE>& D);
  
  //used by operators to guarintee consistancy of implimentation
  inline tensor<DATATYPE,INDEXTYPE>& operate(const tensor<DATATYPE,INDEXTYPE>& operand, DATATYPE (tensor<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&)const) const;
  inline tensor<DATATYPE,INDEXTYPE>& operate(const DATATYPE& operand, DATATYPE (tensor<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&)const) const;
  inline void operateeq(const tensor<DATATYPE,INDEXTYPE>& operand, DATATYPE (tensor<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&)const);
  inline void operateeq(const DATATYPE& operand, DATATYPE (tensor<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&)const);
  inline DATATYPE mul(const DATATYPE& x1, const DATATYPE& x2) const{return x1*x2;}
  inline DATATYPE div(const DATATYPE& x1, const DATATYPE& x2) const{return x1/x2;}
  inline DATATYPE add(const DATATYPE& x1, const DATATYPE& x2) const{return x1+x2;}
  inline DATATYPE sub(const DATATYPE& x1, const DATATYPE& x2) const{return x1-x2;}  
  
public:

  //Empty initializer
  tensor();

  //Initializes the tensor with no guarintee about the values in its indices
  tensor(const std::vector<INDEXTYPE>& D);

  //Initializes the tensor and assigns a value to every position
  tensor(const std::vector<INDEXTYPE>& D, DATATYPE value);

  //Constructor to interpret data which is already in memory
  tensor(const std::vector<INDEXTYPE>& D, DATATYPE* begin);

  //Initializes the tensor with no guarintee about the values in its indices
  tensor(const tensor<INDEXTYPE,INDEXTYPE>& D);

  //Initializes the tensor and assigns a value to every position
  tensor(const tensor<INDEXTYPE,INDEXTYPE>& D, DATATYPE value);

  //Constructor to interpret data which is already in memory
  tensor(const tensor<INDEXTYPE,INDEXTYPE>& D, DATATYPE* begin);

  //Constructor for manual assignment of shape
  tensor(const INDEXTYPE& n, const INDEXTYPE& d1, ...);

  //Destructor, will call delete on origin
  ~tensor();

  //Rewrite the tensor using the different shape, WARNING: his does not do anything fancy, it will make an tensor with the requested shape then write as much as it can from the old tensor
  void Resize(const std::vector<INDEXTYPE>& D);

  //Rewrite the tensor using the different shape, WARNING: his does not do anything fancy, it will make an tensor with the requested shape then write as much as it can from the old tensor. In this version you can write out the new dimensions manually.
  void Resize(const INDEXTYPE& n, const INDEXTYPE& d1, ...);

  //Given an index along the memory, it returns a vector which has the distances along each dimension
  //std::vector<INDEXTYPE> GetIndex(const INDEXTYPE& index) const;
  tensor<INDEXTYPE, INDEXTYPE> GetIndex(const INDEXTYPE& index) const;

  //Given the distances along each dimension, this will return the true index of the desired data point
  //INDEXTYPE GetIndex(const std::vector<INDEXTYPE>& index) const;
  INDEXTYPE GetIndex(const tensor<INDEXTYPE, INDEXTYPE>& index) const;

  //Returns the value at an index location
  //  DATATYPE* GetValueP(const std::vector<INDEXTYPE>& index) const;
  DATATYPE* GetValueP(const tensor<INDEXTYPE, INDEXTYPE>& index) const;

  //Returns the value at the index
  DATATYPE* GetValueP(const INDEXTYPE& index) const;

  //Returns the value at an index location
  DATATYPE GetValue(const tensor<INDEXTYPE, INDEXTYPE>& index) const;

  //Returns the value at the index
  DATATYPE GetValue(const INDEXTYPE& index) const;

  //Returns a pointer to the beginning of the tensor
  DATATYPE* Begin() const;

  //Returns a pointer to the end of the tensor
  DATATYPE* End() const;

  //returns the number of dimensions in the tensor
  INDEXTYPE NDimensions() const;
  
  //Returns the size of the tensor along the given dimension. note1: size(0) gives the size of the whole tensor in memory, note2: size(n) is always 1
  INDEXTYPE Size(const INDEXTYPE& dimension = 0) const;

  //Returns the vector of values describing how long each dimension is
  tensor<INDEXTYPE,INDEXTYPE> Shape() const;

  //Returns the length of one of the dimensions, as specified
  INDEXTYPE Shape(const INDEXTYPE& dimension) const;

  //Returns the beginning and end pointer for each segment of the tensor according to the dimension requested 
  tensor<DATATYPE*,INDEXTYPE> Axis(INDEXTYPE dimension) const;
  
  //Assign all the values from a vector
  void Assign(const std::vector<DATATYPE>& values);

  //Assign all the values from an tensor
  void Assign(const tensor<DATATYPE,INDEXTYPE>& values);
  
  //Assign the values at the pointed to locations
  void Assign(const DATATYPE* databegin, const DATATYPE* dataend);

  //Assign the given value to an index
  void Assign(const std::vector<INDEXTYPE>& index, const DATATYPE& value);

  //Assign the given value to an index
  void Assign(const INDEXTYPE& index, const DATATYPE& value);

  //Assign the given value to every element in the tensor
  void Assign(const DATATYPE& value);

  //Assign a list of values to the tensor, the first argument is the number of values to assign
  void ManualAssign(const INDEXTYPE& n, const DATATYPE& v1, ...);

  //Use a pre-assigned contiguous block of memory as an tensor. WARNING: Two tensors using the same block will cause an error during destruction, it should have expected behaviour during runtime though.
  void UseBlock(const std::vector<INDEXTYPE>& D, DATATYPE* begin);

  // prints the tensor to std::cout
  void Print() const;
  void SetPrintLimit(INDEXTYPE newprintlimit);

  tensor<DATATYPE,INDEXTYPE>& operator*(const tensor<DATATYPE,INDEXTYPE>& operand) const;

  tensor<DATATYPE,INDEXTYPE>& operator*(const DATATYPE& operand) const;
  
  void operator*=(const tensor<DATATYPE,INDEXTYPE>& operand);
  
  void operator*=(const DATATYPE& operand);

  tensor<DATATYPE,INDEXTYPE>& operator+(const tensor<DATATYPE,INDEXTYPE>& operand) const;

  tensor<DATATYPE,INDEXTYPE>& operator+(const DATATYPE& operand) const;

  void operator+=(const tensor<DATATYPE,INDEXTYPE>& operand);

  void operator+=(const DATATYPE& operand);

  tensor<DATATYPE,INDEXTYPE>& operator-(const tensor<DATATYPE,INDEXTYPE>& operand) const;

  tensor<DATATYPE,INDEXTYPE>& operator-(const DATATYPE& operand) const;

  void operator-=(const tensor<DATATYPE,INDEXTYPE>& operand);

  void operator-=(const DATATYPE& operand);

  tensor<DATATYPE,INDEXTYPE>& operator/(const tensor<DATATYPE,INDEXTYPE>& operand) const;

  tensor<DATATYPE,INDEXTYPE>& operator/(const DATATYPE& operand) const;

  void operator/=(const tensor<DATATYPE,INDEXTYPE>& operand);

  void operator/=(const DATATYPE& operand);

  tensor<DATATYPE,INDEXTYPE>& operator()(const INDEXTYPE& index) const;

  DATATYPE& operator[](const INDEXTYPE& index) const;
  
  tensor<DATATYPE,INDEXTYPE>& operator=(const tensor<DATATYPE,INDEXTYPE>& operand);

  bool operator==(const tensor<DATATYPE,INDEXTYPE>& operand) const;
};

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::InitializeSelf(const std::vector<INDEXTYPE>& D){  
  for (INDEXTYPE i = 0; i < D.size(); ++i){
    if (D[i] < (INDEXTYPE) 1){throw std::invalid_argument( "Array cannot have dimension length of < 1" );}
  }
  
  this->ndimensions = D.size();
  this->shape = new INDEXTYPE[this->ndimensions];
  for (INDEXTYPE i = 0; i < this->ndimensions; ++i){this->shape[i] = D[i];}
  this->offsets = new INDEXTYPE[this->ndimensions+1];
  for (INDEXTYPE a = 0; a < this->ndimensions+1; ++a){
    INDEXTYPE axislength = 1;
    for (INDEXTYPE d = a; d < this->ndimensions; ++d){axislength *= this->shape[d];}
    this->offsets[a] = axislength;
  }
  this->printlimit = 150;
}   
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::InitializeSelf(const tensor<INDEXTYPE,INDEXTYPE>& D){
  for (INDEXTYPE i = 0; i < D.Size(); ++i){
    if (D[i] < (INDEXTYPE) 1){throw std::invalid_argument( "Array cannot have dimension length of < 1" );}
  }
  
  this->ndimensions = D.Size();
  this->shape = new INDEXTYPE[this->ndimensions];
  for (INDEXTYPE i = 0; i < this->ndimensions; ++i){this->shape[i] = D[i];}
  this->offsets = new INDEXTYPE[this->ndimensions+1];
  for (INDEXTYPE a = 0; a < this->ndimensions+1; ++a){
    INDEXTYPE axislength = 1;
    for (INDEXTYPE d = a; d < this->ndimensions; ++d){axislength *= this->shape[d];}
    this->offsets[a] = axislength;
  }
  this->printlimit = 150;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE, INDEXTYPE>::tensor(){
  std::vector<INDEXTYPE> D(1,1);
  this->InitializeSelf(D);
  this->origin = new DATATYPE[this->offsets[0]];
  this->usingother = false;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE, INDEXTYPE>::tensor(const std::vector<INDEXTYPE>& D){
  this->InitializeSelf(D);
  this->origin = new DATATYPE[this->offsets[0]];
  this->usingother = false;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE, INDEXTYPE>::tensor(const std::vector<INDEXTYPE>& D, DATATYPE value){
  this->InitializeSelf(D);
  this->origin = new DATATYPE[this->offsets[0]];
  this->usingother = false;
  this->Assign(value);
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE, INDEXTYPE>::tensor(const std::vector<INDEXTYPE>& D, DATATYPE* begin){
  this->InitializeSelf(D);
  this->origin = begin;
  this->usingother = true;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE, INDEXTYPE>::tensor(const tensor<INDEXTYPE,INDEXTYPE>& D){
  this->InitializeSelf(D);
  this->origin = new DATATYPE[this->offsets[0]];
  this->usingother = false;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE, INDEXTYPE>::tensor(const tensor<INDEXTYPE,INDEXTYPE>& D, DATATYPE value){
  this->InitializeSelf(D);
  this->origin = new DATATYPE[this->offsets[0]];
  this->usingother = false;
  this->Assign(value);
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE, INDEXTYPE>::tensor(const tensor<INDEXTYPE,INDEXTYPE>& D, DATATYPE* begin){
  this->InitializeSelf(D);
  this->origin = begin;
  this->usingother = true;
}

template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE, INDEXTYPE>::tensor(const INDEXTYPE& n, const INDEXTYPE& d1, ...){
  std::vector<INDEXTYPE> D(n,d1);
  if (n > 1){
    va_list values;
    va_start(values, d1);
    for (INDEXTYPE i = 1; i <= n; ++i){
      D[i] = va_arg(values, INDEXTYPE);
      D[i] = D[i] < 1 ? 1:D[i];
    }
    va_end(values);
  }
  this->InitializeSelf(D);
  this->origin = new DATATYPE[this->offsets[0]];
  this->usingother = false;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE, INDEXTYPE>::~tensor(){
  if (!this->usingother) {delete[] this->origin;}
  delete[] this->offsets;
  delete[] this->shape;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::Resize(const std::vector<INDEXTYPE>& D){
  bool allequal = true;
  if (D.size() != this->ndimensions){allequal = false;}
  else{for (INDEXTYPE i = 0; i < D.size(); ++i){if (D[i] != this->shape[i]){allequal = false; break;}}}
  
  if (!allequal){
    delete[] this->shape;
    INDEXTYPE oldlength = this->offsets[0];
    delete[] this->offsets;
    DATATYPE* oldorigin = this->origin;
    this->InitializeSelf(D);
    this->origin = new DATATYPE[this->offsets[0]];
    this->Assign(oldorigin, oldorigin + min(oldlength, this->offsets[0]) - 1);
    if (!this->usingother) {delete[] oldorigin;}
    this->usingother = false;
  }
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::Resize(const INDEXTYPE& n, const INDEXTYPE& d1, ...){
  std::vector<INDEXTYPE> D(n,d1);
  if (n > 1){
    va_list values;
    va_start(values, d1);
    for (INDEXTYPE i = 0; i < n; ++i){
      D[i] = va_arg(values, INDEXTYPE);
    }
    va_end(values);
  }
  this->Resize(D);
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<INDEXTYPE, INDEXTYPE> tensor<DATATYPE, INDEXTYPE>::GetIndex(const INDEXTYPE& index) const{
  tensor<INDEXTYPE, INDEXTYPE> newindex(1,this->ndimensions);
  for (INDEXTYPE i = 0; i < this->ndimensions; ++i){
    newindex[i] = (index%this->offsets[i])/this->offsets[i+1];
  }
  return newindex;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
INDEXTYPE tensor<DATATYPE, INDEXTYPE>::GetIndex(const tensor<INDEXTYPE,INDEXTYPE>& index) const{
  INDEXTYPE newindex = 0;
  for (INDEXTYPE i = 0; i < index.Size(); ++i){
    newindex += index[i]*this->offsets[i+1];
  }
  return newindex;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
DATATYPE* tensor<DATATYPE, INDEXTYPE>::GetValueP(const tensor<INDEXTYPE, INDEXTYPE>& index) const{
  return (this->origin + this->GetIndex(index));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
DATATYPE* tensor<DATATYPE, INDEXTYPE>::GetValueP(const INDEXTYPE& index) const{
  return (this->origin + index);
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
DATATYPE tensor<DATATYPE, INDEXTYPE>::GetValue(const tensor<INDEXTYPE, INDEXTYPE>& index) const{
  return *(GetValueP(index));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
DATATYPE tensor<DATATYPE, INDEXTYPE>::GetValue(const INDEXTYPE& index) const{
  return *(GetValueP(index));
}

template <class DATATYPE, class INDEXTYPE>
INDEXTYPE tensor<DATATYPE, INDEXTYPE>::NDimensions() const{
  return this->ndimensions;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
INDEXTYPE tensor<DATATYPE, INDEXTYPE>::Size(const INDEXTYPE& dimension) const{
  return this->offsets[dimension];
}
  
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<INDEXTYPE,INDEXTYPE> tensor<DATATYPE, INDEXTYPE>::Shape() const{
  tensor<INDEXTYPE,INDEXTYPE> returnshape(1,this->ndimensions);
  returnshape.Assign(this->shape, this->shape + this->ndimensions - 1);
  return returnshape;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
INDEXTYPE tensor<DATATYPE, INDEXTYPE>::Shape(const INDEXTYPE& dimension) const{
  return this->shape[dimension];
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
DATATYPE* tensor<DATATYPE, INDEXTYPE>::Begin() const{return this->origin;}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
DATATYPE* tensor<DATATYPE, INDEXTYPE>::End() const{return this->origin+this->offsets[0]-1;}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE*, INDEXTYPE> tensor<DATATYPE, INDEXTYPE>::Axis(INDEXTYPE dimension) const{
  if (dimension > this->ndimensions){throw std::invalid_argument( "Invalid dimension request, greater than ndimensions" );}
  if (dimension < 0){dimension = this->ndimensions - 1;}
  tensor<DATATYPE*, INDEXTYPE> axis(2,this->offsets[0]/this->offsets[dimension], 2);
  for (INDEXTYPE i = 0; i < axis.Shape(0); ++i){
    axis[2*i] = this->origin + i*this->offsets[dimension];
    axis[2*i+1] = this->origin + (i+1)*this->offsets[dimension] - 1;
  }
  return axis;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::Assign(const std::vector<DATATYPE>& values){
  for (INDEXTYPE i = 0; i < values.size(); ++i){
    *(this->origin+i) = values[i];
  }
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::Assign(const tensor<DATATYPE,INDEXTYPE>& values){
  for (INDEXTYPE i = 0; i < values.Size(); ++i){
    *(this->origin+i) = values[i];
  }
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::Assign(const DATATYPE* databegin, const DATATYPE* dataend){
  for (INDEXTYPE i = 0; i < min(dataend-databegin+1, this->offsets[0]); ++i){
    *(this->origin+i) = *(databegin + i);
  }
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::Assign(const std::vector<INDEXTYPE>& index, const DATATYPE& value){
  *(this->origin + this->GetIndex(index)) = value;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::Assign(const INDEXTYPE& index, const DATATYPE& value){
  *(this->origin + index) = value;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::Assign(const DATATYPE& value){
  for (INDEXTYPE i = 0; i < this->offsets[0]; ++i){
    *(this->origin + i) = value;
  }
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::ManualAssign(const INDEXTYPE& n, const DATATYPE& v1, ...){
  *(this->origin) = v1;
  if (n > 1){
    va_list values;
    va_start(values, v1);
    for (INDEXTYPE i = 1; i <= n; ++i){
      *(this->origin+i) = va_arg(values, DATATYPE);
    }
    va_end(values);
  }
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::UseBlock(const std::vector<INDEXTYPE>& D, DATATYPE* begin){
  if ( this->Begin() <=  begin && begin <= this->End() ){std::invalid_argument( " Cannot use self as memory block " );}
  else{
    if (!this->usingother) {delete[] this->origin;}
    this->InitializeSelf(D);
    this->origin = begin;
    this->usingother = true;
  }
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::Print() const{
  if (this->offsets[0] > this->printlimit){std::cout<<"tensor too long: " << this->offsets[0] << " elements" << std::endl; return;}
  
  for (INDEXTYPE i = 0; i < this->offsets[0]; ++i){
    bool comma = true;
    for (INDEXTYPE dim = 0; dim < this->ndimensions; ++dim){
      if (i%this->offsets[dim] == 0){std::cout << "|";comma = false;}
    }
    if (comma){std::cout << ", ";}
    std::cout << this->GetValue(i);
  }
  for (INDEXTYPE dim = 0; dim < this->ndimensions; ++dim){std::cout << "|";}
  std::cout<<std::endl;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::SetPrintLimit(INDEXTYPE newprintlimit){
  this->printlimit = newprintlimit;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE,INDEXTYPE>& tensor<DATATYPE, INDEXTYPE>::operator*(const tensor<DATATYPE,INDEXTYPE>& operand) const{
  return this->operate(operand, (&tensor<DATATYPE, INDEXTYPE>::mul));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE,INDEXTYPE>& tensor<DATATYPE, INDEXTYPE>::operator*(const DATATYPE& operand) const{
  return this->operate(operand, (&tensor<DATATYPE, INDEXTYPE>::mul));
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE,INDEXTYPE>& tensor<DATATYPE, INDEXTYPE>::operator+(const tensor<DATATYPE,INDEXTYPE>& operand) const{
  return this->operate(operand, (&tensor<DATATYPE, INDEXTYPE>::add));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE,INDEXTYPE>& tensor<DATATYPE, INDEXTYPE>::operator+(const DATATYPE& operand) const{
  return this->operate(operand, (&tensor<DATATYPE, INDEXTYPE>::add));
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE,INDEXTYPE>& tensor<DATATYPE, INDEXTYPE>::operator-(const tensor<DATATYPE,INDEXTYPE>& operand) const{
  return this->operate(operand, (&tensor<DATATYPE, INDEXTYPE>::sub));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE,INDEXTYPE>& tensor<DATATYPE, INDEXTYPE>::operator-(const DATATYPE& operand) const{
  return this->operate(operand, (&tensor<DATATYPE, INDEXTYPE>::sub));
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE,INDEXTYPE>& tensor<DATATYPE, INDEXTYPE>::operator/(const tensor<DATATYPE,INDEXTYPE>& operand) const{
  return this->operate(operand, (&tensor<DATATYPE, INDEXTYPE>::div));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE,INDEXTYPE>& tensor<DATATYPE, INDEXTYPE>::operator/(const DATATYPE& operand) const{
  return this->operate(operand, (&tensor<DATATYPE, INDEXTYPE>::div));
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::operator*=(const tensor<DATATYPE,INDEXTYPE>& operand){
  return this->operateeq(operand, (&tensor<DATATYPE, INDEXTYPE>::mul));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::operator*=(const DATATYPE& operand){
  return this->operateeq(operand, (&tensor<DATATYPE, INDEXTYPE>::mul));
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::operator+=(const tensor<DATATYPE,INDEXTYPE>& operand){
  return this->operateeq(operand, (&tensor<DATATYPE, INDEXTYPE>::add));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::operator+=(const DATATYPE& operand){
  return this->operateeq(operand, (&tensor<DATATYPE, INDEXTYPE>::add));
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::operator-=(const tensor<DATATYPE,INDEXTYPE>& operand){
  return this->operateeq(operand, (&tensor<DATATYPE, INDEXTYPE>::sub));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::operator-=(const DATATYPE& operand){
  return this->operateeq(operand, (&tensor<DATATYPE, INDEXTYPE>::sub));
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::operator/=(const tensor<DATATYPE,INDEXTYPE>& operand){
  return this->operateeq(operand, (&tensor<DATATYPE, INDEXTYPE>::div));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void tensor<DATATYPE, INDEXTYPE>::operator/=(const DATATYPE& operand){
  return this->operateeq(operand, (&tensor<DATATYPE, INDEXTYPE>::div));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE,INDEXTYPE>& tensor<DATATYPE, INDEXTYPE>::operator()(const INDEXTYPE& index) const{
  std::vector<INDEXTYPE> newshape;
  if (this->ndimensions > 1){
    for (INDEXTYPE i = 1; i < this->ndimensions; ++i){
      newshape.push_back(*(this->shape+i));
    }
  }else{
    throw std::invalid_argument( "Cannot construct sub-tensor from 1D tensor, (Perhaps you want to use element access [] operator?)" );
  }
  tensor<DATATYPE,INDEXTYPE>* newtensor = new tensor<DATATYPE,INDEXTYPE>(newshape, origin + index*this->offsets[1]);
  return *newtensor;
}


//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
DATATYPE& tensor<DATATYPE, INDEXTYPE>::operator[](const INDEXTYPE& index) const{
  return *this->GetValueP(index);
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
tensor<DATATYPE,INDEXTYPE>& tensor<DATATYPE, INDEXTYPE>::operator=(const tensor<DATATYPE,INDEXTYPE>& operand){
  if (this != &operand){
    if (!this->usingother) {delete[] this->origin;}
    this->InitializeSelf(operand.Shape());
    this->origin = new DATATYPE[this->offsets[0]];
    this->usingother = false;
    this->Assign(operand.Begin(), operand.End());
  }
  return *this;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
bool tensor<DATATYPE, INDEXTYPE>::operator==(const tensor<DATATYPE,INDEXTYPE>& operand) const{
  if (this == &operand){return true;
  }else if (this->ndimensions != operand.NDimensions()){return false;}
  for (INDEXTYPE i = 0; i < this->ndimensions; ++i){if (this->shape[i] != operand.Shape(i)){return false;}}
  for (INDEXTYPE i = 0; i < this->offsets[0]; ++i){if (this->GetValue(i) != operand.GetValue(i)){return false;}}
  return true;
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
inline tensor<DATATYPE,INDEXTYPE>& tensor<DATATYPE, INDEXTYPE>::operate(const tensor<DATATYPE,INDEXTYPE>& operand, DATATYPE (tensor<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&) const) const{
  tensor<DATATYPE,INDEXTYPE>* result = new tensor<DATATYPE,INDEXTYPE>(operand.Size() < this->Size() ? operand.Shape():this->Shape());
      
  DATATYPE *operandbegin = operand.Begin(), *resultbegin = result->Begin();
  for (INDEXTYPE i = 0; i < this->offsets[0]; ++i){
    *(resultbegin + i) = (this->*op)(*(this->origin + i),*(operandbegin+i));
  }
  return *result;
}

template <class DATATYPE, class INDEXTYPE>
inline tensor<DATATYPE,INDEXTYPE>& tensor<DATATYPE, INDEXTYPE>::operate(const DATATYPE& operand, DATATYPE (tensor<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&) const) const{
  tensor<DATATYPE,INDEXTYPE>* result = new tensor(this->Shape());
      
  DATATYPE *resultbegin = result->Begin();
  for (INDEXTYPE i = 0; i < this->offsets[0]; ++i){
    *(resultbegin + i) = (this->*op)(*(this->origin + i), operand);
  }
  return *result;
}

template <class DATATYPE, class INDEXTYPE>
inline void tensor<DATATYPE, INDEXTYPE>::operateeq(const tensor<DATATYPE,INDEXTYPE>& operand, DATATYPE (tensor<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&) const){
  DATATYPE *operandbegin = operand.Begin();
  for (INDEXTYPE i = 0; i < this->offsets[0]; ++i){
    *(this->origin + i) = (this->*op)(*(this->origin + i),*(operandbegin+i));
  }
}

template <class DATATYPE, class INDEXTYPE>
inline void tensor<DATATYPE, INDEXTYPE>::operateeq(const DATATYPE& operand, DATATYPE (tensor<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&) const){
  for (INDEXTYPE i = 0; i < this->offsets[0]; ++i){
    *(this->origin + i) = (this->*op)(*(this->origin + i), operand);
  }
}

#endif
