// $Id$ 

// Purpose: Description (definition) of two dimensional array class

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <a2d.hh> // Two dimensional arrays

#ifndef A2D_HH // Contents have not yet been inserted in current source file
#define A2D_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class

// Standard C headers
#include <cassert> // Assertions

// 3rd party vendors

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()

// Forward declarations (DeD01 p. 500)
template<class val_T> class a2d_cls;
template<class val_T> std::ostream & operator<<(std::ostream &srm_out,const a2d_cls<val_T> &a2d_obj);

// Typedefs

// Define a2d_cls class
template<class val_T> // [obj] Object type
class a2d_cls{ // Two dimensional arrays
public:

  // Friends

  /* Symbol <> in following necessary as per comp.lang.c++ 20020226
     GCC 2.95+ and xlC require '<>' in forward declaration (but forbid it in prototype!) 
     ICC > 5.5 requires '<>' to be removed to compile correctly */

#if defined(__GNUC__) || defined(__xlC__) || defined(__INTEL_COMPILER) || defined(__PATHCC__) || defined(PGI_CXX) // Compiler is GCC, xlC, icpc, pathCC, or pgCC
friend std::ostream & // [srm] Reference to output stream for cascading
operator<< <> // [fnc] Stream insertion operator
(std::ostream &srm_out, // [srm] Output stream
 const a2d_cls &a2d_obj); // [obj] Object to insert in stream
#else // Other compiler (como, SGI CC, Sun CC)
friend std::ostream & // [srm] Reference to output stream for cascading
operator<< // [fnc] Stream insertion operator
(std::ostream &srm_out, // [srm] Output stream
 const a2d_cls &a2d_obj); // [obj] Object to insert in stream
#endif // !__GNUC__ || _AIX || __INTEL_COMPILER
  
  // Static public members
  static int nst_nbr_get(); // [nbr] Number of instantiated class members
  
  // Public member functions
  
  a2d_cls // [fnc] Constructor
  (const size_t &lrv_nbr_arg, // [nbr] Size of least rapidly varying (outer) dimension 
   const size_t &mrv_nbr_arg) // [nbr] Size of most rapidly varying (inner) dimension
  {
    // Purpose: Default constructor for a2d_cls objects
    lrv_nbr=lrv_nbr_arg; // [nbr] Size of least rapidly varying (outer) dimension 
    mrv_nbr=mrv_nbr_arg; // [nbr] Size of most rapidly varying (inner) dimension
    lmn_nbr=lrv_nbr_arg*mrv_nbr_arg; // [nbr] Total size
    if(lrv_nbr == 0 || mrv_nbr == 0) err_prn("a2d_cls::a2d_cls","a2d_cls has 0 size");
    data=new val_T[lmn_nbr]; // [ptr] Pointer to data buffer
    nst_nbr++; // [nbr] Number of instantiated class members
  } // end a2d_cls() prototype

  /* Old constructor as of 20040531---recherche initializers but otherwise OK
  a2d_cls // [fnc] Constructor
  (const size_t &lrv_nbr_arg, // [nbr] Size of least rapidly varying (outer) dimension 
   const size_t &mrv_nbr_arg) : // [nbr] Size of most rapidly varying (inner) dimension
    lrv_nbr(lrv_nbr_arg), // [nbr] Size of least rapidly varying (outer) dimension 
    mrv_nbr(mrv_nbr_arg), // [nbr] Size of most rapidly varying (inner) dimension
    lmn_nbr(lrv_nbr_arg*mrv_nbr_arg), // [nbr] Total size
    data(new val_T[lmn_nbr]) // [ptr] Pointer to data buffer
  {
    // Purpose: Default constructor for a2d_cls objects
    if(lrv_nbr == 0 || mrv_nbr == 0) err_prn("a2d_cls::a2d_cls","a2d_cls has 0 size");
    nst_nbr++; // [nbr] Number of instantiated class members
  } // end a2d_cls() prototype
  */ 

  ~a2d_cls() // [fnc] Destructor
  {
    // Purpose: Destructor for a2d_cls objects
    assert(data != NULL);
    delete[] data; // [ptr] Pointer to data
    nst_nbr--; // [nbr] Number of instantiated class members
    /* 20060816: PGI pgCC does not implement catch/throw
    try{
      if(data == NULL){
	std::cout << "ERROR: a2d_cls::~a2d_cls: Attempting to free memory at location " << data << std::endl;
#define NULL_EXCEPTION 1
	throw NULL_EXCEPTION;
      } // endif
    }catch(...){ // Catch all exceptions
      err_prn("a2d_cls::~a2d_cls","Attempt to free invalid pointer");
    } // end try
    */ // end PGI pgCC
  } // end a2d_cls<val_T>::~a2d_cls()

  const val_T& // [lmn] Read-only reference to a2d_cls element
  operator() // [fnc] Provide read-only reference to a2d_cls element
    (const size_t &lrv_idx, // [idx] Index for least rapidly varying (outer) dimension
     const size_t &mrv_idx) // [idx] Index for most rapidly varying (inner) dimension
    const{ // Read-only access does not change object
    /* Purpose: Return reference to array element for read-only use
       Method is used when element appears on RHS of assignment (rvalue) 
       Use () rather than [] for subscripting because [] only takes one parameter
       whereas () can take multiple parameters (two are needed for 2-D arrays) */
    return data[index(lrv_idx,mrv_idx)]; // [lmn] Requested element
  } // end operator() prototype

  val_T& // [lmn] Read-write reference to a2d_cls element
  operator() // [fnc] Provide read-write reference to a2d_cls element
    (const size_t &lrv_idx, // [idx] Index for least rapidly varying (outer) dimension
     const size_t &mrv_idx) // [idx] Index for most rapidly varying (inner) dimension
  {
    /* Purpose: Return reference to array element for read-write use
       Method is used when element appears on LHS of assignment (lvalue)
       Use () rather than [] for subscripting because [] only takes one parameter
       whereas () can take multiple parameters (two are needed for 2-D arrays) */
    return data[index(lrv_idx,mrv_idx)]; // [lmn] Requested element
  } // end operator() prototype

  const val_T& // [lmn] Read-only reference to a2d_cls element
  operator() // [fnc] Provide read-only reference to a2d_cls element referenced by 1-D subscript
    (const size_t &lmn_idx) // [idx] Index for one-dimensional array
    const{ // Read-only access does not change object
    /* Purpose: Return reference to array element for read-only use
       Method is used when element appears on RHS of assignment (rvalue) */
    return data[index(lmn_idx)]; // [lmn] Requested element
  } // end operator() prototype

  val_T& // [lmn] Read-write reference to a2d_cls element
  operator() // [fnc] Provide read-write reference to a2d_cls element referenced by 1-D subscript
    (const size_t &lmn_idx) // [idx] Index for one-dimensional array
  {
    /* Purpose: Return reference to array element for read-write use
       Method is used when element appears on LHS of assignment (lvalue) */
    return data[index(lmn_idx)]; // [lmn] Requested element
  } // end operator() prototype

  const val_T& // [lmn] Read-only reference to a2d_cls element
  operator[] // [fnc] Provide read-only reference to a2d_cls element referenced by 1-D subscript
  (const size_t &lmn_idx) // [idx] Index for one-dimensional array
    const{ // Read-only access does not change object
    /* Purpose: Return reference to array element for read-only use
       Method is used when element appears on RHS of assignment (rvalue) */
    return data[index(lmn_idx)]; // [lmn] Requested element
  } // end operator[] prototype

  val_T& // [lmn] Read-write reference to a2d_cls element
  operator[] // [fnc] Provide read-write reference to a2d_cls element referenced by 1-D subscript
  (const size_t &lmn_idx) // [idx] Index for one-dimensional array
  {
    /* Purpose: Return reference to array element for read-write use
       Method is used when element appears on LHS of assignment (lvalue) */
    return data[index(lmn_idx)]; // [lmn] Requested element
  } // end operator[] prototype

  a2d_cls<val_T>& // [obj] 
  operator= // [fnc] Assignment operator
  (const a2d_cls<val_T> &rhs) // [obj] Right-hand-side of a2d_cls assignment
  {
    // Purpose: Copy one a2d_cls object into another
    if(this == &rhs) return *this;

    // Assignment operator copies data, not object (ANSI C++)
    if((lrv_nbr == rhs.lrv_nbr) && (mrv_nbr == rhs.mrv_nbr)){
      // Arrays are same size
      for(size_t idx=0;idx<lmn_nbr;idx++) data[idx]=rhs.data[idx];
    }else{ // endif arrays are same size differ
      wrn_prn("a2d_cls::operator=","Assignment of mis-matched sizes");
      // Free old data buffer
      delete[] data; // [ptr] Pointer to data
      // Assign new sizes to old a2d object
      lrv_nbr=rhs.lrv_nbr; // [nbr] Size of least rapidly varying (outer) dimension 
      mrv_nbr=rhs.mrv_nbr; // [nbr] Size of most rapidly varying (inner) dimension
      lmn_nbr=rhs.lmn_nbr; // [nbr] Total size
      // Allocate new data buffer
      data=new val_T[lmn_nbr]; // [ptr] Pointer to data buffer
      // Fill new data buffer with RHS values
      for(size_t idx=0;idx<lmn_nbr;idx++) data[idx]=rhs.data[idx];
    } // endelse arrays sizes differ
    return *this; // [lmn] Requested element
  } // end operator=()

  a2d_cls<val_T>& // [obj] 
  operator= // [fnc] Assign all elements of array to scalar
  (const val_T &rhs) // [frc] Scalar value
  {
    // Purpose: Assign all elements of array to scalar
    for(size_t idx=0;idx<lmn_nbr;idx++) data[idx]=rhs;
    return *this; // [lmn] Requested element
  } // end a2d_cls<val_T>::operator=()

private:
  // Static private members
  static int nst_nbr; // [nbr] Number of instantiated class members

  // Private members

  size_t lrv_nbr; // [nbr] Size of least rapidly varying (outer) dimension 
  size_t mrv_nbr; // [nbr] Size of most rapidly varying (inner) dimension
  size_t lmn_nbr; // [nbr] Total size
  val_T *data; // [ptr] Pointer to data buffer

  size_t // [idx] Validated 1-D offset index
  index // [fnc] Validate and return 1-D offset
  (const size_t &lmn_idx) // [idx] Index for one-dimensional array
    const{ 
    // Purpose: Validate and return 1-D offset
    if(lmn_idx > lmn_nbr-1 || lmn_idx < 0) err_prn("a2d_cls::index","index out of range");
    return lmn_idx;
  } // end index()

  size_t // [idx] Validated 1-D offset index
  index // [fnc] Convert multidimensional indices into validated 1-D offset
  (const size_t &lrv_idx, // [idx] Index for least rapidly varying (outer) dimension
   const size_t &mrv_idx) // [idx] Index for most rapidly varying (inner) dimension
    const{ 
    // Purpose: Convert multidimensional indices into 1-D offset assuming contiguous storage
    const size_t lmn_idx(lrv_idx*mrv_nbr+mrv_idx); // [idx] 1-D offset assuming contiguous storage
    if(lmn_idx > lmn_nbr-1) err_prn("a2d_cls::index","index out of range");
    return lmn_idx;
  } // end index()

}; // end class a2d_cls

// a2d_cls class static member templates
// Initialize static member data of template classes _after_ class definition DeD01 p. 717
template<class val_T> // [obj] Object type
int a2d_cls<val_T>::nst_nbr(0); // [nbr] Number of instantiated class members
template<class val_T> // [obj] Object type
int a2d_cls<val_T>::nst_nbr_get(){return nst_nbr;} // [nbr] Number of instantiated class members

// Prototype global functions with C++ linkages

// NB: Templated friends must appear after template definition in .hh files!
template<class val_T> // [obj] Object type
std::ostream & // [srm] Reference to output stream for cascading
operator<< // [fnc] Stream insertion operator
// NB: Previous line MUST NOT use <> even though forward declaration does!
(std::ostream &srm_out, // [srm] Output stream
 const a2d_cls<val_T> &a2d_obj) // [obj] Object to insert in stream
{
  /* Purpose: Overloaded stream insertion operator for a2d_cls class objects
     Overloaded stream operators discussed on DeD01 p. 529 */
  srm_out << "Public contents of a2d_obj: " << std::endl;
  srm_out << "a2d_obj.nst_nbr_get() = " << a2d_obj.nst_nbr_get() << std::endl;
  srm_out << "Size of a2d elements = sizeof(val_T) = " << sizeof(val_T) << std::endl;
  srm_out << "Private contents of a2d_obj: " << std::endl;
  srm_out << "lrv_nbr = " << a2d_obj.lrv_nbr << std::endl;
  srm_out << "mrv_nbr = " << a2d_obj.mrv_nbr << std::endl;
  srm_out << "lmn_nbr = " << a2d_obj.lmn_nbr << std::endl;
  for(size_t lrv_idx=0;lrv_idx<a2d_obj.lrv_nbr;lrv_idx++){
    for(size_t mrv_idx=0;mrv_idx<a2d_obj.mrv_nbr;mrv_idx++){
      std::cout << "a2d_obj(" << lrv_idx << "," << mrv_idx << ") = " << a2d_obj(lrv_idx,mrv_idx) << std::endl;
    } // end loop over mrv
  } // end loop over lrv
  return srm_out; // [srm] Reference to output stream for cascading
} // end operator<<()

// Define inline'd functions in header so source is visible to calling files

#endif // A2D_HH
