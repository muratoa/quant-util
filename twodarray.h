#ifndef TWO_D_ARRAY_H_INCLUDED
#define TWO_D_ARRAY_H_INCLUDED

/** CPP INCLUDES ***************************************************************/
#include <array>
#include <ostream>

/** APP INCLUDES ***************************************************************/

namespace Util {

template <typename T, size_t D1, size_t D2>
class TwoDArray {
public:

    typedef T value_type;
    typedef value_type* iterator;
    typedef const value_type* const_iterator;
    
    TwoDArray();
    ~TwoDArray();
    
    TwoDArray(const TwoDArray<T,D1,D2>& x);
    TwoDArray& operator=(const TwoDArray<T,D1,D2>& x); 

    size_t nrow() const {return nrow_;}
    size_t ncol() const {return ncol_;}
    size_t size() const {return nrow() * ncol();}
    
    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;

    operator value_type*() {return begin();}
    operator const value_type*() const {return begin();}
    
    inline T& operator()(size_t i);
    inline T& operator ()(size_t i, size_t j);
    inline const T& operator()(size_t i) const;
    inline const T& operator()(size_t i, size_t j) const;

    void dimSwap(); // call after tranpose.
    
    friend std::ostream& operator<<(std::ostream& os,
                                    const TwoDArray<T,D1,D2>& x) {
        x.print(os);
        return os;
    }
    
private:
    size_t getIndex(size_t i, size_t j) const;
    void copyFrom(const TwoDArray<T,D1,D2>& x);
    void print(std::ostream& os) const;
    
private:
    typedef std::array<std::array<T,D2>,D1> data_type;
    
    size_t nrow_;
    size_t ncol_;
    data_type* pData_;
};

} // namespace

#include "twodarray.ch"

/*(+doc.)
  .if
  .SH DESCRIPTION

  .SS

  .SH EXAMPLES

  .SH FILES
  .nf
  util/twodarray.h - class declaration
  .if

  .SH SEE ALSO

  .SH AUTHOR
  Murat Ahmed
  (-doc.)*/

#endif // TWO_D_ARRAY_H included
