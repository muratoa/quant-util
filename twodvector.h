#ifndef TWO_D_VECTOR_H_INCLUDED
#define TWO_D_VECTOR_H_INCLUDED

/** CPP INCLUDES ***************************************************************/
#include <ostream>

/** APP INCLUDES ***************************************************************/

namespace Util {

template <typename T>
class TwoDVector {
public:

    typedef T value_type;
    typedef value_type* iterator;
    typedef const value_type* const_iterator;
    
    TwoDVector(size_t D1, size_t D2);
    TwoDVector(const_iterator first, size_t D1, size_t D2);
    TwoDVector();
    ~TwoDVector();
    
    TwoDVector(const TwoDVector<T>& x);
    TwoDVector& operator=(const TwoDVector<T>& x); 

    void resize(size_t D1, size_t D2);
    void resize(const_iterator first, size_t D1, size_t D2);

    bool empty() const {return size() == 0;}
    
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
    inline T& operator()(size_t i, size_t j);
    inline const T& operator()(size_t i) const;
    inline const T& operator()(size_t i, size_t j) const;

    void dimSwap(); // call after tranpose.
    
    friend std::ostream& operator<<(std::ostream& os,
                                    const TwoDVector<T>& x) {
        x.print(os);
        return os;
    }
    
private:
    size_t getIndex(size_t i, size_t j) const;
    void resizeImpl(size_t D1, size_t D2, const TwoDVector<T>& x);
    void resizeImpl(size_t D1, size_t D2, const_iterator first);
    void copyFrom(const TwoDVector<T>& x);
    void print(std::ostream& os) const;
    
private:
    size_t nrow_;
    size_t ncol_;
    value_type* pData_;
};
    
} // namespace

#include "twodvector.ch"

/*(+doc.)
  .if
  .SH DESCRIPTION

  .SS

  .SH EXAMPLES

  .SH FILES
  .nf
  util/twodvector.h - class declaration
  .if

  .SH SEE ALSO

  .SH AUTHOR
  Murat Ahmed
  (-doc.)*/

#endif // TWO_D_VECTOR_H included
