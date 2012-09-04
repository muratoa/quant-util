#ifndef TWO_D_VECTOR_CH_INCLUDED
#define TWO_D_VECTOR_CH_INCLUDED

/** CPP INCLUDES ***************************************************************/
#include <algorithm>
#include <iterator>

template <typename T>
Util::TwoDVector<T>::TwoDVector(const size_t D1, const size_t D2)
    : rowSize_(0)
    , colSize_(0)
    , pData_(NULL)
{
    resize(D1,D2);
}

template <typename T>
Util::TwoDVector<T>::TwoDVector(const_iterator first,
                                const size_t D1, const size_t D2)
    : rowSize_(0)
    , colSize_(0)
    , pData_(NULL)
{
    resize(first,D1,D2);
}

template <typename T>
Util::TwoDVector<T>::TwoDVector()
    : rowSize_(0)
    , colSize_(0)
    , pData_(NULL)
{
    resize(1,1);
}

template <typename T>
Util::TwoDVector<T>::~TwoDVector()
{
    if (pData_) {
        delete[] pData_;
        pData_ = NULL;
    }
}

template <typename T>
void
Util::TwoDVector<T>::resize(const size_t D1, const size_t D2)
{
    resizeImpl(D1,D2,*this);
}

template <typename T>
void
Util::TwoDVector<T>::resize(const_iterator first,
                            const size_t D1, const size_t D2)
{
    resizeImpl(D1,D2,first);
}

template <typename T>
void
Util::TwoDVector<T>::resizeImpl(const size_t D1, const size_t D2,
                                const Util::TwoDVector<T>& x)
{
    const size_t newSize = D1 * D2;
    const size_t sizeToCopy = newSize <= x.size() ? newSize : x.size();
    
    T* tmp = new T[newSize];
    std::copy(x.begin(),x.begin() + sizeToCopy,tmp);
    
    if (pData_) {
        delete[] pData_;
        pData_ = NULL;
    }
    
    pData_ = tmp;
    rowSize_ = D1;
    colSize_ = D2;
}

template <typename T>
void
Util::TwoDVector<T>::resizeImpl(const size_t D1, const size_t D2,
                                Util::TwoDVector<T>::const_iterator first)
{
    const size_t newSize = D1 * D2;
    
    T* tmp = new T[newSize];
    std::copy(first,first + newSize,tmp);
    
    if (pData_) {
        delete[] pData_;
        pData_ = NULL;
    }
    
    pData_ = tmp;
    rowSize_ = D1;
    colSize_ = D2;
}

template <typename T>
void
Util::TwoDVector<T>::copyFrom(const Util::TwoDVector<T>& x)
{
    resizeImpl(x.rowSize(),x.colSize(),x);
}

template <typename T>
Util::TwoDVector<T>::TwoDVector(const Util::TwoDVector<T>& x)
    : pData_(NULL)
{
    copyFrom(x);
}

template <typename T>
Util::TwoDVector<T>&
Util::TwoDVector<T>::operator=(const Util::TwoDVector<T>& x)
{
    if (this != &x) {
        copyFrom(x);
    }
    
    return *this;
}

template <typename T>
typename Util::TwoDVector<T>::const_iterator
Util::TwoDVector<T>::begin() const
{
    return pData_;
}

template <typename T>
typename Util::TwoDVector<T>::iterator
Util::TwoDVector<T>::begin()
{
    return pData_;
}

template <typename T>
typename Util::TwoDVector<T>::const_iterator
Util::TwoDVector<T>::end() const
{
    return begin() + size();
}

template <typename T>
typename Util::TwoDVector<T>::iterator
Util::TwoDVector<T>::end()
{
    return begin() + size();
}

template <typename T>
T&
Util::TwoDVector<T>::operator()(const size_t i)
{
    return pData_[i];
}

template <typename T>
const T&
Util::TwoDVector<T>::operator()(const size_t i) const
{
    return pData_[i];
}

template <typename T>
T&
Util::TwoDVector<T>::operator()(const size_t i, const size_t j)
{
    return pData_[getIndex(i,j)];
}

template <typename T>
const T&
Util::TwoDVector<T>::operator()(const size_t i, const size_t j) const
{
    return pData_[getIndex(i,j)];
}

template <typename T>
size_t
Util::TwoDVector<T>::getIndex(const size_t i, const size_t j) const
{
    return colSize()*i + j;
}

template <typename T>
void
Util::TwoDVector<T>::dimSwap()
{
    const size_t tmp = rowSize();
    rowSize_ = colSize();
    colSize_ = tmp;
}

template <typename T>
void
Util::TwoDVector<T>::print(std::ostream& os) const
{
    const TwoDVector<T>& x = *this;
    std::ostream_iterator<T> oit(os," ");
    
    for (size_t i = 0; i < x.rowSize(); ++i) {
        
        std::copy(x.begin() + i*x.colSize(),
                  x.begin() + (i+1)*x.colSize(), oit);
        os << std::endl;
    }
}

#endif // TWO_D_VECTOR_CH included
