#include "twodvector.h"

/** CPP INCLUDES ***************************************************************/
#include <algorithm>

template <typename T>
Util::TwoDVector<T>::TwoDVector(const size_t D1, const size_t D2)
    : nrow_(0)
    , ncol_(0)
    , pData_(NULL)
{
    resize(D1,D2);
}

template <typename T>
Util::TwoDVector<T>::TwoDVector()
    : nrow_(0)
    , ncol_(0)
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
    nrow_ = D1;
    ncol_ = D2;
}

template <typename T>
void
Util::TwoDVector<T>::copyFrom(const Util::TwoDVector<T>& x)
{
    resizeImpl(x.nrow(),x.ncol(),x);
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
    return ncol()*i + j;
}

template <typename T>
void
Util::TwoDVector<T>::print(std::ostream& os) const
{
    const TwoDVector<T>& x = *this;
    
    for (size_t i = 0;  i < nrow(); ++i) {
        for (size_t j = 0; j < ncol(); ++j) {
            os << x(i,j) << ' ';
        }
        os << std::endl;
    }
}
