#include "twodarray.h"

/** CPP INCLUDES ***************************************************************/
#include <algorithm>

template <typename T, size_t D1, size_t D2>
Util::TwoDArray<T,D1,D2>::TwoDArray()
    : nrow_(D1)
    , ncol_(D2)
    , pData_(new Util::TwoDArray<T,D1,D2>::data_type())
{
}

template <typename T, size_t D1, size_t D2>
Util::TwoDArray<T,D1,D2>::~TwoDArray()
{
    if (pData_) {
        delete pData_;
        pData_ = NULL;
    }
}

template <typename T, size_t D1, size_t D2>
void
Util::TwoDArray<T,D1,D2>::copyFrom(const Util::TwoDArray<T,D1,D2>& x)
{
    std::copy(x.begin(),x.end(),begin());
}

template <typename T, size_t D1, size_t D2>
Util::TwoDArray<T,D1,D2>::TwoDArray(const Util::TwoDArray<T,D1,D2>& x)
    : pData_(new Util::TwoDArray<T,D1,D2>::data_type())
{
    copyFrom(x);
}

template <typename T, size_t D1, size_t D2>
Util::TwoDArray<T,D1,D2>&
Util::TwoDArray<T,D1,D2>::operator=(const Util::TwoDArray<T,D1,D2>& x)
{
    if (this != &x) {
        copyFrom(x);
    }
    
    return *this;
}

template <typename T, size_t D1, size_t D2>
typename Util::TwoDArray<T,D1,D2>::const_iterator
Util::TwoDArray<T,D1,D2>::begin() const
{
    return &((*pData_)[0][0]);
}

template <typename T, size_t D1, size_t D2>
typename Util::TwoDArray<T,D1,D2>::iterator
Util::TwoDArray<T,D1,D2>::begin()
{
    return &((*pData_)[0][0]);
}

template <typename T, size_t D1, size_t D2>
typename Util::TwoDArray<T,D1,D2>::const_iterator
Util::TwoDArray<T,D1,D2>::end() const
{
    return begin() + size();
}

template <typename T, size_t D1, size_t D2>
typename Util::TwoDArray<T,D1,D2>::iterator
Util::TwoDArray<T,D1,D2>::end()
{
    return begin() + size();
}

template <typename T, size_t D1, size_t D2>
T&
Util::TwoDArray<T,D1,D2>::operator()(const size_t i)
{
    return *(begin() + i);
}

template <typename T, size_t D1, size_t D2>
const T&
Util::TwoDArray<T,D1,D2>::operator()(const size_t i) const
{
    return *(begin() + i);
}

template <typename T, size_t D1, size_t D2>
T&
Util::TwoDArray<T,D1,D2>::operator()(const size_t i, const size_t j)
{
    return *(begin() + getIndex(i,j));
}

template <typename T, size_t D1, size_t D2>
const T&
Util::TwoDArray<T,D1,D2>::operator()(const size_t i, const size_t j) const
{
    return *(begin() + getIndex(i,j));
}

template <typename T, size_t D1, size_t D2>
size_t
Util::TwoDArray<T,D1,D2>::getIndex(const size_t i, const size_t j) const
{
    return ncol()*i + j;
}

template <typename T, size_t D1, size_t D2>
void
Util::TwoDArray<T,D1,D2>::print(std::ostream& os) const
{
    const TwoDArray<T,D1,D2>& x = *this;
    
    for (size_t i = 0;  i < nrow(); ++i) {
        for (size_t j = 0; j < ncol(); ++j) {
            os << x(i,j) << ' ';
        }
        os << std::endl;
    }
}
