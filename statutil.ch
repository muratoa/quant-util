#ifndef STAT_UTIL_CH_INCLUDED
#define STAT_UTIL_CH_INCLUDED

/** CPP INCLUDES ***************************************************************/
#include <algorithm>
#include <cmath>
#include <iterator>

/** APP INCLUDES ***************************************************************/

#define UNUSED(x) ((void)(x))

namespace Util {

/** Impl Classes ***************************************************************/

template <typename T, typename ForwardIterator>
class VarCalcImpl {
public:
    void calc(ForwardIterator first, ForwardIterator last, T& m, T& m2);
};

template <typename T, typename ForwardIterator>
class MeanCalc : VarCalcImpl<T,ForwardIterator> {
public:
    T operator()(ForwardIterator first, ForwardIterator last);
};

template <typename T, typename ForwardIterator>
class StdDevCalc : VarCalcImpl<T,ForwardIterator> {
public:
    T operator()(ForwardIterator first, ForwardIterator last);
};

template <typename T, typename ForwardIterator>
class ScaleCalc : VarCalcImpl<T,ForwardIterator> {
public:
    T operator()(ForwardIterator first, ForwardIterator last, T& mu, T& sd);
};

namespace Detail {

template <bool b> struct assert;
template <> struct assert<true>{};

template <typename T>
struct valid_numeric_type {
    static const bool value = false;
};

template <>
struct valid_numeric_type<float> {
    static const bool value = true;
};

template <>
struct valid_numeric_type<double> {
    static const bool value = true;
};

template <typename T>
int8_t
sgn(T x) {
    return (x > T(0)) - (x < T(0));
}

template <class InputIterator, typename T>
void
scaleImpl(InputIterator first, InputIterator last, T& center, T& scale)
{
    ScaleCalc<T,InputIterator> scalecalc;
    scalecalc(first,last,center,scale);
}

template <typename T>
class Standardize {
public:
    Standardize(const T mu, const T sd) : mu_(mu), sd_(sd) {};
    T operator()(const T x) {return (x - mu_) / sd_;}
    
private:
    const T mu_;
    const T sd_;
};

} // namespace

} // namespace

template <typename T>
void
Util::transpose(const T* __restrict__ a, const int m, const int n,
                T* __restrict__ b, const int nb)
{
    for (int i = 0; i < m; i += nb) {
        for (int j = 0; j < n; j += nb) {
            
            for (int k = i; (k < i + nb) && (k < m); ++k) {
                for (int l = j; (l < j + nb) && (l < n); ++l) {
                    
                    b[k + l*m] = a[k*n + l];
                }
            }
        }
    }
}

template <typename RandomAccessIterator>
void
Util::transpose(RandomAccessIterator first,
                RandomAccessIterator last,
                const size_t colSize)
{
    const size_t sz = (last - first - 1);
    const size_t rowSize = (last - first) / colSize;
    
    std::vector<bool> visited(last - first);
    
    RandomAccessIterator cycle = first;
    while (++cycle != last) {
        
        if (visited[cycle - first])
            continue;
        
        int idx = cycle - first;
        
        do {
            idx = (idx == sz) ? sz : (rowSize*idx) % sz;
            std::swap(*(first + idx),*cycle);
            visited[idx] = true;
            
        } while ((first + idx) != cycle);
    }
}

template <typename T, typename ForwardIterator>
void
Util::VarCalcImpl<T,ForwardIterator>::calc(ForwardIterator first,
                                           ForwardIterator last,
                                           T& m, T& m2)
{
    typedef typename std::iterator_traits<ForwardIterator>::value_type Tp;
    m = T(0);
    m2 = T(0);
    T delta = T(0);

    int n = 1;
    while (first != last) {

        Tp val = *first;
        delta = val - m;
        m += delta / n;
        m2 += delta * (val - m);

        ++first;
        ++n;
    }
}

template <typename T, typename ForwardIterator>
T
Util::MeanCalc<T,ForwardIterator>::operator()(ForwardIterator first,
                                              ForwardIterator last)
{
    T m;
    T m2;
    calc(first,last,m,m2);

    return m;
}

template <typename T, typename ForwardIterator>
T
Util::StdDevCalc<T,ForwardIterator>::operator()(ForwardIterator first,
                                                ForwardIterator last)
{
    T m;
    T m2;
    calc(first,last,m,m2);
    const size_t n = std::distance(first,last);
    
    return (n > 1) ? std::sqrt(m2 / (n - 1)) : T(-1);
}

template <typename T, typename ForwardIterator>
T
Util::ScaleCalc<T,ForwardIterator>::operator()(ForwardIterator first,
                                               ForwardIterator last,
                                               T& mu, T& sd)
{
    T m;
    T m2;
    calc(first,last,m,m2);
    const size_t n = std::distance(first,last);
    
    mu = m;
    sd = (n > 1) ? std::sqrt(m2 / (n - 1)) : T(-1);
}

template <typename T, typename ForwardIterator>
T
Util::mean(ForwardIterator first, ForwardIterator last)
{
    MeanCalc<T,ForwardIterator> mean;

    return mean(first,last);
}

template <typename T, typename ForwardIterator>
T
Util::stddev(ForwardIterator first, ForwardIterator last)
{
    StdDevCalc<T,ForwardIterator> stddev;

    return stddev(first,last);
}

template <typename T, typename ForwardIterator>
T
Util::sumOfSquares(ForwardIterator first, ForwardIterator last)
{
    T ss = T(0);
    for (; first!=last; ++first)
        ss += (*first) * (*first);
    
    return ss;
}

template <typename T, typename ForwardIterator>
T
Util::meanAbsDiff(ForwardIterator first, ForwardIterator last)
{
    const size_t n = std::distance(first,last);
    if (n <= 1)
        return T(-1);

    T mad = T(0);
    while (first != last) {

        ForwardIterator curr = first;
        ++curr;
        while (curr != last) {

            mad += std::abs(*first - *curr);
            ++curr;
        }
        ++first;
    }
    
    return 2 * mad / (n*(n-1));
}

template <typename T, typename ForwardIterator, class F>
std::vector<T>
Util::rowApply(ForwardIterator first,
               ForwardIterator last,
               const size_t rowSize, F f)
{
    const size_t colSize = (last - first) / rowSize;   
    std::vector<T> res(rowSize);
    
    for (int i = 0; i < rowSize; i++) {
        res[i] = f(first + i*colSize,
                   first + (i+1)*colSize);
    }
    
    return res;
}

template <typename T, typename ForwardIterator, class F>
void
Util::rowApply(ForwardIterator first,
               ForwardIterator last,
               const size_t rowSize, F f,
               std::vector<T>& res)
{
    const size_t colSize = (last - first) / rowSize;   
    res.resize(rowSize);
    
    for (int i = 0; i < rowSize; i++) {
        res[i] = f(first + i*colSize,
                   first + (i+1)*colSize);
    }
}

template <typename T, typename ForwardIterator, class F>
std::vector<T>
Util::colApply(ForwardIterator first,
               ForwardIterator last,
               const size_t colSize, F f)
{
    std::vector<T> res;
    Util::colApply(first,last,colSize,f,res);
    
    return res;
}

template <typename T, typename ForwardIterator, class F>
void
Util::colApply(ForwardIterator first,
               ForwardIterator last,
               const size_t colSize, F f,
               std::vector<T>& res)
{
    typedef typename std::iterator_traits<ForwardIterator>::value_type Tp;
    const size_t rowSize = (last - first) / colSize;
    
    Util::transpose(const_cast<Tp*>(first),
                    const_cast<Tp*>(last),colSize);
    
    Util::rowApply(first,last,colSize,f,res);
    
    Util::transpose(const_cast<Tp*>(first),
                    const_cast<Tp*>(last),rowSize);
}

template <typename T, typename ForwardIterator>
std::vector<T>
Util::rowMeans(ForwardIterator first,
               ForwardIterator last,
               const size_t rowSize)  
{
    MeanCalc<T,ForwardIterator> mean;
    
    return rowApply<T>(first,last,rowSize,mean);
}

template <typename T, typename ForwardIterator>
void
Util::rowMeans(ForwardIterator first,
               ForwardIterator last,
               const size_t rowSize,
               std::vector<T>& res)  
{
    MeanCalc<T,ForwardIterator> mean;
    rowApply(first,last,rowSize,mean,res);
}

template <typename T, typename ForwardIterator>
std::vector<T>
Util::colMeans(ForwardIterator first,
               ForwardIterator last,
               const size_t colSize)
{
    MeanCalc<T,ForwardIterator> mean;
    
    return colApply<T>(first,last,colSize,mean);
}

template <typename T, typename ForwardIterator>
void
Util::colMeans(ForwardIterator first,
               ForwardIterator last,
               const size_t colSize,
               std::vector<T>& res)
{
    MeanCalc<T,ForwardIterator> mean;
    colApply(first,last,colSize,mean,res);
}

template <typename T, typename ForwardIterator, class F>
std::vector<T>
Util::rowStdDevs(ForwardIterator first,
                 ForwardIterator last,
                 const size_t rowSize)
{
    StdDevCalc<T,ForwardIterator> stddev;
    
    return rowApply(first,last,rowSize,stddev);
}

template <typename T, typename ForwardIterator>
void
Util::rowStdDevs(ForwardIterator first,
                 ForwardIterator last,
                 const size_t rowSize,
                 std::vector<T>& res)
{
    StdDevCalc<T,ForwardIterator> stddev;
    rowApply(first,last,rowSize,stddev,res);
}

template <typename T, typename ForwardIterator>
std::vector<T>
Util::colStdDevs(ForwardIterator first,
                 ForwardIterator last,
                 const size_t colSize)
{
    StdDevCalc<T,ForwardIterator> stddev;
    
    return colApply(first,last,colSize,stddev);
}

template <typename T, typename ForwardIterator>
void
Util::colStdDevs(ForwardIterator first,
                 ForwardIterator last,
                 const size_t colSize,
                 std::vector<T>& res)
{
    StdDevCalc<T,ForwardIterator> stddev;
    colApply(first,last,colSize,stddev,res);
}

template <class InputIterator, class OutputIterator>
OutputIterator
Util::standardize(InputIterator first,
                  InputIterator last,
                  OutputIterator result,
                  bool center=true, bool scale=true)
{
    typedef typename std::iterator_traits<OutputIterator>::value_type T;
    Detail::assert<Detail::valid_numeric_type<T>::value> invalid_type;
    UNUSED(invalid_type);

    T mu;
    T sd;
    Detail::scaleImpl(first,last,mu,sd);

    mu = center ? mu : T(0);
    sd = scale ? (sd > T(0) ? sd : T(1)) : T(1);
    Detail::Standardize<T> op(mu,sd);

    return std::transform(first,last,result,op);
}

template <class InputIterator, class OutputIterator, typename T>
OutputIterator
Util::standardize(InputIterator first,
                  InputIterator last,
                  OutputIterator result,
                  const T center, T scale)
{
    typedef typename std::iterator_traits<OutputIterator>::value_type Tp;
    Detail::assert<Detail::valid_numeric_type<Tp>::value> invalid_type;
    UNUSED(invalid_type);

    scale = (scale > T(0)) ? scale : T(1);
    Detail::Standardize<T> op(center,scale);
    
    return std::transform(first,last,result,op);
}

template <class InputIterator>
void
Util::standardizeRows(InputIterator first,
                      InputIterator last,
                      const size_t rowSize,
                      bool center=true, bool scale=true)
{
    typedef typename std::iterator_traits<InputIterator>::value_type T;
    Detail::assert<Detail::valid_numeric_type<T>::value> invalid_type;
    UNUSED(invalid_type);
    
    const size_t colSize = (last - first) / rowSize;
    
    for (int i = 0; i < rowSize; i++) {
        standardize(first + i*colSize,
                    first + (i+1)*colSize,
                    first + i*colSize,
                    center,scale);
    }
}

template <class InputIterator>
void
Util::standardizeCols(InputIterator first,
                      InputIterator last,
                      const size_t colSize,
                      bool center=true, bool scale=true)
{
    const size_t rowSize = (last - first) / colSize;
    
    Util::transpose(first,last,colSize);
    Util::standardizeRows(first,last,colSize,center,scale);
    Util::transpose(first,last,rowSize);
}

template <class InputIterator, typename T>
void
Util::standardizeRows(InputIterator first,
                      InputIterator last,
                      const size_t rowSize,
                      const std::vector<T>& means,
                      const std::vector<T>& stddevs)
{
    typedef typename std::iterator_traits<InputIterator>::value_type Tp;
    Detail::assert<Detail::valid_numeric_type<Tp>::value> invalid_type;
    UNUSED(invalid_type);
    
    const size_t colSize = (last - first) / rowSize;
    
    for (int i = 0; i < rowSize; i++) {
        standardize(first + i*colSize,
                    first + (i+1)*colSize,
                    first + i*colSize,
                    means[i],stddevs[i]);
    }
}

template <class InputIterator, typename T>
void
Util::standardizeCols(InputIterator first,
                      InputIterator last,
                      const size_t colSize,
                      const std::vector<T>& means,
                      const std::vector<T>& stddevs)
{
    const size_t rowSize = (last - first) / colSize;
    
    Util::transpose(first,last,colSize);
    Util::standardizeRows(first,last,colSize,means,stddevs);
    Util::transpose(first,last,rowSize);
}

#endif // STAT_UTIL_CH included
