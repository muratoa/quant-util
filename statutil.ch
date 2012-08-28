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
                RandomAccessIterator last, const int m)
{
    const int sz = (last - first - 1);
    const int n = (last - first) / m;
    
    std::vector<bool> visited(last - first);

    RandomAccessIterator cycle = first;
    while (++cycle != last) {

        if (visited[cycle - first])
            continue;

        int idx = cycle - first;

        do { 
            idx = (idx == sz) ? sz : (n*idx) % sz;
            std::swap(*(first + idx), *cycle);
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

template <typename RetT, typename T, class F>
std::vector<RetT>
Util::rowApply(const T* a, int m, int n, F f)
{
    std::vector<RetT> res(m);
    for (int i = 0; i < m; i++) {
        res[i] = f(a[i*n], a[(i+1)*n - 1]);
    }

    return res;
}

template <typename RetT, typename T, class F>
void
Util::rowApply(const T* a, int m, int n, F f, std::vector<RetT>& res)
{
    res.resize(m);
    for (int i = 0; i < m; i++) {
        res[i] = f(a + i*n, a + (i+1)*n);
    }
}

template <typename RetT, typename T, class F>
std::vector<RetT>
Util::colApply(const T* a, const int m, const int n, F f)
{
    T b[m*n];
    transpose(a,m,n,b);
    
    std::vector<RetT> res(n);
    for (int i = 0; i < n; i++) {
        res[i] = f(b + i*m, b + (i+1)*m);
    }

    return res;
}

template <typename RetT, typename T, class F>
void
Util::colApply(const T* a, const int m, const int n, F f, std::vector<RetT>& res)
{
    T b[m*n];
    transpose(a,m,n,b);
    
    res.resize(n);
    for (int i = 0; i < n; i++) {
        res[i] = f(b + i*m, b + (i+1)*m);
    }
}

template <typename RetT, typename T>
std::vector<RetT>
Util::rowMeans(const T* a, const int m, const int n)
{
    MeanCalc<RetT,const T*> mean;
    
    return rowApply(a,m,n,mean);
}

template <typename RetT, typename T>
void
Util::rowMeans(const T* a, const int m, const int n, std::vector<RetT>& res)
{
    MeanCalc<RetT,const T*> mean;
    rowApply(a,m,n,mean,res);
}

template <typename RetT, typename T>
std::vector<RetT>
Util::colMeans(const T* a, const int m, const int n)
{
    MeanCalc<RetT,const T*> mean;

    return colApply(a,m,n,mean);
}

template <typename RetT, typename T>
void
Util::colMeans(const T* a, const int m, const int n, std::vector<RetT>& res)
{
    MeanCalc<RetT,const T*> mean;
    colApply(a,m,n,mean,res);
}

template <typename RetT, typename T>
std::vector<RetT>
Util::rowStdDevs(const T* a, const int m, const int n)
{
    StdDevCalc<RetT,const T*> stddev;
    
    return rowApply(a,m,n,stddev);
}

template <typename RetT, typename T>
void
Util::rowStdDevs(const T* a, const int m, const int n, std::vector<RetT>& res)
{
    StdDevCalc<RetT,const T*> stddev;
    rowApply(a,m,n,stddev,res);
}

template <typename RetT, typename T>
std::vector<RetT>
Util::colStdDevs(const T* a, const int m, const int n)
{
    StdDevCalc<RetT,const T*> stddev;

    return colApply(a,m,n,stddev);
}

template <typename RetT, typename T>
void
Util::colStdDevs(const T* a, const int m, const int n, std::vector<RetT>& res)
{
    StdDevCalc<RetT,const T*> stddev;
    colApply(a,m,n,stddev,res);
}

template <class InputIterator, class OutputIterator>
OutputIterator
Util::standardize(InputIterator first, InputIterator last, OutputIterator result,
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
Util::standardize(InputIterator first, InputIterator last, OutputIterator result,
                  T center, T scale)
{
    typedef typename std::iterator_traits<OutputIterator>::value_type Tp;
    Detail::assert<Detail::valid_numeric_type<Tp>::value> invalid_type;
    UNUSED(invalid_type);

    scale = (scale > T(0)) ? scale : T(1);
    Detail::Standardize<T> op(center,scale);
    
    return std::transform(first,last,result,op);
}

template <typename T>
void
Util::standardizeRows(T* a, const int m, const int n,
                      bool center=true, bool scale=true)
{
    Detail::assert<Detail::valid_numeric_type<T>::value> invalid_type;
    UNUSED(invalid_type);
    
    for (int i = 0; i < m; i++) {
        standardize(a + i*n,
                    a + (i+1)*n,
                    a + i*n,
                    center, scale);
    }
}

template <typename T>
void
Util::standardizeCols(T* a, const int m, const int n,
                      bool center=true, bool scale=true)
{
    Detail::assert<Detail::valid_numeric_type<T>::value> invalid_type;
    UNUSED(invalid_type);

    T b[m*n];
    transpose(a,m,n,b);
    
    for (int i = 0; i < n; i++) {
        standardize(b + i*m,
                    b + (i+1)*m,
                    b + i*m,
                    center,scale);
    }
    
    transpose(b,n,m,a);
}

template <typename T, typename RetT>
void
Util::standardizeRows(T* a, const int m, const int n,
                      const std::vector<RetT>& means,
                      const std::vector<RetT>& stddevs)
{
    Detail::assert<Detail::valid_numeric_type<RetT>::value> invalid_type;
    UNUSED(invalid_type);

    for (int i = 0; i < m; i++) {
        standardize(a + i*n,
                    a + (i+1)*n,
                    a + i*n,
                    means[i],stddevs[i]);
    }
}

template <typename T, typename RetT>
void
Util::standardizeCols(T* a, const int m, const int n,
                      const std::vector<RetT>& means,
                      const std::vector<RetT>& stddevs)
{
    Detail::assert<Detail::valid_numeric_type<RetT>::value> invalid_type;
    UNUSED(invalid_type);

    T b[m*n];
    transpose(a,m,n,b);
    
    for (int i = 0; i < n; i++) {
        standardize(b + i*m,
                    b + (i+1)*m,
                    b + i*m,
                    means[i],stddevs[i]);
    }
    
    transpose(b,n,m,a);
}

#endif // STAT_UTIL_CH included
