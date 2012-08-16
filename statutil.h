#ifndef STAT_UTIL_H_INCLUDED
#define STAT_UTIL_H_INCLUDED

/** CPP INCLUDES ***************************************************************/
#include <vector>

/** APP INCLUDES ***************************************************************/

namespace Util {

/** Utility Fns ****************************************************************/
    
template <typename T>
void transpose(const T* a, int m, int n, T* b, int nb=64);

template <typename RandomAccessIterator>
void transpose(RandomAccessIterator first, RandomAccessIterator last, int m);
    
/** Descriptive Fns ************************************************************/

template <typename T, typename ForwardIterator>
T mean(ForwardIterator first, ForwardIterator last);

template <typename T, typename ForwardIterator>
T stddev(ForwardIterator first, ForwardIterator last);

template <typename T, typename ForwardIterator>
T sumOfSquares(ForwardIterator first, ForwardIterator last);

template <typename T, typename ForwardIterator>
T meanAbsDiff(ForwardIterator first, ForwardIterator last);    

/** Matrix Apply Fns ***********************************************************/

template <typename RetT, typename T, class F>
std::vector<RetT> rowApply(const T* a, int m, int n, F f);

template <typename RetT, typename T, class F>
std::vector<RetT> colApply(const T* a, int m, int n, F f);

template <typename RetT, typename T>
std::vector<RetT> rowMeans(const T* a, int m, int n);    

template <typename RetT, typename T>
std::vector<RetT> colMeans(const T* a, int m, int n);

template <typename RetT, typename T>
std::vector<RetT> rowStdDevs(const T* a, int m, int n);

template <typename RetT, typename T>
std::vector<RetT> colStdDevs(const T* a, int m, int n);

template <typename RetT, typename T, class F>
void rowApply(const T* a, int m, int n, F f, std::vector<RetT>& res);

template <typename RetT, typename T, class F>
void colApply(const T* a, int m, int n, F f, std::vector<RetT>& res);    

template <typename RetT, typename T>
void rowMeans(const T* a, int m, int n, std::vector<RetT>& res);    

template <typename RetT, typename T>
void colMeans(const T* a, int m, int n, std::vector<RetT>& res);

template <typename RetT, typename T>
void rowStdDevs(const T* a, int m, int n, std::vector<RetT>& res);

template <typename RetT, typename T>
void colStdDevs(const T* a, int m, int n, std::vector<RetT>& res);    

/** Scale Fns ******************************************************************/

template <class InputIterator, class OutputIterator>
OutputIterator
standardize(InputIterator first, InputIterator last, OutputIterator result,
            bool center=true, bool scale=true);

template <class InputIterator, class OutputIterator, typename T>
OutputIterator
standardize(InputIterator first, InputIterator last, OutputIterator result,
            T center, T scale);

template <typename T>
void
standardizeRows(T* a, int m, int n, bool center=true, bool scale=true);

template <typename T>
void
standardizeCols(T* a, int m, int n, bool center=true, bool scale=true);
    
template <typename T, typename T2>
void
standardizeRows(T* a, int m, int n,
                const std::vector<T2>& means,
                const std::vector<T2>& stddevs);    
    
template <typename T, typename T2>
void
standardizeCols(T* a, int m, int n,
                const std::vector<T2>& means,
                const std::vector<T2>& stddevs);
    
} // namespace

#include "statutil.ch"

/*(+doc.)
  .if
  .SH DESCRIPTION

  .SS

  .SH EaAMPLES

  .SH FILES
  .nf
  util/statutil.h - class declaration
  .if

  .SH SEE ALSO

  .SH AUTHOR
  Murat Ahmed
  (-doc.)*/


#endif // STAT_UTIL_H included
