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
void transpose(RandomAccessIterator first, RandomAccessIterator last,
               size_t colSize);
    
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

template <typename T, typename ForwardIterator, class F>
std::vector<T> rowApply(ForwardIterator first,
                        ForwardIterator last,
                        size_t rowSize, F f);

template <typename T, typename ForwardIterator, class F>
std::vector<T> colApply(ForwardIterator first,
                        ForwardIterator last,
                        size_t colSize, F f);

template <typename T, typename ForwardIterator>
std::vector<T> rowMeans(ForwardIterator first,
                        ForwardIterator last,
                        size_t rowSize);

template <typename T, typename ForwardIterator>
std::vector<T> colMeans(ForwardIterator first,
                        ForwardIterator last,
                        size_t colSize);

template <typename T, typename ForwardIterator>
std::vector<T> rowStdDevs(ForwardIterator first,
                          ForwardIterator last,
                          size_t rowSize);

template <typename T, typename ForwardIterator>
std::vector<T> colStdDevs(ForwardIterator first,
                          ForwardIterator last,
                          size_t colSize);

template <typename T, typename ForwardIterator, class F>
void rowApply(ForwardIterator first,
              ForwardIterator last,
              size_t rowSize,
              F f, std::vector<T>& res);

template <typename T, typename ForwardIterator, class F>
void colApply(ForwardIterator first,
              ForwardIterator last,
              size_t colSize,
              F f, std::vector<T>& res);

template <typename T, typename ForwardIterator>
void rowMeans(ForwardIterator first,
              ForwardIterator last,
              size_t rowSize,
              std::vector<T>& res);

template <typename T, typename ForwardIterator>
void colMeans(ForwardIterator first,
              ForwardIterator last,
              size_t colSize,
              std::vector<T>& res);

template <typename T, typename ForwardIterator>
void rowStdDevs(ForwardIterator first,
                ForwardIterator last,
                size_t rowSize,
                std::vector<T>& res);
    
template <typename T, typename ForwardIterator>
void colStdDevs(ForwardIterator first,
                ForwardIterator last,
                size_t colSize,
                std::vector<T>& res);

/** Scale Fns ******************************************************************/

template <class InputIterator, class OutputIterator>
OutputIterator
standardize(InputIterator first,
            InputIterator last,
            OutputIterator result,
            bool center=true, bool scale=true);

template <class InputIterator, class OutputIterator, typename T>
OutputIterator
standardize(InputIterator first,
            InputIterator last,
            OutputIterator result,
            T center, T scale);

template <class InputIterator>
void
standardizeRows(InputIterator first,
                InputIterator last,
                size_t rowSize,
                bool center=true, bool scale=true);

template <class InputIterator>
void
standardizeCols(InputIterator first,
                InputIterator last,
                size_t colSize,
                bool center=true, bool scale=true);
    
template <class InputIterator, typename T>
void
standardizeRows(InputIterator first,
                InputIterator last,
                size_t rowSize,
                const std::vector<T>& means,
                const std::vector<T>& stddevs);
    
template <class InputIterator, typename T>
void
standardizeCols(InputIterator first,
                InputIterator last,
                size_t rowSize,
                const std::vector<T>& means,
                const std::vector<T>& stddevs);
    
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
