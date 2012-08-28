#ifndef MULT_REGRESS_H_INCLUDED
#define MULT_REGRESS_H_INCLUDED

/** CPP INCLUDES ***************************************************************/
#include <ostream>
#include <vector>

/** APP INCLUDES ***************************************************************/

namespace Util {

template <class MatrixT>    
class MultRegress {
public:

    MultRegress(const std::vector<double>& y,
                MatrixT& X,
                const bool calcInferenceStats=false,
                const bool calcInputStats=false);

    MultRegress(const std::vector<double>& y,
                const MatrixT& X,
                const bool calcInferenceStats=false,
                const bool calcInputStats=false);

    MultRegress(const std::vector<double>& y,
                MatrixT& X,
                const std::vector<std::string>& varNames,
                const bool calcInferenceStats=false,
                const bool calcInputStats=false);
    
    MultRegress(const std::vector<double>& y,
                const MatrixT& X,
                const std::vector<std::string>& varNames,
                const bool calcInferenceStats=false,
                const bool calcInputStats=false);

    MultRegress();

    const std::vector<double>& beta() const {return bp_;}
    const std::vector<double>& stderr() const {return se_;}
    const std::vector<double>& resid() const {return resid_;}
    const std::vector<double>& regressorMeans() const {return mu_;}
    const std::vector<double>& regressorStdDevs() const {return sd_;}
    
    double R2() const;
    double R2Adj() const;
    double stdDevResid() const;
    double predict(const double* first) const;
    
    void predict(const double* first, double& yhat, double& yerr) const;
    void rmseCalc(const std::vector<double>& y,
                  const MatrixT& X,
                  double& nullRmse,
                  double& predRmse) const;
    
    friend std::ostream& operator<<(std::ostream& os,
                                    const MultRegress<MatrixT>& x) {
        x.print(os);
        return os;
    }
    
private:
    void init(const std::vector<double>& y,
              MatrixT& X,
              bool calcInferenceStats,
              bool calcInputStats);
    
    void init(const std::vector<double>& y,
              const MatrixT& X,
              bool calcInferenceStats,
              bool calcInputStats);
    
    void coeffCalc(const std::vector<double>& y,
                   MatrixT& X,
                   bool calcInferenceStats);
    void coeffCalcImpl(const std::vector<double>& y, MatrixT& X);
    void coeffCalcInferenceImpl(const std::vector<double>& y, const MatrixT& X);
    void inferenceCalc(const std::vector<double>& y, const MatrixT& X);
    void inpStatsCalc(const MatrixT& X);
    void print(std::ostream& os) const;
    void printModelImpl(std::ostream& os) const;

private:
    const int ip_;
    const int ldx_;
    const std::vector<std::string> varNames_;

    std::vector<double> mu_;
    std::vector<double> sd_;

    double tss_;
    double rss_;
    std::vector<double> bp_;
    std::vector<double> se_;
    std::vector<double> resid_;
    MatrixT cov_;
};
    
} // namespace

#include "multregress.ch"

/*(+doc.)
  .if
  .SH DESCRIPTION

  .SS

  .SH EXAMPLES

  .SH FILES
  .nf
  util/multregress.h - class declaration
  .if

  .SH SEE ALSO

  .SH AUTHOR
  Murat Ahmed
  (-doc.)*/

#endif // MULT_REGRESS_H_INCLUDED
