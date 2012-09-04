#ifndef MULT_REGRESS_CH_INCLUDED
#define MULT_REGRESS_CH_INCLUDED

/** CPP INCLUDES ***************************************************************/
#include <iomanip>
#include <numeric>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>

/** APP INCLUDES ***************************************************************/
#include "statutil.h"

template <class MatrixT>
Util::MultRegress<MatrixT>::MultRegress()
    : ip_(0)
    , ldx_(0)
    , varNames_()
    , mu_()
    , sd_()
    , tss_(-1.)
    , rss_(-1.)
    , bp_()
    , se_()
    , resid_()
    , cov_()
{
}

template <class MatrixT>
Util::MultRegress<MatrixT>::MultRegress(const std::vector<double>& y,
                                        MatrixT& X,
                                        const bool calcInferenceStats,
                                        const bool calcInputStats)
    : ip_(X.colSize())
    , ldx_(X.rowSize())
    , varNames_()
    , mu_()
    , sd_()
    , tss_(-1.)
    , rss_(-1.)
    , bp_(ip_)
    , se_()
    , resid_()
    , cov_()
{
    init(y,X,
         calcInferenceStats,
         calcInputStats);
}

template <class MatrixT>
Util::MultRegress<MatrixT>::MultRegress(const std::vector<double>& y,
                                        const MatrixT& X,
                                        const bool calcInferenceStats,
                                        const bool calcInputStats)
    : ip_(X.colSize())
    , ldx_(X.rowSize())
    , varNames_()
    , mu_()
    , sd_()
    , tss_(-1.)
    , rss_(-1.)
    , bp_(ip_)
    , se_()
    , resid_()
    , cov_()
{
    init(y,X,
         calcInferenceStats,
         calcInputStats);
}

template <class MatrixT>
Util::MultRegress<MatrixT>::MultRegress(const std::vector<double>& y,
                                        MatrixT& X,
                                        const std::vector<std::string>& varNames,
                                        const bool calcInferenceStats,
                                        const bool calcInputStats)
    : ip_(X.colSize())
    , ldx_(X.rowSize())
    , varNames_(varNames)
    , mu_()
    , sd_()
    , tss_(-1.)
    , rss_(-1.)
    , bp_(ip_)
    , se_()
    , resid_()
    , cov_()
{
    init(y,X,
         calcInferenceStats,
         calcInputStats);
}

template <class MatrixT>
Util::MultRegress<MatrixT>::MultRegress(const std::vector<double>& y,
                                        const MatrixT& X,
                                        const std::vector<std::string>& varNames,
                                        const bool calcInferenceStats,
                                        const bool calcInputStats)
    : ip_(X.colSize())
    , ldx_(X.rowSize())
    , varNames_(varNames)
    , mu_()
    , sd_()
    , tss_(-1.)
    , rss_(-1.)
    , bp_(ip_)
    , se_()
    , resid_()
    , cov_()
{
    init(y,X,
         calcInferenceStats,
         calcInputStats);
}

template <class MatrixT>
void
Util::MultRegress<MatrixT>::init(const std::vector<double>& y,
                                 MatrixT& X,
                                 const bool calcInferenceStats,
                                 const bool calcInputStats)
{
    coeffCalc(y,X,calcInferenceStats);

    if (calcInputStats)
        inpStatsCalc(X);
}

template <class MatrixT>
void
Util::MultRegress<MatrixT>::init(const std::vector<double>& y,
                                 const MatrixT& X,
                                 const bool calcInferenceStats,
                                 const bool calcInputStats)
{
    MatrixT XCopy = X;
    coeffCalc(y,XCopy,calcInferenceStats);
    
    if (calcInputStats)
        inpStatsCalc(X);
}

template <class MatrixT>
void
Util::MultRegress<MatrixT>::coeffCalc(const std::vector<double>& y,
                                      MatrixT& X,
                                      const bool calcInferenceStats)
{
    if (calcInferenceStats)
        coeffCalcInferenceImpl(y,X);
    else
        coeffCalcImpl(y,X);
}

template <class MatrixT>
void
Util::MultRegress<MatrixT>::coeffCalcImpl(const std::vector<double>& y,
                                          MatrixT& X)
{
    const size_t m = X.rowSize();
    const size_t n = X.colSize();

    gsl_matrix_view a = gsl_matrix_view_array(X,m,n);
    gsl_vector_view x = gsl_vector_view_array(&bp_[0],bp_.size());
    gsl_vector_const_view b = gsl_vector_const_view_array(&y[0],y.size());

    {
        gsl_vector* resid = gsl_vector_alloc(m);
        gsl_vector* tau = gsl_vector_alloc(std::min(m,n));
        
        gsl_linalg_QR_decomp(&a.matrix,tau);
        gsl_linalg_QR_lssolve(&a.matrix,tau,&b.vector,&x.vector,resid);
        
        gsl_vector_free(tau);
        gsl_vector_free(resid);
    }
}

template <class MatrixT>
void
Util::MultRegress<MatrixT>::coeffCalcInferenceImpl(const std::vector<double>& y,
                                                   const MatrixT& X)
{
    const size_t m = X.rowSize();
    const size_t n = X.colSize();

    cov_.resize(n,n);
    
    gsl_matrix_const_view a = gsl_matrix_const_view_array(X,m,n);
    gsl_vector_const_view b = gsl_vector_const_view_array(&y[0],y.size());
    gsl_vector_view x = gsl_vector_view_array(&bp_[0],bp_.size());
    gsl_matrix_view cov = gsl_matrix_view_array(cov_,n,n);
    
    {
        gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(m,n);
        gsl_multifit_linear(&a.matrix,&b.vector,&x.vector,&cov.matrix,
                            &rss_,work);
        gsl_multifit_linear_free(work);
    }
    
    inferenceCalc(y,X);
}

template <class MatrixT>
void
Util::MultRegress<MatrixT>::inferenceCalc(const std::vector<double>& y,
                                          const MatrixT& X)
{
    const double sd = Util::stddev<double>(y.begin(),y.end());
    tss_ = (y.size() - 1) * sd * sd;
    resid_.resize(ldx_);
    
    gsl_vector_view r = gsl_vector_view_array(&resid_[0],resid_.size());
    gsl_vector_const_view x = gsl_vector_const_view_array(&bp_[0],bp_.size());
    gsl_vector_const_view b = gsl_vector_const_view_array(&y[0],y.size());
    gsl_matrix_const_view a = gsl_matrix_const_view_array(X,X.rowSize(),X.colSize());

    {
        gsl_multifit_linear_residuals(&a.matrix,&b.vector,&x.vector,&r.vector);
    }

    se_.reserve(ip_);
    for (int i = 0; i < ip_; i++) {    
        se_.push_back(std::sqrt(cov_(i,i)));
    }
}

template <class MatrixT>
void
Util::MultRegress<MatrixT>::inpStatsCalc(const MatrixT& X)
{   
    Util::colMeans(X.begin(),X.end(),X.colSize(),mu_);
    Util::colStdDevs(X.begin(),X.end(),X.colSize(),sd_);
}

template <class MatrixT>
void
Util::MultRegress<MatrixT>::predict(const double* first,
                                    double& yhat,
                                    double& yerr) const
{
    if (cov_.empty()) {
        
        yhat = predict(first);
        return;
    }
    
    const size_t n = bp_.size();
    {
        gsl_vector_const_view x = gsl_vector_const_view_array(first,n);
        gsl_vector_const_view c = gsl_vector_const_view_array(&bp_[0],n);
        gsl_matrix_const_view cov = gsl_matrix_const_view_array(cov_,n,n);       
        gsl_multifit_linear_est(&x.vector,&c.vector,&cov.matrix,&yhat,&yerr);
    }
}

template <class MatrixT>
double
Util::MultRegress<MatrixT>::predict(const double* first) const
{
    double yhat = 0.;
    
    return std::inner_product(bp_.begin(),bp_.end(),first,yhat);
}

template <class MatrixT>
void
Util::MultRegress<MatrixT>::rmseCalc(const std::vector<double>& y,
                                     const MatrixT& X,
                                     double& nullRmse,
                                     double& predRmse) const
{
    const size_t sz = y.size();
    if ((X.rowSize() != sz) || (X.colSize() != ip_))
        return;

    double ssq = 0.;
    double nssq = 0.;
    
    for (int i = 0; i < sz; i++) {

        const typename MatrixT::const_iterator first = &X(i,0);
        const double val = y[i];
        const double err = val - predict(first);
        ssq += err * err;
        nssq += val * val;
    }

    nullRmse = std::sqrt(nssq / sz);
    predRmse = std::sqrt(ssq / sz);
}

template <class MatrixT>
double
Util::MultRegress<MatrixT>::R2() const
{
    return 1. - rss_ / tss_;
}

template <class MatrixT>
double
Util::MultRegress<MatrixT>::R2Adj() const
{
    return 1. - rss_ / tss_ * (ldx_ - 1) / (ldx_ - ip_);
}

template <class MatrixT>
double
Util::MultRegress<MatrixT>::stdDevResid() const
{
    return (ip_ > 0) ? std::sqrt(rss_ / (ldx_ - ip_)) : -1.;
}

template <class MatrixT>
void
Util::MultRegress<MatrixT>::print(std::ostream& os) const
{
    const int pr = os.precision(4);
    
    os << std::setw(25) << '\t'
       << std::setw(10) << "Value" << '\t'
       << std::setw(10) << "Std. Error" << '\t'
       << std::setw(10) << "t value" << '\t'
       << std::setw(10) << "Pr(>|t|)" << std::endl;

    if ((bp_.size() == ip_) && (se_.size() == ip_))
        printModelImpl(os);
    
    const double fstat = (ldx_ - ip_) * (tss_ / rss_ - 1) / (ip_ - 1);
    const double fpval = 1 - gsl_cdf_fdist_P(fstat,ip_ - 1,ldx_ - ip_);

    os << std::endl;
    os << "Residual standard error: " << stdDevResid()
       << " on " << ldx_ - ip_ << " degrees of freedom\n"
       << "Multiple R-Squared: " << R2() << ",     "
       << "Adjusted R-Squared: " << R2Adj() << '\n'
       << "F-statistic: " << fstat
       << " on " << ip_ - 1
       << " and " << ldx_ - ip_
       << " DF,  p-value: " << fpval
       << std::endl;

    os.precision(pr);
}

template <class MatrixT>
void
Util::MultRegress<MatrixT>::printModelImpl(std::ostream& os) const
{
    for (int i = 0; i < ip_; i++) {

        //int ifail = 1;
        const double bp = bp_[i];
        const double se = se_[i];
        const double tpval = 2. * (1 - gsl_cdf_tdist_P(std::abs(bp / se),
                                                       ldx_ - ip_));
        
        std::stringstream sstr;
        sstr << "#VAR_";
        if (i == 0)
            sstr << "(Intercept)";
        else {
            if (varNames_.empty())
                sstr << "x" << i;
            else
                sstr << varNames_[i - 1];
        }
        
        os << std::setw(30) << sstr.str() << '\t'
           << std::setw(10) << bp << '\t'
           << std::setw(10) << se << '\t'
           << std::setw(10) << bp / se << '\t'
           << std::setw(10) << tpval << std::endl;
    }
}

#endif // MULT_REGRESS_CH included
