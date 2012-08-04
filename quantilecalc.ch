#ifndef QUANTILE_CALC_CH_INCLUDED
#define QUANTILE_CALC_CH_INCLUDED

/** CPP INCLUDES ***************************************************************/
#include <cmath>
#include <limits>

/** APP INCLUDES ***************************************************************/

namespace Util {

namespace Detail {
    
    double
    erf(double x)
    {
        const double a1 = 0.254829592;
        const double a2 = -0.284496736;
        const double a3 = 1.421413741;
        const double a4 = -1.453152027;
        const double a5 = 1.061405429;
        const double p = 0.3275911;

        const int sg = (x < 0.) ? -1 : 1;
        x = std::abs(x);
        
        const double t = 1. / (1. + p*x);
        const double z = (1. - ((((a5*t + a4)*t + a3)*t + a2)*t + a1) *
                          t * std::exp(-x*x));

        return sg * z;
    }
    
    template <typename T>
    int8_t
    sgn(const T val)
    {
        return (T(0) < val) - (val < T(0));
    }

    template <typename T>
    void
    boundVal(T& val, const T minVal, const T maxVal)
    {
        if (val < minVal) {
            val = minVal;
        }
        else if (val > maxVal) {
            val = maxVal;
        }
    }
    
    template <class NodeT>
    float
    upperBound(NodeT* pleaf)
    {
        float val = std::numeric_limits<float>::max();
        NodeT* pcurr = pleaf->right();
        
        while (pcurr) {
            val = pcurr->value();
            pcurr = pcurr->left();
        }
        
        return val;
    }

    template <class NodeT>
    float
    lowerBound(NodeT* pleaf)
    {
        float val = -std::numeric_limits<float>::max();
        NodeT* pcurr = pleaf->left();

        while (pcurr) {
            val = pcurr->value();
            pcurr = pcurr->right();
        }

        return val;
    }
    
} // namespace
    
} // namespace

template <typename T>
Util::EfficientMovingMedian<T>::EfficientMovingMedian(const float tau,
                                                      const float initVal)
    : tau_(tau)
    , stddev_(sqrt(tau / (2.f - tau)))
    , numObs_(1)
    , value_(initVal)
    , signEma_(0.f)
    , scale_(initVal)
{
    assert(tau < 1.f && tau > 0.f);
}

template <typename T>
void
Util::EfficientMovingMedian<T>::newData(const T value,
                                        const float minVal,
                                        const float maxVal)
{
    scale_ += tau_ * (std::abs(value) - scale_);
    value_ += tau_ * scale_ * signedZScore(value);
    Detail::boundVal(value_,minVal,maxVal);
    
    ++numObs_;
}

template <typename T>
float
Util::EfficientMovingMedian<T>::signedZScore(const T value)
{
    const int8_t sg = Detail::sgn(value - value_);
    const float corrctn = sqrt(1.f - std::pow(1 - tau_,2*numObs_));
    
    signEma_ += tau_ * (sg - signEma_);
    
    const float z = std::abs(signEma_) / (stddev_ * corrctn);
    
    return sg * Detail::erf(z / sqrt(2));
}

template <typename T>
Util::MovingQuantileTree<T>::node::node(const float tau, const T val)
    : left_(NULL)
    , right_(NULL)
    , median_(tau,val)
{
}

template <typename T>
void
Util::MovingQuantileTree<T>::node::newData(const T val)
{
    median_.newData(val,
                    Detail::lowerBound(this),
                    Detail::upperBound(this));
}

template <typename T>
Util::MovingQuantileTree<T>::MovingQuantileTree(const int depth,
                                                const float lowerBound,
                                                const float upperBound,
                                                const float tau)
    : depth_(depth)
    , tau_(tau)
    , head_(NULL)
    , quantiles_()
{
    init(lowerBound,upperBound);
}

template <typename T>
Util::MovingQuantileTree<T>::~MovingQuantileTree()
{
    destroy(head_);
}

template <typename T>
void
Util::MovingQuantileTree<T>::init(const float lowerBound, const float upperBound)
{
    const float val = (lowerBound + upperBound) / 2.f;
    insert(val);
    create(val,lowerBound,upperBound,0);
}

template <typename T>
void
Util::MovingQuantileTree<T>::create(const float value,
                                    const float lowerBound,
                                    const float upperBound,
                                    const int depth)
{
    if (depth > depth_)
        return;

    const float offset = (upperBound - lowerBound) * std::pow(2.f,-(depth + 2));
    const float lval = value - offset;
    const float rval = value + offset;

    insert(lval);
    insert(rval);
    
    create(lval,lowerBound,upperBound,depth+1);
    create(rval,lowerBound,upperBound,depth+1);
}

template <typename T>
void
Util::MovingQuantileTree<T>::insert(const float val)
{
    if (head_)
        insert(val,head_);
    else
        head_ = new node(tau_,val);
}

template <typename T>
void
Util::MovingQuantileTree<T>::insert(const float value, node* pleaf)
{
    if (value < pleaf->value()) {
        
        if (pleaf->left())
            insert(value,pleaf->left());
        else
            pleaf->setLeft(new node(tau_,value));
    }
    else {

        if (pleaf->right())
            insert(value,pleaf->right());
        else
            pleaf->setRight(new node(tau_,value));
    }
}

template <typename T>
void
Util::MovingQuantileTree<T>::newData(const T val)
{
    newData(val,head_);
}

template <typename T>
void
Util::MovingQuantileTree<T>::newData(const T val, node* pleaf)
{
    if (pleaf) {
        
        pleaf->newData(val);
        
        if (val < pleaf->value())
            newData(val,pleaf->left());
        else
            newData(val,pleaf->right());
    }
}

template <typename T>
void
Util::MovingQuantileTree<T>::destroy(node* pleaf)
{
    if (pleaf) {

        destroy(pleaf->left());
        destroy(pleaf->right());
        delete pleaf;
    }
}

template <typename T>
void
Util::MovingQuantileTree<T>::getValues() const
{
    quantiles_.clear();
    quantiles_.reserve(1 << (depth_ + 1) - 1);

    getValuesImpl(head_);
}

template <typename T>
void
Util::MovingQuantileTree<T>::getValuesImpl(node* pcurr) const
{
    if (!pcurr)
        return;
    
    if (pcurr->left())
        getValuesImpl(pcurr->left());
    
    quantiles_.push_back(pcurr->value());
    
    if (pcurr->right())
        getValuesImpl(pcurr->right());   
}

template <typename T>
const std::vector<float>&
Util::MovingQuantileTree<T>::quantiles() const
{
    getValues();
    
    return quantiles_;
}

template <typename T>
float
Util::MovingQuantileTree<T>::symmetricMedian() const
{
    const int n = 1 << depth_;
    const std::vector<float>& q = quantiles();

    std::vector<float>::const_iterator head = q.begin();
    std::vector<float>::const_iterator tail = q.end();
    --tail;

    int cnt = 0;
    float median = 0.5f * q[n-1];
    
    while (head != tail) {

        median += std::pow(2.f,cnt - n - 1) * (*head + *tail);
        ++head;
        --tail;
        ++cnt;
    }

    median /= 1 - std::pow(2.f,-n);

    return median;
}

template <typename T>
void
Util:: MovingQuantileTree<T>::print(std::ostream& os) const
{
    const std::vector<float>& q = this->quantiles();
    
    for (std::vector<float>::const_iterator it = q.begin();
         it != q.end(); ++it) {
        
        os << *it << ' ';
    }
    os << std::endl;
}

#endif // QUANTILE_CALC_CH included
