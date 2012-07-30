/** CPP INCLUDES ***************************************************************/
#include <algorithm>
#include <cmath>

template <typename T, size_t N>
Util::ExpDecayValue<T,N>::ExpDecayValue(const T* vals)
    : pValue_(new ExpDecayObj(vals))
{
}

template <typename T, size_t N>
Util::ExpDecayValue<T,N>::~ExpDecayValue()
{
    // Delete when ref cnt hits 1 as init vals are stored in static map.
    if (pValue_->refCount() == 2)
        pValue_->delRef();
}

template <typename T, size_t N>
const T&
Util::ExpDecayValue<T,N>::operator[](int index) const
{
    return pValue_->data()[index];
}

template <typename T, size_t N>
bool
Util::ExpDecayValue<T,N>::isNull() const
{
    return pValue_->refCount() == 0;
}

template <typename T, size_t N>
Util::ExpDecayValue<T,N>::ExpDecayObj::ExpDecayObj(const T* vals)
    : data_(NULL)
{
    init(vals);
}

template <typename T, size_t N>
Util::ExpDecayValue<T,N>::ExpDecayObj::ExpDecayObj(const ExpDecayObj& rhs)
{
    init(rhs.data_);
}

template <typename T, size_t N>
Util::ExpDecayValue<T,N>::ExpDecayObj::~ExpDecayObj()
{
    delete[] data_;
}

template <typename T, size_t N>
void
Util::ExpDecayValue<T,N>::ExpDecayObj::init(const T* val)
{
    data_ = new T[N];
    std::copy(val,val + N,data_);
}

template <size_t N>
Util::ExpDecayMap<N>::ExpDecayMap()
    : data_()
{
}

template <size_t N>
Util::ExpDecayValue<float,N>
Util::ExpDecayMap<N>::operator[](const unsigned int tau)
{
    typename map_type::iterator it = data_.find(tau);

    if (it != data_.end()) {
        
        if (!elemIsNull(it))
            return it->second;
        else
            data_.erase(it);
    }
    
    float vals[N];
    float ftau = static_cast<float>(tau);
    
    for (int i = 0; i < N; i++) {
        vals[i] = 1.f - std::exp(-(i+1) / ftau);
    }
    
    ExpDecayValue<float,N> edv(vals);
    data_.insert(std::make_pair(tau,edv));
    
    return edv;
}

template <size_t N>
bool
Util::ExpDecayMap<N>::elemIsNull(const typename map_type::const_iterator& it) const
{
    return it->second.isNull();
}

template <class TimeT, class ConvT, size_t N>
Util::ExpDecayMap<N> Util::DecayConvOp<TimeT,ConvT,N>::expDecayMap_;

template <class TimeT, class ConvT, size_t N>
Util::DecayConvOp<TimeT,ConvT,N>::DecayConvOp(const float initVal,
                                              const unsigned int tau)
    : expDecayVals_(expDecayMap_[tau])
    , value_(initVal)
    , prevTime_()
{
}

template <class TimeT, class ConvT, size_t N>
void
Util::DecayConvOp<TimeT,ConvT,N>::newData(float val)
{
    Util::newDataImpl(*this,val);
}

template <class TimeT, class ConvT, size_t N>
void
Util::DecayConvOp<TimeT,ConvT,N>::resetTau(const unsigned int tau)
{
    expDecayVals_ = expDecayMap_[tau];
}

template <class TimeT, class ConvT, size_t N>
int
Util::DecayConvOp<TimeT,ConvT,N>::decayIdx(const TimeT& currTime) const
{
    return Util::decayIdxImpl(*this,currTime);
}

template <class TimeT, class ConvT, size_t N>
float
Util::DecayConvOp<TimeT,ConvT,N>::value() const
{
    return Util::valueImpl(*this);
}

template <class ConvT, size_t N>
void
Util::newDataImpl(DecayConvOp<Util::SimuTime,ConvT,N>& x,
                  const float val)
{
    const Util::SimuTime currTime = getSimuTime();
    const float decay = x.expDecayVals_[x.decayIdx(currTime)];
    
    x.valueUpdate(x.value_,val,decay);
    x.prevTime_ = currTime;
}

template <class ConvT, size_t N>
void
Util::newDataImpl(DecayConvOp<Util::SimuTimeval,ConvT,N>& x,
                  const float val)
{
    const Util::SimuTimeval& currTime = getSimuTimeval();
    const float decay = x.expDecayVals_[x.decayIdx(currTime)];
    
    x.valueUpdate(x.value_,val,decay);
    x.prevTime_ = currTime;
}

template <class TimeT, size_t N>
float
Util::valueImpl(const DecayConvOp<TimeT,Util::EMA,N>& x)
{
    return x.value_;
}

template <size_t N>
float
Util::valueImpl(const DecayConvOp<Util::SimuTime,Util::EDA,N>& x)
{
    const Util::SimuTime currTime = getSimuTime();
    const float decay = x.expDecayVals_[x.decayIdx(currTime)];
    
    return x.value_ * (1 - decay);
}

template <size_t N>
float
Util::valueImpl(const DecayConvOp<Util::SimuTimeval,Util::EDA,N>& x)
{
    const Util::SimuTimeval& currTime = getSimuTimeval();
    const float decay = x.expDecayVals_[x.decayIdx(currTime)];

    return x.value_ * (1 - decay);
}

template <class ConvT, size_t N>
int
Util::decayIdxImpl(const DecayConvOp<Util::SimuTime,ConvT,N>& x,
                   const SimuTime& currTime)
{
    const int idx = currTime - x.prevTime_;
    
    if (idx < 0)
        return 0;
    
    return (idx < N) ? idx : (N - 1);
}

template <class ConvT, size_t N>
int
Util::decayIdxImpl(const DecayConvOp<Util::SimuTimeval,ConvT,N>& x,
                   const Util::SimuTimeval& currTime)
{
    const double res = currTime - x.prevTime_;
    
    if (res < 0.)
        return 0;
    
    const int idx = static_cast<int>(res * 1000);

    return (idx < N) ? idx : (N - 1);
}

void
Util::EMA::valueUpdate(float& value,
                       const float newValue, const float decayValue)
{
    value += decayValue * (newValue - value);
}

void
Util::EDA::valueUpdate(float& value,
                       const float newValue, const float decayValue)
{
    value += newValue - decayValue * value;
}
