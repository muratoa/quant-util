#ifndef CONVOLUTION_OPERATORS_H_INCLUDED
#define CONVOLUTION_OPERATORS_H_INCLUDED

/** CPP INCLUDES ***************************************************************/
#include <map>

/** APP INCLUDES ***************************************************************/
#include "util/refcounter.h"
#include "util/simutimeval.h"

namespace Util {

/******************************** ExpDecayValue ********************************/

template <typename T, size_t N>
class ExpDecayValue {
public:
    ExpDecayValue(const T* vals);
    ~ExpDecayValue();
    
    const T& operator[](int index) const;
    bool isNull() const;
    
private:
    struct ExpDecayObj : public Util::RefCntObj {
        
        ExpDecayObj(const T* vals);
        ExpDecayObj(const ExpDecayObj& rhs);
        ~ExpDecayObj();
        
        void init(const T* val);
        const T* data() const {return data_;}
        
    private:
        T* data_;
    };

private:
    Util::RefCntPtr<ExpDecayObj> pValue_;
};

    
/******************************** DecayConvOp **********************************/

template <class TimeT, class ConvT, size_t N> class DecayConvOp;

/******************************** ConvTypes ************************************/
    
class EMA {
public:
    void valueUpdate(float& value, float newValue, float decayValue);
};

class EDA {
public:
    void valueUpdate(float& value, float newValue, float decayValue);
};
    
/******************************** Friend Helpers *******************************/

template <class ConvT, size_t N>
void newDataImpl(DecayConvOp<SimuTime,ConvT,N>& x, float val);

template <class ConvT, size_t N>
void newDataImpl(DecayConvOp<SimuTimeval,ConvT,N>& x, float val);    
    
template <class TimeT, size_t N>
float valueImpl(const DecayConvOp<TimeT,EMA,N>& x);
    
template <size_t N>
float valueImpl(const DecayConvOp<SimuTime,EDA,N>& x);

template <size_t N>
float valueImpl(const DecayConvOp<SimuTimeval,EDA,N>& x);

template <class ConvT, size_t N>
int decayIdxImpl(const DecayConvOp<SimuTime,ConvT,N>& x,
                 const SimuTime& currTime);

template <class ConvT, size_t N>
int decayIdxImpl(const DecayConvOp<SimuTimeval,ConvT,N>& x,
                 const SimuTimeval& currTime);  
    
    
/******************************** ExpDecayMap **********************************/

template <size_t N>
class ExpDecayMap {
public:
    typedef std::map<unsigned int, ExpDecayValue<float,N> > map_type;
    
    ExpDecayMap();
    ExpDecayValue<float,N> operator[](unsigned int tau);

private:
    bool elemIsNull(const typename map_type::const_iterator& it) const;
    
private:
    map_type data_;
};
    
/******************************** DecayConvOp **********************************/    
    
template <class TimeT, class ConvT, size_t N>
class DecayConvOp : private ConvT {
public:
    DecayConvOp(float initVal, unsigned int tau);
    
    void newData(float val);
    void resetTau(unsigned int tau);

    float value() const;
    
public:
    friend void newDataImpl <>(DecayConvOp& x, float val);
    friend float valueImpl <>(const DecayConvOp& x);
    friend int decayIdxImpl <>(const DecayConvOp& x, const TimeT& currTime);
    
private:
    int decayIdx(const TimeT& currTime) const;
    
private:
    static ExpDecayMap<N> expDecayMap_;
    ExpDecayValue<float,N> expDecayVals_;
    
    float value_;
    TimeT prevTime_;
};

/******************************** typedefs *************************************/

typedef DecayConvOp<SimuTime,EMA,1000> ExpMvAvg;
typedef DecayConvOp<SimuTimeval,EMA,2000> ExpMvAvgMs;

typedef DecayConvOp<SimuTime,EDA,1000> ExpDecayAvg;
typedef DecayConvOp<SimuTimeval,EDA,2000> ExpDecayAvgMs;
    
} // namespace

#include "convolutionoperators.ch"

/*(+doc.)
  .if
  .SH DESCRIPTION

  .SS

  .SH EXAMPLES

  .SH FILES
  .nf
  util/convolutionoperators.h - class declaration
  .if

  .SH SEE ALSO

  .SH AUTHOR
  Murat Ahmed
  (-doc.)*/

#endif // CONVOLUTION_OPERATORS_H included
