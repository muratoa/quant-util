#ifndef SIMU_TIMEVAL_H_INCLUDED
#define SIMU_TIMEVAL_H_INCLUDED

/** CPP INCLUDES ***************************************************************/
#include <sys/time.h>

/** APP INCLUDES ***************************************************************/

namespace Util {
    
/******************************** SimuTimeval **********************************/
    
class SimuTimeval {
public:

    SimuTimeval();
    SimuTimeval(const struct timeval& t);
    SimuTimeval(unsigned long sec);
    SimuTimeval(unsigned long sec, long usec);
    
    explicit SimuTimeval(double ts);
    operator double() const;
    
    SimuTimeval(const SimuTimeval& rhs);
    ~SimuTimeval();
    
    SimuTimeval& operator=(const SimuTimeval& rhs);
    SimuTimeval& operator+=(const SimuTimeval& rhs);
    SimuTimeval& operator-=(const SimuTimeval& rhs);
    
    bool operator< (const SimuTimeval& rhs) const;
    bool operator> (const SimuTimeval& rhs) const;
    bool operator==(const SimuTimeval& rhs) const;
    bool operator!=(const SimuTimeval& rhs) const;
    
    void reset(long sec=0, long usec=0);

private:
    static void standardize(struct timeval& tv);
    
private:
    struct timeval currTV_;
};

/******************************** SimuTime *************************************/

class SimuTime {
public:

    SimuTime();
    explicit SimuTime(unsigned long sec);
    explicit SimuTime(const SimuTimeval& t);
    
    operator long() const;
    
    SimuTime(const SimuTime& rhs);
    ~SimuTime();

    SimuTime& operator=(const SimuTime& rhs);
    SimuTime& operator+=(const SimuTime& rhs);
    SimuTime& operator-=(const SimuTime& rhs);
    
    bool operator< (const SimuTime& rhs) const;
    bool operator> (const SimuTime& rhs) const;
    bool operator==(const SimuTime& rhs) const;
    bool operator!=(const SimuTime& rhs) const;

    void reset(long sec=0);
    
private:
    long sec_;
};

/******************************** Free Fns *************************************/    
    
void setSimuTimeval(const SimuTimeval& t);

const SimuTime getSimuTime();
const SimuTimeval& getSimuTimeval();

const SimuTime operator+(const SimuTime& lhs, const SimuTime& rhs);
const SimuTime operator-(const SimuTime& lhs, const SimuTime& rhs);
const SimuTimeval operator+(const SimuTimeval& lhs, const SimuTimeval& rhs);
const SimuTimeval operator-(const SimuTimeval& lhs, const SimuTimeval& rhs);

} // namespace

#endif // SIMU_TIME_VAL_H included
