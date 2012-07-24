#include "simutimeval.h"

/** CPP INCLUDES ***************************************************************/

/** APP INCLUDES ***************************************************************/

// SimuTimeval declaration and fns for extern declaration.
Util::SimuTimeval simuTimeval_;

const Util::SimuTimeval&
Util::getSimuTimeval()
{
    return simuTimeval_;
}

const Util::SimuTime
Util::getSimuTime()
{
    return SimuTime(simuTimeval_);
}

void
Util::setSimuTimeval(const SimuTimeval& t)
{
    simuTimeval_ = t;
}

Util::SimuTimeval::SimuTimeval()
    : currTV_()
{
    reset();
}

Util::SimuTimeval::SimuTimeval(const struct timeval& t)
    : currTV_()
{
    reset(t.tv_sec,t.tv_usec);
}

Util::SimuTimeval::SimuTimeval(const unsigned long sec)
    : currTV_()
{
    reset(sec);
}

Util::SimuTimeval::SimuTimeval(const unsigned long sec, const long usec)
    : currTV_()
{
    reset(sec,usec);
}

Util::SimuTimeval::SimuTimeval(const double ts)
    : currTV_()
{
    const long sec = long(ts);
    const long usec = long((ts - sec) * 1000000);
    reset(sec,usec);
}

Util::SimuTimeval::SimuTimeval(const SimuTimeval& rhs)
    : currTV_(rhs.currTV_)
{
}

Util::SimuTimeval::~SimuTimeval()
{
}

Util::SimuTimeval::operator double() const
{
    const double sec = double(currTV_.tv_sec);

    return sec + currTV_.tv_usec / 1000000.;
}

Util::SimuTimeval&
Util::SimuTimeval::operator=(const SimuTimeval& rhs)
{
    if (this != &rhs) {
        currTV_ = rhs.currTV_;
    }

    return *this;
}

Util::SimuTimeval&
Util::SimuTimeval::operator+=(const SimuTimeval& rhs)
{
    currTV_.tv_sec += rhs.currTV_.tv_sec;
    currTV_.tv_usec += rhs.currTV_.tv_usec;
    standardize(currTV_);

    return *this;
}

Util::SimuTimeval&
Util::SimuTimeval::operator-=(const SimuTimeval& rhs)
{
    currTV_.tv_sec -= rhs.currTV_.tv_sec;
    currTV_.tv_usec -= rhs.currTV_.tv_usec;
    standardize(currTV_);

    return *this;
}

bool
Util::SimuTimeval::operator<(const SimuTimeval& rhs) const
{
    return ((currTV_.tv_sec < rhs.currTV_.tv_sec) ||
            (currTV_.tv_sec == rhs.currTV_.tv_sec &&
             currTV_.tv_usec < rhs.currTV_.tv_usec));
}

bool
Util::SimuTimeval::operator>(const SimuTimeval& rhs) const
{
    return ((currTV_.tv_sec > rhs.currTV_.tv_sec) ||
            (currTV_.tv_sec == rhs.currTV_.tv_sec &&
             currTV_.tv_usec > rhs.currTV_.tv_usec));
}

bool
Util::SimuTimeval::operator==(const SimuTimeval& rhs) const
{
    return ((currTV_.tv_sec == rhs.currTV_.tv_sec) &&
            (currTV_.tv_usec == rhs.currTV_.tv_usec));
}

bool
Util::SimuTimeval::operator!=(const SimuTimeval& rhs) const
{
    return ((currTV_.tv_sec != rhs.currTV_.tv_sec) ||
            (currTV_.tv_usec != rhs.currTV_.tv_usec));
}

void
Util::SimuTimeval::reset(const long sec, const long usec)
{
    currTV_.tv_sec = sec;
    currTV_.tv_usec = usec;
    standardize(currTV_);
}

void
Util::SimuTimeval::standardize(struct timeval& t)
{
    long sec = t.tv_usec / 1000000;

    if (sec != 0) {
        t.tv_sec += sec;
        t.tv_usec -= (sec * 1000000);
    }
}

Util::SimuTime::SimuTime()
    : sec_(0)
{
}

Util::SimuTime::SimuTime(const unsigned long sec)
    : sec_(sec)
{
}

Util::SimuTime::SimuTime(const SimuTimeval& t)
    : sec_(static_cast<long>(t))
{
}

Util::SimuTime::operator long() const
{
    return sec_;
}

Util::SimuTime::SimuTime(const SimuTime& rhs)
    : sec_(rhs.sec_)
{
}

Util::SimuTime::~SimuTime()
{
}

Util::SimuTime&
Util::SimuTime::operator=(const SimuTime& rhs)
{
    if (this != &rhs) {
        sec_ = rhs.sec_;
    }

    return *this;
}

Util::SimuTime&
Util::SimuTime::operator+=(const SimuTime& rhs)
{
    sec_ += rhs.sec_;

    return *this;
}

Util::SimuTime&
Util::SimuTime::operator-=(const SimuTime& rhs)
{
    sec_ -= rhs.sec_;
    
    return *this;
}

bool
Util::SimuTime::operator<(const SimuTime& rhs) const
{
    return sec_ < rhs.sec_;
}

bool
Util::SimuTime::operator>(const SimuTime& rhs) const
{
    return sec_ > rhs.sec_;
}

bool
Util::SimuTime::operator==(const SimuTime& rhs) const
{
    return sec_ == rhs.sec_;
}

bool
Util::SimuTime::operator!=(const SimuTime& rhs) const
{
    return sec_ == rhs.sec_;
}

void
Util::SimuTime::reset(const long sec)
{
    sec_ = sec;
}

const Util::SimuTimeval
Util::operator+(const SimuTimeval& lhs, const SimuTimeval& rhs)
{
    return SimuTimeval(lhs) += rhs;
}

const Util::SimuTimeval
Util::operator-(const SimuTimeval& lhs, const SimuTimeval& rhs)
{
    return SimuTimeval(lhs) -= rhs;
}

const Util::SimuTime
Util::operator+(const SimuTime& lhs, const SimuTime& rhs)
{
    return SimuTime(lhs) += rhs;
}

const Util::SimuTime
Util::operator-(const SimuTime& lhs, const SimuTime& rhs)
{
    return SimuTime(lhs) -= rhs;
}
