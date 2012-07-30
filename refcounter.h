#ifndef REF_COUNTER_H_INCLUDED
#define REF_COUNTER_H_INCLUDED

/** CPP INCLUDES ***************************************************************/

/** APP INCLUDES ***************************************************************/

namespace Util {

template <class T>
class RefCntPtr {
public:
    RefCntPtr(T* p=0);
    RefCntPtr(const RefCntPtr& rhs);
    ~RefCntPtr();

    RefCntPtr& operator=(const RefCntPtr<T>& rhs);

    T* operator->() const;
    T& operator*() const;

private:
    void init();
    
private:
    T* pT_;
};
    
class RefCntObj {
public:
    RefCntObj() : shareable_(true), refCount_(0) {};
    
    RefCntObj(const RefCntObj& rhs) : shareable_(true), refCount_(0) {};
    RefCntObj& operator=(const RefCntObj& rhs) {return *this;}

    virtual ~RefCntObj() = 0;
    
    void addRef();
    void delRef();
    size_t refCount() const {return refCount_;}
    
    bool isShared() const;
    bool shareable() const {return shareable_;}
    void setUnshareable() {shareable_ = false;}
    
private:
    bool shareable_;
    size_t refCount_;
};
    
} // namespace

#include "refcounter.ch"

/*(+doc.)
  .if
  .SH DESCRIPTION
  Adaptation of Scott Meyers More Effective C++ Item 29.
  
  .SS

  .SH EXAMPLES

  .SH FILES
  .nf
  util/refcounter.h - class declaration
  .if

  .SH SEE ALSO

  .SH AUTHOR
  (-doc.)*/

#endif // REF_COUNTER_H included
