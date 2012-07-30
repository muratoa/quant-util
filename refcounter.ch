template <class T>
Util::RefCntPtr<T>::RefCntPtr(T* p)
    : pT_(p)
{
    init();
}

template <class T>
Util::RefCntPtr<T>::RefCntPtr(const Util::RefCntPtr<T>& rhs)
    : pT_(rhs.pT_)
{
    init();
}

template <class T>
Util::RefCntPtr<T>::~RefCntPtr()
{
    if (pT_)
        pT_->delRef();
}

template <class T>
Util::RefCntPtr<T>&
Util::RefCntPtr<T>::operator=(const Util::RefCntPtr<T>& rhs)
{
    if (pT_ != rhs.pT_) {
        
        T* ptmp = pT_;
        pT_ = rhs.pT_;
        init();

        if (ptmp)
            ptmp->delRef();
    }

    return *this;
}

template <class T>
void
Util::RefCntPtr<T>::init()
{
    if (!pT_)
        return;

    if (!pT_->shareable()) {
        pT_ = new T(*pT_);
    }

    pT_->addRef();
}

template <class T>
T*
Util::RefCntPtr<T>::operator->() const
{
    return pT_;
}

template <class T>
T&
Util::RefCntPtr<T>::operator*() const
{
    return *pT_;
}

Util::RefCntObj::~RefCntObj()
{
}

void
Util::RefCntObj::addRef()
{
    ++refCount_;
}

void
Util::RefCntObj::delRef()
{
    if (--refCount_ == 0)
        delete this;
}

bool
Util::RefCntObj::isShared() const
{
    return refCount_ > 1;
}
