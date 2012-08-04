#ifndef QUANTILE_CALC_H_INCLUDED
#define QUANTILE_CALC_H_INCLUDED

/** CPP INCLUDES ***************************************************************/
#include <vector>

/** APP INCLUDES ***************************************************************/

namespace Util {

/**************************** EfficientMovingMedian ****************************/

template <typename T>
class EfficientMovingMedian {
public:
    EfficientMovingMedian(float tau=0.005f, float initVal=0.f);

    void newData(T val, float minVal, float maxVal);
    float value() const {return value_;}

private:
    float signedZScore(T value);

private:
    const float tau_;
    const float stddev_;

    int numObs_;
    float value_;
    float signEma_;
    float scale_;
};

/***************************** MovingQuantileTree ******************************/

template <typename T>
class MovingQuantileTree {
public:
    MovingQuantileTree(int depth, float lowerBound, float upperBound,
                       float tau=0.005f);
    ~MovingQuantileTree();
    
    void newData(T val);
    float symmetricMedian() const;  // for symmetric distributions
    const std::vector<float>& quantiles() const;

    friend std::ostream& operator<<(std::ostream& os,
                                    const MovingQuantileTree<T>& x) {
        x.print(os);
        return os;
    }
    
private:
    class node;
    
    void init(float lowerBound, float upperBound);
    void create(float val, float lowerBound, float upperBound, int depth);
    void destroy(node* pleaf);
    void insert(float val);
    void insert(float val, node* pleaf);
    void newData(T val, node* pleaf);
    void getValues() const;
    void getValuesImpl(node* pcurr) const;
    void print(std::ostream& os) const;
    
private:
    class node {
    public:
        node(float tau, T val);

        void newData(T val);
        float value() const {return median_.value();}

        node* const left() const {return left_;}
        node* const right() const {return right_;}

        void setLeft(node* p) {left_ = p;}
        void setRight(node* p) {right_ = p;}

    private:
        node* left_;
        node* right_;
        EfficientMovingMedian<T> median_;
    };
    
private:
    const int depth_;
    const float tau_;

    node* head_;
    mutable std::vector<float> quantiles_;
};

} // namespace

#include "quantilecalc.ch"

/*(+doc.)
  .if
  .SH DESCRIPTION

  .SS

  .SH EXAMPLES

  .SH FILES
  .nf
  util/quantilecalc.h - class declaration
  .if

  .SH SEE ALSO

  .SH AUTHOR
  Murat Ahmed
  (-doc.)*/

#endif // QUANTILE_CALC_H included
