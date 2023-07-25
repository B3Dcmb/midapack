#ifndef MAPPRAISER_RUNNING_SUM_H
#define MAPPRAISER_RUNNING_SUM_H

#include <cassert>
#include <cstddef>
#include <iostream>

template<typename ValueType = double, typename ValidType = bool>
class RunningSum {
public:
    explicit RunningSum(int w0);

    __attribute__((unused)) RunningSum(int w0, double thresh_low);

    void process_index(int nval, int idx, const ValueType *values, const ValidType *valid);

    double baseline_value() const;

    void print(std::ostream &out) const;

private:
    void internal_process_index(int nval, int idx, const ValueType *values, const ValidType *valid);

    void add(int pos, ValueType value, ValidType valid);

    void rmv(int pos, ValueType value, ValidType valid);

    double m_acc;           // current accumulated sum
    int    m_count;         // current (actual) number of samples accumulated in the sum
    int    m_w0;            // half size of the original window
    int    m_w;             // half size of the current window being used
    double m_threshold_low; // threshold to increase window size
    int    m_first_elem;    // smallest index currently part of the sum (whether valid or not)
    int    m_last_elem;     // largest index currently part of the sum (whether valid or not)
};
template<typename ValueType, typename ValidType>
void RunningSum<ValueType, ValidType>::print(std::ostream &out) const {
    out << "count " << m_count << " (first " << m_first_elem << ", last " << m_last_elem << ")";
}

template<typename ValueType, typename ValidType>
RunningSum<ValueType, ValidType>::RunningSum(int w0)
    : m_acc(0), m_count(0), m_w0(w0), m_w(w0), m_threshold_low(0.2), m_first_elem(1), m_last_elem(0) {}

template<typename ValueType, typename ValidType>
__attribute__((unused)) RunningSum<ValueType, ValidType>::RunningSum(int w0, double thresh_low)
    : m_acc(0), m_count(0), m_w0(w0), m_w(w0), m_threshold_low(thresh_low), m_first_elem(1), m_last_elem(0) {}

template<typename ValueType, typename ValidType>
void RunningSum<ValueType, ValidType>::process_index(int nval, int idx, const ValueType *values,
                                                     const ValidType *valid) {
    // process current index with current window
    // the current window depends on what happened while processing previous indices
    internal_process_index(nval, idx, values, valid);

    // nbr of samples in the average when there are no gaps
    int nbr_samples_ref = 2 * m_w0 + 1;

    // check if we have enough samples
    while (m_count < nbr_samples_ref * m_threshold_low) {
        // less than acceptable threshold

        // std::cout << "idx " << idx << " ratio " << static_cast<double>(m_count) / static_cast<double>
        // (nbr_samples_ref) << " < " << m_threshold_low << ", increasing window size" << std::endl;

        // increase window by 10 percent the original size
        m_w += 0.1 * m_w0;
        internal_process_index(nval, idx, values, valid);

        // if we have some samples and current window is larger than 4 times w0, stop
        if (m_count > 0 && m_w > 4 * m_w0) break;
    }

    while (m_w > m_w0) {
        if (m_count > nbr_samples_ref * 2 * m_threshold_low) {
            // std::cout << "idx " << idx << " above 2*thresh, decreasing window size" << std::endl;
            // try to go back to smaller window
            m_w -= 0.1 * m_w0;
            internal_process_index(nval, idx, values, valid);
        } else break;
    }
}

template<typename ValueType, typename ValidType>
void RunningSum<ValueType, ValidType>::internal_process_index(int nval, int idx, const ValueType *values,
                                                              const ValidType *valid) {
    int start = std::max(0, idx - m_w);
    int end   = std::min(nval - 1, idx + m_w - 1);

    // add elements between start and m_first_elem
    for (int i = m_first_elem - 1; i >= start; --i) {
        assert(0 <= i && i < nval);
        add(i, values[i], valid[i]);
    }

    // remove elements between m_first_elem and start
    for (int i = m_first_elem; i < start; ++i) {
        assert(0 <= i && i < nval);
        rmv(i, values[i], valid[i]);
    }

    // add elements between m_last_elem and end
    for (int i = m_last_elem + 1; i < end + 1; ++i) {
        assert(0 <= i && i < nval);
        add(i, values[i], valid[i]);
    }

    // remove elements between end and m_last_elem
    for (int i = m_last_elem; i > end; --i) {
        assert(0 <= i && i < nval);
        rmv(i, values[i], valid[i]);
    }
}

template<typename ValueType, typename ValidType>
void RunningSum<ValueType, ValidType>::add(int pos, ValueType value, ValidType valid) {
    if (valid) {
        ++m_count;
        m_acc += value;
    }
    // just added a value, update the first and last elements
    if (pos > m_last_elem) m_last_elem = pos;
    if (pos < m_first_elem) m_first_elem = pos;
}

template<typename ValueType, typename ValidType>
void RunningSum<ValueType, ValidType>::rmv(int pos, ValueType value, ValidType valid) {
    if (valid) {
        --m_count;
        m_acc -= value;
    }
    // just removed a value, update the first and last elements
    if (pos == m_last_elem) m_last_elem = pos - 1;
    if (pos == m_first_elem) m_first_elem = pos + 1;
}

template<typename ValueType, typename ValidType>
double RunningSum<ValueType, ValidType>::baseline_value() const {
    return m_acc / static_cast<ValueType>(m_count);
}


#endif // MAPPRAISER_RUNNING_SUM_H
