#ifndef MAPPRAISER_RUNNING_SUM_H
#define MAPPRAISER_RUNNING_SUM_H

#include <cstddef>
#include <iostream>

template<typename ValueType = double, typename ValidType = bool>
class RunningSum {
public:
    explicit RunningSum(int window);

    void process_index(int nval, int idx, const ValueType *values, const ValidType *valid);

    void add(int pos, ValueType value, ValidType valid);

    void rmv(int pos, ValueType value, ValidType valid);

    double baseline_value() const;

private:
    double m_acc;        // current accumulated sum
    int    m_count;      // current (actual) number of samples accumulated in the sum
    int    m_w0;         // half size of the window
    int    m_first_elem; // smallest index currently part of the sum (whether valid or not)
    int    m_last_elem;  // largest index currently part of the sum (whether valid or not)
};

template<typename ValueType, typename ValidType>
RunningSum<ValueType, ValidType>::RunningSum(int window)
    : m_acc(0), m_count(0), m_w0(window / 2), m_first_elem(1), m_last_elem(0) {}

template<typename ValueType, typename ValidType>
void RunningSum<ValueType, ValidType>::process_index(int nval, int idx, const ValueType *values,
                                                     const ValidType *valid) {
    int start = idx - m_w0;
    int end   = idx + m_w0 - 1;

    // add elements between start and m_first_elem
    for (int i = start; i < m_first_elem; ++i) {
        if (0 <= i && i < nval) { add(i, values[i], valid[i]); }
    }

    // remove elements between m_first_elem and start
    for (int i = m_first_elem; i < start; ++i) {
        if (0 <= i && i < nval) { rmv(i, values[i], valid[i]); }
    }

    // add elements between m_last_elem and end
    for (int i = m_last_elem + 1; i < end + 1; ++i) {
        if (0 <= i && i < nval) { add(i, values[i], valid[i]); }
    }

    // remove elements between end and m_last_elem
    for (int i = end + 1; i < m_last_elem + 1; ++i) {
        if (0 <= i && i < nval) { rmv(i, values[i], valid[i]); }
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
