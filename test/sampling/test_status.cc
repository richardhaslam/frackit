#include <stdexcept>
#include <iostream>

#include <frackit/common/id.hh>
#include <frackit/sampling/status.hh>

//! test sampling status
int main()
{
    using namespace Frackit;
    SamplingStatus status;

    status.setTargetCount(Id(0), 2);
    status.increaseCounter(Id(0));
    status.increaseCounter(Id(0));

    if (!status.finished())
        throw std::runtime_error("Test 1 failed");
    if (!status.finished(Id(0)))
        throw std::runtime_error("Test 2 failed");

    status.resetCounter(Id(0));
    if (status.finished())
        throw std::runtime_error("Test 3 failed");
    if (status.finished(Id(0)))
        throw std::runtime_error("Test 4 failed");

    status.setTargetCount(Id(0), 2);
    status.setTargetCount(Id(1), 2);

    status.increaseCounter(Id(0));
    status.increaseCounter(Id(0));
    if (!status.finished(Id(0)))
        throw std::runtime_error("Test 5 failed");
    if (status.finished(Id(1)))
        throw std::runtime_error("Test 6 failed");
    if (status.finished())
        throw std::runtime_error("Test 7 failed");

    status.resetCounters();
    if (status.finished(Id(0)))
        throw std::runtime_error("Test 8 failed");
    if (status.finished(Id(1)))
        throw std::runtime_error("Test 9 failed");
    if (status.finished())
        throw std::runtime_error("Test 10 failed");

    // this should print a warning
    status.increaseCounter(Id(0));
    status.increaseCounter(Id(0));
    std::cout << std::endl;
    std::cout << "If the test was successful, a warning should be printed below" << std::endl;
    status.setTargetCount(Id(0), 1);
    std::cout << std::endl;

    std::cout << "All tests passed" << std::endl;
    return 0;
}
