#ifndef FisherRepo_H
#define FisherRepo_H

#include <vector>

class FisherRepo {
  private:
    std::vector<std::vector<double>> _ps;
  public:
    FisherRepo(const std::size_t allocationsPerArm);
    ~FisherRepo(){}
    const double getP(const std::size_t a, const std::size_t c); // should be unsigned int, but R does not have a corresponding type, so avoiding casting here
};

#endif
