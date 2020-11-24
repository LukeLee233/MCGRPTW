#include "RNG.h"

unsigned long RNG::Randint(int low, int high)
{
    return (unsigned long) (gsl_rng_uniform(gr) * (high - low) + low
        + 0.5); // +0.5 para que al truncar tanto low como high est��n incluidos
}

double RNG::Randfloat(double low, double high)
{
    return (gsl_rng_uniform(gr) * (high - low) + low);
}

void RNG::change(long int seed)
{
    gsl_rng_set(gr, seed);
    _seed = seed;
}

namespace sample{

random_device rd_;
mt19937 gen_(rd_());
uniform_real_distribution<double> dis_(0.0, 1.0);


int uniform_sample(const vector<double>& probs, mt19937& gen, uniform_real_distribution<double>& dis){
    double p = dis(gen);
    int i = 0;
    while( i < probs.size() && (p -= probs[i]) > 0)
        ++i;
    return i;
}

}