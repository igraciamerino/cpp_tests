#include <iostream>
#include <string>
#include <random>
#include <chrono>

template <size_t N1, size_t N2, size_t d>
double cluster_fcn(double (&W1)[N1][d], double (&W2)[N2][d], double (&weight1)[N1], bool &postprocess);

template <size_t N1>
bool check_assignment_fcn(int (&J_assigned)[N1], size_t N2);

template <size_t N1, size_t N2, size_t d>
double lloyd_fcn(double (&W1)[N1][d], double (&W2)[N2][d], double (&weights1)[N1], size_t &N_iterations_lloyd);

template <size_t d>
double dist2_sq_fcn(double (&x)[d], double (&y)[d]);

int main()
{
    const int N1 = 100000;
    const int N2 = 200;
    const int d = 4;
    size_t N_iterations_lloyd = 10;
    double std_w = 0.1;
    bool postprocess = true;
    int N_clustering_steps = 20;
    static double W1[N1][d] = {};
    double W2[N2][d] = {};
    double W_dist = 0.0;
    static double weights1[N1];

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, std_w);

    // Generate sets of samples and give uniform weights to W1
    for (size_t i = 0; i < N1; i++)
    {
        W1[i][0] = distribution(generator);
        W1[i][1] = distribution(generator);
        weights1[i] = 1.0 / N1;
    }

    for (size_t i = 0; i < N2; i++)
    {
        W2[i][0] = distribution(generator);
        W2[i][1] = distribution(generator);
    }

    auto start = std::chrono::system_clock::now();

    W_dist = lloyd_fcn(W1, W2, weights1, N_iterations_lloyd);

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << elapsed.count() << '\n';

    return 0;
}

template <size_t N1, size_t N2, size_t d>
double cluster_fcn(double (&W1)[N1][d], double (&W2)[N2][d], double (&weights1)[N1], bool &postprocess)
{
    double W2_recentered[N2][d] = {};
    double x1[d];
    double x2[d];
    double W_dist = 0.0;
    double dist_max_sq = 1000.0;
    double weights2[N2] = {};
    double dist_sq;
    double dist_sq_current;
    static int J_assigned[N1] = {};
    bool check_assignment = true;

    for (size_t i = 0; i < N2; i++)
    {
        W2_recentered[i][0] = 0.0;
        W2_recentered[i][1] = 0.0;
    }

    for (size_t j = 0; j < N2; j++)
    {
        weights2[j] = 0.0;
    }

    // Compute W_dist
    for (size_t i = 0; i < N1; i++)
    {
        x1[0] = W1[i][0];
        x1[1] = W1[i][1];
        dist_sq = dist_max_sq;
        for (size_t j = 0; j < N2; j++)
        {
            x2[0] = W2[j][0];
            x2[1] = W2[j][1];
            dist_sq_current = dist2_sq_fcn(x1, x2);
            if (dist_sq_current < dist_sq)
            {
                dist_sq = dist_sq_current;
                J_assigned[i] = j;
            }
        }
        weights2[J_assigned[i]] += weights1[i];
        if (postprocess)
        {
            W2_recentered[J_assigned[i]][0] += weights1[i] * x1[0];
            W2_recentered[J_assigned[i]][1] += weights1[i] * x1[1];
        }
        else
        {
            W_dist += weights1[i] * dist_sq;
        }
    }

    // Check if all j have an i assigned to them
    check_assignment = check_assignment_fcn(J_assigned, N2);
    if (!check_assignment)
    {
        std::cout << "Error: not all j have an i assigned to them..." << std::endl;
    }

    // Re-center W2
    if (postprocess)
    {
        for (size_t j = 0; j < N2; j++)
        {
            W2[j][0] = W2_recentered[j][0] / weights2[j];
            W2[j][1] = W2_recentered[j][1] / weights2[j];
        }
        for (size_t i = 0; i < N1; i++)
        {
            x1[0] = W1[i][0];
            x1[1] = W1[i][1];
            x2[0] = W2[J_assigned[i]][0];
            x2[1] = W2[J_assigned[i]][1];
            dist_sq = dist2_sq_fcn(x1, x2);
            W_dist += weights1[i] * dist_sq;
        }
    }

    W_dist = sqrt(W_dist);

    return W_dist;
}

template <size_t N1, size_t N2, size_t d>
double lloyd_fcn(double (&W1)[N1][d], double (&W2)[N2][d], double (&weights1)[N1], size_t &N_iterations_lloyd)
{
    double W_dist = 0.0;
    bool postprocess = true;
    std::cout << "Starting Lloyd's algorithm..." << W_dist << std::endl;
    for (size_t i = 0; i < N_iterations_lloyd; i++)
    {
        W_dist = cluster_fcn(W1, W2, weights1, postprocess);
        std::cout << "Iteration " << i << " of " << N_iterations_lloyd << ". Wass dist = " << W_dist << std::endl;
    }

    return W_dist;
}

template <size_t N1>
bool check_assignment_fcn(int (&J_assigned)[N1], size_t N2)
{
    bool has_some_i_assigned[N2];
    size_t counter = 0;
    size_t j;

    for (size_t j = 0; j < N2; j++)
    {
        has_some_i_assigned[j] = false;
    }

    for (size_t i = 0; i < N1; i++)
    {
        j = J_assigned[i];
        if (!has_some_i_assigned[j])
        {
            has_some_i_assigned[j] = true;
            counter++;
        }
    }
    
    if (counter == N2)
    {
        return true;
    }
    else
    {
        return false;
    }
}

template <size_t d>
double dist2_sq_fcn(double (&x)[d], double (&y)[d])
{   
    if (d == 2)
    {
        return (x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]);
    }
    else if (d == 3)
    {
        return (x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]) + (x[2] - y[2])*(x[2] - y[2]);
    }
    else if (d == 4)
    {
        return (x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]) + (x[2] - y[2])*(x[2] - y[2]) + (x[3] - y[3])*(x[3] - y[3]);
    }
    return 0.0;
}
