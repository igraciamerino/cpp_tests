#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <chrono>
#include <unordered_map>

/**
 * Class for representing a point. coordinate_type must be a numeric type.
 */
template <typename coordinate_type, size_t dimensions>
class point
{
public:
    point(std::array<coordinate_type, dimensions> c) : coords_(c) {}
    point(std::initializer_list<coordinate_type> list)
    {
        size_t n = std::min(dimensions, list.size());
        std::copy_n(list.begin(), n, coords_.begin());
    }
    /**
     * Returns the coordinate in the given dimension.
     *
     * @param index dimension index (zero based)
     * @return coordinate in the given dimension
     */
    coordinate_type get(size_t index) const
    {
        return coords_[index];
    }

    void set(std::array<coordinate_type, dimensions> coord_values)
    {
        coords_ = coord_values;
    }
    /**
     * Returns the distance squared from this point to another
     * point.
     *
     * @param pt another point
     * @return distance squared from this point to the other point
     */
    double distance(const point &pt) const
    {
        double dist = 0;
        for (size_t i = 0; i < dimensions; ++i)
        {
            double d = get(i) - pt.get(i);
            dist += d * d;
        }
        return dist;
    }

private:
    std::array<coordinate_type, dimensions> coords_;
};

template <typename coordinate_type, size_t dimensions>
std::ostream &operator<<(std::ostream &out, const point<coordinate_type, dimensions> &pt)
{
    out << '(';
    for (size_t i = 0; i < dimensions; ++i)
    {
        if (i > 0)
            out << ", ";
        out << pt.get(i);
    }
    out << ')';
    return out;
}

/**
 * C++ k-d tree implementation, based on the C version at rosettacode.org.
 */
template <typename coordinate_type, size_t dimensions>
class kdtree
{
public:
    typedef point<coordinate_type, dimensions> point_type;

private:
    struct node
    {
        node(const point_type &pt) : point_(pt), left_(nullptr), right_(nullptr) {}
        coordinate_type get(size_t index) const
        {
            return point_.get(index);
        }
        double distance(const point_type &pt) const
        {
            return point_.distance(pt);
        }
        point_type point_;
        size_t node_index_;
        node *left_;
        node *right_;
    };
    node *root_ = nullptr;
    node *best_ = nullptr;
    double best_dist_sq = 0;
    size_t visited_ = 0;
    std::vector<node> nodes_;

    struct node_cmp
    {
        node_cmp(size_t index) : index_(index) {}
        bool operator()(const node &n1, const node &n2) const
        {
            return n1.point_.get(index_) < n2.point_.get(index_);
        }
        size_t index_;
    };

    node *make_tree(size_t begin, size_t end, size_t index)
    {
        if (end <= begin)
            return nullptr;
        size_t n = begin + (end - begin) / 2;
        auto i = nodes_.begin();
        std::nth_element(i + begin, i + n, i + end, node_cmp(index));
        index = (index + 1) % dimensions;
        nodes_[n].left_ = make_tree(begin, n, index);
        nodes_[n].right_ = make_tree(n + 1, end, index);
        // std::cout << n << std::endl;
        return &nodes_[n];
    }

    void nearest(node *root, const point_type &point, size_t index)
    {
        if (root == nullptr)
            return;
        ++visited_;
        double d = root->distance(point);
        if (best_ == nullptr || d < best_dist_sq)
        {
            best_dist_sq = d;
            best_ = root;
        }
        if (best_dist_sq == 0)
            return;
        double dx = root->get(index) - point.get(index);
        index = (index + 1) % dimensions;
        nearest(dx > 0 ? root->left_ : root->right_, point, index);
        if (dx * dx >= best_dist_sq)
            return;
        nearest(dx > 0 ? root->right_ : root->left_, point, index);
    }

public:
    kdtree(const kdtree &) = delete;
    kdtree &operator=(const kdtree &) = delete;
    /**
     * Constructor taking a pair of iterators. Adds each
     * point in the range [begin, end) to the tree.
     *
     * @param begin start of range
     * @param end end of range
     */
    template <typename iterator>
    kdtree(iterator begin, iterator end) : nodes_(begin, end)
    {
        root_ = make_tree(0, nodes_.size(), 0);

        // My functionality to add node indices
        for (size_t i = 0; i < nodes_.size(); i++)
        {
            nodes_.at(i).node_index_ = i;
        }
    }

    /**
     * Constructor taking a function object that generates
     * points. The function object will be called n times
     * to populate the tree.
     *
     * @param f function that returns a point
     * @param n number of points to add
     */
    template <typename func>
    kdtree(func &&f, size_t n)
    {
        nodes_.reserve(n);
        for (size_t i = 0; i < n; ++i)
            nodes_.push_back(f());
        root_ = make_tree(0, nodes_.size(), 0);
    }

    /**
     * Returns true if the tree is empty, false otherwise.
     */
    bool empty() const { return nodes_.empty(); }

    /**
     * Returns the number of nodes visited by the last call
     * to nearest().
     */
    size_t visited() const { return visited_; }

    /**
     * Returns the distance between the input point and return value
     * from the last call to nearest().
     */
    double distance() const { return std::sqrt(best_dist_sq); }

    double distance_sq() const { return best_dist_sq; }

    /**
     * Finds the nearest point in the tree to the given point.
     * It is not valid to call this function if the tree is empty.
     *
     * @param pt a point
     * @return the nearest point in the tree to the given point
     */

    // Old nearest fcn
    const point_type &nearest(const point_type &pt)
    {
        if (root_ == nullptr)
            throw std::logic_error("tree is empty");
        best_ = nullptr;
        visited_ = 0;
        best_dist_sq = 0;
        nearest(root_, pt, 0);
        return best_->point_;
    }

    // My nearest fcn
    const std::pair<point_type, size_t> nearest_with_index(const point_type &pt)
    {
        if (root_ == nullptr)
            throw std::logic_error("tree is empty");
        best_ = nullptr;
        visited_ = 0;
        best_dist_sq = 0;
        nearest(root_, pt, 0);
        return std::make_pair(best_->point_, best_->node_index_);
    }

    const point_type get_point_from_ind(size_t ind)
    {
        return nodes_.at(ind).point_;
    }
};

void test_wikipedia()
{
    typedef point<int, 2> point2d;
    typedef kdtree<int, 2> tree2d;

    point2d points[] = {{2, 3}, {5, 4}, {9, 6}, {4, 7}, {8, 1}, {7, 2}};

    tree2d tree(std::begin(points), std::end(points));
    point2d n = tree.nearest({9, 2});

    std::cout << "Wikipedia example data:\n";
    std::cout << "nearest point: " << n << '\n';
    std::cout << "distance: " << tree.distance() << '\n';
    std::cout << "nodes visited: " << tree.visited() << '\n';
}

typedef point<double, 3> point3d;
typedef point<double, 2> point2d;
typedef kdtree<double, 3> tree3d;
typedef kdtree<double, 2> tree2d;

struct random_point_generator
{
    random_point_generator(double min, double max)
        : engine_(std::random_device()()), distribution_(min, max) {}

    point3d operator()()
    {
        double x = distribution_(engine_);
        double y = distribution_(engine_);
        double z = distribution_(engine_);
        return point3d({x, y, z});
    }

    std::mt19937 engine_;
    std::uniform_real_distribution<double> distribution_;
};

struct random_point_generator_gaussian
{
    random_point_generator_gaussian(double mu, double std)
        : engine_(std::random_device()()), distribution_(mu, std) {}

    point2d operator()()
    {
        double x = distribution_(engine_);
        double y = distribution_(engine_);
        return point2d({x, y});
    }

    std::default_random_engine engine_;
    std::normal_distribution<double> distribution_;
};

void test_random(size_t count)
{
    random_point_generator rpg(0, 1);
    tree3d tree(rpg, count);
    point3d pt(rpg());
    point3d n = tree.nearest(pt);

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < 1'000'000; i++)
    {
        point3d pt(rpg());
        point3d n = tree.nearest(pt);
    }

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Search time = " << elapsed.count() << '\n';

    std::cout << "Random data (" << count << " points):\n";
    std::cout << "point: " << pt << '\n';
    std::cout << "nearest point: " << n << '\n';
    std::cout << "distance: " << tree.distance() << '\n';
    std::cout << "nodes visited: " << tree.visited() << '\n';
}
struct MyArrayHasher
{
    std::size_t operator()(const std::array<double, 2> &a) const
    {
        std::size_t h = 0;

        for (auto e : a)
        {
            h ^= std::hash<int>{}(e) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

std::pair<bool, size_t> check_assignment_fcn(std::vector<size_t> &J_assigned, size_t N2)
{
    std::vector<bool> has_some_i_assigned(N2, false);
    size_t counter = 0;
    size_t j;

    for (size_t i = 0; i < J_assigned.size(); i++)
    {
        if (!has_some_i_assigned.at(J_assigned.at(i)))
        {
            has_some_i_assigned.at(J_assigned.at(i)) = true;
            counter++;
        }
    }

    if (counter == N2)
    {
        return std::pair(true, N2);
    }
    else
    {
        return std::pair(false, counter);
    }
}

double dist2_sq_2D_fcn(std::array<double, 2> x, std::array<double, 2> y)
{
    return (x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]);
}

/*
double cluster_fcn_slow(std::vector<point2d> &W1, tree2d &tree, std::vector<std::array<double, 2>> &W2_recentered, bool &postprocess)
{
    size_t N2 = W2_recentered.size();
    point2d x1 = {0.0, 0.0};
    point2d x2 = {0.0, 0.0};
    std::array<double, 2> x1_arr = {0.0, 0.0};
    std::array<double, 2> x2_arr = {0.0, 0.0};
    double W_dist = 0.0;
    double dist_max_sq = 1'000'000.0;
    std::vector<double> weights2(N2, 0.0);
    std::unordered_map<std::array<double, 2>, size_t, MyArrayHasher> J_index;
    size_t j_index;
    double dist_sq;
    double dist_sq_current;
    std::unordered_map<std::array<double, 2>, std::array<double, 2>, MyArrayHasher> J_assigned;
    std::pair<bool, size_t> check_assignment(true, N2);

    // Compute W_dist
    for (size_t i = 0; i < W1.size(); i++)
    {
        x1 = W1.at(i);
        x1_arr = {x1.get(0), x1.get(1)};
        x2 = tree.nearest(W1.at(i));
        x2_arr = {x2.get(0), x2.get(1)};
        J_assigned[x1_arr] = x2_arr;
        j_index = J_index[x2_arr];
        weights2.at(j_index) += 1.0;

        if (postprocess)
        {
            W2_recentered.at(j_index)[0] += x1_arr[0];
            W2_recentered.at(j_index)[1] += x1_arr[1];
        }
        else
        {
            W_dist += tree.distance_sq();
        }
    }

    // Check if all j have an i assigned to them
    check_assignment = check_assignment_fcn(J_index, N2);
    if (!check_assignment.first)
    {
        std::cout << "Error: not all j have an i assigned to them..." << std::endl;
        std::cout << "Assigned: " << check_assignment.second << " of " << N2 << std::endl;
    }

    // Re-center W2
    if (postprocess)
    {
        for (size_t j = 0; j < N2; j++)
        {
            W2_recentered.at(j)[0] = W2_recentered.at(j)[0] / weights2.at(j);
            W2_recentered.at(j)[1] = W2_recentered.at(j)[1] / weights2.at(j);
        }
        for (size_t i = 0; i < W1.size(); i++)
        {
            x1 = W1.at(i);
            x2_arr = W2_recentered.at(J_index[{x1.get(0), x1.get(1)}]);
            x2.set(x2_arr);
            dist_sq = x1.distance(x2);
            W_dist += dist_sq;
        }
    }

    W_dist = sqrt(W_dist / W1.size());

    return W_dist;
}
*/

double cluster_fcn(std::vector<std::array<double, 2>> &W1, std::vector<std::array<double, 2>> &W2, std::vector<double> &weights2, size_t p, bool &postprocess)
{
    size_t N1 = W1.size();
    size_t N2 = W2.size();
    point2d x2 = {0.0, 0.0};
    std::array<double, 2> x1 = {0.0, 0.0};
    std::array<double, 2> x2_arr = {0.0, 0.0};
    double W_dist = 0.0;
    size_t j_index;
    double dist_sq;
    std::vector<size_t> J_assigned(N1, 0);
    std::pair<bool, size_t> check_assignment(true, N2);

    // Initialize weights2 to zero
    std::fill(weights2.begin(), weights2.end(), 0.0);

    // Generate tree for W2
    tree2d tree(std::begin(W2), std::end(W2));
    std::cout << "Tree constructed... " << '\n';

    // Compute W_dist
    for (size_t i = 0; i < N1; i++)
    {
        x1 = W1.at(i);
        j_index = tree.nearest_with_index(W1.at(i)).second;
        // x2_arr = {x2.get(0), x2.get(1)};
        J_assigned.at(i) = j_index;
        weights2.at(j_index) += 1.0;

        if (postprocess)
        {
            W2.at(j_index)[0] += x1[0];
            W2.at(j_index)[1] += x1[1];
        }
        else
        {
            W_dist += tree.distance_sq();
        }
    }

    // Check if all j have an i assigned to them
    check_assignment = check_assignment_fcn(J_assigned, N2);
    if (!check_assignment.first)
    {
        std::cout << "Error: not all j have an i assigned to them..." << std::endl;
        std::cout << "Assigned: " << check_assignment.second << " of " << N2 << std::endl;
    }

    // Re-center W2
    if (postprocess)
    {
        for (size_t j = 0; j < N2; j++)
        {
            W2.at(j)[0] = W2.at(j)[0] / weights2.at(j);
            W2.at(j)[1] = W2.at(j)[1] / weights2.at(j);
            weights2.at(j) = weights2.at(j) / N1;
        }
        for (size_t i = 0; i < N1; i++)
        {
            x1 = W1.at(i);
            x2_arr = W2.at(J_assigned.at(i));
            dist_sq = dist2_sq_2D_fcn(x1, x2_arr);
            if (p == 2)
            {
                W_dist += dist_sq;
            }
            else if (p == 1)
            {
                W_dist += sqrt(dist_sq);
            }
        }
    }
    if (p == 2)
    {
        W_dist = sqrt(W_dist / N1);
    }
    else if (p == 1)
    {
        W_dist = W_dist / N1;
    }

    return W_dist;
}

double lloyd_fcn(std::vector<std::array<double, 2>> &W1, std::vector<std::array<double, 2>> &W2, std::vector<double> &weights2, size_t &p, size_t &N_iterations_lloyd)
{
    bool postprocess = true;

    double W_dist;

    std::cout << "Starting Lloyd's algorithm..." << std::endl;
    for (size_t i = 0; i < N_iterations_lloyd; i++)
    {
        W_dist = cluster_fcn(W1, W2, weights2, p, postprocess);
        std::cout << "Iteration " << i << " of " << N_iterations_lloyd << ". Wass dist = " << W_dist << std::endl;
    }

    return W_dist;
}

std::array<double, 2> simulate_fcn(double (&A)[4][4], double (&D)[4][2], size_t &T, std::default_random_engine &generator, std::normal_distribution<double> &distribution)
{
    std::array<double, 4> x = {0.0, 0.0, 0.0, 0.0};
    // std::array<double, 2> x_proj;
    std::array<double, 2> w;

    for (size_t t = 0; t < T - 1; t++)
    {
        w[0] = distribution(generator);
        w[1] = distribution(generator);
        for (size_t i = 0; i < 4; i++)
        {
            x[i] = A[i][0] * x[0] + A[i][1] * x[1] + A[i][2] * x[2] + A[i][3] * x[3] + D[i][0] * w[0] + D[i][1] * w[1];
        }
    }
    // x_proj[0] = x[0];
    // x_proj[1] = x[1];

    return {x[0], x[1]}; // x_proj;
}

void generate_samples(std::vector<std::array<double, 2>> &W1, std::vector<std::array<double, 2>> &W2, size_t &N1, size_t &N2, double (&A)[4][4], double (&D)[4][2], size_t &T, std::default_random_engine &generator, std::normal_distribution<double> &distribution)
{
    std::array<double, 2> x;

    // Generate sets of samples
    for (size_t i = 0; i < N1; i++)
    {
        x = simulate_fcn(A, D, T, generator, distribution);
        W1.push_back(x);
    }

    for (size_t j = 0; j < N2; j++)
    {
        x = simulate_fcn(A, D, T, generator, distribution);
        W2.push_back(x);
    }
    std::cout << "W1, W2 defined... " << '\n';
}

int main()
{
    try
    {
        // Input params
        size_t N1 = 10'000'000;
        size_t N2 = 10'000;
        size_t p = 1;
        size_t N_iterations_lloyd = 10;
        size_t N = 1'000'000'000;
        size_t T = 40;
        double mu = 0.0;
        double std = 1.0;
        

        std::array<double, 2> x;
        std::array<double, 2> y;
        std::default_random_engine generator;
        std::normal_distribution<double> distribution(mu, std);
        double W_dist;
        std::pair<point2d, size_t> nearest({0.0, 0.0}, 0);
        std::vector<std::array<double, 2>> W1;
        std::vector<std::array<double, 2>> W2;
        
        std::vector<double> weights2(N2, 0.0);
        
        double dist_sq;
        double A[4][4] = {{0.945883, 0, 0.0765802, 0},
                      {0, 0.945883, 0, 0.0765802},
                      {-1.082331, 0, 0.531604, 0},
                      {0, -1.082331, 0, 0.531604}};
        // double D[4][2] = {{0, 0}, {0, 0}, {0.04472136, 0.031622776}, {0.031622776, 0.04472136}}; // TS
        double D[4][2] = {{0, 0}, {0, 0}, {0.031622777, 0.0}, {0.0, 0.031622777}}; // Gaussian

        // Generate the data
        generate_samples(W1, W2, N1, N2, A, D, T, generator, distribution);

        // Call Lloyd's algorithm
        auto start = std::chrono::system_clock::now();
        W_dist = lloyd_fcn(W1, W2, weights2, p, N_iterations_lloyd);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cout << "Cluster time = " << elapsed.count() << " ms." << '\n';
        std::cout << "W_dist = " << W_dist << '\n';

        // Construct the tree
        tree2d tree(std::begin(W2), std::end(W2));
        std::cout << "Tree constructed... " << '\n';

        // Clustering N
        // Initialize weights2 to zero
        std::fill(weights2.begin(), weights2.end(), 0.0);
        std::cout << "Clustering N... " << '\n';
        auto start1 = std::chrono::system_clock::now();
        for (size_t i = 0; i < N; i++)
        {
            x = simulate_fcn(A, D, T, generator, distribution);
            nearest = tree.nearest_with_index(x);
            weights2.at(nearest.second) += 1.0;
            y[0] = nearest.first.get(0);
            y[1] = nearest.first.get(1);
            dist_sq = dist2_sq_2D_fcn(x, y);
            if (p == 2)
            {
                W_dist += dist_sq;
            }
            else if (p == 1)
            {
                W_dist += sqrt(dist_sq);
            }
        }

        for (size_t j = 0; j < N2; j++)
        {
            weights2.at(j) = weights2.at(j) / N;
        }

        if (p == 2)
        {
            W_dist = sqrt(W_dist / N);
        }
        else if (p == 1)
        {
            W_dist = W_dist / N;
        }

        auto end1 = std::chrono::system_clock::now();
        auto elapsed1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1);

        std::cout << "Cluster time = " << elapsed1.count() << " ms." << '\n';
        std::cout << "W_dist = " << W_dist << '\n';
        
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }
    return 0;
}
