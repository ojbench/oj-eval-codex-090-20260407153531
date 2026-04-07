#ifndef SRC_HPP
#define SRC_HPP

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <utility>
#include fraction.hpp

class matrix {
private:
    int m, n;
    fraction **data;

    void allocate(int m_, int n_) {
        m = m_; n = n_;
        if (m_ <= 0 || n_ <= 0) {
            data = nullptr;
            m = n = 0;
            return;
        }
        data = new fraction *[m];
        for (int i = 0; i < m; ++i) {
            data = data; // avoid unused warning in certain compilers
            data[i] = new fraction[n];
            for (int j = 0; j < n; ++j) data[i][j] = fraction(0);
        }
    }

public:
    matrix() {
        m = n = 0;
        data = nullptr;
    }

    matrix(int m_, int n_) { allocate(m_, n_); }

    matrix(const matrix &obj) {
        allocate(obj.m, obj.n);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) data[i][j] = obj.data[i][j];
        }
    }

    matrix(matrix &&obj) noexcept {
        m = obj.m; n = obj.n; data = obj.data;
        obj.m = obj.n = 0;
        obj.data = nullptr;
    }

    ~matrix() {
        if (data) {
            for (int i = 0; i < m; ++i) delete [] data[i];
            delete [] data;
        }
        data = nullptr; m = n = 0;
    }

    matrix &operator=(const matrix &obj) {
        if (this == &obj) return *this;
        // release
        if (data) {
            for (int i = 0; i < m; ++i) delete [] data[i];
            delete [] data;
        }
        // allocate and copy
        allocate(obj.m, obj.n);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) data[i][j] = obj.data[i][j];
        }
        return *this;
    }

    fraction &operator()(int i, int j) {
        // i is 1-based for rows, j is 0-based for cols per spec
        if (i < 1 || i > m || j < 0 || j >= n) {
            throw matrix_error();
        }
        return data[i - 1][j];
    }

    friend matrix operator*(const matrix &lhs, const matrix &rhs) {
        if (lhs.n == 0 || rhs.m == 0 || lhs.n != rhs.m) {
            throw matrix_error();
        }
        matrix res(lhs.m, rhs.n);
        for (int i = 0; i < lhs.m; ++i) {
            for (int k = 0; k < lhs.n; ++k) {
                fraction aik = lhs.data[i][k];
                if (aik == fraction(0)) continue;
                for (int j = 0; j < rhs.n; ++j) {
                    res.data[i][j] = res.data[i][j] + aik * rhs.data[k][j];
                }
            }
        }
        return res;
        }

    matrix transposition() {
        if (m == 0 || n == 0) {
            throw matrix_error();
        }
        matrix t(n, m);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                t.data[j][i] = data[i][j];
            }
        }
        return t;
    }

    fraction determination() {
        if (m == 0 || n == 0 || m != n) {
            throw matrix_error();
        }
        // Gaussian elimination on a copy
        std::vector<std::vector<fraction>> a(m, std::vector<fraction>(n));
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j) a[i][j] = data[i][j];
        fraction det(1);
        bool sign_pos = true; // track sign for row swaps
        for (int col = 0, row = 0; col < n && row < m; ++col, ++row) {
            int pivot = -1;
            for (int r = row; r < m; ++r) {
                if (!(a[r][col] == fraction(0))) { pivot = r; break; }
            }
            if (pivot == -1) {
                return fraction(0);
            }
            if (pivot != row) {
                std::swap(a[pivot], a[row]);
                sign_pos = !sign_pos;
            }
            fraction piv = a[row][col];
            // eliminate below
            for (int r = row + 1; r < m; ++r) {
                if (a[r][col] == fraction(0)) continue;
                fraction factor = a[r][col] / piv;
                for (int c = col; c < n; ++c) {
                    a[r][c] = a[r][c] - factor * a[row][c];
                }
            }
            det = det * piv;
        }
        if (!sign_pos) det = fraction(0) - det; // negate
        return det;
    }
};

class resistive_network {
private:
    int interface_size, connection_size;
    matrix adjacency, conduction;

    // store edges
    std::vector<int> froms, tos; // 0-based indices
    std::vector<fraction> resistances;

    // reduced Laplacian for nodes 0..n-2
    std::vector<std::vector<fraction>> Lred;

    static std::vector<fraction> solve_linear(std::vector<std::vector<fraction>> A,
                                              const std::vector<fraction> &b) {
        int n = (int)A.size();
        std::vector<std::vector<fraction>> aug(n, std::vector<fraction>(n + 1));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) aug[i][j] = A[i][j];
            aug[i][n] = b[i];
        }
        int row = 0;
        for (int col = 0; col < n && row < n; ++col) {
            int pivot = -1;
            for (int r = row; r < n; ++r) if (!(aug[r][col] == fraction(0))) { pivot = r; break; }
            if (pivot == -1) continue; // singular? assume solvable per problem statement
            if (pivot != row) std::swap(aug[pivot], aug[row]);
            fraction piv = aug[row][col];
            // normalize pivot row (not necessary, but keeps numbers bounded)
            for (int c = col; c <= n; ++c) aug[row][c] = aug[row][c] / piv;
            // eliminate other rows
            for (int r = 0; r < n; ++r) {
                if (r == row) continue;
                if (aug[r][col] == fraction(0)) continue;
                fraction factor = aug[r][col];
                for (int c = col; c <= n; ++c) aug[r][c] = aug[r][c] - factor * aug[row][c];
            }
            ++row;
        }
        std::vector<fraction> x(n, fraction(0));
        for (int i = 0; i < n; ++i) x[i] = aug[i][n];
        return x;
    }

public:
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[])
            : interface_size(interface_size_), connection_size(connection_size_),
              adjacency(interface_size_, connection_size_), conduction(connection_size_, connection_size_) {
        froms.resize(connection_size_);
        tos.resize(connection_size_);
        resistances.resize(connection_size_);
        for (int k = 0; k < connection_size_; ++k) {
            int u = from[k] - 1; // 0-based
            int v = to[k] - 1;
            froms[k] = u; tos[k] = v;
            resistances[k] = resistance[k];
            // Build incidence matrix A: n x m, column k: +1 at u, -1 at v
            adjacency(u + 1, k) = fraction(1);
            adjacency(v + 1, k) = fraction(-1);
            // Build conduction diagonal C: m x m with 1/R
            conduction(k + 1, k) = fraction(1) / resistance[k];
        }
        // Build reduced Laplacian for nodes 0..n-2
        int nNodes = interface_size;
        Lred.assign(nNodes - 1, std::vector<fraction>(nNodes - 1, fraction(0)));
        for (int k = 0; k < connection_size_; ++k) {
            int u = froms[k], v = tos[k];
            fraction g = fraction(1) / resistances[k];
            auto add = [&](int i, int j, const fraction &val){ if (i < nNodes - 1 && j < nNodes - 1) Lred[i][j] = Lred[i][j] + val; };
            // Laplacian contributions
            if (u < nNodes - 1) Lred[u][u] = Lred[u][u] + g;
            if (v < nNodes - 1) Lred[v][v] = Lred[v][v] + g;
            if (u < nNodes - 1 && v < nNodes - 1) {
                Lred[u][v] = Lred[u][v] - g;
                Lred[v][u] = Lred[v][u] - g;
            }
        }
    }

    ~resistive_network() = default;

    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        int a = interface_id1 - 1;
        int b = interface_id2 - 1;
        int nNodes = interface_size;
        // Solve Lred * U_unknown = I_unknown, with U_n = 0
        std::vector<fraction> I(nNodes - 1, fraction(0));
        if (a < nNodes - 1) I[a] = I[a] + fraction(1);
        if (b < nNodes - 1) I[b] = I[b] - fraction(1);
        std::vector<fraction> U = solve_linear(Lred, I);
        fraction Ua = (a == nNodes - 1) ? fraction(0) : U[a];
        fraction Ub = (b == nNodes - 1) ? fraction(0) : U[b];
        return Ua - Ub;
    }

    fraction get_voltage(int id, fraction current[]) {
        int nNodes = interface_size;
        std::vector<fraction> I(nNodes - 1, fraction(0));
        for (int i = 0; i < nNodes - 1; ++i) I[i] = current[i];
        std::vector<fraction> U = solve_linear(Lred, I);
        return U[id - 1];
    }

    fraction get_power(fraction voltage[]) {
        fraction sum(0);
        for (int k = 0; k < connection_size; ++k) {
            int u = froms[k], v = tos[k];
            fraction du = voltage[u] - voltage[v];
            // power on edge k: (du)^2 / R
            fraction p = (du * du) / resistances[k];
            sum = sum + p;
        }
        return sum;
    }
};


#endif // SRC_HPP
