#ifndef SRC_HPP
#define SRC_HPP

#include <vector>
#include <algorithm>
#include <utility>
#define STRIFY(x) #x
#include STRIFY(fraction.hpp)

class matrix {
private:
    int m, n;
    fraction **data;

    void allocate(int m_, int n_) {
        m = m_; n = n_;
        if (m_ <= 0 || n_ <= 0) { m = n = 0; data = nullptr; return; }
        data = new fraction*[m];
        for (int i = 0; i < m; ++i) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; ++j) data[i][j] = fraction(0);
        }
    }

public:
    matrix(): m(0), n(0), data(nullptr) {}
    matrix(int m_, int n_) { allocate(m_, n_); }
    matrix(const matrix &obj) { allocate(obj.m, obj.n); for (int i=0;i<m;++i) for(int j=0;j<n;++j) data[i][j]=obj.data[i][j]; }
    matrix(matrix &&obj) noexcept : m(obj.m), n(obj.n), data(obj.data) { obj.m=obj.n=0; obj.data=nullptr; }
    ~matrix(){ if (data){ for(int i=0;i<m;++i) delete[] data[i]; delete[] data; } }

    matrix &operator=(const matrix &obj){
        if (this==&obj) return *this;
        if (data){ for(int i=0;i<m;++i) delete[] data[i]; delete[] data; }
        allocate(obj.m, obj.n);
        for (int i=0;i<m;++i) for(int j=0;j<n;++j) data[i][j]=obj.data[i][j];
        return *this;
    }

    fraction &operator()(int i, int j){
        if (i < 1 || i > m || j < 0 || j >= n) throw matrix_error();
        return data[i-1][j];
    }

    friend matrix operator*(const matrix &lhs, const matrix &rhs){
        if (lhs.n == 0 || rhs.m == 0 || lhs.n != rhs.m) throw matrix_error();
        matrix res(lhs.m, rhs.n);
        for (int i=0;i<lhs.m;++i){
            for (int k=0;k<lhs.n;++k){
                fraction aik = lhs.data[i][k];
                if (aik == fraction(0)) continue;
                for (int j=0;j<rhs.n;++j){ res.data[i][j] = res.data[i][j] + aik * rhs.data[k][j]; }
            }
        }
        return res;
    }

    matrix transposition(){
        if (m==0 || n==0) throw matrix_error();
        matrix t(n, m);
        for (int i=0;i<m;++i) for(int j=0;j<n;++j) t.data[j][i] = data[i][j];
        return t;
    }

    fraction determination(){
        if (m==0 || n==0 || m!=n) throw matrix_error();
        std::vector<std::vector<fraction>> a(m, std::vector<fraction>(n));
        for (int i=0;i<m;++i) for(int j=0;j<n;++j) a[i][j]=data[i][j];
        fraction det(1);
        bool sign_pos = true;
        for (int col=0, row=0; col<n && row<m; ++col, ++row){
            int pivot=-1; for (int r=row; r<m; ++r){ if (!(a[r][col]==fraction(0))){ pivot=r; break; } }
            if (pivot==-1) return fraction(0);
            if (pivot!=row){ std::swap(a[pivot], a[row]); sign_pos = !sign_pos; }
            fraction piv = a[row][col];
            for (int r=row+1; r<m; ++r){ if (a[r][col]==fraction(0)) continue; fraction factor = a[r][col] / piv; for (int c=col; c<n; ++c) a[r][c] = a[r][c] - factor * a[row][c]; }
            det = det * piv;
        }
        if (!sign_pos) det = fraction(0) - det;
        return det;
    }
};

class resistive_network {
private:
    int interface_size, connection_size;
    matrix adjacency, conduction;
    std::vector<int> froms, tos;
    std::vector<fraction> resistances;
    std::vector<std::vector<fraction>> Lred; // reduced Laplacian (exclude last node)

    static std::vector<fraction> solve_linear(std::vector<std::vector<fraction>> A, const std::vector<fraction> &b){
        int n = (int)A.size();
        std::vector<std::vector<fraction>> aug(n, std::vector<fraction>(n+1));
        for (int i=0;i<n;++i){ for (int j=0;j<n;++j) aug[i][j]=A[i][j]; aug[i][n]=b[i]; }
        int row=0;
        for (int col=0; col<n && row<n; ++col){
            int pivot=-1; for (int r=row; r<n; ++r) if (!(aug[r][col]==fraction(0))){ pivot=r; break; }
            if (pivot==-1) continue; // assume solvable
            if (pivot!=row) std::swap(aug[pivot], aug[row]);
            fraction piv = aug[row][col];
            for (int c=col; c<=n; ++c) aug[row][c] = aug[row][c] / piv;
            for (int r=0; r<n; ++r){ if (r==row) continue; if (aug[r][col]==fraction(0)) continue; fraction factor = aug[r][col]; for (int c=col; c<=n; ++c) aug[r][c] = aug[r][c] - factor * aug[row][c]; }
            ++row;
        }
        std::vector<fraction> x(n, fraction(0)); for (int i=0;i<n;++i) x[i]=aug[i][n]; return x;
    }

public:
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[])
        : interface_size(interface_size_), connection_size(connection_size_), adjacency(interface_size_, connection_size_), conduction(connection_size_, connection_size_){
        froms.resize(connection_size_); tos.resize(connection_size_); resistances.resize(connection_size_);
        for (int k=0;k<connection_size_;++k){
            int u = from[k]-1; int v = to[k]-1; froms[k]=u; tos[k]=v; resistances[k]=resistance[k];
            adjacency(u+1, k) = fraction(1); adjacency(v+1, k) = fraction(-1);
            conduction(k+1, k) = fraction(1) / resistance[k];
        }
        int n = interface_size;
        Lred.assign(n-1, std::vector<fraction>(n-1, fraction(0)));
        for (int k=0;k<connection_size_;++k){
            int u = froms[k], v = tos[k]; fraction g = fraction(1) / resistances[k];
            if (u < n-1) Lred[u][u] = Lred[u][u] + g;
            if (v < n-1) Lred[v][v] = Lred[v][v] + g;
            if (u < n-1 && v < n-1){ Lred[u][v] = Lred[u][v] - g; Lred[v][u] = Lred[v][u] - g; }
        }
    }

    ~resistive_network() = default;

    fraction get_equivalent_resistance(int interface_id1, int interface_id2){
        int a = interface_id1-1, b = interface_id2-1; int n = interface_size;
        std::vector<fraction> I(n-1, fraction(0));
        if (a < n-1) I[a] = I[a] + fraction(1);
        if (b < n-1) I[b] = I[b] - fraction(1);
        std::vector<fraction> U = solve_linear(Lred, I);
        fraction Ua = (a == n-1) ? fraction(0) : U[a];
        fraction Ub = (b == n-1) ? fraction(0) : U[b];
        return Ua - Ub; // total injected current is 1
    }

    fraction get_voltage(int id, fraction current[]){
        int n = interface_size; std::vector<fraction> I(n-1, fraction(0));
        for (int i=0;i<n-1;++i) I[i] = current[i];
        std::vector<fraction> U = solve_linear(Lred, I);
        return U[id-1];
    }

    fraction get_power(fraction voltage[]){
        fraction sum(0);
        for (int k=0;k<connection_size;++k){
            int u = froms[k], v = tos[k]; fraction du = voltage[u] - voltage[v];
            sum = sum + (du * du) / resistances[k];
        }
        return sum;
    }
};

#endif // SRC_HPP
