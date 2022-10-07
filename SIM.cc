#include <map>
#include <cstdio>
#include <vector>
#include <random>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;

template<typename T>
void pvv(const char *desc, ostream &out, const vector<vector<T>> &vv) {
    out << desc << "\n";
    out << "---------------------\n";
    for (const auto& v: vv) {
        for (auto i: v) out << i << "\t";
        out << "\n";
    }
    out << "---------------------\n\n";
}

template<typename T>
void pv(const vector<T> &v) {
    cout << "{ ";
    for (auto i: v) cout << i << " ";
    cout << "}";
}

struct SimInfo {
    SimInfo(int m, int n) : m(m), n(n), N(), B(), Y() {
        N.resize(m);
        B.resize(m);
        Y.resize(m);

        for (auto& v: N) v.resize(n, 0);
        for (auto& v: Y) v.resize(n, -1);
    }

    SimInfo(const SimInfo &rhs)
        : N(rhs.N), B(rhs.B), Y(rhs.Y) { }

    SimInfo(SimInfo &&rhs) noexcept
        : N(move(rhs.N)), B(move(rhs.B)), Y(move(rhs.Y)) { }

    friend ostream& operator<<(ostream& out, const SimInfo& sim_info) {
        out << "m = " << sim_info.m << ", n = " << sim_info.n << "\n\n";
        pvv("SimInfo#N", out, sim_info.N);
        pvv("SimInfo#B", out, sim_info.B);
        pvv("SimInfo#Y", out, sim_info.Y);
        return out;
    }

    int m;
    int n;
    vector<vector<int>> N;
    vector<vector<int>> B;
    vector<vector<double>> Y;
};

struct RNG {
    RNG(int seed) : mt(seed), u1(0, 1) { }

    double random(double y) {
        return u1(mt) * y;
    }

    mt19937 mt;
    uniform_real_distribution<double> u1;
};

struct SimAlgo {
    SimAlgo(const SimInfo &sim_info, int runs, int seed)
        : sim_info(sim_info), runs(runs), rng(seed), cache(), p() {

        Y.resize(sim_info.m);
        for (auto &v: Y) v.resize(sim_info.n);
    }

    double run() {
        const int j = sim_info.m - 1;
        p.clear();

        for (int i = 0; i < runs; ++i) {
            realize();
            cache.clear();

            vector<int> pre(sim_info.m, -1);
            T(j, pre);

            update_result(pre);
        }


        return 0;
    }

    void update_result(const vector<int> &pre) {
        vector<int> pi;
        for (int x = sim_info.m-1; x != -1; x = pre[x]) {
            pi.push_back(x+1);
        }
        ++p[pi];
    }

    void print_result() {
        for (const auto& [pi, cnt]: p) {
            cout << ":";
            for (int i=(int)pi.size()-1; i>0; --i) {
                cout << "a" << pi[i] << "/" << pi[i-1];
                cout << (i > 0 ? ":" : ",");
            }
            cout << " " << setprecision(5) << scientific << (double)cnt / runs << "\n";
        }
    }

    void realize() {
        for (int i = 0; i < sim_info.m; ++i) {
            for (int j = 0; j < sim_info.m; ++j) {
                if (sim_info.Y[i][j] > 0) {
                    Y[i][j] = rng.random(sim_info.Y[i][j]);
                }
            }
        }
    }

    double T(int j, vector<int> &pre) {
        if (cache.count(j)) {
            return cache[j];
        }

        int k = 0; // 0-based
        int l = 0;
        double tmax = 0.;

        while (l < (int)sim_info.B[j].size()) {
            if (sim_info.N[j][k] == -1) {
                int i = 0; // 0-based
                while (sim_info.N[i][k] != 1) {
                    i++;
                }
                double t = T(i, pre) + Y[i][j];
                if (t >= tmax) {
                    tmax = t;
                    pre[j] = i;
                }
                l++;
            }
            k++;
        }

        return cache[j] = tmax;
    }

    const SimInfo& sim_info;
    int runs;
    RNG rng;
    map<int, double> cache;
    vector<vector<double>> Y;
    map<vector<int>, int> p;
};

pair<int, int> determine_m_and_n(const char *net_file) {
    ifstream inf(net_file);
    int a, b;
    double c;
    int m = 0, n = 0;

    while (inf >> a >> b >> c) {
        m = max(a, m);
        m = max(b, m);
        ++n;
    }

    return { m, n };
}

SimInfo read_into_sim_info(const char *net_file, int m, int n) {
    SimInfo sim_info(m, n);

    ifstream inf(net_file);
    int a, b;
    double c;

    int j = 0;
    while (inf >> a >> b >> c) {
        --a, --b;
        sim_info.N[a][j] = 1;
        sim_info.N[b][j] = -1;
        sim_info.Y[a][b] = c;
        sim_info.B[b].push_back(a);
        ++j;
    }

    return sim_info;
}

SimInfo init_net_from_file(const char *net_file) {
    auto [m, n] = determine_m_and_n(net_file);
    return read_into_sim_info(net_file, m, n);
}

void run(const char *net_file, int runs, int seed) {
    SimInfo sim_info = init_net_from_file(net_file);
    SimAlgo sim_algo(sim_info, runs, seed);

    sim_algo.run();
    sim_algo.print_result();
}

int main(int argc, char **argv) {
    if (argc != 4) {
        cout << "Usage: " << *argv << " NET RUNS SEED\n";
        return 0;
    }

    // run(/*net_file=*/argv[1], /*runs=*/atoi(argv[2]), /*seed=*/atoi(argv[2]));
    run(/*net_file=*/argv[3], /*runs=*/atoi(argv[2]), /*seed=*/atoi(argv[1]));

    return 0;
}
