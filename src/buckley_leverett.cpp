// 1D Buckley-Leverett waterflood, sequential implicit, linear algebra by AMGCL.
//
// Dimensionless incompressible two-phase displacement on x in [0, 1]:
//
//   dS/dt + d f(S)/dx = 0
//   d/dx ( k(x) lambda_t(S) dp/dx ) = 0,   -k lambda_t dp/dx = 1,   p(1) = 0
//
// t is pore volumes injected and the total velocity is fixed by the inlet
// rate, so the saturation profile is the classical Buckley-Leverett solution
// (Welge). Heterogeneous k changes the pressure, not the saturation.
//
// Each time step:
//   1. Implicit upwind finite volume for S. Newton linear systems are
//      nonsymmetric and are solved with BiCGStab, preconditioned by
//      smoothed-aggregation AMG with ILU(0) relaxation.
//   2. TPFA pressure system. It is SPD and is solved with CG, preconditioned
//      by smoothed-aggregation AMG with SPAI(0) relaxation.
//
// Build: see CMakeLists.txt. Run: buckley [cells] [pore_volumes] [cfl]

#ifndef AMGCL_NO_BOOST
#define AMGCL_NO_BOOST
#endif

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/cg.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

using Backend = amgcl::backend::builtin<double>;

// Elliptic pressure block: SPD, so CG + AMG(SPAI0).
using PressureSolver = amgcl::make_solver<
    amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::spai0>,
    amgcl::solver::cg<Backend>>;

// Transport Newton matrix: nonsymmetric lower bidiagonal, so BiCGStab + AMG(ILU0).
using TransportSolver = amgcl::make_solver<
    amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::ilu0>,
    amgcl::solver::bicgstab<Backend>>;

struct Fluid {
    double Swc = 0.2;
    double Sor = 0.2;
    double nw = 2.0;
    double no = 2.0;
    double krw_end = 1.0;
    double kro_end = 1.0;
    double mu_w = 1.0;  // cP, relative
    double mu_o = 2.0;

    double span() const { return 1.0 - Swc - Sor; }
    double Smax() const { return 1.0 - Sor; }

    // Endpoint mobility ratio M = (krw_end/mu_w) / (kro_end/mu_o).
    double mobility_ratio() const {
        return (krw_end / mu_w) / (kro_end / mu_o);
    }

    void mobilities(double S, double& lw, double& lo, double& dlw, double& dlo) const {
        const double width = span();
        const double s = std::clamp((S - Swc) / width, 0.0, 1.0);
        const double ds_dS = 1.0 / width;
        const double Aw = krw_end / mu_w;
        const double Ao = kro_end / mu_o;

        if (s <= 0.0) {
            lw = 0.0;
            lo = Ao;
            dlw = (nw == 1.0) ? Aw * ds_dS : 0.0;
            dlo = 0.0;
            return;
        }
        if (s >= 1.0) {
            lw = Aw;
            lo = 0.0;
            dlw = 0.0;
            dlo = (no == 1.0) ? -Ao * ds_dS : 0.0;
            return;
        }

        lw = Aw * std::pow(s, nw);
        lo = Ao * std::pow(1.0 - s, no);
        dlw = Aw * nw * std::pow(s, nw - 1.0) * ds_dS;
        dlo = -Ao * no * std::pow(1.0 - s, no - 1.0) * ds_dS;
    }

    void flow(double S, double& f, double& df) const {
        double lw, lo, dlw, dlo;
        mobilities(S, lw, lo, dlw, dlo);
        const double lt = lw + lo;
        f = lw / lt;
        df = (dlw * lo - lw * dlo) / (lt * lt);
    }

    double fractional_flow(double S) const {
        double f, df;
        flow(S, f, df);
        return f;
    }

    double total_mobility(double S) const {
        double lw, lo, dlw, dlo;
        mobilities(S, lw, lo, dlw, dlo);
        return lw + lo;
    }
};

// Welge shock: the chord from (Swc, 0) that is tangent to f(S).
struct Welge {
    double shock_saturation = 0.0;
    double shock_speed = 0.0;
    double tangent_residual = 0.0;

    explicit Welge(const Fluid& fluid) {
        const double Swc = fluid.Swc;
        const double Smax = fluid.Smax();
        auto sigma = [&](double S) { return fluid.fractional_flow(S) / (S - Swc); };

        // Golden-section search for the maximum of the chord slope.
        const double phi = 0.5 * (std::sqrt(5.0) - 1.0);
        double a = Swc + 1e-8;
        double b = Smax - 1e-8;
        double c = b - phi * (b - a);
        double d = a + phi * (b - a);
        for (int i = 0; i < 80; ++i) {
            if (sigma(c) > sigma(d)) {
                b = d;
                d = c;
                c = b - phi * (b - a);
            } else {
                a = c;
                c = d;
                d = a + phi * (b - a);
            }
        }
        double S = 0.5 * (a + b);

        // Polish the root of g(S) = f'(S) (S - Swc) - f(S) = 0.
        for (int it = 0; it < 25; ++it) {
            double fv, dfv;
            fluid.flow(S, fv, dfv);
            const double g = dfv * (S - Swc) - fv;
            const double eps = 1e-7;
            double fp, dfp, fm, dfm;
            fluid.flow(S + eps, fp, dfp);
            fluid.flow(S - eps, fm, dfm);
            const double gp = dfp * ((S + eps) - Swc) - fp;
            const double gm = dfm * ((S - eps) - Swc) - fm;
            const double dg = (gp - gm) / (2.0 * eps);
            if (!(std::abs(dg) > 0.0) || !std::isfinite(dg)) break;
            const double step = g / dg;
            const double Sn = std::clamp(S - step, Swc + 1e-6, Smax - 1e-6);
            if (std::abs(Sn - S) < 1e-14) {
                S = Sn;
                break;
            }
            S = Sn;
        }

        double fv, dfv;
        fluid.flow(S, fv, dfv);
        shock_saturation = S;
        shock_speed = fv / (S - Swc);
        tangent_residual = dfv - shock_speed;
    }

    double saturation(const Fluid& fluid, double x_over_t) const {
        if (!(x_over_t < shock_speed)) return fluid.Swc;
        const double Smax = fluid.Smax();
        if (Smax - shock_saturation < 1e-8) return Smax;

        // f' decreases from the shock speed at S_f to about 0 at Smax.
        double lo = shock_saturation;
        double hi = Smax;
        for (int i = 0; i < 70; ++i) {
            const double mid = 0.5 * (lo + hi);
            double f, df;
            fluid.flow(mid, f, df);
            if (df > x_over_t) lo = mid;
            else hi = mid;
        }
        return 0.5 * (lo + hi);
    }
};

double permeability(double x) {
    // Dimensionless k/k_ref. The tight streak steepens dp/dx and does not
    // change the saturation when the inlet rate is prescribed.
    if (x < 0.45) return 1.0;
    if (x < 0.55) return 0.1;
    return 1.0;
}

struct Crs {
    int n = 0;
    std::vector<int> ptr;
    std::vector<int> col;
    std::vector<double> val;
};

struct SolveStats {
    size_t calls = 0;
    size_t iterations = 0;
    double max_residual = 0.0;
    std::string hierarchy;
    bool show_next = false;
};

template <class Solver>
void solve_linear(typename Solver::params const& prm, Crs& A, std::vector<double> const& rhs,
                  std::vector<double>& x, SolveStats& stats) {
    int n = A.n;
    Solver solver(std::tie(n, A.ptr, A.col, A.val), prm);
    if (stats.show_next) {
        std::ostringstream os;
        os << solver << "\n";
        stats.hierarchy = os.str();
        stats.show_next = false;
    }

    size_t iters = 0;
    double error = 0.0;
    std::tie(iters, error) = solver(rhs, x);
    if (!std::isfinite(error) || error > 1e-8) {
        std::ostringstream os;
        os << "AMGCL did not converge: iterations=" << iters << " residual=" << error;
        throw std::runtime_error(os.str());
    }
    ++stats.calls;
    stats.iterations += iters;
    stats.max_residual = std::max(stats.max_residual, error);
}

double derivative_check(Fluid const& fluid) {
    double worst = 0.0;
    for (int i = 1; i < 40; ++i) {
        const double s = static_cast<double>(i) / 40.0;
        const double S = fluid.Swc + s * fluid.span();
        double f, df;
        fluid.flow(S, f, df);
        const double eps = 1e-6;
        const double fd =
            (fluid.fractional_flow(S + eps) - fluid.fractional_flow(S - eps)) / (2.0 * eps);
        worst = std::max(worst, std::abs(df - fd));
    }
    return worst;
}

struct Grid {
    int n = 0;
    double dx = 0.0;
    std::vector<double> x;
    std::vector<double> k;
    std::vector<double> S;
    std::vector<double> p;
};

Grid make_grid(int n, Fluid const& fluid) {
    Grid g;
    g.n = n;
    g.dx = 1.0 / static_cast<double>(n);
    g.x.resize(n);
    g.k.resize(n);
    g.S.assign(n, fluid.Swc);
    g.p.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        g.x[i] = (static_cast<double>(i) + 0.5) * g.dx;
        g.k[i] = permeability(g.x[i]);
    }
    return g;
}

double half_transmissibility(Grid const& g, Fluid const& fluid, int i) {
    // Center-to-face distance is dx/2.
    return 2.0 * g.k[i] * fluid.total_mobility(g.S[i]) / g.dx;
}

double face_transmissibility(Grid const& g, Fluid const& fluid, int left) {
    const double a = half_transmissibility(g, fluid, left);
    const double b = half_transmissibility(g, fluid, left + 1);
    return (a * b) / (a + b);
}

Crs assemble_pressure(Grid const& g, Fluid const& fluid, std::vector<double>& rhs) {
    const int n = g.n;
    Crs A;
    A.n = n;
    A.ptr.reserve(n + 1);
    A.col.reserve(3 * n);
    A.val.reserve(3 * n);
    rhs.assign(n, 0.0);
    A.ptr.push_back(0);

    for (int i = 0; i < n; ++i) {
        if (i == 0) {
            const double T = face_transmissibility(g, fluid, 0);
            A.col.push_back(0);
            A.val.push_back(T);
            A.col.push_back(1);
            A.val.push_back(-T);
            rhs[0] = 1.0;  // prescribed total rate
        } else if (i == n - 1) {
            const double T = face_transmissibility(g, fluid, n - 2);
            const double Tout = half_transmissibility(g, fluid, n - 1);  // p(face=1) = 0
            A.col.push_back(n - 2);
            A.val.push_back(-T);
            A.col.push_back(n - 1);
            A.val.push_back(T + Tout);
        } else {
            const double Tl = face_transmissibility(g, fluid, i - 1);
            const double Tr = face_transmissibility(g, fluid, i);
            A.col.push_back(i - 1);
            A.val.push_back(-Tl);
            A.col.push_back(i);
            A.val.push_back(Tl + Tr);
            A.col.push_back(i + 1);
            A.val.push_back(-Tr);
        }
        A.ptr.push_back(static_cast<int>(A.col.size()));
    }
    return A;
}

std::vector<double> exact_pressure(Grid const& g, Fluid const& fluid) {
    // Exact solution of the discrete TPFA system: the face flux equals 1,
    // and p at the outlet face is 0. Piecewise-constant k*lambda makes this
    // the same as integrating dp = -dx / (k lambda) from cell to cell.
    const int n = g.n;
    std::vector<double> p(n);
    p[n - 1] = 1.0 / half_transmissibility(g, fluid, n - 1);
    for (int i = n - 2; i >= 0; --i) {
        p[i] = p[i + 1] + 1.0 / face_transmissibility(g, fluid, i);
    }
    return p;
}

double max_flux_error(Grid const& g, Fluid const& fluid) {
    double worst = 0.0;
    const int n = g.n;
    for (int i = 0; i < n - 1; ++i) {
        const double flux = face_transmissibility(g, fluid, i) * (g.p[i] - g.p[i + 1]);
        worst = std::max(worst, std::abs(flux - 1.0));
    }
    const double out = half_transmissibility(g, fluid, n - 1) * g.p[n - 1];
    return std::max(worst, std::abs(out - 1.0));
}

double residual_linf(Grid const& g, Fluid const& fluid, std::vector<double> const& Sold, double dt,
                     std::vector<double> const& S) {
    const double f_inj = 1.0;  // kro(Smax) = 0
    double worst = 0.0;
    for (int i = 0; i < g.n; ++i) {
        const double up = (i == 0) ? f_inj : fluid.fractional_flow(S[i - 1]);
        const double R = (g.dx / dt) * (S[i] - Sold[i]) + fluid.fractional_flow(S[i]) - up;
        worst = std::max(worst, std::abs(R));
    }
    return worst;
}

Crs assemble_transport(Grid const& g, Fluid const& fluid, std::vector<double> const& Sold, double dt,
                       std::vector<double> const& S, std::vector<double>& rhs) {
    const int n = g.n;
    const double f_inj = 1.0;
    const double accum = g.dx / dt;
    Crs A;
    A.n = n;
    A.ptr.resize(n + 1);
    A.col.reserve(2 * n);
    A.val.reserve(2 * n);
    rhs.resize(n);

    std::vector<double> f(n), df(n);
    for (int i = 0; i < n; ++i) fluid.flow(S[i], f[i], df[i]);

    int nnz = 0;
    for (int i = 0; i < n; ++i) {
        A.ptr[i] = nnz;
        if (i > 0) {
            double df_up;
            double dummy;
            fluid.flow(S[i - 1], dummy, df_up);
            A.col.push_back(i - 1);
            A.val.push_back(-df_up);
            ++nnz;
        }
        A.col.push_back(i);
        A.val.push_back(accum + df[i]);
        ++nnz;

        const double up = (i == 0) ? f_inj : f[i - 1];
        const double R = accum * (S[i] - Sold[i]) + f[i] - up;
        rhs[i] = -R;
    }
    A.ptr[n] = nnz;
    return A;
}

void newton_saturation(Grid& g, Fluid const& fluid, std::vector<double> const& Sold, double dt,
                       TransportSolver::params const& prm, SolveStats& stats) {
    const double tol = 1e-9;
    const int max_it = 15;
    for (int it = 0; it < max_it; ++it) {
        const double rn = residual_linf(g, fluid, Sold, dt, g.S);
        if (rn < tol) return;

        std::vector<double> rhs;
        Crs A = assemble_transport(g, fluid, Sold, dt, g.S, rhs);
        std::vector<double> dS(g.n, 0.0);
        solve_linear<TransportSolver>(prm, A, rhs, dS, stats);

        double omega = 1.0;
        bool accepted = false;
        std::vector<double> trial(g.n);
        for (int ls = 0; ls < 12; ++ls) {
            for (int i = 0; i < g.n; ++i) {
                trial[i] = std::clamp(g.S[i] + omega * dS[i], fluid.Swc, fluid.Smax());
            }
            const double rt = residual_linf(g, fluid, Sold, dt, trial);
            if (rt <= rn * (1.0 + 1e-8)) {
                g.S.swap(trial);
                accepted = true;
                break;
            }
            omega *= 0.5;
        }
        if (!accepted) {
            throw std::runtime_error("Newton line search rejected the saturation update");
        }
    }
    throw std::runtime_error("Newton did not converge");
}

void solve_pressure(Grid& g, Fluid const& fluid, PressureSolver::params const& prm, SolveStats& stats) {
    std::vector<double> rhs;
    Crs A = assemble_pressure(g, fluid, rhs);
    solve_linear<PressureSolver>(prm, A, rhs, g.p, stats);
}

double mean_saturation_gain(Grid const& g, Fluid const& fluid) {
    double sum = 0.0;
    for (double S : g.S) sum += S - fluid.Swc;
    return sum * g.dx;
}

struct Errors {
    double saturation_l1 = 0.0;
    double saturation_linf = 0.0;
    double pressure = 0.0;
    double flux = 0.0;
    double conservation = 0.0;
    double shock_x = 0.0;
};

Errors measure(Grid const& g, Fluid const& fluid, Welge const& welge, double time, double produced) {
    Errors e;
    e.flux = max_flux_error(g, fluid);
    const auto p_exact = exact_pressure(g, fluid);
    for (int i = 0; i < g.n; ++i) {
        const double Se = welge.saturation(fluid, g.x[i] / time);
        e.saturation_l1 += std::abs(g.S[i] - Se);
        e.saturation_linf = std::max(e.saturation_linf, std::abs(g.S[i] - Se));
        e.pressure = std::max(e.pressure, std::abs(g.p[i] - p_exact[i]));
    }
    e.saturation_l1 *= g.dx;
    // Backward Euler finite volume is conservative: gained water = injected - produced.
    e.conservation = std::abs(mean_saturation_gain(g, fluid) - (time - produced));

    const double level = 0.5 * (welge.shock_saturation + fluid.Swc);
    e.shock_x = 1.0;
    for (int i = 0; i < g.n; ++i) {
        if (g.S[i] < level) {
            e.shock_x = g.x[i];
            break;
        }
    }
    return e;
}

void write_outputs(Grid const& g, Fluid const& fluid, Welge const& welge, double time,
                   std::filesystem::path const& dir) {
    std::filesystem::create_directories(dir);
    const auto p_exact = exact_pressure(g, fluid);

    {
        std::ofstream out(dir / "profile.csv");
        out << std::setprecision(16);
        out << "x,S,S_exact,p,p_exact,k,lambda_t,f\n";
        for (int i = 0; i < g.n; ++i) {
            const double Se = welge.saturation(fluid, g.x[i] / time);
            out << g.x[i] << ',' << g.S[i] << ',' << Se << ',' << g.p[i] << ',' << p_exact[i] << ','
                << g.k[i] << ',' << fluid.total_mobility(g.S[i]) << ',' << fluid.fractional_flow(g.S[i])
                << '\n';
        }
    }

    {
        std::ofstream out(dir / "fractional_flow.csv");
        out << std::setprecision(16);
        out << "S,krw,kro,f,convex_hull\n";
        const int m = 400;
        for (int i = 0; i <= m; ++i) {
            const double S = fluid.Swc + fluid.span() * (static_cast<double>(i) / m);
            double lw, lo, dlw, dlo;
            fluid.mobilities(S, lw, lo, dlw, dlo);
            const double f = fluid.fractional_flow(S);
            const double hull = (S <= welge.shock_saturation) ? welge.shock_speed * (S - fluid.Swc) : f;
            const double sn = (S - fluid.Swc) / fluid.span();
            out << S << ',' << fluid.krw_end * std::pow(sn, fluid.nw) << ','
                << fluid.kro_end * std::pow(1.0 - sn, fluid.no) << ',' << f << ',' << hull << '\n';
        }
    }

    // Two-panel SVG so the shock and the pressure drop can be opened directly.
    {
        const int width = 880;
        const int height = 720;
        const double left = 70, right = 840, top1 = 40, bot1 = 330, top2 = 390, bot2 = 680;
        auto sx = [&](double x) { return left + x * (right - left); };

        double pmax = 0.0;
        for (double v : g.p) pmax = std::max(pmax, v);
        for (double v : p_exact) pmax = std::max(pmax, v);
        pmax *= 1.05;
        auto syS = [&](double S) {
            const double z = (S - fluid.Swc) / fluid.span();
            return bot1 - z * (bot1 - top1);
        };
        auto syP = [&](double p) { return bot2 - (p / pmax) * (bot2 - top2); };

        std::ostringstream svg;
        svg << std::setprecision(8);
        svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\"" << height
            << "\" font-family=\"Consolas, monospace\" font-size=\"14\">\n";
        svg << "<rect width=\"100%\" height=\"100%\" fill=\"#f7f4ee\"/>\n";
        svg << "<text x=\"70\" y=\"24\">Saturation at t = " << time
            << " PVI. Blue: implicit upwind. Red: Welge.</text>\n";
        svg << "<rect x=\"" << left << "\" y=\"" << top1 << "\" width=\"" << (right - left) << "\" height=\""
            << (bot1 - top1) << "\" fill=\"none\" stroke=\"#222\"/>\n";
        svg << "<path d=\"";
        for (int i = 0; i < g.n; ++i) {
            svg << (i ? "L" : "M") << sx(g.x[i]) << ' ' << syS(g.S[i]) << ' ';
        }
        svg << "\" fill=\"none\" stroke=\"#1d4e89\" stroke-width=\"2\"/>\n<path d=\"";
        for (int i = 0; i < g.n; ++i) {
            const double Se = welge.saturation(fluid, g.x[i] / time);
            svg << (i ? "L" : "M") << sx(g.x[i]) << ' ' << syS(Se) << ' ';
        }
        svg << "\" fill=\"none\" stroke=\"#b33c2e\" stroke-width=\"1.5\" stroke-dasharray=\"5 4\"/>\n";
        svg << "<text x=\"70\" y=\"360\">Pressure. Blue: AMGCL (CG + AMG). The drop steepens in the tight "
               "streak x in [0.45, 0.55].</text>\n";
        svg << "<rect x=\"" << left << "\" y=\"" << top2 << "\" width=\"" << (right - left) << "\" height=\""
            << (bot2 - top2) << "\" fill=\"none\" stroke=\"#222\"/>\n<path d=\"";
        for (int i = 0; i < g.n; ++i) {
            svg << (i ? "L" : "M") << sx(g.x[i]) << ' ' << syP(g.p[i]) << ' ';
        }
        svg << "\" fill=\"none\" stroke=\"#1d4e89\" stroke-width=\"2\"/>\n";
        svg << "<text x=\"70\" y=\"708\">x (dimensionless)</text>\n</svg>\n";
        std::ofstream out(dir / "profile.svg");
        out << svg.str();
    }
}

std::string ascii_saturation(Grid const& g, Fluid const& fluid) {
    const int cols = 72;
    const int rows = 12;
    std::vector<std::string> canvas(rows, std::string(cols, ' '));
    for (int c = 0; c < cols; ++c) {
        const int i = std::min(g.n - 1, static_cast<int>(std::floor((c + 0.5) / cols * g.n)));
        const double z = std::clamp((g.S[i] - fluid.Swc) / fluid.span(), 0.0, 1.0);
        const int r = rows - 1 - std::min(rows - 1, static_cast<int>(std::lround(z * (rows - 1))));
        canvas[r][c] = '#';
    }
    std::ostringstream os;
    os << "Smax |" << canvas[0] << "\n";
    for (int r = 1; r < rows - 1; ++r) os << "     |" << canvas[r] << "\n";
    os << "Swc  |" << canvas[rows - 1] << "\n";
    os << "     +" << std::string(cols, '-') << "\n";
    os << "      0" << std::string(cols - 8, ' ') << "x=1\n";
    return os.str();
}

template <class Solver>
typename Solver::params solver_params() {
    typename Solver::params prm;
    prm.solver.tol = 1e-10;
    prm.solver.maxiter = 200;
    // Stop coarsening while the coarse problem is still a 1D chain, so the
    // printed hierarchy is actually multilevel on this grid.
    prm.precond.coarse_enough = 60;
    prm.precond.npre = 1;
    prm.precond.npost = 1;
    return prm;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        int n = 800;
        double cfl = 1.0;
        bool t_end_set = false;
        double t_end = 0.0;
        if (argc > 1) n = std::stoi(argv[1]);
        if (argc > 2) {
            t_end = std::stod(argv[2]);
            t_end_set = true;
        }
        if (argc > 3) cfl = std::stod(argv[3]);
        if (n < 16 || !(cfl > 0.0) || (t_end_set && !(t_end > 0.0))) {
            throw std::runtime_error("Usage: buckley [cells>=16] [pore_volumes] [cfl>0]");
        }

        Fluid fluid;
        if (!(fluid.span() > 0.0)) throw std::runtime_error("Swc + Sor must be < 1");
        Welge welge(fluid);
        if (!t_end_set) t_end = 0.55 / welge.shock_speed;

        const double df_error = derivative_check(fluid);
        Grid grid = make_grid(n, fluid);
        const double dt_target = cfl * grid.dx / welge.shock_speed;

        auto pressure_prm = solver_params<PressureSolver>();
        auto transport_prm = solver_params<TransportSolver>();
        SolveStats pressure_stats;
        SolveStats transport_stats;
        pressure_stats.show_next = true;

        std::ostringstream report;
        report << std::setprecision(6);
        report << "1D Buckley-Leverett, sequential implicit, AMGCL\n";
        report << "Corey nw=" << fluid.nw << " no=" << fluid.no << " Swc=" << fluid.Swc
               << " Sor=" << fluid.Sor << " mu_o/mu_w=" << (fluid.mu_o / fluid.mu_w)
               << " endpoint M=" << fluid.mobility_ratio() << "\n";
        report << "Shock saturation Sf=" << welge.shock_saturation
               << " speed=" << welge.shock_speed
               << " breakthrough PVI=" << (1.0 / welge.shock_speed)
               << " tangent residual=" << welge.tangent_residual << "\n";
        report << "Grid n=" << n << " dx=" << grid.dx << " t_end=" << t_end
               << " PVI  dt=" << dt_target << " (CFL=" << cfl << " on the shock)\n";
        report << "df/dS finite-difference check: " << df_error << "\n\n";

        solve_pressure(grid, fluid, pressure_prm, pressure_stats);
        const double p_inlet_initial = grid.p.front();

        double time = 0.0;
        double produced = 0.0;
        int steps = 0;
        int cuts = 0;
        bool transport_hierarchy_done = false;
        const double hierarchy_time = 0.45 * t_end;

        while (time < t_end - 1e-14) {
            double dt = std::min(dt_target, t_end - time);
            const std::vector<double> Sold = grid.S;
            const std::vector<double> pold = grid.p;
            bool advanced = false;
            for (int attempt = 0; attempt < 12; ++attempt) {
                try {
                    if (!transport_hierarchy_done && time >= hierarchy_time) {
                        transport_stats.show_next = true;
                    }
                    newton_saturation(grid, fluid, Sold, dt, transport_prm, transport_stats);
                    solve_pressure(grid, fluid, pressure_prm, pressure_stats);
                    const Errors step_err = measure(grid, fluid, welge, time + dt, produced + fluid.fractional_flow(grid.S.back()) * dt);
                    if (step_err.pressure > 1e-6 || step_err.flux > 1e-6) {
                        throw std::runtime_error("pressure solution does not match the TPFA flux");
                    }
                    advanced = true;
                    if (transport_stats.show_next == false && !transport_hierarchy_done && time >= hierarchy_time) {
                        transport_hierarchy_done = true;
                    }
                    break;
                } catch (std::exception const& err) {
                    grid.S = Sold;
                    grid.p = pold;
                    dt *= 0.5;
                    ++cuts;
                    std::cerr << "cutting dt to " << dt << " after: " << err.what() << "\n";
                    if (attempt == 11 || !(dt > 1e-10)) {
                        throw std::runtime_error(std::string("time step failed: ") + err.what());
                    }
                }
            }
            if (!advanced) throw std::runtime_error("time step did not advance");
            produced += fluid.fractional_flow(grid.S.back()) * dt;
            time += dt;
            ++steps;
        }

        const Errors err = measure(grid, fluid, welge, time, produced);
        const double p_inlet_final = grid.p.front();
        const double pressure_iters =
            pressure_stats.calls ? static_cast<double>(pressure_stats.iterations) / pressure_stats.calls : 0.0;
        const double transport_iters =
            transport_stats.calls ? static_cast<double>(transport_stats.iterations) / transport_stats.calls : 0.0;

        report << "Pressure AMG (CG + smoothed aggregation + SPAI0), first setup:\n";
        report << pressure_stats.hierarchy << "\n";
        report << "Transport AMG (BiCGStab + smoothed aggregation + ILU0), mid-flood Jacobian:\n";
        report << (transport_stats.hierarchy.empty() ? "(hierarchy was not captured)\n" : transport_stats.hierarchy)
               << "\n";
        report << "Steps=" << steps << " dt cuts=" << cuts << "\n";
        report << "Pressure solves=" << pressure_stats.calls
               << " average iterations=" << pressure_iters
               << " max residual=" << pressure_stats.max_residual << "\n";
        report << "Transport solves=" << transport_stats.calls
               << " average iterations=" << transport_iters
               << " max residual=" << transport_stats.max_residual << "\n";
        report << "Inlet pressure t=0: " << p_inlet_initial << "   t=end: " << p_inlet_final << "\n";
        report << "Inlet pressure rises because lambda_t has a minimum between Swc and Smax\n";
        report << "(two-phase interference), so the same rate needs a larger drawdown.\n";
        report << "Transport iterations stay at 1: the Newton matrix is bidiagonal and ILU(0)\n";
        report << "factors that pattern exactly. The tridiagonal pressure system is the one\n";
        report << "that actually exercises the AMG V-cycle (CG iterations above).\n";
        report << "Saturation L1=" << err.saturation_l1 << " Linf=" << err.saturation_linf << "\n";
        report << "Shock x numeric=" << err.shock_x << " exact=" << (welge.shock_speed * time) << "\n";
        report << "Pressure max-abs=" << err.pressure << " face-flux max-abs=" << err.flux << "\n";
        report << "Conservation |gained - (injected-produced)|=" << err.conservation << "\n\n";
        report << ascii_saturation(grid, fluid);

        const bool ok = df_error < 1e-5 && std::abs(welge.tangent_residual) < 1e-6 &&
                        err.pressure < 1e-6 && err.flux < 1e-6 && err.conservation < 1e-8 &&
                        err.saturation_l1 < 0.03 && pressure_stats.calls > 0 && transport_stats.calls > 0;

        report << (ok ? "CHECKS PASSED\n" : "CHECKS FAILED\n");
        std::cout << report.str();

        const std::filesystem::path dir = "output";
        write_outputs(grid, fluid, welge, time, dir);
        std::ofstream(dir / "report.txt") << report.str();
        std::cout << "Wrote " << (dir / "profile.csv").string() << ", "
                  << (dir / "fractional_flow.csv").string() << ", "
                  << (dir / "profile.svg").string() << " and "
                  << (dir / "report.txt").string() << "\n";
        return ok ? 0 : 1;
    } catch (std::exception const& err) {
        std::cerr << "error: " << err.what() << "\n";
        return 1;
    }
}
