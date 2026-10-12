// apps/checks/test_math/main.cpp -- known-answer self-checks for the math core (math.hpp / src/math/math.cpp).
// Each case prints "[TEST] <name> expected=<> actual=<> <PASS|FAIL>"; exit code 1 on any FAIL, 2 on a thrown exception.
// Units: lengths mm, angles rad unless a case name or local comment says otherwise.
#include "sgkquadric/math.hpp"

#include <Foundation/init.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace sggk;
using namespace sgkq;

namespace
{

// Angles in degrees; conversion and trig use double-width pi.
const double kPi = 3.14159265358979323846;

int g_total = 0;
int g_failed = 0;

// Prints one [TEST] line and accumulates the totals.
void report(const std::string& name, double expected, double actual, bool pass)
{
    ++g_total;
    if (!pass) ++g_failed;
    std::cout << "[TEST] " << name
              << " expected=" << expected
              << " actual=" << actual
              << (pass ? " PASS" : " FAIL") << std::endl;
}

// Relative compare: |actual - expected| <= tol * max(1, |expected|); the floor of 1 keeps a near-zero expected from over-sensitizing.
void checkClose(const std::string& name, double expected, double actual, double tol)
{
    const double scale = std::max(1.0, std::fabs(expected));
    report(name, expected, actual, std::fabs(actual - expected) <= tol * scale);
}

// Sign-insensitive axis angle in degrees (|cos|); a degenerate zero vector returns 180 to force a FAIL.
double axisAngleDeg(const Vector3D& a, const Vector3D& b)
{
    const double la = a.Length();
    const double lb = b.Length();
    if (la <= 0.0 || lb <= 0.0) return 180.0;
    double c = std::fabs(a.Dot(b) / (la * lb));
    if (c > 1.0) c = 1.0;
    return std::acos(c) * 180.0 / kPi;
}

// 1) symEig3: A = Q diag(2,1,3) Q^T with a known rotation Q; d is deliberately unsorted to test sorting.
void testSymEig3()
{
    Vector3D axis(1.0, 2.0, 2.0);
    axis.Normalize();
    const double ang = 0.7;
    const double c = std::cos(ang), s = std::sin(ang), t = 1.0 - c;
    const double x = axis.X(), y = axis.Y(), z = axis.Z();

    const double Q[3][3] = {
        { t * x * x + c,     t * x * y - s * z, t * x * z + s * y },
        { t * x * y + s * z, t * y * y + c,     t * y * z - s * x },
        { t * x * z - s * y, t * y * z + s * x, t * z * z + c     }
    };
    const double d[3] = { 2.0, 1.0, 3.0 };
    double A[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            double v = 0.0;
            for (int k = 0; k < 3; ++k) v += Q[i][k] * d[k] * Q[j][k];
            A[i][j] = v;
        }

    const SymEig3 e = symEig3(A);
    checkClose("symEig3.val0", 1.0, e.val[0], 1e-12);
    checkClose("symEig3.val1", 2.0, e.val[1], 1e-12);
    checkClose("symEig3.val2", 3.0, e.val[2], 1e-12);

    // eigvec for lambda 1/2/3 is column 2/1/3 of Q; values checked to 1e-12, directions to 1e-9 deg.
    const Vector3D q0(Q[0][1], Q[1][1], Q[2][1]);
    const Vector3D q1(Q[0][0], Q[1][0], Q[2][0]);
    const Vector3D q2(Q[0][2], Q[1][2], Q[2][2]);
    checkClose("symEig3.vec0", 0.0, axisAngleDeg(e.vec[0], q0), 1e-9);
    checkClose("symEig3.vec1", 0.0, axisAngleDeg(e.vec[1], q1), 1e-9);
    checkClose("symEig3.vec2", 0.0, axisAngleDeg(e.vec[2], q2), 1e-9);
}

// 2) Cone axis, half-angle 60 deg: the centered-covariance min eigenvector is the true axis, while the
// uncentered second moment's LARGEST eigenvector is perpendicular to it (the classic cone-axis trap).
void testConeAxis()
{
    const double alpha = 60.0 * kPi / 180.0;
    const double sa = std::sin(alpha), ca = std::cos(alpha);

    std::vector<Vector3D> gen;
    for (int i = 0; i < 36; ++i)
    {
        const double th = 2.0 * kPi * i / 36.0;
        gen.push_back(Vector3D(sa * std::cos(th), sa * std::sin(th), ca));
    }

    const Vector3D a = minVarDirection(gen);
    const double ang = axisAngleDeg(a, Vector3D(0.0, 0.0, 1.0));
    report("cone_axis_covariance", 0.0, ang, ang < 1e-10);

    double M[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
    for (size_t i = 0; i < gen.size(); ++i)
    {
        const double u[3] = { gen[i].X(), gen[i].Y(), gen[i].Z() };
        for (int r = 0; r < 3; ++r)
            for (int cc = 0; cc < 3; ++cc)
                M[r][cc] += u[r] * u[cc];
    }
    const SymEig3 em = symEig3(M);
    const double angMax = axisAngleDeg(em.vec[2], Vector3D(0.0, 0.0, 1.0));
    report("cone_axis_second_moment_trap", 90.0, angMax, std::fabs(angMax - 90.0) < 1e-6);
}

// 3) Cylinder-normal axis recovery: full circle and a 30 deg trimmed arc; the axis must stay z even for the partial arc.
void testCylinderAxis()
{
    std::vector<Vector3D> full;
    for (int i = 0; i < 36; ++i)
    {
        const double th = 2.0 * kPi * i / 36.0;
        full.push_back(Vector3D(std::cos(th), std::sin(th), 0.0));
    }
    const double ang1 = axisAngleDeg(minVarDirection(full), Vector3D(0.0, 0.0, 1.0));
    report("cylinder_axis_full", 0.0, ang1, ang1 < 1e-10);

    // 36 points over a 30 deg arc: partial coverage must still pin the axis.
    std::vector<Vector3D> arc;
    for (int i = 0; i < 36; ++i)
    {
        const double th = (30.0 * kPi / 180.0) * i / 35.0;
        arc.push_back(Vector3D(std::cos(th), std::sin(th), 0.0));
    }
    const double ang2 = axisAngleDeg(minVarDirection(arc), Vector3D(0.0, 0.0, 1.0));
    report("cylinder_axis_arc30", 0.0, ang2, ang2 < 1e-10);
}

// 4a) solveLinear3: well-conditioned solve, plus a singular matrix (row2 = 2*row1) that must be rejected.
void testSolveLinear3()
{
    // SPD A: any stable solver copes with it; the singular S below is the real probe.
    const double A[3][3] = { { 2.0, 1.0, 1.0 }, { 1.0, 3.0, 2.0 }, { 1.0, 2.0, 4.0 } };
    const double xt[3] = { 1.0, 2.0, 3.0 };
    double b[3];
    for (int i = 0; i < 3; ++i)
    {
        double v = 0.0;
        for (int j = 0; j < 3; ++j) v += A[i][j] * xt[j];
        b[i] = v;
    }
    double x[3];
    const bool ok = solveLinear3(A, b, x);
    double err = 0.0;
    if (ok)
        for (int i = 0; i < 3; ++i) err = std::max(err, std::fabs(x[i] - xt[i]));
    report("solveLinear3", 0.0, ok ? err : 1e30, ok && err < 1e-12);

    const double S[3][3] = { { 1.0, 2.0, 3.0 }, { 2.0, 4.0, 6.0 }, { 1.0, 1.0, 1.0 } };
    const double sb[3] = { 1.0, 2.0, 3.0 };
    double sx[3];
    const bool rejected = !solveLinear3(S, sb, sx);
    report("solveLinear3_singular", 1.0, rejected ? 1.0 : 0.0, rejected);
}

// 4b) solveLineConcurrency: 20 random lines all passing through C.
void testLineConcurrency()
{
    const Point3D C(10.0, 20.0, 30.0);
    // Fixed seed makes the random line set reproducible.
    std::mt19937 rng(20261009u);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);

    std::vector<Point3D> pts;
    std::vector<Vector3D> dirs;
    while (dirs.size() < 20)
    {
        Vector3D d(uni(rng), uni(rng), uni(rng));
        if (d.Length() < 0.1) continue;
        d.Normalize();
        const double t = 1.0 + 9.0 * std::fabs(uni(rng));
        dirs.push_back(d);
        pts.push_back(C + d * t);
    }

    Point3D X;
    const bool ok = solveLineConcurrency(pts, dirs, X);
    const double err = ok ? Point3D::Distance(X, C) : 1e30;
    report("solveLineConcurrency", 0.0, err, ok && err < 1e-10);
}

// 5) solvePluckerAxis: torus samples + true normals, closed form
// p = C + (R + r cos v)(cos u, sin u, 0) + r sin v z, n = cos v (cos u, sin u, 0) + sin v z.
void testPluckerAxis()
{
    const double Cx = 0.0, Cy = 0.0, Cz = -5.0;
    const double R = 100.0, r = 10.0;

    std::vector<Point3D> pts;
    std::vector<Vector3D> nrm;
    for (int iu = 0; iu < 24; ++iu)
    {
        const double u = 2.0 * kPi * iu / 24.0;
        const double cu = std::cos(u), su = std::sin(u);
        for (int iv = 0; iv < 12; ++iv)
        {
            const double v = 2.0 * kPi * iv / 12.0;
            const double cv = std::cos(v), sv = std::sin(v);
            pts.push_back(Point3D(Cx + (R + r * cv) * cu,
                                  Cy + (R + r * cv) * su,
                                  Cz + r * sv));
            nrm.push_back(Vector3D(cv * cu, cv * su, sv));
        }
    }

    Vector3D axisDir;
    Point3D axisPt;   // returned but not asserted; only the direction is checked.
    const bool ok = solvePluckerAxis(pts, nrm, axisDir, axisPt);
    const double ang = ok ? axisAngleDeg(axisDir, Vector3D(0.0, 0.0, 1.0)) : 180.0;
    report("solvePluckerAxis.dir", 0.0, ang, ok && ang < 1e-8);
}

// 6) fitCircle2D: circle center (3,-4) r=5, 20 samples with 1e-9 noise; the result must beat 1e-6.
void testFitCircle2D()
{
    // 1e-9 noise is far below the 1e-6 acceptance, so this probes fit stability rather than noise robustness.
    const double tcx = 3.0, tcy = -4.0, tr = 5.0;
    std::mt19937 rng(1234u);
    std::uniform_real_distribution<double> noise(-1e-9, 1e-9);

    std::vector<std::array<double, 2>> pts;
    for (int i = 0; i < 20; ++i)
    {
        const double th = 2.0 * kPi * i / 20.0;
        pts.push_back(std::array<double, 2>{
            tcx + tr * std::cos(th) + noise(rng),
            tcy + tr * std::sin(th) + noise(rng) });
    }

    double cx = 0.0, cy = 0.0, r = 0.0;
    const bool ok = fitCircle2D(pts, cx, cy, r);
    double err = 1e30;
    if (ok)
        err = std::max(std::max(std::fabs(cx - tcx), std::fabs(cy - tcy)), std::fabs(r - tr));
    report("fitCircle2D", 0.0, err, ok && err < 1e-6);
}

// 7) Five closed-form orthogonal distances: on-surface points < 1e-12; a known 0.5 mm offset (relative tol < 1e-9).
void testDistances()
{
    // On-surface cases use an absolute tol of 1e-12; offset cases use a relative tol of 1e-9.
    const Point3D pp0(1.0, 2.0, 3.0);
    const Vector3D pn(0.0, 0.0, 1.0);
    checkClose("distPlane.on", 0.0, distPlane(Point3D(5.0, -7.0, 3.0), pp0, pn), 1e-12);
    checkClose("distPlane.off0.5", 0.5, distPlane(Point3D(5.0, -7.0, 3.5), pp0, pn), 1e-9);

    const Point3D cpt(0.0, 0.0, 0.0);
    const Vector3D cax(0.0, 0.0, 1.0);
    const double cr = 7.0;
    checkClose("distCylinder.on", 0.0, distCylinder(Point3D(cr, 0.0, 3.0), cpt, cax, cr), 1e-12);
    checkClose("distCylinder.off0.5", 0.5, distCylinder(Point3D(cr + 0.5, 0.0, 3.0), cpt, cax, cr), 1e-9);
    checkClose("distCylinder.off0.5b", 0.5, distCylinder(Point3D(-(cr + 0.5), 0.0, -4.0), cpt, cax, cr), 1e-9);

    const Point3D sc(1.0, 1.0, 1.0);
    const double sr = 10.0;
    checkClose("distSphere.on", 0.0, distSphere(Point3D(11.0, 1.0, 1.0), sc, sr), 1e-12);
    checkClose("distSphere.off0.5", 0.5, distSphere(Point3D(11.5, 1.0, 1.0), sc, sr), 1e-9);

    // Cone half-angle 60 deg; outward normal is (cos a, 0, -sin a).
    const Point3D apex(0.0, 0.0, 0.0);
    const Vector3D cax2(0.0, 0.0, 1.0);
    const double alpha = 60.0 * kPi / 180.0;
    const double tanA = std::tan(alpha);
    const double nOx = std::cos(alpha), nOz = -std::sin(alpha);

    // Near apex (h=1) and far from it (h=10): both must give the same offset behavior.
    {
        const double h = 1.0, rho = h * tanA;
        checkClose("distCone.nearApex.on", 0.0,
                   distCone(Point3D(rho, 0.0, h), apex, cax2, alpha), 1e-12);
        const Point3D q(rho + 0.5 * nOx, 0.0, h + 0.5 * nOz);
        checkClose("distCone.nearApex.off0.5", 0.5, distCone(q, apex, cax2, alpha), 1e-9);
    }
    {
        const double h = 10.0, rho = h * tanA;
        checkClose("distCone.farApex.on", 0.0,
                   distCone(Point3D(rho, 0.0, h), apex, cax2, alpha), 1e-12);
        const Point3D q(rho + 0.5 * nOx, 0.0, h + 0.5 * nOz);
        checkClose("distCone.farApex.off0.5", 0.5, distCone(q, apex, cax2, alpha), 1e-9);
    }

    // Torus: outer and inner equator.
    const Point3D tc(0.0, 0.0, 0.0);
    const Vector3D tax(0.0, 0.0, 1.0);
    const double R = 100.0, rr = 10.0;
    checkClose("distTorus.outer.on", 0.0, distTorus(Point3D(R + rr, 0.0, 0.0), tc, tax, R, rr), 1e-12);
    checkClose("distTorus.outer.off0.5", 0.5,
               distTorus(Point3D(R + rr + 0.5, 0.0, 0.0), tc, tax, R, rr), 1e-9);
    checkClose("distTorus.inner.on", 0.0, distTorus(Point3D(R - rr, 0.0, 0.0), tc, tax, R, rr), 1e-12);
    checkClose("distTorus.inner.off0.5", 0.5,
               distTorus(Point3D(R - rr - 0.5, 0.0, 0.0), tc, tax, R, rr), 1e-9);
}

} // namespace

// Entry: runs the self-checks in order; a thrown exception (incl. init/fini failure) yields exit code 2.
int main()
{
    try
    {
        init();

        testSymEig3();
        testConeAxis();
        testCylinderAxis();
        testSolveLinear3();
        testLineConcurrency();
        testPluckerAxis();
        testFitCircle2D();
        testDistances();

        fini();
    }
    catch (const std::exception& e)
    {
        std::cerr << "[FAIL] exception: " << e.what() << std::endl;
        return 2;
    }
    catch (...)
    {
        std::cerr << "[FAIL] unknown exception" << std::endl;
        return 2;
    }

    std::cout << "[SUMMARY] total=" << g_total << " failed=" << g_failed << std::endl;
    return (g_failed == 0) ? 0 : 1;
}
