// sgkqcore: closed-form numerical core (symmetric eigen-solvers, small linear systems,
// axis/circle fitting, point-to-quadric distances). Touches no SGGK evaluation/topology
// types, so it is unit-testable on synthetic point sets.
//
// Units: angles in radians, lengths in model units (mm); no conversion is done here.
// Fitted directions are returned normalized; solvers report success via bool and leave
// outputs undefined on failure, so callers must check the return before reading outputs.
#include "sgkquadric/math.hpp"

#include <algorithm>
#include <cmath>

namespace sgkq
{

namespace
{

// Relative threshold for "singular/ill-conditioned" pivots, scaled by max |A| over the system.
const double kPivotEps = 1e-13;

// Cyclic Jacobi eigen-decomposition of a symmetric matrix; the single eigen routine in this file.
//
// @param a [in,out] n x n row-major symmetric matrix, DESTROYED in place: on return its diagonal
//        holds the eigenvalues and the off-diagonals are driven to zero (rounding residue only).
//        Not symmetric-checked; asymmetric input silently produces wrong results.
// @param v [out] n x n row-major; column j is the j-th eigenvector (orthogonal, unit columns).
// @pre `a` must be exactly symmetric (callers only feed symmetric accumulations like sum w w^T).
// O(n^3) time, O(1) extra space, no allocation, no exceptions.
//
// Convergence per sweep (off = sum of squared upper off-diagonals, diag = sum of squared diagonals):
//   off <= 1e-300        absolute underflow guard (elements ~1e-150, below double range).
//   off <= 1e-30 * diag  relative test giving ~1e-15 relative eigenvalue perturbation; scale-free
//                        so it holds for covariance matrices (length^2) at any model size.
// maxSweeps=100 is an unreachable guard: 3x3 needs ~3-5 sweeps, 6x6 ~6-9.
// The rotation uses the smaller-magnitude root of t^2 + 2*theta*t - 1 = 0 so |phi| <= pi/4,
// which is numerically the best-conditioned choice (large rotations amplify rounding).
void jacobiEigen(double* a, double* v, int n)
{
    // v starts as the identity; accumulated rotations make its columns the eigenvectors.
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            v[i * n + j] = (i == j) ? 1.0 : 0.0;

    const int maxSweeps = 100;
    for (int sweep = 0; sweep < maxSweeps; ++sweep)
    {
        double off = 0.0;
        double diag = 0.0;
        for (int p = 0; p < n; ++p)
        {
            diag += a[p * n + p] * a[p * n + p];
            for (int q = p + 1; q < n; ++q)
                off += a[p * n + q] * a[p * n + q];
        }
        if (off <= 1e-300) break;
        if (off <= 1e-30 * diag) break;

        for (int p = 0; p < n; ++p)
        {
            for (int q = p + 1; q < n; ++q)
            {
                const double apq = a[p * n + q];
                if (std::fabs(apq) <= 1e-300) continue;

                const double app = a[p * n + p];
                const double aqq = a[q * n + q];
                // Smaller-magnitude root of t^2 + 2*theta*t - 1 = 0 annihilates a[p][q].
                const double theta = (aqq - app) / (2.0 * apq);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                                 (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0);
                const double s = t * c;

                for (int i = 0; i < n; ++i)
                {
                    if (i == p || i == q) continue;
                    const double aip = a[i * n + p];
                    const double aiq = a[i * n + q];
                    a[i * n + p] = c * aip - s * aiq;
                    a[p * n + i] = a[i * n + p];
                    a[i * n + q] = s * aip + c * aiq;
                    a[q * n + i] = a[i * n + q];
                }
                a[p * n + p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
                a[q * n + q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
                // The pair is zeroed exactly by construction; write 0 to avoid residue build-up.
                a[p * n + q] = 0.0;
                a[q * n + p] = 0.0;

                for (int i = 0; i < n; ++i)
                {
                    const double vip = v[i * n + p];
                    const double viq = v[i * n + q];
                    v[i * n + p] = c * vip - s * viq;
                    v[i * n + q] = s * vip + c * viq;
                }
            }
        }
    }
}

} // namespace

//! @brief Eigen-decomposition of a symmetric 3x3 matrix; eigenvalues ascending.
//! @param [in] m symmetric, row-major; copied internally and left untouched.
//! @return val[k] ascending eigenvalues, vec[k] the aligned unit eigenvectors. val[0]/vec[0] is
//!         the smallest pair, which minVarDirection relies on.
//! @pre m must be symmetric (not checked); callers pass accumulations such as sum of outer products.
//! Eigenvector sign is arbitrary (only direction is used downstream). Degenerate input (all-zero
//! matrix) yields equal zero eigenvalues with orthogonal but arbitrary directions; the caller
//! judges validity by geometry. Eigenvalues carry the input units squared, eigenvectors are unitless.
//! Implementation: copy, jacobiEigen, sort indices by diagonal, then back-fill val/vec together.
SymEig3 symEig3(const double m[3][3])
{
    double a[9];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            a[i * 3 + j] = m[i][j];

    double v[9];
    jacobiEigen(a, v, 3);

    int idx[3] = { 0, 1, 2 };
    std::sort(idx, idx + 3, [&](int p, int q) { return a[p * 3 + p] < a[q * 3 + q]; });

    SymEig3 r;
    for (int k = 0; k < 3; ++k)
    {
        const int j = idx[k];
        r.val[k] = a[j * 3 + j];
        r.vec[k] = Vector3D(v[0 * 3 + j], v[1 * 3 + j], v[2 * 3 + j]);
    }
    return r;
}

//! @brief Eigen-decomposition of a symmetric 6x6 matrix; eigenvalues ascending.
//! @param [in] m symmetric, row-major; copied internally and left untouched.
//! @return val[k] ascending eigenvalues, vec[k] the aligned unit eigenvectors. Component i of the
//!         k-th eigenvector is vec[k][i] = v[i*6 + j] (i indexes components, not eigenpairs).
//! Used by solvePluckerAxis: the smallest eigenvector is the common axis (a, m).
//! @pre m must be symmetric (not checked). Degenerate input (rank < 6, e.g. concurrent normals)
//! gives a near-zero smallest eigenvalue but a non-unique eigenvector; solvePluckerAxis screens
//! those out by the |a| ~ 0 test. O(6^3), converges within ~9 sweeps, no allocation.
SymEig6 symEig6(const double m[6][6])
{
    double a[36];
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            a[i * 6 + j] = m[i][j];

    double v[36];
    jacobiEigen(a, v, 6);

    int idx[6] = { 0, 1, 2, 3, 4, 5 };
    std::sort(idx, idx + 6, [&](int p, int q) { return a[p * 6 + p] < a[q * 6 + q]; });

    SymEig6 r;
    for (int k = 0; k < 6; ++k)
    {
        const int j = idx[k];
        r.val[k] = a[j * 6 + j];
        for (int i = 0; i < 6; ++i)
            r.vec[k][i] = v[i * 6 + j];
    }
    return r;
}

//! @brief Smallest-eigenvector of the covariance of `v` = the constant-projection direction,
//! i.e. the axis of a cylinder or cone whose unit normals all lie in the plane normal to it.
//! @return Unit direction, or the zero vector for empty input (caller tests Length()>0).
//!
//! Covariance MUST be built from mean-subtracted samples. Using the raw second moment sum v v^T
//! instead makes the LARGEST eigenvector the axis for cones with half-angle > 54.7356 deg, which
//! rejects the whole cone class. This is locked by a regression test; do not "simplify" it away.
//! The covariance is translation/scale covariant, so no input normalization is needed.
Vector3D minVarDirection(const std::vector<Vector3D>& v)
{
    const size_t n = v.size();
    if (n == 0) return Vector3D(0.0, 0.0, 0.0);

    // First pass: mean (the covariance must be centered on it).
    double mx = 0.0, my = 0.0, mz = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        mx += v[i].X();
        my += v[i].Y();
        mz += v[i].Z();
    }
    const double inv = 1.0 / static_cast<double>(n);
    mx *= inv; my *= inv; mz *= inv;

    // C = (1/N) sum (v_i - mean)(v_i - mean)^T.
    double C[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
    for (size_t i = 0; i < n; ++i)
    {
        const double d[3] = { v[i].X() - mx, v[i].Y() - my, v[i].Z() - mz };
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                C[r][c] += d[r] * d[c];
    }
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            C[r][c] *= inv;

    // The smallest eigenvalue (ascending index 0) gives the direction of least variance.
    const SymEig3 e = symEig3(C);
    return e.vec[0];
}

//! @brief Solve the 3x3 system A x = b by partial-pivoted Gaussian elimination with back-substitution.
//! @param [out] x valid only when true is returned; undefined otherwise.
//! @return false when max|A|,|b| <= 0 (all-zero input) or a pivot falls to kPivotEps*maxAbs
//!         (singular / ill-conditioned); the threshold scales with the input magnitude so it is
//!         dimensionally consistent for mm- or m-scale models. A and b are not modified.
//! Preferred over Cramer's rule, whose explicit formulas are unstable near singularity and give
//! no natural degeneracy test. Consumers: fitCircle2D (Kasa normal equations) and the concurrency
//! solve used for the sphere center and cone apex.
bool solveLinear3(const double A[3][3], const double b[3], double x[3])
{
    double a[3][4];
    double maxAbs = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            a[i][j] = A[i][j];
            maxAbs = std::max(maxAbs, std::fabs(A[i][j]));
        }
        a[i][3] = b[i];
        maxAbs = std::max(maxAbs, std::fabs(b[i]));
    }
    if (maxAbs <= 0.0) return false;

    // Pivot threshold tracked to the largest input element, so its units match A and b.
    const double thresh = kPivotEps * maxAbs;

    // Forward elimination with partial pivoting.
    for (int k = 0; k < 3; ++k)
    {
        // Largest-magnitude row in the column limits rounding growth.
        int piv = k;
        for (int i = k + 1; i < 3; ++i)
            if (std::fabs(a[i][k]) > std::fabs(a[piv][k])) piv = i;

        if (std::fabs(a[piv][k]) <= thresh) return false; // singular / ill-conditioned
        if (piv != k)
            for (int j = 0; j < 4; ++j) std::swap(a[k][j], a[piv][j]);

        for (int i = k + 1; i < 3; ++i)
        {
            const double f = a[i][k] / a[k][k];
            for (int j = k; j < 4; ++j)
                a[i][j] -= f * a[k][j];
        }
    }

    // Back-substitution.
    for (int i = 2; i >= 0; --i)
    {
        double s = a[i][3];
        for (int j = i + 1; j < 3; ++j) s -= a[i][j] * x[j];
        x[i] = s / a[i][i];
    }
    return true;
}

//! @brief Least-squares concurrency point: minimizer of the summed squared distances to the
//! lines {p_i, n_i}; used for the sphere center and cone apex (all their normals are concurrent).
//! @param [in] p a point on each line; @param [in] n the line direction (normalized internally).
//! @return false on size mismatch, fewer than 3 points, a zero direction, or a singular system.
//!
//! Normal equations: minimize sum ||(I - n_i n_i^T)(X - p_i)||^2 gives A X = b with
//! A = sum M_i, b = sum M_i p_i, M_i = I - n_i n_i^T. M_i projects onto the plane normal to n_i,
//! so A is symmetric (positive definite unless degenerate); the solve is delegated to solveLinear3.
//! Not applicable to the torus axis, whose normals share no common point; use solvePluckerAxis.
bool solveLineConcurrency(const std::vector<Point3D>& p, const std::vector<Vector3D>& n, Point3D& out)
{
    if (p.size() != n.size() || p.size() < 3) return false;

    double A[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
    double b[3] = { 0.0, 0.0, 0.0 };

    for (size_t i = 0; i < p.size(); ++i)
    {
        Vector3D ni = n[i];
        if (!ni.Normalize()) return false; // zero direction

        const double nx = ni.X(), ny = ni.Y(), nz = ni.Z();
        // M_i = I - n_i n_i^T; accumulate M_i into A and M_i p_i into b.
        const double M[3][3] = {
            { 1.0 - nx * nx, -nx * ny,       -nx * nz },
            { -nx * ny,      1.0 - ny * ny,  -ny * nz },
            { -nx * nz,      -ny * nz,       1.0 - nz * nz }
        };
        const double pi[3] = { p[i].X(), p[i].Y(), p[i].Z() };

        for (int r = 0; r < 3; ++r)
        {
            for (int c = 0; c < 3; ++c) A[r][c] += M[r][c];
            double v = 0.0;
            for (int c = 0; c < 3; ++c) v += M[r][c] * pi[c];
            b[r] += v;
        }
    }

    double x[3];
    if (!solveLinear3(A, b, x)) return false;
    out = Point3D(x[0], x[1], x[2]);
    return true;
}

//! @brief Common axis intersecting all lines {p_i, n_i} (least squares); used for the torus axis.
//! @param [in] p a point on each line; @param [in] n the line direction (normalized internally).
//! @param [out] axisDir unit axis (sign arbitrary); @param [out] axisPoint a representative point.
//! @return false on size mismatch, fewer than 3 points, a zero direction, or |a| ~ 0 (no common
//!         axis / deficient null space).
//!
//! Derivation. With unit axis a and axis point Q (Q perp a), the Plucker moment is m = Q x a and
//! Q = a x m. Writing each row as w_i = [(p_i x n_i)^T | n_i^T], the incidence condition
//! n_i . m + (p_i x n_i) . a = 0 becomes w_i . (a, m) = 0, so (a, m) is the smallest eigenvector
//! of sum w_i w_i^T (symEig6 gives it at index 0). The eigenvector carries an unknown scalar
//! lambda with a = lambda a_true, m = lambda m_true, hence a x m = lambda^2 Q, i.e. Q = a x m / |a|^2.
//! The final axis point is Q shifted along the axis to the mean axial projection of the samples,
//! which centers the axial origin for a trimmed torus (detect aligns h to the tube-section plane).
//! Unlike solveLineConcurrency this requires only a common axis, not a common point.
bool solvePluckerAxis(const std::vector<Point3D>& p, const std::vector<Vector3D>& n,
                      Vector3D& axisDir, Point3D& axisPoint)
{
    if (p.size() != n.size() || p.size() < 3) return false;

    double M[6][6];
    for (int r = 0; r < 6; ++r)
        for (int c = 0; c < 6; ++c)
            M[r][c] = 0.0;

    for (size_t i = 0; i < p.size(); ++i)
    {
        Vector3D ni = n[i];
        if (!ni.Normalize()) return false;

        const Vector3D pi(p[i].X(), p[i].Y(), p[i].Z());
        const Vector3D mi = pi.Cross(ni);
        // Plucker row w_i = [m_i^T | n_i^T], with m_i = p_i x n_i.
        const double w[6] = { mi.X(), mi.Y(), mi.Z(), ni.X(), ni.Y(), ni.Z() };

        for (int r = 0; r < 6; ++r)
            for (int c = 0; c < 6; ++c)
                M[r][c] += w[r] * w[c];
    }

    const SymEig6 e = symEig6(M);
    const double* v = e.vec[0]; // smallest eigenvector = (a, m)

    const Vector3D av(v[0], v[1], v[2]);
    const Vector3D mv(v[3], v[4], v[5]);

    const double la = av.Length();
    if (la <= 1e-300) return false;

    const Vector3D aHat = av / la;
    const Vector3D q0 = av.Cross(mv) / (la * la); // nearest point on the axis to the origin

    double tsum = 0.0;
    for (size_t i = 0; i < p.size(); ++i)
    {
        const Vector3D d = p[i] - q0; // Point3D - Point3D -> Vector3D
        tsum += d.Dot(aHat);
    }
    const double tmean = tsum / static_cast<double>(p.size());

    axisDir = aHat;
    const Vector3D axisPt = q0 + aHat * tmean;
    axisPoint = Point3D(axisPt.X(), axisPt.Y(), axisPt.Z());
    return true;
}

//! @brief Algebraic (Kasa) least-squares circle fit in the plane: x^2+y^2+D x+E y+F=0.
//! @param [in] pts planar points {(x_i, y_i)}, in mm (or axis-frame rho, h, also mm).
//! @param [out] cx, cy center; @param [out] r radius (> 0). Outputs valid only on true.
//! @return false on fewer than 3 points, a singular normal system (collinear points), or r^2 <= 0.
//!
//! Center = (-D/2, -E/2), radius = sqrt(cx^2+cy^2-F). Algebraic least squares is used instead of an
//! iterative geometric fit because exact on-circle data (the analytic inputs have no noise) is
//! recovered to machine precision with no iteration, matching the closed-form, deterministic policy.
//! It is weighted by (x^2+y^2), so a poorly distributed arc leans slightly toward larger radii; for
//! the "circle on a genuine torus" usage the bias is far below tolerance. Free-form (rho, h) scatter
//! instead leaves a residual that is a significant fraction of its spread, which is how the torus is
//! separated from a general free-form surface downstream. The normal equations scale with the
//! coordinates, so callers normalize the residual by the bounding-box diagonal rather than here.
bool fitCircle2D(const std::vector<std::array<double, 2>>& pts, double& cx, double& cy, double& r)
{
    const double n = static_cast<double>(pts.size());
    if (pts.size() < 3) return false;

    double Sxx = 0.0, Sxy = 0.0, Syy = 0.0, Sx = 0.0, Sy = 0.0, Sz = 0.0, Sxz = 0.0, Syz = 0.0;
    for (size_t i = 0; i < pts.size(); ++i)
    {
        const double x = pts[i][0];
        const double y = pts[i][1];
        const double z = x * x + y * y;
        Sxx += x * x;
        Sxy += x * y;
        Syy += y * y;
        Sx += x;
        Sy += y;
        Sz += z;
        Sxz += x * z;
        Syz += y * z;
    }

    // Normal system: [Sxx Sxy Sx; Sxy Syy Sy; Sx Sy n] [D E F]^T = -[Sxz Syz Sz]^T
    const double A[3][3] = {
        { Sxx, Sxy, Sx },
        { Sxy, Syy, Sy },
        { Sx,  Sy,  n  }
    };
    const double b[3] = { -Sxz, -Syz, -Sz };

    double sol[3];
    if (!solveLinear3(A, b, sol)) return false; // degenerate (collinear, etc.)

    const double D = sol[0], E = sol[1], F = sol[2];
    const double ccx = -D / 2.0;
    const double ccy = -E / 2.0;
    const double rr = ccx * ccx + ccy * ccy - F;
    if (rr <= 0.0) return false;

    cx = ccx;
    cy = ccy;
    r = std::sqrt(rr);
    return true;
}

//! @brief Absolute distance from a point to an infinite plane: |(q - p0) . n_hat|.
//! Zero normal (degenerate) returns 0.0; the absolute value drops the side information, which this
//! consistency check does not need. Used by the plane-precision verification, not by classification.
double distPlane(const Point3D& q, const Point3D& p0, const Vector3D& n)
{
    Vector3D nn = n;
    if (!nn.Normalize()) return 0.0; // degenerate normal
    return std::fabs((q - p0).Dot(nn));
}

//! @brief Orthogonal distance from a point to an infinite cylinder of radius r: |rho - r|,
//! with rho the perpendicular distance to the axis. No end caps (trimming is a sampling concern).
//! Zero axis direction returns 0.0. No apex degeneracy: the formula is continuous everywhere.
double distCylinder(const Point3D& q, const Point3D& axisPt, const Vector3D& axisDir, double r)
{
    Vector3D a = axisDir;
    if (!a.Normalize()) return 0.0;

    // Split (q - axisPt) into axial h and radial parts; distance = | |radial| - r |.
    const Vector3D d = q - axisPt;
    const double h = d.Dot(a);
    const Vector3D radial = d - a * h;
    return std::fabs(radial.Length() - r);
}

//! @brief Orthogonal distance from a point to a single-napped infinite cone (apex handled).
//! @param halfAngle half-angle alpha in radians, in the open interval (0, pi/2).
//!
//! In the axial half-plane (rho, h) the generator direction is (sin a, cos a) and the foot
//! parameter is t = rho*sin a + h*cos a. For t >= 0 the foot lies on the generator, giving
//! |rho*cos a - h*sin a|; for t < 0 the foot falls behind the apex and the nearest point on the
//! single-napped cone degenerates to the apex, giving |q - apex|. Without that branch the result
//! would be the distance to the opposite nappe. Zero axis direction returns 0.0.
double distCone(const Point3D& q, const Point3D& apex, const Vector3D& axisDir, double halfAngle)
{
    Vector3D a = axisDir;
    if (!a.Normalize()) return 0.0;

    const Vector3D d = q - apex;
    const double h = d.Dot(a);
    const double rho = (d - a * h).Length();

    const double sa = std::sin(halfAngle);
    const double ca = std::cos(halfAngle);

    const double t = rho * sa + h * ca;
    if (t < 0.0) return std::sqrt(rho * rho + h * h);
    return std::fabs(rho * ca - h * sa);
}

//! @brief Orthogonal distance from a point to a sphere: | |q - center| - r |.
//! Smooth everywhere; no branches or edge cases. At the center (norm 0) it returns r, the natural
//! limit, without a special case.
double distSphere(const Point3D& q, const Point3D& center, double r)
{
    const Vector3D d = q - center;
    return std::fabs(d.Length() - r);
}

//! @brief Orthogonal distance from a point to a torus: |sqrt((rho-R)^2 + h^2) - r|,
//! with rho the distance to the axis and h the axial coordinate. The torus is the distance locus
//! of the tube-section circle of center (R, 0) and radius r in the axial (rho, h) frame, so the
//! 3D distance reduces to a 2D circle distance, which is more stable than a componentwise solve.
//! @param R major radius, r tube radius, both > 0.
//! h MUST be measured from the TUBE-section center plane: for a trimmed torus callers must shift it
//! (H = h - cy) or the surface is wrongly rejected. Zero axis direction returns 0.0.
double distTorus(const Point3D& q, const Point3D& center, const Vector3D& axisDir, double R, double r)
{
    Vector3D a = axisDir;
    if (!a.Normalize()) return 0.0;

    const Vector3D d = q - center;
    const double h = d.Dot(a);
    const double rho = (d - a * h).Length();
    const double drho = rho - R;
    return std::fabs(std::sqrt(drho * drho + h * h) - r);
}

} // namespace sgkq
