// sgkquadric parameter fitting: per-class closed-form initial estimates, Gauss-Newton
// refinement, and closed-form distance verification. Inputs are pre-sampled point/normal
// pairs (FitSample); no SDK geometry objects are handled here. Units are mm and radians.
#include "sgkquadric/fit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <sstream>
#include <vector>

namespace sgkq
{

namespace
{

const double kPi = 3.14159265358979323846;

//! Cone half-angle valid range (rad): sin(alpha) to 0 (near line) and alpha to pi/2 (near plane)
//! are both ill-conditioned degeneracies and must be rejected.
const double kMinHalfAngle = 1e-6;
const double kMaxHalfAngle = kPi * 0.5 - 1e-6;

//! Constant-projection constant c = mean(n.a) = +/-sin(alpha); |c| below this threshold means the
//! cone has degenerated to a line/cylinder (indistinguishable from a cylinder) and is rejected.
const double kMinSinAlpha = 1e-6;

//! Cone self-consistency tolerance (rad): the angle reconstructed from generator lines,
//! acos(mean(u_i.a)), must match alpha = asin(|mean(n_i.a)|) to within this margin.
const double kConeSelfConsistencyTol = 1e-6;

//! Constant-projection relative test: lambda_min(Cov) <= kConstProjEps * lambda_max(Cov) (1e-3 is a wide margin).
const double kConstProjEps = 1e-3;

//! Absolute lower bound on the normal spread within the plane orthogonal to the axis.
const double kMinNormalSpread = 1e-10;

//! Relative degeneracy threshold, multiplied by the point-set scale (PointScale).
const double kMinRelSize = 1e-9;

//! Minimum identifiable sample count per type (closed-form identifiability floor).
const size_t kMinSamplesPlane = 3;
const size_t kMinSamplesCylinder = 6;
const size_t kMinSamplesCone = 8;
const size_t kMinSamplesSphere = 6;
const size_t kMinSamplesTorus = 12;


inline Vector3D AsVec(const Point3D& p) { return Vector3D(p.X(), p.Y(), p.Z()); }
inline Point3D  AsPnt(const Vector3D& v) { return Point3D(v.X(), v.Y(), v.Z()); }

std::string ToStr(double v)
{
    std::ostringstream os;
    os << v;
    return os.str();
}

inline bool Unit(const Vector3D& v, Vector3D& out)
{
    out = v;
    return out.Normalize();
}

//! @brief Build a right-handed frame (e1, e2, axis) orthonormal to `axis`.
inline void OrthoFrame(const Vector3D& axis, Vector3D& e1, Vector3D& e2)
{
    Vector3D a;
    if (!Unit(axis, a)) a = Vector3D(0.0, 0.0, 1.0);   // zero axis -> fall back to +Z

    Vector3D ref(0.0, 0.0, 1.0);   // reference axis kept well off a (|ref . a| <= 0.9) to avoid a degenerate cross product
    if (std::fabs(a.Dot(ref)) > 0.9) ref = Vector3D(1.0, 0.0, 0.0);

    e1 = ref.Cross(a); e1.Normalize();
    e2 = a.Cross(e1);  e2.Normalize();
}

//! @brief Point-set scale: max distance from the samples to their centroid (base of relative thresholds).
double PointScale(const std::vector<FitSample>& samples)
{
    if (samples.empty()) return 0.0;

    Vector3D mean(0.0, 0.0, 0.0);
    for (const FitSample& s : samples) mean += AsVec(s.p);
    mean /= static_cast<double>(samples.size());

    double r = 0.0;
    for (const FitSample& s : samples)
        r = std::max(r, (AsVec(s.p) - mean).Length());
    return r;
}

//! @brief Collect valid unit normals, skipping zero normals; needs at least 3.
bool CollectNormals(const std::vector<FitSample>& samples, std::vector<Vector3D>& out)
{
    out.clear();
    out.reserve(samples.size());
    for (const FitSample& s : samples)
    {
        Vector3D n;
        if (Unit(s.n, n)) out.push_back(n);
    }
    return out.size() >= 3;
}

//! @brief Collect paired points + unit normals (the concurrency/colinearity solvers require equal lengths).
void CollectPN(const std::vector<FitSample>& samples,
               std::vector<Point3D>& p, std::vector<Vector3D>& n)
{
    p.clear();
    n.clear();
    p.reserve(samples.size());
    n.reserve(samples.size());
    for (const FitSample& s : samples)
    {
        Vector3D u;
        if (!Unit(s.n, u)) continue;
        p.push_back(s.p);
        n.push_back(u);
    }
}

//! @brief Covariance eigen-decomposition (ascending) of a unit-vector set; used only for degeneracy tests.
bool CovEig3(const std::vector<Vector3D>& v, SymEig3& e)
{
    const size_t m = v.size();
    if (m == 0) return false;

    Vector3D mean(0.0, 0.0, 0.0);
    for (const Vector3D& x : v) mean += x;
    mean /= static_cast<double>(m);

    double C[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
    for (const Vector3D& x : v)
    {
        const double d[3] = { x.X() - mean.X(), x.Y() - mean.Y(), x.Z() - mean.Z() };
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                C[r][c] += d[r] * d[c] / static_cast<double>(m);
    }
    e = symEig3(C);
    return true;
}

Vector3D Centroid(const std::vector<Point3D>& p)
{
    Vector3D c(0.0, 0.0, 0.0);
    if (p.empty()) return c;
    for (const Point3D& q : p) c += AsVec(q);
    c /= static_cast<double>(p.size());
    return c;
}


//! Plane: normal = normalized mean(n_i); d = mean(p_i . n).
bool FitPlane(const std::vector<FitSample>& samples, SurfaceParams& out, std::string& reason)
{
    std::vector<Vector3D> n;
    if (!CollectNormals(samples, n) || samples.size() < kMinSamplesPlane)
    {
        reason = "plane fit: insufficient valid samples or normals";
        return false;
    }

    Vector3D sum(0.0, 0.0, 0.0);
    for (const Vector3D& x : n) sum += x;

    Vector3D nrm;
    const Vector3D meanN = sum / static_cast<double>(n.size());
    if (!Unit(meanN, nrm) || meanN.Length() < 0.5)
    {
        reason = "plane fit: mean normal length = " + ToStr(meanN.Length()) +
                 " < 0.5 (normal spread too large / mutually canceling), samples not from a single plane";
        return false;
    }

    double d = 0.0;
    for (const FitSample& s : samples) d += AsVec(s.p).Dot(nrm);
    d /= static_cast<double>(samples.size());

    out.type = Canonical::Plane;
    out.axisDir = nrm;
    out.axisPoint = AsPnt(nrm * d);
    out.radius = 0.0;
    out.minorRadius = 0.0;
    out.halfAngle = 0.0;
    out.apex = out.axisPoint;
    reason.clear();
    return true;
}

//! Cylinder: axis = smallest eigenvector of Cov(n); axis position and radius from a 2D circle fit of
//! the projected points. Never use the projected-point mean as the center: for a trimmed arc it is not the circle center.
bool FitCylinder(const std::vector<FitSample>& samples, SurfaceParams& out, std::string& reason)
{
    if (samples.size() < kMinSamplesCylinder)
    {
        reason = "cylinder fit: insufficient samples (need ≥ 6)";
        return false;
    }

    std::vector<Vector3D> n;
    if (!CollectNormals(samples, n))
    {
        reason = "cylinder fit: fewer than 3 valid normals";
        return false;
    }

    Vector3D axis;
    if (!Unit(minVarDirection(n), axis))
    {
        reason = "cylinder fit: degenerate axis direction (zero normal covariance)";
        return false;
    }

    // Degeneracy: (1) the constant-projection condition n.a = 0 must hold; (2) the normals must have
    // a non-degenerate spread in the plane orthogonal to a, else the axis is unidentifiable.
    SymEig3 e;
    if (!CovEig3(n, e))
    {
        reason = "cylinder fit: normal statistics failed";
        return false;
    }
    if (!(e.val[0] <= kConstProjEps * e.val[2]))
    {
        reason = "cylinder fit: normals violate the constant-projection condition n·a ≡ 0, not a cylinder";
        return false;
    }
    if (e.val[1] <= kMinNormalSpread)
    {
        reason = "cylinder fit: insufficient normal spread (near plane or near line), axis direction not identifiable";
        return false;
    }

    // Project the points onto the plane orthogonal to a and fit a 2D circle: center gives the axis point, radius gives r.
    Vector3D e1, e2;
    OrthoFrame(axis, e1, e2);

    std::vector<std::array<double, 2>> uv;
    uv.reserve(samples.size());
    for (const FitSample& s : samples)
    {
        const Vector3D p = AsVec(s.p);
        uv.push_back(std::array<double, 2>{ p.Dot(e1), p.Dot(e2) });
    }

    double cx = 0.0, cy = 0.0, r = 0.0;
    if (!fitCircle2D(uv, cx, cy, r))
    {
        reason = "cylinder fit: degenerate 2D circle fit of projected points (collinear/single position), radius → 0";
        return false;
    }

    const double scale = PointScale(samples);
    if (r <= kMinRelSize * std::max(scale, 1e-300))
    {
        reason = "cylinder fit: degenerate radius (≈ 0), input is nearly a straight line";
        return false;
    }

    out.type = Canonical::Cylinder;
    out.axisDir = axis;
    out.axisPoint = AsPnt(e1 * cx + e2 * cy);
    out.radius = r;
    out.minorRadius = 0.0;
    out.halfAngle = 0.0;
    out.apex = out.axisPoint;
    reason.clear();
    return true;
}

// Cone. The axis MUST be the SMALLEST eigenvector of Cov(n): the largest eigenvector is
// perpendicular to the true axis once half-angle > 54.7356 deg (e.g. a 60 deg countersink).
// Half-angle comes straight from the normals: sin(alpha) = |mean(n . a)|.
// Apex from generator lines, not normal lines: the normal lines are NOT concurrent - they meet
// the axis at h_i / cos^2(alpha). Only generator lines (apex to sample) pass through the apex.
// Generator direction g_i = cos(alpha)*a + sin(alpha)*rhat_i. Take X0 = concurrence solution of
// the generators, then apex = X0 + median_i(s_i)*a, with s_i the generators axial intersection.
bool FitCone(const std::vector<FitSample>& samples, SurfaceParams& out, std::string& reason)
{
    std::vector<Point3D> p;
    std::vector<Vector3D> n;
    CollectPN(samples, p, n);
    if (p.size() < kMinSamplesCone)
    {
        reason = "cone fit: insufficient valid samples (need ≥ 8)";
        return false;
    }

    Vector3D a;
    if (!Unit(minVarDirection(n), a))
    {
        reason = "cone fit: degenerate axis direction (zero normal covariance)";
        return false;
    }

    // Constant projection n.a = c must hold (variance along a ~ 0), with a non-degenerate normal
    // spread orthogonal to a, else (alpha to pi/2 or thin sampling) the axis is unidentifiable.
    SymEig3 e;
    if (!CovEig3(n, e))
    {
        reason = "cone fit: normal statistics failed";
        return false;
    }
    if (!(e.val[0] <= kConstProjEps * e.val[2]) || e.val[1] <= kMinNormalSpread)
    {
        reason = "cone fit: normals violate the constant-projection condition or have insufficient spread, not a cone (ill-conditioned)";
        return false;
    }

    double nDotSum = 0.0;
    for (const Vector3D& x : n) nDotSum += x.Dot(a);
    const double cSigned = nDotSum / static_cast<double>(n.size());

    const double sinAlpha = std::fabs(cSigned);

    // c to 0 (near line / indistinguishable from a cylinder) is ill-conditioned; reject.
    if (sinAlpha <= kMinSinAlpha)
    {
        reason = "cone fit: sin α = |mean(n·a)| = " + ToStr(sinAlpha) +
                 " ≈ 0, cone degenerates to a line-cone/cylinder (ill-conditioned), rejected";
        return false;
    }

    // Orient a so mean(n.a) < 0: this keeps the generator directions on the sample side (no flip).
    if (cSigned > 0.0) a = -a;

    const double alpha = std::asin(std::min(1.0, sinAlpha));
    if (!(alpha > kMinHalfAngle && alpha < kMaxHalfAngle))
    {
        reason = "cone fit: degenerate half-angle α = " + ToStr(alpha) +
                 " rad (α → 0 approaches a line-cone / α → π/2 approaches a plane, ill-conditioned)";
        return false;
    }
    const double sinA = std::sin(alpha);
    const double cosA = std::cos(alpha);

    std::vector<Point3D> pc;
    std::vector<Vector3D> g;
    pc.reserve(p.size());
    g.reserve(p.size());
    for (size_t i = 0; i < p.size(); ++i)
    {
        Vector3D rad = n[i] - a * n[i].Dot(a);
        if (!rad.Normalize()) continue;
        g.push_back(a * cosA + rad * sinA);
        pc.push_back(p[i]);
    }
    if (pc.size() < kMinSamplesCone)
    {
        reason = "cone fit: failed to rebuild generator-line directions (insufficient usable samples)";
        return false;
    }

    Point3D X0;
    if (!solveLineConcurrency(pc, g, X0))
    {
        reason = "cone fit: singular generator concurrency solve (apex not identifiable), not a cone";
        return false;
    }
    const Vector3D x0v = AsVec(X0);

    // Axial coordinate where generator i meets the axis: s_i = (E - cos(alpha)*D) / sin^2(alpha),
    // with D = g_i.(p_i - X0), E = a.(p_i - X0) and g_i.a = cos(alpha) for every i.
    std::vector<double> sv;
    sv.reserve(pc.size());
    for (size_t i = 0; i < pc.size(); ++i)
    {
        const Vector3D w = AsVec(pc[i]) - x0v;
        const double D = g[i].Dot(w);
        const double E = a.Dot(w);
        sv.push_back((E - cosA * D) / (sinA * sinA));
    }
    std::sort(sv.begin(), sv.end());
    const double mAxial = (sv.size() % 2 == 1)
        ? sv[sv.size() / 2]
        : 0.5 * (sv[sv.size() / 2 - 1] + sv[sv.size() / 2]);
    Point3D apex = AsPnt(x0v + a * mAxial);

    // Orient a from the apex toward the samples (mean((p_i-apex).a) > 0) to match the distance convention.
    const Vector3D apexV = AsVec(apex);
    double orientSum = 0.0;
    for (const Point3D& q : pc) orientSum += (AsVec(q) - apexV).Dot(a);
    if (orientSum < 0.0) a = -a;

    double cosSum = 0.0;
    size_t used = 0;
    for (const Point3D& q : pc)
    {
        Vector3D u = AsVec(q) - apexV;
        if (!u.Normalize()) continue;
        cosSum += u.Dot(a);
        ++used;
    }
    if (used < 3)
    {
        reason = "cone fit: degenerate generators (all samples fall on the apex)";
        return false;
    }
    const double caGen = std::min(1.0, std::max(-1.0, cosSum / static_cast<double>(used)));
    const double alphaGen = std::acos(caGen);
    if (std::fabs(alphaGen - alpha) > kConeSelfConsistencyTol)
    {
        reason = "cone fit: generator/normal self-consistency check failed (α_normal = " + ToStr(alpha) +
                 " rad inconsistent with α_generator = " + ToStr(alphaGen) + " rad)";
        return false;
    }

    out.type = Canonical::Cone;
    out.axisDir = a;
    out.apex = apex;
    out.axisPoint = apex;
    out.halfAngle = alpha;
    out.radius = 0.0;
    out.minorRadius = 0.0;
    reason.clear();
    return true;
}

//! Sphere: center = normal-line concurrency solution; r = mean(|p_i - C|).
bool FitSphere(const std::vector<FitSample>& samples, SurfaceParams& out, std::string& reason)
{
    std::vector<Point3D> p;
    std::vector<Vector3D> n;
    CollectPN(samples, p, n);
    if (p.size() < kMinSamplesSphere)
    {
        reason = "sphere fit: insufficient valid samples (need ≥ 6)";
        return false;
    }

    Point3D center;
    if (!solveLineConcurrency(p, n, center))
    {
        reason = "sphere fit: normal lines not concurrent (singular center solve), not a sphere (may be a plane/cylinder)";
        return false;
    }

    const Vector3D c = AsVec(center);
    double r = 0.0;
    for (const Point3D& q : p) r += (AsVec(q) - c).Length();
    r /= static_cast<double>(p.size());

    const double scale = PointScale(samples);
    if (r <= kMinRelSize * std::max(scale, 1e-300))
    {
        reason = "sphere fit: degenerate radius (≈ 0), all sample points collapse to one point";
        return false;
    }

    out.type = Canonical::Sphere;
    out.axisPoint = center;
    out.axisDir = Vector3D(0.0, 0.0, 1.0);
    out.radius = r;
    out.minorRadius = 0.0;
    out.halfAngle = 0.0;
    out.apex = center;
    reason.clear();
    return true;
}

//! Torus: axis (direction + position) = Plucker 6x6 null-space solution. The normal lines meet the
//! axis but the intersection point varies, so the concurrency form does not apply. R and r come from
//! a 2D circle fit of (rho_i, h_i) in the axis frame.
bool FitTorus(const std::vector<FitSample>& samples, SurfaceParams& out, std::string& reason)
{
    std::vector<Point3D> p;
    std::vector<Vector3D> n;
    CollectPN(samples, p, n);
    if (p.size() < kMinSamplesTorus)
    {
        reason = "torus fit: insufficient valid samples (need ≥ 12)";
        return false;
    }

    // Center by translation: the Plucker moment m_i = p_i x n_i grows with the distance to the origin,
    // so centering markedly improves the conditioning of the 6x6 eigen-decomposition.
    const Vector3D c0 = Centroid(p);
    std::vector<Point3D> pc;
    pc.reserve(p.size());
    for (const Point3D& q : p) pc.push_back(AsPnt(AsVec(q) - c0));

    Vector3D axis;
    Point3D axisPtLocal;
    if (!solvePluckerAxis(pc, n, axis, axisPtLocal))
    {
        reason = "torus fit: Plücker collinearity solve failed (too few points or degenerate), not a torus";
        return false;
    }
    const Point3D axisPt = AsPnt(AsVec(axisPtLocal) + c0);

    std::vector<std::array<double, 2>> rh;
    rh.reserve(p.size());
    for (const Point3D& q : p)
    {
        const Vector3D d = AsVec(q) - AsVec(axisPt);
        const double h = d.Dot(axis);
        const double rho = (d - axis * h).Length();
        rh.push_back(std::array<double, 2>{ rho, h });
    }

    double cx = 0.0, cy = 0.0, r = 0.0;
    if (!fitCircle2D(rh, cx, cy, r))
    {
        reason = "torus fit: degenerate 2D circle fit of (ρ,h) in the axis frame (collinear), not a torus";
        return false;
    }

    const double scale = PointScale(samples);
    const double relMin = kMinRelSize * std::max(scale, 1e-300);
    if (r <= relMin || std::fabs(cx) <= relMin)
    {
        reason = "torus fit: degenerate major/minor radii (≈ 0), input is not a torus";
        return false;
    }

    out.type = Canonical::Torus;
    out.axisDir = axis;
    // Move the axis point to the center plane (h = 0); R = cx.
    out.axisPoint = AsPnt(AsVec(axisPt) + axis * cy);
    out.radius = cx;
    out.minorRadius = r;
    out.halfAngle = 0.0;
    out.apex = out.axisPoint;

    const char* shape = "Doughnut (R>r)";
    if (std::fabs(cx - r) <= 1e-12 * std::max(1.0, std::fabs(r))) shape = "Apple (R=r)";
    else if (cx < r) shape = "Lemon (R<r)";
    reason = std::string("ok; torus shape=") + shape;
    return true;
}


//! Unsigned closed-form orthogonal distance from q to the fitted surface (delegates to math.hpp).
double ClosedDistance(const SurfaceParams& s, const Point3D& q)
{
    switch (s.type)
    {
    case Canonical::Plane:    return distPlane(q, s.axisPoint, s.axisDir);
    case Canonical::Cylinder: return distCylinder(q, s.axisPoint, s.axisDir, s.radius);
    case Canonical::Cone:     return distCone(q, s.apex, s.axisDir, s.halfAngle);
    case Canonical::Sphere:   return distSphere(q, s.axisPoint, s.radius);
    case Canonical::Torus:    return distTorus(q, s.axisPoint, s.axisDir, s.radius, s.minorRadius);
    default:                  return 0.0;
    }
}

//! Signed versions of the same formulas, for Gauss-Newton (the |.| kink of the unsigned form breaks a
//! numerical Jacobian): plane (p-P0).n; cylinder rho-r; sphere |p-C|-r; torus sqrt((rho-R)^2+h^2)-r; cone rho*cos(alpha)-h*sin(alpha) (|p-apex| if the foot is behind the apex).
double SignedResidual(const SurfaceParams& s, const Point3D& q)
{
    const Vector3D P = AsVec(q);

    switch (s.type)
    {
    case Canonical::Plane:
    {
        Vector3D n(0.0, 0.0, 1.0);
        Unit(s.axisDir, n);
        return (P - AsVec(s.axisPoint)).Dot(n);
    }
    case Canonical::Cylinder:
    {
        Vector3D a(0.0, 0.0, 1.0);
        Unit(s.axisDir, a);
        const Vector3D d = P - AsVec(s.axisPoint);
        const double h = d.Dot(a);
        return (d - a * h).Length() - s.radius;
    }
    case Canonical::Cone:
    {
        Vector3D a(0.0, 0.0, 1.0);
        Unit(s.axisDir, a);
        const Vector3D d = P - AsVec(s.apex);
        const double h = d.Dot(a);
        const double rho = (d - a * h).Length();
        const double sa = std::sin(s.halfAngle);
        const double ca = std::cos(s.halfAngle);
        const double t = rho * sa + h * ca;
        if (t < 0.0) return std::sqrt(rho * rho + h * h);
        return rho * ca - h * sa;
    }
    case Canonical::Sphere:
        return (P - AsVec(s.axisPoint)).Length() - s.radius;

    case Canonical::Torus:
    {
        Vector3D a(0.0, 0.0, 1.0);
        Unit(s.axisDir, a);
        const Vector3D d = P - AsVec(s.axisPoint);
        const double h = d.Dot(a);
        const double rho = (d - a * h).Length();
        const double dr = rho - s.radius;
        return std::sqrt(dr * dr + h * h) - s.minorRadius;
    }
    default:
        return 0.0;
    }
}

//! Residual vector and its squared sum (the Gauss-Newton cost).
double Residuals(const SurfaceParams& s, const std::vector<FitSample>& samples, std::vector<double>& r)
{
    r.resize(samples.size());
    double cost = 0.0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        const double d = SignedResidual(s, samples[i].p);
        r[i] = d;
        cost += d * d;
    }
    return cost;
}


//! Refinement parameterization (minimal DOF, no redundancy); one entry per type:
//!   plane    [t1,t2,d]              3 DOF: normal 2 + offset 1 (in-plane sliding is invariant).
//!   cylinder [s1,s2,d1,d2,r]         5 DOF: axis point 2 (orthogonal component only; sliding along the
//!                                    axis leaves every rho_i unchanged) + direction 2 + radius 1.
//!   cone     [t1,t2,t3,d1,d2,alpha]  6 DOF: apex 3 (full 3D; sliding changes rho = h*tan(alpha)) +
//!                                    direction 2 + half-angle 1.
//!   torus    [s1,s2,s3,d1,d2,R,r]    7 DOF: axis point 3 + direction 2 + R 1 + r 1 (the axial part of
//!                                    the axis point is the center-plane location, a real DOF).
//! Direction uses a local chart a(d) = normalize(aRef + d1*e1 + d2*e2); d stays in the chart, avoiding the +/-z singularity.
struct Chart
{
    Canonical type = Canonical::Unknown;
    Vector3D  aRef;
    Vector3D  e1, e2;
    Point3D   pRef;
    int       dim = 0;
    double    scale = 1.0;
    bool      isAngle[8] = { false, false, false, false, false, false, false, false };
    double    step[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
};

inline Vector3D DecodeDir(const Chart& c, double d1, double d2)
{
    Vector3D a = c.aRef + c.e1 * d1 + c.e2 * d2;
    if (!a.Normalize()) a = c.aRef;
    return a;
}

void Decode(const Chart& c, const double* x, SurfaceParams& out)
{
    out = SurfaceParams();
    out.type = c.type;

    switch (c.type)
    {
    case Canonical::Plane:
    {
        const Vector3D n = DecodeDir(c, x[0], x[1]);
        out.axisDir = n;
        out.axisPoint = AsPnt(n * x[2]);
        break;
    }
    case Canonical::Cylinder:
    {
        out.axisDir = DecodeDir(c, x[2], x[3]);
        out.axisPoint = AsPnt(AsVec(c.pRef) + c.e1 * x[0] + c.e2 * x[1]);
        out.radius = x[4];
        break;
    }
    case Canonical::Cone:
    {
        out.axisDir = DecodeDir(c, x[3], x[4]);
        out.apex = AsPnt(AsVec(c.pRef) + c.e1 * x[0] + c.e2 * x[1] + c.aRef * x[2]);
        out.axisPoint = out.apex;
        out.halfAngle = x[5];
        break;
    }
    case Canonical::Sphere:
    {
        out.axisPoint = AsPnt(AsVec(c.pRef) + Vector3D(x[0], x[1], x[2]));
        out.axisDir = Vector3D(0.0, 0.0, 1.0);
        out.radius = x[3];
        break;
    }
    case Canonical::Torus:
    {
        out.axisDir = DecodeDir(c, x[3], x[4]);
        out.axisPoint = AsPnt(AsVec(c.pRef) + c.e1 * x[0] + c.e2 * x[1] + c.aRef * x[2]);
        out.radius = x[5];
        out.minorRadius = x[6];
        break;
    }
    default:
        break;
    }
}

bool BuildChart(const SurfaceParams& io, const std::vector<FitSample>& samples, Chart& c, double x[8])
{
    for (int j = 0; j < 8; ++j) x[j] = 0.0;

    c.type = io.type;
    c.scale = std::max(1.0, PointScale(samples));

    switch (io.type)
    {
    case Canonical::Plane:
    {
        Vector3D n;
        if (!Unit(io.axisDir, n)) return false;
        c.dim = 3;
        c.aRef = n;
        OrthoFrame(n, c.e1, c.e2);
        c.pRef = io.axisPoint;
        x[2] = AsVec(io.axisPoint).Dot(n);
        c.isAngle[0] = c.isAngle[1] = true;
        break;
    }
    case Canonical::Cylinder:
    {
        Vector3D a;
        if (!Unit(io.axisDir, a)) return false;
        c.dim = 5;
        c.aRef = a;
        OrthoFrame(a, c.e1, c.e2);
        c.pRef = io.axisPoint;
        x[4] = io.radius;
        c.isAngle[2] = c.isAngle[3] = true;
        break;
    }
    case Canonical::Cone:
    {
        Vector3D a;
        if (!Unit(io.axisDir, a)) return false;
        c.dim = 6;
        c.aRef = a;
        OrthoFrame(a, c.e1, c.e2);
        c.pRef = io.apex;
        x[5] = io.halfAngle;
        c.isAngle[3] = c.isAngle[4] = true;
        break;
    }
    case Canonical::Sphere:
    {
        c.dim = 4;
        c.aRef = Vector3D(0.0, 0.0, 1.0);
        OrthoFrame(c.aRef, c.e1, c.e2);
        c.pRef = io.axisPoint;
        x[3] = io.radius;
        break;
    }
    case Canonical::Torus:
    {
        Vector3D a;
        if (!Unit(io.axisDir, a)) return false;
        c.dim = 7;
        c.aRef = a;
        OrthoFrame(a, c.e1, c.e2);
        c.pRef = io.axisPoint;
        x[5] = io.radius;
        x[6] = io.minorRadius;
        c.isAngle[3] = c.isAngle[4] = true;
        break;
    }
    default:
        return false;
    }

    // Difference step: angle 1e-7 rad; length 1e-6 * max(1, |x|, scale).
    for (int j = 0; j < c.dim; ++j)
        c.step[j] = c.isAngle[j] ? 1e-7 : 1e-6 * std::max(1.0, std::max(std::fabs(x[j]), c.scale));
    return true;
}

//! Dense m x m solve by partial-pivoted Gaussian elimination (A and b destroyed in place).
bool SolveDense(std::vector<double>& A, std::vector<double>& b, int m, std::vector<double>& x)
{
    double maxAbs = 0.0;
    for (int i = 0; i < m; ++i)
    {
        maxAbs = std::max(maxAbs, std::fabs(b[i]));
        for (int j = 0; j < m; ++j) maxAbs = std::max(maxAbs, std::fabs(A[(size_t)i * m + j]));
    }
    if (maxAbs <= 0.0) return false;
    const double thresh = 1e-13 * maxAbs;

    for (int k = 0; k < m; ++k)
    {
        int piv = k;
        for (int i = k + 1; i < m; ++i)
            if (std::fabs(A[(size_t)i * m + k]) > std::fabs(A[(size_t)piv * m + k])) piv = i;
        if (std::fabs(A[(size_t)piv * m + k]) <= thresh) return false;

        if (piv != k)
        {
            for (int j = 0; j < m; ++j) std::swap(A[(size_t)k * m + j], A[(size_t)piv * m + j]);
            std::swap(b[k], b[piv]);
        }
        for (int i = k + 1; i < m; ++i)
        {
            const double f = A[(size_t)i * m + k] / A[(size_t)k * m + k];
            if (f == 0.0) continue;
            for (int j = k; j < m; ++j) A[(size_t)i * m + j] -= f * A[(size_t)k * m + j];
            b[i] -= f * b[k];
        }
    }

    x.assign(m, 0.0);
    for (int i = m - 1; i >= 0; --i)
    {
        double s = b[i];
        for (int j = i + 1; j < m; ++j) s -= A[(size_t)i * m + j] * x[j];
        x[i] = s / A[(size_t)i * m + i];
    }
    return true;
}

}


bool FitSurface(Canonical type,
                const std::vector<FitSample>& samples,
                SurfaceParams& out,
                std::string& reason)
{
    out = SurfaceParams();
    out.type = type;
    reason.clear();

    switch (type)
    {
    case Canonical::Plane:    return FitPlane(samples, out, reason);
    case Canonical::Cylinder: return FitCylinder(samples, out, reason);
    case Canonical::Cone:     return FitCone(samples, out, reason);
    case Canonical::Sphere:   return FitSphere(samples, out, reason);
    case Canonical::Torus:    return FitTorus(samples, out, reason);
    case Canonical::Unknown:
    default:
        reason = "fit: target type is Unknown, cannot choose a parameterization";
        return false;
    }
}

bool RefineSurface(SurfaceParams& io,
                   const std::vector<FitSample>& samples,
                   int maxIter,
                   double tol,
                   std::string& reason)
{
    reason.clear();
    if (io.type == Canonical::Unknown)
    {
        reason = "refine: params.type == Unknown, cannot choose a parameterization";
        return false;
    }
    if (maxIter <= 0)
    {
        reason = "refine: maxIter ≤ 0 is invalid";
        return false;
    }
    if (!(tol > 0.0))
    {
        reason = "refine: tol ≤ 0 is invalid";
        return false;
    }
    if (samples.size() < 6)
    {
        reason = "refine: insufficient valid samples (need ≥ 6)";
        return false;
    }

    Chart ch;
    double x[8];
    if (!BuildChart(io, samples, ch, x))
    {
        reason = "refine: invalid initial guess (e.g. zero axis direction), cannot build a parameterization";
        return false;
    }

    const int m = ch.dim;
    const int N = static_cast<int>(samples.size());

    std::vector<double> r, rTry;
    SurfaceParams cur, cand;
    Decode(ch, x, cur);
    double cost = Residuals(cur, samples, r);

    std::vector<double> J(static_cast<size_t>(N) * m, 0.0);
    std::vector<double> G(static_cast<size_t>(m) * m, 0.0);
    std::vector<double> g(m, 0.0), b(m, 0.0), delta(m, 0.0);

    bool converged = false;
    for (int it = 0; it < maxIter; ++it)
    {
        for (int j = 0; j < m; ++j)
        {
            const double h = ch.step[j];
            double xp[8], xm[8];
            for (int k = 0; k < 8; ++k) { xp[k] = 0.0; xm[k] = 0.0; }
            for (int k = 0; k < m; ++k) { xp[k] = x[k]; xm[k] = x[k]; }
            xp[j] += h;
            xm[j] -= h;

            SurfaceParams sp, sm;
            Decode(ch, xp, sp);
            Decode(ch, xm, sm);
            for (int i = 0; i < N; ++i)
            {
                const double rp = SignedResidual(sp, samples[i].p);
                const double rm = SignedResidual(sm, samples[i].p);
                J[static_cast<size_t>(i) * m + j] = (rp - rm) / (2.0 * h);
            }
        }

        for (int a = 0; a < m; ++a)
        {
            g[a] = 0.0;
            for (int c = 0; c < m; ++c) G[static_cast<size_t>(a) * m + c] = 0.0;
        }
        for (int i = 0; i < N; ++i)
        {
            const double* Ji = &J[static_cast<size_t>(i) * m];
            for (int a = 0; a < m; ++a)
            {
                g[a] += Ji[a] * r[i];
                for (int c = 0; c < m; ++c) G[static_cast<size_t>(a) * m + c] += Ji[a] * Ji[c];
            }
        }
        for (int a = 0; a < m; ++a) b[a] = -g[a];

        if (!SolveDense(G, b, m, delta))
        {
            reason = "refine: singular normal equations (parameters not identifiable or samples insufficient for this type)";
            return false;
        }

        double stepScale = 1.0;
        double newCost = cost;
        bool accepted = false;
        for (int k = 0; k < 24; ++k)
        {
            double xt[8];
            for (int q = 0; q < 8; ++q) xt[q] = 0.0;
            for (int q = 0; q < m; ++q) xt[q] = x[q] + stepScale * delta[q];

            Decode(ch, xt, cand);
            const double c2 = Residuals(cand, samples, rTry);
            if (c2 <= cost)
            {
                newCost = c2;
                accepted = true;
                r.swap(rTry);
                break;
            }
            stepScale *= 0.5;
        }
        if (!accepted)
        {
            reason = "refine: line search failed (residual cannot decrease), final residual Σd² = " + ToStr(cost);
            return false;
        }

        double maxStep = 0.0;
        for (int q = 0; q < m; ++q)
        {
            x[q] += stepScale * delta[q];
            maxStep = std::max(maxStep, std::fabs(stepScale * delta[q]));
        }
        const double dCost = cost - newCost;
        cost = newCost;
        cur = cand;

        // Converged when the parameter step is small or the relative cost change is small.
        if (maxStep <= tol || dCost <= tol * std::max(1.0, cost))
        {
            converged = true;
            break;
        }
    }

    if (!converged)
    {
        reason = "refine: reached iteration limit maxIter = " + ToStr(static_cast<double>(maxIter)) +
                 " without meeting the convergence criterion, final residual Σd² = " + ToStr(cost);
        return false;
    }

    io = cur;
    reason.clear();
    return true;
}

AccuracyReport MeasureAccuracy(const SurfaceParams& params,
                               const std::vector<FitSample>& samples,
                               const Thresholds& th)
{
    AccuracyReport rep;

    if (samples.empty() || params.type == Canonical::Unknown)
    {
        rep.accepted = false;
        return rep;
    }

    double mx = 0.0;
    double sum = 0.0;
    for (const FitSample& s : samples)
    {
        const double d = ClosedDistance(params, s.p);
        mx = std::max(mx, d);
        sum += d;
    }

    rep.maxDist = mx;
    rep.avgDist = sum / static_cast<double>(samples.size());
    rep.accepted = (rep.maxDist <= th.maxDist);
    // boundaryCurveHausdorff / g1Before / g1After are filled in by later stages (boundary, G1).
    return rep;
}

}
