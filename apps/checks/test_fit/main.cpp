// apps/checks/test_fit/main.cpp - known-answer self-test for the parameter-fit module (spec 11-fit.md).
// Cases: five-type parameter recovery (analytic + ToBSpline sampling paths), 60-degree clipped arc,
// refinement under 1e-4 noise, closed-form maxDist < 1e-12, and rejection of degenerate inputs.
// Prints "[TEST] <name> expected=<> actual=<> <PASS|FAIL>"; any FAIL -> exit 1 (2 on exception).
#include "sgkquadric/fit.hpp"                              // under test: FitSurface/RefineSurface/MeasureAccuracy/FitSample

#include <Foundation/init.h>                               // init()/fini(): SDK global lifecycle
#include <GeomBase/Dir3D.h>                                // unit direction vector
#include <GeomBase/Interval.h>                             // 1D closed interval (parameter domain)
#include <GeomBase/Number.h>                               // basic numeric types
#include <GeomBase/Ucs3D.h>                                // local coordinate system (axis/origin source)
#include <GeomBase/UVRange.h>                              // U/V parameter domain
#include <Geometry/3D/Surface/BSplineSurface.h>            // ToBSpline product (path b)
#include <Geometry/3D/Surface/ConeSurface.h>               // analytic cone
#include <Geometry/3D/Surface/CylinderSurface.h>           // analytic cylinder
#include <Geometry/3D/Surface/Plane.h>                     // analytic plane
#include <Geometry/3D/Surface/SphereSurface.h>             // analytic sphere
#include <Geometry/3D/Surface/TorusSurface.h>              // analytic torus

#include <algorithm>                                       // std::max
#include <cmath>                                           // trig / isfinite / fabs
#include <exception>                                       // std::exception
#include <iostream>                                        // result printing
#include <random>                                          // uniform noise and fixed seed
#include <string>                                          // case names / failure reasons
#include <vector>                                          // sample container

using namespace sggk;   // SDK geometry primitives/surfaces
using namespace sgkq;   // this project's fit core

namespace
{

// pi at full double width: for angle conversion and numeric assertions.
const double kPi = 3.14159265358979323846;   // pi (radians)

int g_total = 0;    //!< cases judged
int g_failed = 0;   //!< FAIL cases; >0 makes main return 1

//! Print one test result and update the counters.
//! @param expected/actual share the same unit (mm or rad, per case).
//! @param pass computed by the caller from the case's criterion.
void report(const std::string& name, double expected, double actual, bool pass)   // result printer
{
    ++g_total;                 // count every case
    if (!pass) ++g_failed;     // count failures
    std::cout << "[TEST] " << name          // case name
              << " expected=" << expected   // expected
              << " actual=" << actual       // actual
              << (pass ? " PASS" : " FAIL") << std::endl;   // verdict
}

//! Angle between two directions (sign ignored), radians.
//! @return [0, pi/2]; sentinel 1e30 if either vector is degenerate (forces a FAIL).
double AngleRad(const Vector3D& a, const Vector3D& b)   // direction angle (rad)
{
    const double la = a.Length();   // |a|
    const double lb = b.Length();   // |b|
    if (la <= 0.0 || lb <= 0.0) return 1e30;   // degenerate: sentinel
    double c = std::fabs(a.Dot(b) / (la * lb));   // |cos(theta)|
    if (c > 1.0) c = 1.0;                          // clamp float overshoot
    return std::acos(c);                           // [0, pi/2] radians
}

//! Perpendicular distance rho from point q to line {pt, dir}.
//! @param dir normalized internally; sentinel 1e30 if it is a zero vector.
double DistToLine(const Point3D& q, const Point3D& pt, const Vector3D& dir)   // point-line distance (mm)
{
    Vector3D a = dir;                  // copy, normalize locally (leaves the argument alone)
    if (!a.Normalize()) return 1e30;   // degenerate direction: sentinel
    const Vector3D d(q.X() - pt.X(), q.Y() - pt.Y(), q.Z() - pt.Z());   // q - pt
    const double h = d.Dot(a);         // projection along the axis (mm)
    return (d - a * h).Length();       // perpendicular component length (mm)
}

//! Signed axial offset h of point q along line {pt, dir}.
//! @return signed axial projection (mm); sentinel 1e30 if dir is zero.
double OffsetAlong(const Point3D& q, const Point3D& pt, const Vector3D& dir)   // signed axial offset (mm)
{
    Vector3D a = dir;                  // copy, normalize
    if (!a.Normalize()) return 1e30;   // degenerate direction: sentinel
    const Vector3D d(q.X() - pt.X(), q.Y() - pt.Y(), q.Z() - pt.Z());   // q - pt
    return d.Dot(a);                   // signed axial component (mm)
}

// Geometry truth; cylinder/cone/torus axis and plane normal are the local z axis.
struct Truth
{
    Vector3D axis;        //!< axis direction / plane normal
    Point3D  refPt;       //!< a point on the axis / sphere center / cone apex / a point on the plane
    double   R = 0.0;     //!< cylinder / sphere / torus major radius (mm)
    double   r = 0.0;     //!< torus minor radius (mm)
    double   half = 0.0;  //!< cone half-angle (radians)
};

//! Test local frame: z = (0, 0.6, 0.8), deliberately not world-aligned (and z perpendicular to x).
//! @return Ucs3D with origin (12,-7,5), x axis +X, z axis (0,0.6,0.8).
Ucs3D TestUcs()   // shared local frame
{
    return Ucs3D(Point3D(12.0, -7.0, 5.0), Dir3D(1.0, 0.0, 0.0), Dir3D(0.0, 0.6, 0.8));   // origin + x + z
}

//! Build the expected Truth from a type and frame; sizes match the spec (and the surface objects).
//! @param type Unknown/Plane leaves the size fields at 0.
Truth TruthOf(Canonical type, const Ucs3D& ucs)   // build truth
{
    Truth tr;                       // truth container (all-zero default)
    const Dir3D z = ucs.UcsZ();     // local z = revolution axis / plane normal
    tr.axis = Vector3D(z.X(), z.Y(), z.Z());   // axis direction
    tr.refPt = ucs.Origin();        // point on axis = frame origin
    switch (type)                   // fill sizes by type
    {
    case Canonical::Cylinder: tr.R = 200.0; break;                     // radius 200 mm
    case Canonical::Sphere:   tr.R = 100.0; break;                     // radius 100 mm
    case Canonical::Cone:     tr.half = 60.0 * kPi / 180.0; break;     // half-angle 60 deg (radians)
    case Canonical::Torus:    tr.R = 100.0; tr.r = 10.0; break;        // major 100 / minor 10 mm
    default: break;                                                    // Plane/Unknown: no size
    }
    return tr;   // return truth
}

bool SaneParam(double x) { return std::isfinite(x) && std::fabs(x) < 1e20; }   // finite and reasonably scaled

//! U sampling interval: use the domain when bounded and reasonably scaled, else the fallback.
//! @note CalcPoint still returns an on-surface point outside the domain, so a fallback is safe.
void PickU(const Surface& srf, double fb0, double fb1, double& lo, double& hi)   // pick U interval
{
    const Interval d = srf.DomainU();                     // U domain
    const double a = d.Min(), b = d.Max();                // domain endpoints
    // Domain-size cap 100: beyond it (or infinite) use the fallback, so the span cannot blow up absolute-error criteria.
    if (SaneParam(a) && SaneParam(b) && (b - a) > 1e-9 && (b - a) <= 100.0) { lo = a; hi = b; }   // use actual domain
    else { lo = fb0; hi = fb1; }                          // fall back to given range
}

//! V sampling interval (same rule as PickU).
void PickV(const Surface& srf, double fb0, double fb1, double& lo, double& hi)   // pick V interval
{
    const Interval d = srf.DomainV();                     // V domain
    const double a = d.Min(), b = d.Max();                // domain endpoints
    if (SaneParam(a) && SaneParam(b) && (b - a) > 1e-9 && (b - a) <= 100.0) { lo = a; hi = b; }   // use actual domain
    else { lo = fb0; hi = fb1; }                          // fall back to given range
}

//! Sample an nu x nv grid (cell centers, avoiding domain edges) over [u0,u1]x[v0,v1] via CalcNormal;
//! degenerate points (isDegen) and near-zero normals are dropped.
//! @return points + unit normals; throwing/degenerate samples are silently skipped.
std::vector<FitSample> SampleGrid(const Surface& srf,     // grid sampler
                                  double u0, double u1, int nu,   // u interval and cell count
                                  double v0, double v1, int nv)   // v interval and cell count
{
    std::vector<FitSample> out;                       // output set
    out.reserve(static_cast<size_t>(nu) * nv);   // reserve upper bound

    for (int i = 0; i < nu; ++i)                      // over u
    {
        const double u = u0 + (u1 - u0) * (i + 0.5) / nu;   // cell center: offset by 0.5
        for (int j = 0; j < nv; ++j)                  // over v
        {
            const double v = v0 + (v1 - v0) * (j + 0.5) / nv;   // cell center

            Point3D pt;         // out: surface point
            bool degen = true;  // out: degenerate flag (init true, CalcNormal overwrites)
            Dir3D nrm;          // out: normal
            try                 // CalcNormal may throw on out-of-domain params
            {
                nrm = srf.CalcNormal(u, v, pt, degen);
            }
            catch (...)         // any exception
            {
                continue;   // out of range etc: skip this sample
            }
            if (degen) continue;   // skip poles/apex etc.

            const Vector3D n(nrm.X(), nrm.Y(), nrm.Z());   // normal to vector
            if (n.Length() < 0.5) continue;   // drop zero/near-zero normals (unit should be ~1)

            FitSample s;        // assemble sample
            s.p = pt;   // point (mm)
            s.n = n;    // normal (unit)
            out.push_back(s);   // append
        }
    }
    return out;   // result
}

// Item 1: five-type parameter recovery (analytic sampling + ToBSpline sampling).

//! Closed-form fit a sample set, then assert each parameter against the truth.
//! Per-type assertions: plane normal+point, cylinder axis/radius, sphere center/radius,
//! cone axis+apex+half-angle, torus axis+radii.
//! @param samples fewer than 20 samples fails sampleCount outright.
void CheckRecovery(const std::string& tag, Canonical type, const Truth& tr,   // single-path recovery check
                   const std::vector<FitSample>& samples)
{
    report(tag + ".sampleCount", 20.0, static_cast<double>(samples.size()), samples.size() >= 20);   // sample lower bound 20

    SurfaceParams params;   // fit output
    std::string reason;     // failure/diagnostic text
    const bool ok = FitSurface(type, samples, params, reason);   // closed-form fit entry

    report(tag + ".FitSurface", 1.0, ok ? 1.0 : 0.0, ok);   // fit must succeed
    if (!ok)                // fit failed
    {
        std::cout << "        reason: " << reason << std::endl;   // print failure reason
        return;             // end this case early
    }
    if (!reason.empty())    // non-fatal note
        std::cout << "        note: " << reason << std::endl;     // print it

    switch (type)           // assert parameters per type
    {
    case Canonical::Plane:
    {
        const double ang = AngleRad(params.axisDir, tr.axis);              // normal angle (rad)
        const double dp = distPlane(params.axisPoint, tr.refPt, tr.axis);  // plane-point distance (mm)
        report(tag + ".axisAngleRad", 0.0, ang, ang < 1e-8);               // normal angle < 1e-8 rad
        report(tag + ".planePointDist", 0.0, dp, dp < 1e-6);               // point distance < 1e-6 mm
        break;              // end plane branch
    }
    case Canonical::Cylinder:
    {
        const double ang = AngleRad(params.axisDir, tr.axis);              // axis angle (rad)
        const double rho = DistToLine(params.axisPoint, tr.refPt, tr.axis); // axis-point distance (mm)
        const double rel = std::fabs(params.radius - tr.R) / tr.R;         // relative radius error
        report(tag + ".axisAngleRad", 0.0, ang, ang < 1e-8);               // axis angle < 1e-8 rad
        report(tag + ".axisPointDist", 0.0, rho, rho < 1e-6);              // axis-point distance < 1e-6 mm
        report(tag + ".radiusRelErr", 0.0, rel, rel < 1e-8);               // relative radius error < 1e-8
        break;              // end cylinder branch
    }
    case Canonical::Sphere:
    {
        const double dc = Point3D::Distance(params.axisPoint, tr.refPt);   // center distance (mm)
        const double rel = std::fabs(params.radius - tr.R) / tr.R;         // relative radius error
        report(tag + ".centerDist", 0.0, dc, dc < 1e-6);                   // center distance < 1e-6 mm
        report(tag + ".radiusRelErr", 0.0, rel, rel < 1e-8);               // relative radius error < 1e-8
        break;              // end sphere branch
    }
    case Canonical::Cone:
    {
        const double ang = AngleRad(params.axisDir, tr.axis);              // axis angle (rad)
        const double da = Point3D::Distance(params.apex, tr.refPt);        // apex distance (mm)
        const double dh = std::fabs(params.halfAngle - tr.half);           // half-angle absolute error (rad)
        report(tag + ".axisAngleRad", 0.0, ang, ang < 1e-8);               // axis angle < 1e-8 rad
        report(tag + ".apexDist", 0.0, da, da < 1e-6);                     // apex distance < 1e-6 mm
        report(tag + ".halfAngleErr", 0.0, dh, dh < 1e-8);                 // half-angle error < 1e-8 rad
        break;              // end cone branch
    }
    case Canonical::Torus:
    {
        const double ang = AngleRad(params.axisDir, tr.axis);              // axis angle (rad)
        const double rho = DistToLine(params.axisPoint, tr.refPt, tr.axis); // axis-point distance (mm)
        const double hh = std::fabs(OffsetAlong(params.axisPoint, tr.refPt, tr.axis));   // axial offset (mm)
        const double relR = std::fabs(params.radius - tr.R) / tr.R;        // relative major-radius error
        const double relr = std::fabs(params.minorRadius - tr.r) / tr.r;   // relative minor-radius error
        report(tag + ".axisAngleRad", 0.0, ang, ang < 1e-8);               // axis angle < 1e-8 rad
        report(tag + ".axisPointRho", 0.0, rho, rho < 1e-6);               // distance < 1e-6 mm
        report(tag + ".centerPlaneOffset", 0.0, hh, hh < 1e-6);            // axial offset < 1e-6 mm (ring center on axis)
        report(tag + ".majorRadiusRelErr", 0.0, relR, relR < 1e-8);        // relative major-radius error < 1e-8
        report(tag + ".minorRadiusRelErr", 0.0, relr, relr < 1e-8);        // relative minor-radius error < 1e-8
        break;              // end torus branch
    }
    default:
        report(tag + ".unknownType", 1.0, 0.0, false);   // should never happen
        break;              // end default branch
    }
}

//! Two sampling paths for one surface type: the analytic domain + the ToBSpline own domain.
//! fb* are fallbacks when the domain is unbounded or unreasonably scaled.
void RunRecovery(const Surface& srf, Canonical type, const std::string& name,   // two-path recovery
                 const Ucs3D& ucs, int nu, int nv,
                 double fbU0, double fbU1, double fbV0, double fbV1)
{
    const Truth tr = TruthOf(type, ucs);   // truth

    // (a) analytic surface: u/v from the domain, fallback otherwise
    double u0 = fbU0, u1 = fbU1, v0 = fbV0, v1 = fbV1;   // init to fallback
    PickU(srf, fbU0, fbU1, u0, u1);      // try analytic U domain
    PickV(srf, fbV0, fbV1, v0, v1);      // try analytic V domain
    const std::vector<FitSample> analytic = SampleGrid(srf, u0, u1, nu, v0, v1, nv);   // analytic sampling
    CheckRecovery(name + "_analytic", type, tr, analytic);   // path (a)

    // (b) ToBSpline(analytic domain) -> sample its own domain (U is reparameterized)
    SurfacePtr spline;                   // spline surface
    try                                  // ToBSpline may throw
    {
        spline = srf.ToBSpline(UVRange(Interval(u0, u1), Interval(v0, v1)));   // pass analytic domain
    }
    catch (...)                          // any exception
    {
        spline = nullptr;   // failed conversion -> null
    }
    if (!spline)                         // conversion failed
    {
        report(name + "_bspline.ToBSpline", 1.0, 0.0, false);   // must convert
        return;                          // end this path early
    }

    double su0 = 0.0, su1 = 1.0, sv0 = v0, sv1 = v1;   // spline defaults [0,1] (U reparameterized)
    PickU(*spline, 0.0, 1.0, su0, su1);  // try spline U domain
    PickV(*spline, v0, v1, sv0, sv1);    // try spline V domain
    const std::vector<FitSample> splineSamples = SampleGrid(*spline, su0, su1, nu, sv0, sv1, nv);   // spline sampling
    CheckRecovery(name + "_bspline", type, tr, splineSamples);   // path (b)
}

// Item 2: clipped arc -- only a 60-degree arc of the cylinder.
void testClippedArc(const Surface& cyl)   // clipped-arc test
{
    const Ucs3D ucs = TestUcs();                          // local frame
    const Truth tr = TruthOf(Canonical::Cylinder, ucs);   // truth: axis/origin/radius

    double v0 = 0.0, v1 = 10.0;      // v init (mm)
    PickV(cyl, 0.0, 10.0, v0, v1);   // v interval = axial height (mm)

    const std::vector<FitSample> arc = SampleGrid(cyl, 0.0, kPi / 3.0, 24, v0, v1, 5);   // 60-deg arc, 24x5
    report("arc.sampleCount", 60.0, static_cast<double>(arc.size()), arc.size() >= 60);   // expect >= 60 samples

    SurfaceParams params;   // fit output
    std::string reason;     // failure reason
    const bool ok = FitSurface(Canonical::Cylinder, arc, params, reason);   // fit the 1/6 arc only
    report("arc.FitSurface", 1.0, ok ? 1.0 : 0.0, ok);   // must succeed
    if (!ok)                // failure
    {
        std::cout << "        reason: " << reason << std::endl;   // print reason
        return;             // end early
    }

    const double ang = AngleRad(params.axisDir, tr.axis);               // axis angle (rad)
    const double rho = DistToLine(params.axisPoint, tr.refPt, tr.axis);  // axis-point distance (mm)
    const double rel = std::fabs(params.radius - tr.R) / tr.R;          // relative radius error
    report("arc.axisAngleRad", 0.0, ang, ang < 1e-8);                   // axis angle < 1e-8 rad
    report("arc.axisPointDist", 0.0, rho, rho < 1e-6);                  // axis-point distance < 1e-6 mm
    report("arc.radiusRelErr", 0.0, rel, rel < 1e-8);                   // relative radius error < 1e-8

    // Regression guard: the axis must not be the mean of projected points (an arc mean is not the center).
    // Correct axis-to-centroid distance is about r*|mean(radial unit vectors)| ~= 0.955*r; taking the mean as the center gives 0.
    Vector3D centroid(0.0, 0.0, 0.0);   // arc centroid (accumulator)
    for (const FitSample& s : arc)                   // over samples
        centroid += Vector3D(s.p.X(), s.p.Y(), s.p.Z());   // accumulate position vectors
    centroid /= static_cast<double>(arc.size());   // average
    const Point3D cpt(centroid.X(), centroid.Y(), centroid.Z());   // centroid point
    const double cDist = DistToLine(cpt, params.axisPoint, params.axisDir);   // centroid-to-axis distance (mm)
    report("arc.axisNotAtArcCentroid", 0.9 * tr.R, cDist, cDist > 0.5 * tr.R);   // must be clearly >0 (not the centroid)
}

// Item 3: refinement (1e-4 noise: coarse fit on clipped arc -> full-turn refinement).

//! Cylinder sampling with 1e-4 uniform noise on both points and normals; the fixed seed makes it reproducible.
//! @param rng caller-supplied engine so the sequence is reproducible.
//! @return noisy samples; degenerate/throwing points are dropped.
std::vector<FitSample> NoisyCylinderSamples(const Surface& srf,   // noisy cylinder sampling
                                            double u0, double u1, int nu,   // u interval and cell count
                                            double v0, double v1, int nv,   // v interval and cell count
                                            std::mt19937& rng)              // random engine
{
    std::uniform_real_distribution<double> noise(-1e-4, 1e-4);   // +/-1e-4 uniform noise (mm for points, unitless for normals)
    std::vector<FitSample> out;   // output set

    for (int i = 0; i < nu; ++i)                  // over u
    {
        const double u = u0 + (u1 - u0) * (i + 0.5) / nu;   // u cell center
        for (int j = 0; j < nv; ++j)              // over v
        {
            const double v = v0 + (v1 - v0) * (j + 0.5) / nv;   // v cell center

            Point3D pt;         // surface point
            bool degen = true;  // degenerate flag
            Dir3D nrm;          // normal
            try                 // CalcNormal may throw
            {
                nrm = srf.CalcNormal(u, v, pt, degen);
            }
            catch (...)         // any exception
            {
                continue;   // skip out of range
            }
            if (degen || nrm.Length() < 0.5) continue;   // skip degenerate/zero-normal

            FitSample s;    // assemble noisy sample
            s.p = Point3D(pt.X() + noise(rng), pt.Y() + noise(rng), pt.Z() + noise(rng));   // noisy point (mm)
            Vector3D n(nrm.X() + noise(rng), nrm.Y() + noise(rng), nrm.Z() + noise(rng));   // noisy normal
            n.Normalize();   // renormalize
            s.n = n;    // store normal
            out.push_back(s);   // append
        }
    }
    return out;   // result
}

//! Scalar cylinder parameter error: |dr| + axis-point-to-true-axis distance + r * (axis angle).
double ParamError(const SurfaceParams& p, const Truth& tr)   // cylinder parameter error (mm)
{
    return std::fabs(p.radius - tr.R)                      // absolute radius error (mm)
         + DistToLine(p.axisPoint, tr.refPt, tr.axis)      // axis-point distance (mm)
         + tr.R * AngleRad(p.axisDir, tr.axis);            // angle (rad) * R -> arc length (mm)
}

//! Closed-form sum of squared orthogonal distances (via public distCylinder, same units as refinement).
double CostCylinder(const SurfaceParams& p, const std::vector<FitSample>& samples)   // residual sum of squares (mm^2)
{
    double c = 0.0;   // residual sum of squares (mm^2)
    for (const FitSample& s : samples)               // over samples
    {
        const double d = distCylinder(s.p, p.axisPoint, p.axisDir, p.radius);   // per-point orthogonal distance (mm)
        c += d * d;   // accumulate square
    }
    return c;   // cost
}

void testRefine(const Surface& cyl, const Ucs3D& ucs)   // refinement test
{
    const Truth tr = TruthOf(Canonical::Cylinder, ucs);   // truth

    double v0 = 0.0, v1 = 10.0;      // v init (mm)
    PickV(cyl, 0.0, 10.0, v0, v1);   // axial interval (mm)

    std::mt19937 rng(20261009u);   // fixed seed: reproducible noise sequence

    // full-turn noisy samples (for refinement): all u, 24x6
    const std::vector<FitSample> all =
        NoisyCylinderSamples(cyl, 0.0, 2.0 * kPi, 24, v0, v1, 6, rng);   // 144 samples

    // clipped-arc noisy samples (coarse fit): only a 40-deg arc, 4x3 -> params clearly off the truth
    const std::vector<FitSample> few =
        NoisyCylinderSamples(cyl, 0.0, 40.0 * kPi / 180.0, 4, v0, v1, 3, rng);   // 12 samples

    report("refine.sampleCount", 120.0, static_cast<double>(all.size()), all.size() >= 120);   // full-turn lower bound

    SurfaceParams coarse;   // coarse fit output
    std::string reason;     // failure reason
    const bool okCoarse = FitSurface(Canonical::Cylinder, few, coarse, reason);   // 12 samples only
    report("refine.coarseFit", 1.0, okCoarse ? 1.0 : 0.0, okCoarse);   // coarse fit must succeed
    if (!okCoarse)          // failure
    {
        std::cout << "        reason: " << reason << std::endl;   // print reason
        return;             // end early
    }

    const double costCoarse = CostCylinder(coarse, all);   // coarse-solution cost over full turn (mm^2)
    const double errCoarse = ParamError(coarse, tr);       // coarse parameter error (mm)

    SurfaceParams refined = coarse;   // refine starting from the coarse solution
    const bool okRefine = RefineSurface(refined, all, 100, 1e-13, reason);   // 100 iterations, tolerance 1e-13
    report("refine.RefineSurface", 1.0, okRefine ? 1.0 : 0.0, okRefine);     // refinement must succeed
    if (!okRefine)          // failure
        std::cout << "        reason: " << reason << std::endl;   // print reason

    const double costRefined = CostCylinder(refined, all);   // refined cost (mm^2)
    const double errRefined = ParamError(refined, tr);       // refined parameter error (mm)

    report("refine.costDecrease", costCoarse, costRefined, costRefined < costCoarse);            // cost must drop
    report("refine.paramErrorDecrease", errCoarse, errRefined, errRefined < errCoarse);          // parameter error must drop
    report("refine.paramErrorHalved", 0.5 * errCoarse, errRefined, errRefined < 0.5 * errCoarse); // and at least halve
}

// Item 4: closed-form distance (MeasureAccuracy).
void testClosedDistance(const Surface& plane, const Surface& cyl, const Surface& sph,   // closed-form distance test
                        const Surface& cone, const Surface& tor)
{
    const Thresholds th;   // maxDist = 1e-5 mm (spec criterion)

    struct Case
    {
        Canonical   type;       //!< surface type
        const char* name;       //!< case-name prefix
        const Surface* srf;     //!< surface pointer (owned by the caller)
        int    nu, nv;          //!< sampling grid
        double fbU0, fbU1, fbV0, fbV1;   //!< fallback parameter ranges
    };
    const Case cases[] = {
        { Canonical::Plane,    "acc_plane",    &plane, 7,  7,  -10.0, 10.0, -10.0, 10.0 },        //!< plane square window 20x20 mm
        { Canonical::Cylinder, "acc_cylinder", &cyl,   24, 5,  0.0,  2.0 * kPi, 0.0, 10.0 },     //!< cylinder full turn, height 10 mm
        { Canonical::Sphere,   "acc_sphere",   &sph,   24, 8,  0.0,  2.0 * kPi, -1.4, 1.4 },     //!< sphere band (avoid poles, |v|<=1.4 rad)
        { Canonical::Cone,     "acc_cone",     &cone,  24, 6,  0.0,  2.0 * kPi, 0.05, 0.55 },    //!< cone band (v avoids apex)
        { Canonical::Torus,    "acc_torus",    &tor,   24, 12, 0.0,  2.0 * kPi, 0.0, 2.0 * kPi },//!< torus full turn
    };

    for (const Case& c : cases)     // per type
    {
        double u0 = c.fbU0, u1 = c.fbU1, v0 = c.fbV0, v1 = c.fbV1;   // init to fallback
        PickU(*c.srf, c.fbU0, c.fbU1, u0, u1);   // actual U domain
        PickV(*c.srf, c.fbV0, c.fbV1, v0, v1);   // actual V domain

        const std::vector<FitSample> samples = SampleGrid(*c.srf, u0, u1, c.nu, v0, v1, c.nv);   // sample
        report(std::string(c.name) + ".sampleCount", 30.0,
               static_cast<double>(samples.size()), samples.size() >= 30);   // lower bound 30

        SurfaceParams params;   // fit output
        std::string reason;     // failure reason
        const bool ok = FitSurface(c.type, samples, params, reason);   // closed-form fit
        if (!ok)                // failure
        {
            report(std::string(c.name) + ".FitSurface", 1.0, 0.0, false);   // record FAIL
            std::cout << "        reason: " << reason << std::endl;         // print reason
            continue;           // next type
        }

        // Closed-form initial guess: the 6x6 Pluecker eigen-decomposition (torus axis) has a ~1e-10 rad
        // direction floor (see test_math's 1e-8 deg assertion for solvePluckerAxis), ~1e-8 mm in distance;
        // so only loosely verify "the guess lies on the surface"; the strict 1e-12 check comes after refinement.
        const AccuracyReport closed = MeasureAccuracy(params, samples, th);   // initial-guess accuracy
        report(std::string(c.name) + ".closedForm", 0.0, closed.maxDist, closed.maxDist < 1e-5);   // initial: maxDist < 1e-5 mm

        // refine to machine precision, then re-check (spec item 4: maxDist < 1e-12)
        SurfaceParams polished = params;   // copy to refine
        const bool okRefine = RefineSurface(polished, samples, 100, 1e-13, reason);   // 100 iterations, tolerance 1e-13
        report(std::string(c.name) + ".RefineSurface", 1.0, okRefine ? 1.0 : 0.0, okRefine);   // refinement must succeed
        if (!okRefine)          // failure
        {
            std::cout << "        reason: " << reason << std::endl;   // print reason
            continue;           // next type
        }

        const AccuracyReport rep = MeasureAccuracy(polished, samples, th);   // post-refinement accuracy
        report(std::string(c.name) + ".maxDist", 0.0, rep.maxDist, rep.maxDist < 1e-12);   // max distance < 1e-12 mm
        report(std::string(c.name) + ".avgDist", 0.0, rep.avgDist, rep.avgDist < 1e-12);   // average distance < 1e-12 mm
        report(std::string(c.name) + ".accepted", 1.0, rep.accepted ? 1.0 : 0.0, rep.accepted);   // must be accepted
    }
}

// Item 5: degenerate inputs must be rejected explicitly (no crash).
void ExpectReject(const std::string& name, Canonical type, const std::vector<FitSample>& samples)   // expect rejection
{
    SurfaceParams params;   // fit output (ignored)
    std::string reason;     // failure reason
    bool ok = true;         // whether the fit succeeded
    try                     // degenerate input should not throw
    {
        ok = FitSurface(type, samples, params, reason);   // expect false, not throw
    }
    catch (const std::exception& e)   // standard exception
    {
        report(name, 0.0, 2.0, false);       // throwing violates "must not crash"
        std::cout << "        exception: " << e.what() << std::endl;   // print exception info
        return;             // done
    }

    report(name, 0.0, ok ? 1.0 : 0.0, !ok && !reason.empty());   // criterion: rejected with a non-empty reason
    if (!ok && !reason.empty())   // rejected with a reason
        std::cout << "        reason: " << reason << std::endl;   // print reason
}

void testDegenerate()   // degenerate-input test
{
    // (a) near-line cone alpha = 1e-9: points/normals generated analytically (no reliance on the SDK building such a surface)
    {
        const double alpha = 1e-9;                    // half-angle -> 0 (radians)
        const double sa = std::sin(alpha), ca = std::cos(alpha);   // sin(alpha)~1e-9, cos(alpha)~1
        std::vector<FitSample> samples;               // sample set
        for (int i = 0; i < 40; ++i)                  // 40 samples
        {
            const double th = 2.0 * kPi * i / 40.0;   // uniform angle
            const double t = 10.0 + 0.5 * i;          // distance along the generatrix, t in [10,29.5] mm
            FitSample s;                              // assemble sample
            s.p = Point3D(t * sa * std::cos(th), t * sa * std::sin(th), t * ca);   // point near the z axis
            s.n = Vector3D(ca * std::cos(th), ca * std::sin(th), -sa);             // near-radial normal
            samples.push_back(s);                     // append
        }
        ExpectReject("degenerate.cone_alpha_to_0", Canonical::Cone, samples);   // must be rejected
    }

    // (b) r -> 0 cylinder: all sample points lie on one line (degenerate projection -> circle fit fails)
    {
        std::vector<FitSample> samples;   // sample set
        for (int i = 0; i < 40; ++i)      // 40 samples
        {
            const double th = 2.0 * kPi * i / 40.0;   // spreads the normal
            FitSample s;                              // assemble sample
            s.p = Point3D(0.3, -0.7, -5.0 + 0.25 * i);          // a straight line (radius -> 0)
            s.n = Vector3D(std::cos(th), std::sin(th), 0.0);     // perpendicular to axis, spread
            samples.push_back(s);     // append
        }
        ExpectReject("degenerate.cylinder_r_to_0", Canonical::Cylinder, samples);   // must be rejected
    }

    // (c) too few samples
    {
        std::vector<FitSample> samples;   // sample set
        for (int i = 0; i < 2; ++i)       // only 2
        {
            FitSample s;                  // assemble sample
            s.p = Point3D(1.0 * i, 2.0 * i, 3.0 * i);   // two collinear points
            s.n = Vector3D(0.0, 0.0, 1.0);              // normal
            samples.push_back(s);   // append
        }
        ExpectReject("degenerate.tooFewSamples", Canonical::Cylinder, samples);   // too few samples must be rejected
    }

    // (d) unknown type
    {
        std::vector<FitSample> samples;   // sample set
        for (int i = 0; i < 10; ++i)      // 10 samples
        {
            FitSample s;                  // assemble sample
            s.p = Point3D(1.0 * i, 2.0, 3.0);   // 10 points on one line
            s.n = Vector3D(0.0, 0.0, 1.0);      // normal
            samples.push_back(s);   // append
        }
        ExpectReject("degenerate.unknownType", Canonical::Unknown, samples);   // unknown type must be rejected
    }

    // (e) illegal refinement arguments (maxIter / tol)
    {
        SurfaceParams p;                          // valid initial guess
        p.type = Canonical::Cylinder;             // cylinder
        p.axisDir = Vector3D(0.0, 0.0, 1.0);      // axis +z
        p.axisPoint = Point3D(0.0, 0.0, 0.0);     // point on axis
        p.radius = 10.0;                          // radius 10 mm

        std::vector<FitSample> samples;           // valid samples
        for (int i = 0; i < 40; ++i)              // 40 samples
        {
            const double th = 2.0 * kPi * i / 40.0;   // uniform angle
            FitSample s;                              // assemble sample
            s.p = Point3D(10.0 * std::cos(th), 10.0 * std::sin(th), 0.5 * i);   // cylinder-surface points
            s.n = Vector3D(std::cos(th), std::sin(th), 0.0);                    // outward unit normal
            samples.push_back(s);   // append
        }

        std::string reason;   // rejection reason
        const bool okIter = RefineSurface(p, samples, 0, 1e-12, reason);   // maxIter=0 illegal
        report("degenerate.refine_maxIter0", 0.0, okIter ? 1.0 : 0.0, !okIter && !reason.empty());   // must reject with reason

        const bool okTol = RefineSurface(p, samples, 10, 0.0, reason);   // tol=0 illegal
        report("degenerate.refine_tol0", 0.0, okTol ? 1.0 : 0.0, !okTol && !reason.empty());   // must reject with reason
    }
}

} // namespace

// Entry: five-type recovery (two paths) -> clipped arc -> refinement -> closed-form distance -> degenerate rejection.
int main()   // test entry
{
    try                                     // catch everything
    {
        init();   // initialize SDK global state

        const Ucs3D ucs = TestUcs();   // shared local frame (z=(0,0.6,0.8))

        // analytic surfaces (sizes from the spec: cylinder r=200, sphere r=100, cone half-angle 60 deg, torus R=100/r=10)
        const Plane         plane(ucs);                                          // plane: normal = ucs.z
        const CylinderSurface cyl(ucs, 200.0);                                   // cylinder r=200 mm
        const SphereSurface sphere(ucs, 100.0);                                  // sphere r=100 mm
        const ConeSurface   cone = ConeSurface::CreateFromHalfAngle(ucs, 60.0 * kPi / 180.0);   // cone half-angle 60 deg
        const TorusSurface  torus(ucs, 100.0, 10.0);                             // torus R=100, r=10 mm

        // 1) five-type recovery: analytic + ToBSpline, two paths
        RunRecovery(plane,  Canonical::Plane,    "fit_plane",    ucs, 7,  7,  -10.0, 10.0, -10.0, 10.0);   // plane
        RunRecovery(cyl,    Canonical::Cylinder, "fit_cylinder", ucs, 24, 5,  0.0, 2.0 * kPi, 0.0, 10.0);   // cylinder
        RunRecovery(sphere, Canonical::Sphere,   "fit_sphere",   ucs, 24, 8,  0.0, 2.0 * kPi, -1.4, 1.4);   // sphere
        RunRecovery(cone,   Canonical::Cone,     "fit_cone",     ucs, 24, 6,  0.0, 2.0 * kPi, 0.05, 0.55);   // cone
        RunRecovery(torus,  Canonical::Torus,    "fit_torus",    ucs, 24, 12, 0.0, 2.0 * kPi, 0.0, 2.0 * kPi);   // torus

        // 2) clipped arc / 3) refinement / 4) closed-form distance / 5) degenerate
        testClippedArc(cyl);                          // 2) 60-deg clipped arc
        testRefine(cyl, ucs);                         // 3) noisy refinement
        testClosedDistance(plane, cyl, sphere, cone, torus);   // 4) closed-form distance
        testDegenerate();                             // 5) degenerate rejection

        fini();   // release SDK global state
    }
    catch (const std::exception& e)   // standard exception
    {
        std::cerr << "[FAIL] exception: " << e.what() << std::endl;   // print exception
        return 2;   // environment/interface failure
    }
    catch (...)   // unknown exception
    {
        std::cerr << "[FAIL] unknown exception" << std::endl;   // print message
        return 2;   // same as 2
    }

    std::cout << "[SUMMARY] total=" << g_total << " failed=" << g_failed << std::endl;   // summary
    return (g_failed == 0) ? 0 : 1;   // exit code decided only by FAIL count
}
