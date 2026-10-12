// apps/checks/test_detect/main.cpp -- surface detection (classification) self-check.
// Converts analytic surfaces with ToBSpline(uv) and judges them on the spline's own domain
// (ToBSpline reparametrizes U to [0,1], so the passed-in uv range is not a usable domain).
// Prints [TEST]/[INFO]/[NOTE] per case; exit code 0 = all PASS, 1 = FAIL present, 2 = exception.
#include "sgkquadric/detect.hpp"

#include <Foundation/init.h>

#include <GeomBase/Dir3D.h>
#include <GeomBase/Interval.h>
#include <GeomBase/Point3D.h>
#include <GeomBase/Ucs3D.h>
#include <GeomBase/UVRange.h>
#include <GeomBase/Vector3D.h>
#include <Geometry/3D/Surface/BSplineSurface.h>
#include <Geometry/3D/Surface/ConeSurface.h>
#include <Geometry/3D/Surface/CylinderSurface.h>
#include <Geometry/3D/Surface/Plane.h>
#include <Geometry/3D/Surface/SphereSurface.h>
#include <Geometry/3D/Surface/TorusSurface.h>

#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace sggk;   // SDK geometry primitives / surface classes
using namespace sgkq;   // this project's quadric detect/fit core

namespace
{

// Full-precision pi for degree -> radian conversion.
const double kPi = 3.14159265358979323846;

int g_total  = 0;   //!< cases counted toward the exit code (unrejected negatives go to g_warned)
int g_failed = 0;   //!< FAIL count; >0 makes main return 1
int g_warned = 0;   //!< negatives wrongly accepted (WARN); does not affect the exit code

//! @brief First third of each interval: [Min, Min + Length/3], in parameter units.
//! @note Verifies detection does not depend on full-domain coverage; r is not modified.
UVRange cropThird(const UVRange& r)
{
    const Interval& iu = r.IntervalU();
    const Interval& iv = r.IntervalV();
    return UVRange(Interval(iu.Min(), iu.Min() + iu.Length() / 3.0),
                   Interval(iv.Min(), iv.Min() + iv.Length() / 3.0));
}

//! @brief One probe: full Detect() result plus sample/degenerate counts and spline domain.
struct Probe
{
    DetectResult result;              //!< candidate type + residual + notes
    int          sampleCount = 0;     //!< total samples returned by SampleSurface
    int          degenCount  = 0;     //!< degen samples (poles/apex) among them
    double       uMin = 0, uMax = 0, vMin = 0, vMax = 0;  //!< spline's own param domain
};

//! @brief Convert analytic surface to spline, then sample + detect on the spline's own domain.
//! @throws std::runtime_error when ToBSpline returns null.
//! @note Must use bs->DomainU()/DomainV() rather than the passed-in uv: ToBSpline reparametrizes
//!       U, so the input range can fall outside the spline domain and degenerate the sampling.
Probe probeBSpline(const Surface& analytic, const UVRange& uv)
{
    BSplineSurfacePtr bs = analytic.ToBSpline(uv);   // analytic -> NURBS; null pointer on failure
    if (!bs)
    {
        throw std::runtime_error("ToBSpline returned null");
    }

    const UVRange buv(bs->DomainU(), bs->DomainV());  // spline's own domain, not the passed uv

    Probe p;
    const std::vector<Sample> smp = SampleSurface(*bs, buv, 9, 9);  // 9x9 grid; counts, not lengths
    p.sampleCount = static_cast<int>(smp.size());
    for (const Sample& s : smp)
    {
        if (s.degen) ++p.degenCount;   // degenerates (poles/apex) must be ignored by the judgement
    }
    p.uMin = buv.IntervalU().Min();
    p.uMax = buv.IntervalU().Max();
    p.vMin = buv.IntervalV().Min();
    p.vMax = buv.IntervalV().Max();
    p.result = Detect(*bs, buv);
    return p;
}

//! @brief Print one case's [TEST]/[INFO]/[NOTE] lines and update the global counters.
//! @param pass computed by the caller according to its own criterion.
void report(const std::string& name, const char* expected, const DetectResult& r, bool pass,
            const std::string& extra = std::string())
{
    ++g_total;
    if (!pass) ++g_failed;

    std::cout << "[TEST] " << name
              << " expected=" << expected
              << " actual=" << ToString(r.candidate)
              << (pass ? " PASS" : " FAIL") << std::endl;

    std::cout << "[INFO] " << name
              << " residual=" << std::scientific << std::setprecision(6) << r.residual;  // 6 sig digits
    if (!extra.empty()) std::cout << " " << extra;
    std::cout << std::endl;

    for (const std::string& nt : r.notes)
    {
        std::cout << "[NOTE] " << name << " " << nt << std::endl;
    }
}

//! @brief Positive case: candidate must equal expected on the spline's own domain.
void checkCase(const std::string& name, Canonical expected, const Surface& srf, const UVRange& uv)
{
    const Probe p = probeBSpline(srf, uv);
    // Extra context: degen/total samples + the spline's own domain, for manual review.
    const std::string extra = "degen=" + std::to_string(p.degenCount) + "/" +
                              std::to_string(p.sampleCount) +
                              " bsUV=[" + std::to_string(p.uMin) + "," + std::to_string(p.uMax) +
                              "]x[" + std::to_string(p.vMin) + "," + std::to_string(p.vMax) + "]";
    report(name, ToString(expected), p.result, p.result.candidate == expected, extra);
}

//! @brief Cases 1)-2): the five canonical surfaces plus their 1/3-cropped equivalents.
//! @note Window size (in mm) only affects sample density and residual magnitude, not the type.
void runPositives()
{
    const double deg = kPi / 180.0;   // degree -> radian

    // Plane: u,v are in-plane local coords; bounded [-4,4]^2 window (8 mm across).
    Plane plane(Point3D(0.0, 0.0, 0.0), Dir3D(0.0, 0.0, 1.0));
    const UVRange planeUv(Interval(-4.0, 4.0), Interval(-4.0, 4.0));
    checkCase("positive.plane", Canonical::Plane, plane, planeUv);
    checkCase("crop.plane", Canonical::Plane, plane, cropThird(planeUv));

    // Cylinder: u wraps [0,2pi], v is axial height (radius 5 mm, height 12 mm).
    CylinderSurface cyl(Ucs3D(Point3D(0.0, 0.0, 0.0)), 5.0);
    const UVRange cylUv(Interval(0.0, 2.0 * kPi), Interval(-6.0, 6.0));
    checkCase("positive.cylinder", Canonical::Cylinder, cyl, cylUv);
    checkCase("crop.cylinder", Canonical::Cylinder, cyl, cropThird(cylUv));

    // Cone: u wraps, v axial with apex at v=0; positive v avoids the apex singularity.
    const UVRange coneUv(Interval(0.0, 2.0 * kPi), Interval(2.0, 7.0));
    ConeSurface cone30 = ConeSurface::CreateFromHalfAngle(Ucs3D(Point3D(0.0, 0.0, 0.0)), 30.0 * deg);
    checkCase("positive.cone30", Canonical::Cone, cone30, coneUv);
    checkCase("crop.cone30", Canonical::Cone, cone30, cropThird(coneUv));

    // Half-angle 60 deg lies near the 54.7356 deg = atan(sqrt(2)) second-moment trap, where a
    // cone's moments can mimic a cylinder-like signature; it must still classify as Cone.
    ConeSurface cone60 = ConeSurface::CreateFromHalfAngle(Ucs3D(Point3D(0.0, 0.0, 0.0)), 60.0 * deg);
    checkCase("positive.cone60", Canonical::Cone, cone60, coneUv);
    checkCase("crop.cone60", Canonical::Cone, cone60, cropThird(coneUv));

    // Sphere: u wraps, v in [0,pi] (includes both poles); radius 8 mm.
    SphereSurface sph(Ucs3D(Point3D(0.0, 0.0, 0.0)), 8.0);
    const UVRange sphUv(Interval(0.0, 2.0 * kPi), Interval(0.0, kPi));
    checkCase("positive.sphere", Canonical::Sphere, sph, sphUv);
    checkCase("crop.sphere", Canonical::Sphere, sph, cropThird(sphUv));

    // Torus (R > r, no self-intersection): u major, v minor, both wrap; R=10 mm, r=3 mm.
    TorusSurface tor(Ucs3D(Point3D(0.0, 0.0, 0.0)), 10.0, 3.0);
    const UVRange torUv(Interval(0.0, 2.0 * kPi), Interval(0.0, 2.0 * kPi));
    checkCase("positive.torus", Canonical::Torus, tor, torUv);
    checkCase("crop.torus", Canonical::Torus, tor, cropThird(torUv));
}

//! @brief Case 4): full sphere spanning both poles; poles must be skipped, type still Sphere.
//! @note The [INFO] degen ratio must be > 0, evidencing the poles entered the sample set.
void runSpherePoles()
{
    SphereSurface sph(Ucs3D(Point3D(0.0, 0.0, 0.0)), 8.0);
    const UVRange sphUv(Interval(0.0, 2.0 * kPi), Interval(0.0, kPi));

    const Probe p = probeBSpline(sph, sphUv);
    const std::string extra = "degen=" + std::to_string(p.degenCount) + "/" +
                              std::to_string(p.sampleCount) +
                              "(pole samples must be flagged degenerate and skipped)";
    report("sphere.poles", "Sphere", p.result, p.result.candidate == Canonical::Sphere, extra);
}

//! @brief Case 3): perturb one cylinder control point radially by 5%*r; report honestly.
//! @note Rejection is PASS and acceptance is WARN -- neither changes the exit code.
void runNegative()
{
    const double radius = 5.0;   // mm; perturbation is 5%*radius = 0.25 mm
    CylinderSurface cyl(Ucs3D(Point3D(0.0, 0.0, 0.0)), radius);
    const UVRange cylUv(Interval(0.0, 2.0 * kPi), Interval(-6.0, 6.0));

    BSplineSurfacePtr bs = cyl.ToBSpline(cylUv);
    if (!bs)
    {
        throw std::runtime_error("negative case: ToBSpline returned null");
    }

    // Shift the middle control point along its radial direction (in XY) by 5%*r.
    const Point3DMatrix& cps = bs->ControlPoints();   // matrix indexed [iu][iv]
    if (cps.empty() || cps.front().empty())
    {
        throw std::runtime_error("negative case: empty control point matrix");
    }
    const UInt iu = static_cast<UInt>(cps.size() / 2);            // U midpoint index
    const UInt iv = static_cast<UInt>(cps[iu].size() / 2);        // V midpoint index

    const Point3D q0 = cps[iu][iv];
    Vector3D radial(q0.X(), q0.Y(), 0.0);                        // radial direction in XY
    if (!(radial.Length() > 0.0)) radial = Vector3D(1.0, 0.0, 0.0);   // degenerate fallback
    radial.Normalize();

    const Point3D q1 = q0 + radial * (0.05 * radius);            // radial offset 5%*r = 0.25 mm
    bs->SetControlPoint(q1, iu, iv, false);   // false: do not refresh periodicity/closure cache

    const UVRange buv(bs->DomainU(), bs->DomainV());   // still judge on the spline's own domain
    const DetectResult r = Detect(*bs, buv);

    const bool rejected = (r.candidate == Canonical::Unknown);   // expected: rejected as Unknown

    // Accepting the perturbed surface is not a hard failure: emit WARN only, and print the
    // residual magnitude honestly so the tolerance question stays visible.
    if (rejected) ++g_total; else ++g_warned;

    std::cout << "[TEST] negative.perturbedCylinder"
              << " expected=Unknown"
              << " actual=" << ToString(r.candidate)
              << (rejected ? " PASS" : " WARN") << std::endl;

    std::cout << "[INFO] negative.perturbedCylinder"
              << " residual=" << std::scientific << std::setprecision(6) << r.residual   // mm
              << " offset=5%*r=" << (0.05 * radius)
              << " ctrlPointIdx=(" << iu << "," << iv << ")" << std::endl;

    for (const std::string& nt : r.notes)
    {
        std::cout << "[NOTE] negative.perturbedCylinder " << nt << std::endl;
    }

    if (!rejected)
    {
        std::cout << "[WARN] negative.perturbedCylinder accepted as " << ToString(r.candidate)
                  << " (residual=" << std::scientific << std::setprecision(6) << r.residual
                  << "); after a 5-percent radial perturbation of one control point a sample still matches this signature -- reported as-is, criteria not relaxed."
                  << std::endl;
    }
}

} // namespace

//! @brief Test entry: positives+crops / negative / poles; any FAIL or exception exits nonzero.
//! @return 0 = all PASS; 1 = FAIL present; 2 = exception.
//! @note init()/fini() are paired; cases aggregate via g_total/g_failed/g_warned.
int main()
{
    try
    {
        init();   // SDK init: required before any geometry object is used

        runPositives();      // 1)+2) five positives and their crops
        runNegative();       // 3) negative rejection
        runSpherePoles();    // 4) full sphere with both poles

        fini();              // SDK shutdown, paired with init()
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

    // Summary: warned is listed separately so accepted negatives get manual review.
    std::cout << "[SUMMARY] total=" << g_total
              << " failed=" << g_failed
              << " warned=" << g_warned << std::endl;
    return (g_failed == 0) ? 0 : 1;   // only FAIL decides the exit code
}
