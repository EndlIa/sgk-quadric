// apps/checks/test_pipeline/main.cpp -- end-to-end self-check: ProcessFace (decide -> fit -> verify -> write back).
// Builds an all-BSpline solid (cylinder - wedge, each face ToBSpline'd), runs ProcessFace per face, asserts
// decision/accuracy/G1 bounds, then confirms a perturbed face is rejected and left unwritten.
// Lengths in mm (kAccTol=1e-5 mm); G1 angles in rad. Exit: 0 pass / 1 fail / 2 exception.

#include "sgkquadric/pipeline.hpp"
#include "sgkquadric/types.hpp"

#include <Foundation/init.h>
#include <Foundation/RuntimeClass.h>

#include <GeomBase/Dir3D.h>
#include <GeomBase/Number.h>
#include <GeomBase/Point3D.h>
#include <GeomBase/Ucs3D.h>
#include <GeomBase/UVRange.h>
#include <GeomBase/Vector3D.h>

#include <Geometry/SurfaceType.h>
#include <Geometry/3D/Surface/BSplineSurface.h>
#include <Geometry/3D/Surface/Surface.h>

#include <ModelingBase/ModelingRet.h>
#include <ModelingPrim/API.h>
#include <ModelingPrim/Primitives/PrimitivesRet.h>
#include <Boolean/API.h>
#include <Boolean/BooleanOpts.h>
#include <Boolean/BooleanRet.h>

#include <Topology/Brep/Body.h>
#include <Topology/Brep/Face.h>
#include <Topology/Geom/ModelSurface.h>
#include <Topology/Tools/TopoBuilder.h>

#include <algorithm>
#include <cmath>   // std::isfinite
#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace sggk;
using namespace sgkq;

namespace
{

// Non-atomic counters, single-threaded; consumed by the [SUMMARY] line and exit code.
int g_total       = 0;
int g_failed      = 0;   // >0 -> exit 1
int g_converted   = 0;
int g_rejected    = 0;
int g_failedState = 0;

// Write-back accuracy threshold [mm]: maxDist / boundaryCurveHausdorff / boundaryOnSurface must all be <= this.
constexpr double kAccTol = 1e-5;

// Print one assertion and update counters; format "[TEST] name expected=.. actual=.. PASS|FAIL".
void check(const std::string& name, const std::string& expected,
           const std::string& actual, bool pass)
{
    ++g_total;
    if (!pass) ++g_failed;
    std::cout << "[TEST] " << name
              << " expected=" << expected
              << " actual=" << actual
              << (pass ? " PASS" : " FAIL") << std::endl;
}

// Decision enum -> stable text label; unlisted values yield "?".
const char* DecisionStr(Decision d)
{
    switch (d)
    {
    case Decision::Converted: return "Converted";
    case Decision::Rejected:  return "Rejected";        // accuracy/G1 below bar
    case Decision::Failed:    return "Failed";
    default:                  return "?";
    }
}

// SurfaceType enum -> stable text label; full enum coverage, unknown -> "Unknown".
const char* SurfaceTypeName(SurfaceType t)
{
    switch (t)
    {
    case SurfaceType::Plane:    return "Plane";
    case SurfaceType::Cylinder: return "Cylinder";
    case SurfaceType::Cone:     return "Cone";
    case SurfaceType::Sphere:   return "Sphere";
    case SurfaceType::Torus:    return "Torus";
    case SurfaceType::BSpline:  return "BSpline";
    case SurfaceType::Bezier:   return "Bezier";
    case SurfaceType::RectTrimmed: return "RectTrimmed";
    case SurfaceType::Offset:   return "Offset";
    case SurfaceType::Extrude:  return "Extrude";
    default:                    return "Unknown";
    }
}

// True for the five analytic surface types (the write-back target set).
bool IsAnalyticType(SurfaceType t)
{
    return t == SurfaceType::Plane || t == SurfaceType::Cylinder || t == SurfaceType::Cone
        || t == SurfaceType::Sphere || t == SurfaceType::Torus;
}

// Format a value in scientific notation with 6 fractional digits (unit depends on caller: mm or rad).
std::string Fmt(double v)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(6) << v;
    return os.str();
}

// Build an all-BSpline solid: cylinder minus wedge, then ToBSpline every face; throws on any failure.
BodyPtr MakeNurbsSolid()
{
    auto rc = api_make_solid_cylinder(Ucs3D(), 200.0, 500.0, PI2, true);   // r=200mm, h=500mm
    if (!rc || !rc->Succeeded())
        throw std::runtime_error("api_make_solid_cylinder failed");
    BodyPtr cyl = rc->ResultBody();

    // Wedge is a box x in [100,300], y in [100,300], z in [-100,600]; non-tangent with the cylinder
    // avoids boolean degeneracy (its top/bottom are not coplanar with cylinder z=0 / z=500).
    auto rw = api_make_solid_wedge(Point3D(100.0, 100.0, -100.0),
                                   Point3D(300.0, 300.0, 600.0));
    if (!rw || !rw->Succeeded())
        throw std::runtime_error("api_make_solid_wedge failed");
    BodyPtr wedge = rw->ResultBody();

    BooleanOpts opts;
    opts.SetBooleanType(BooleanType::SUBTRACTION);                         // cyl - wedge
    auto br = api_boolean(cyl, wedge, opts);
    if (!br || !br->Succeeded())
        throw std::runtime_error("api_boolean(SUBTRACTION) failed");

    // A boolean may yield several bodies; keep the first, falling back to the cylinder if empty.
    const BodyList& results = br->ResultBodies();
    BodyPtr body = results.empty() ? cyl : results.front();
    if (!body)
        throw std::runtime_error("boolean result body is null");

    const FaceList faces = body->QueryFaces();
    // Convert every remaining face to an equivalent spline surface.
    for (const FacePtr& f : faces)
    {
        if (!f) continue;
        SurfacePtr srf = f->GeomSurface();
        if (!srf) continue;
        if (srf->IsKindOf<BSplineSurface>()) continue;                     // already a spline

        // CalcUVBound() supplies the UV domain used as the ToBSpline reparameterization range.
        const UVRange uv = f->CalcUVBound();
        BSplineSurfacePtr bs = srf->ToBSpline(uv);
        if (!bs)
            throw std::runtime_error("ToBSpline returned null during NURBS-ification");

        // Re-attach as ModelSurface (false = no geometry copy) keeping the original Sense() so the
        // face's normal convention is preserved; the later G1 checks rely on this.
        const ModelSurfacePtr ms = TopoBuilder::MakeModelSurface(bs, false);
        TopoBuilder::FaceAttachSurface(f, ms, f->Sense());
    }
    return body;
}

// End-to-end flow: build all-spline solid -> ProcessFace per face -> accuracy/G1 checks -> perturbed
// negative case. Always returns 0; pass/fail rides on the check() counters. Throws on setup failure.
int Run()
{
    std::cout << std::setprecision(6);

    // Units: lengths in mm (kAccTol=1e-5 mm), G1 angles in rad (th.g1AngleRad = 0.1 deg).
    // 1) NURBS-only solid: every face must be a BSplineSurface.
    BodyPtr body = MakeNurbsSolid();
    // Fetch faces once; the per-face pipeline pass below reuses this list.
    const FaceList faces = body->QueryFaces();
    std::cout << "[INFO] step1: nurbs-only solid face count = " << faces.size() << std::endl;

    int    nonSpline = 0;   // non-spline faces seen; must stay 0
    std::size_t fi = 0;
    for (const FacePtr& f : faces)
    {
        const SurfacePtr srf = f ? f->GeomSurface() : SurfacePtr();
        const SurfaceType st = srf ? srf->SurfaceType() : SurfaceType::Unknown;
        const bool isBs = srf && srf->IsKindOf<BSplineSurface>();
        std::cout << "[INFO] step1 face#" << fi << " SurfaceType=" << SurfaceTypeName(st)
                  << " IsBSpline=" << (isBs ? "true" : "false") << std::endl;
        if (!isBs) ++nonSpline;
        ++fi;
    }
    // A non-spline face would invalidate steps 2-4, so assert the input precondition first.
    check("step1.allFacesBSpline", "0 non-spline faces",
          std::to_string(nonSpline) + " non-spline", nonSpline == 0);
    if (faces.empty())
        throw std::runtime_error("no faces produced; aborting");  // later steps are meaningless without faces

    // 2) Every face should convert; the written-back type must be one of the five analytic ones.
    // Defaults: accuracy tolerances and the G1 angle threshold (see Thresholds).
    const Thresholds      th;
    const PipelineOptions opt;

    // Keep per-face reports for the accuracy/G1 passes so ProcessFace is not called twice per face.
    std::vector<FaceReport> reports;
    reports.reserve(faces.size());

    std::size_t idx = 0;
    // Run the pipeline on every face; classify each decision for the [SUMMARY] counters.
    for (const FacePtr& f : faces)
    {
        const FaceReport r = ProcessFace(f, th, opt);       // interface under test: decide -> fit -> verify -> write back
        reports.push_back(r);                               // keep for the accuracy/G1 passes below

        switch (r.decision)
        {
        case Decision::Converted: ++g_converted;   break;
        case Decision::Rejected:  ++g_rejected;    break;
        default:                  ++g_failedState; break;   // unexpected path
        }

        // [REPORT]: all verification quantities for this face (lengths mm, angles rad).
        std::cout << "[REPORT] face#" << idx
                  << " id=" << r.faceId
                  << " src=" << ToString(r.sourceType)
                  << " decision=" << DecisionStr(r.decision)
                  << " target=" << ToString(r.params.type)
                  << " maxDist=" << Fmt(r.accuracy.maxDist)                // [mm]
                  << " avgDist=" << Fmt(r.accuracy.avgDist)                // [mm]
                  << " boundaryCurveHausdorff=" << Fmt(r.accuracy.boundaryCurveHausdorff)  // two boundary curves, bidirectional Hausdorff [mm]
                  << " boundaryOnSurface=" << Fmt(r.accuracy.boundaryOnSurface)  // 3D boundary edge -> new surface, max orthogonal dist [mm]
                  << " g1Before=" << Fmt(r.accuracy.g1Before)              // [rad]
                  << " g1After=" << Fmt(r.accuracy.g1After)                // [rad]
                  << " err=" << r.errorCode << std::endl;                  // 0 = none
        std::cout << "[MSG] face#" << idx << " " << r.message << std::endl;

        // This synthetic solid is exactly analytic, so every face must convert.
        check("step2.face" + std::to_string(idx) + ".converted",
              "Converted", DecisionStr(r.decision), r.decision == Decision::Converted);

        // Read back the type: must be one of the five analytic types, proving a real replacement.
        const SurfacePtr s2 = f ? f->GeomSurface() : SurfacePtr();
        const SurfaceType t2 = s2 ? s2->SurfaceType() : SurfaceType::Unknown;
        const bool analytic = IsAnalyticType(t2);
        // A successful decision alone is not enough; the face geometry itself must have changed.
        check("step2.face" + std::to_string(idx) + ".analyticAfter",
              "one of {Plane,Cylinder,Cone,Sphere,Torus}",
              SurfaceTypeName(t2), analytic);

        ++idx;
    }

    // 3) Accuracy on converted faces: maxDist / boundaryCurveHausdorff / boundaryOnSurface all <= kAccTol [mm].
    // Only converted faces carry meaningful accuracy numbers.
    for (std::size_t i = 0; i < reports.size(); ++i)
    {
        const FaceReport& r = reports[i];
        if (r.decision != Decision::Converted)
            continue;                                       // rejected/failed faces are unwritten

        // maxDist: max point deviation inside the face.
        // boundaryCurveHausdorff: bidirectional Hausdorff between the two boundary curves (input 3D boundary
        //   curve vs. boundary curve generated on the new analytic surface) -- the PDF accuracy target.
        // boundaryOnSurface: max orthogonal distance from the original 3D boundary edge to the new surface.
        // All three deviations are lengths [mm], each compared against kAccTol.
        const bool ok = (r.accuracy.maxDist <= kAccTol)
                     && (r.accuracy.boundaryCurveHausdorff <= kAccTol)
                     && (r.accuracy.boundaryOnSurface <= kAccTol);
        check("step3.face" + std::to_string(i) + ".accuracy",
              "maxDist<=1e-5 && boundaryCurveHausdorff<=1e-5 && boundaryOnSurface<=1e-5",
              "maxDist=" + Fmt(r.accuracy.maxDist)
                  + ",boundaryCurveHausdorff=" + Fmt(r.accuracy.boundaryCurveHausdorff)
                  + ",boundaryOnSurface=" + Fmt(r.accuracy.boundaryOnSurface),
              ok);

        // Separate line exposing the boundaryCurveHausdorff magnitude for reports; same tolerance as above.
        check("step3.face" + std::to_string(i) + ".newFieldBoundaryCurveHausdorff",
              "printed for report reference (<=1e-5 mm)",
              "boundaryCurveHausdorff=" + Fmt(r.accuracy.boundaryCurveHausdorff) + " mm",
              std::isfinite(r.accuracy.boundaryCurveHausdorff) && (r.accuracy.boundaryCurveHausdorff <= kAccTol));
    }

    // 4) G1 must not degrade: g1After <= max(g1Before, th.g1AngleRad); the threshold leaves headroom
    //    for seams that were already poor before replacement.
    // Only converted (written-back) faces count.
    for (std::size_t i = 0; i < reports.size(); ++i)
    {
        const FaceReport& r = reports[i];
        if (r.decision != Decision::Converted)
            continue;

        // Bound = max(before, threshold): an already-poor seam may stay poor, but must not worsen.
        const double bound = (std::max)(r.accuracy.g1Before, th.g1AngleRad);  // allowed G1 angle bound [rad]
        const bool   ok    = r.accuracy.g1After <= bound;
        check("step4.face" + std::to_string(i) + ".g1NoWorse",
              "g1After<=max(g1Before,th.g1AngleRad)=" + Fmt(bound),
              "g1Before=" + Fmt(r.accuracy.g1Before) + ",g1After=" + Fmt(r.accuracy.g1After),
              ok);
    }

    // 5) Negative case: perturb one face's spline control points -> must be Rejected/Failed, face unchanged.
    {
        // Separate cylinder (r=200mm, h=500mm) so steps 1-4 are not affected.
        const double radius = 200.0;
        auto rc = api_make_solid_cylinder(Ucs3D(), radius, 500.0, PI2, true);
        if (!rc || !rc->Succeeded())
            throw std::runtime_error("negative: api_make_solid_cylinder failed");
        BodyPtr negBody = rc->ResultBody();

        // The lateral face is the unique Cylinder face; it is the perturbation target.
        FacePtr lat;
        for (const FacePtr& f : negBody->QueryFaces())
        {
            if (!f || !f->GeomSurface()) continue;
            if (f->GeomSurface()->SurfaceType() == SurfaceType::Cylinder) { lat = f; break; }
        }
        if (!lat)
            throw std::runtime_error("negative: lateral cylinder face not found");

        // Lateral face -> spline face, same as MakeNurbsSolid (keep sense, no geometry copy).
        {
            const UVRange uv = lat->CalcUVBound();
            BSplineSurfacePtr bs = lat->GeomSurface()->ToBSpline(uv);
            if (!bs)
                throw std::runtime_error("negative: ToBSpline returned null");
            const ModelSurfacePtr ms = TopoBuilder::MakeModelSurface(bs, false);
            TopoBuilder::FaceAttachSurface(lat, ms, lat->Sense());
        }

        SurfacePtr negSrf = lat->GeomSurface();
        BSplineSurfacePtr bs = sg_cast<BSplineSurface>(negSrf);
        if (!bs)
            throw std::runtime_error("negative: face is not BSpline after conversion");

        // Perturb a mid control point radially by 5% of the radius: on a cylinder the radial offset most
        // directly breaks the constant-radius constraint, pushing the fit beyond tolerance.
        const Point3DMatrix& cps = bs->ControlPoints();
        if (cps.empty() || cps.front().empty())
            throw std::runtime_error("negative: empty control point matrix");
        const UInt iu = static_cast<UInt>(cps.size() / 2);                  // mid U index
        const UInt iv = static_cast<UInt>(cps[iu].size() / 2);              // mid V index

        const Point3D q0 = cps[iu][iv];
        // Radial direction in the z=0 plane from the axis to q0; degenerate on-axis -> +X.
        Vector3D radial(q0.X(), q0.Y(), 0.0);
        if (!(radial.Length() > 0.0)) radial = Vector3D(1.0, 0.0, 0.0);
        radial.Normalize();
        const Point3D q1 = q0 + radial * (0.05 * radius);                   // 0.05*200mm = 10 mm
        bs->SetControlPoint(q1, iu, iv, false);                            // false = do not refresh periodic/closed cache

        // ProcessFace again: expect rejection due to out-of-tolerance accuracy/G1, no write-back.
        const FaceReport r = ProcessFace(lat, th, opt);
        // Classify the negative-case decision so [SUMMARY] stays consistent.
        switch (r.decision)
        {
        case Decision::Converted: ++g_converted;   break;   // should not happen
        case Decision::Rejected:  ++g_rejected;    break;   // expected
        default:                  ++g_failedState; break;
        }

        // Report the rejection residual and error code for diagnosis.
        std::cout << "[REPORT] negative id=" << r.faceId
                  << " decision=" << DecisionStr(r.decision)               // expect Rejected/Failed
                  << " target=" << ToString(r.params.type)
                  << " maxDist=" << Fmt(r.accuracy.maxDist)               // out-of-tolerance residual [mm]
                  << " boundaryCurveHausdorff=" << Fmt(r.accuracy.boundaryCurveHausdorff)  // [mm]
                  << " err=" << r.errorCode << std::endl;
        std::cout << "[MSG] negative " << r.message << std::endl;

        // Must be Rejected or Failed -- never silently write back out-of-tolerance geometry.
        const bool notConverted = (r.decision == Decision::Rejected || r.decision == Decision::Failed);
        check("step5.perturbedNotConverted", "Rejected or Failed",
              DecisionStr(r.decision), notConverted);

        // Rejection must leave the face untouched (still the original spline surface).
        const SurfacePtr s5 = lat->GeomSurface();
        const bool stillSpline = s5 && s5->IsKindOf<BSplineSurface>();
        check("step5.perturbedStillSpline", "true",
              std::string(stillSpline ? "true" : "false"), stillSpline);
    }

    return 0;
}

} // namespace

// Entry point: init runtime -> Run() -> summary/exit code. Exceptions count as failures (exit 2).
int main(int /*argc*/, char* /*argv*/[])
{
    try
    {
        init();
        Run();
        fini();
    }
    catch (const std::exception& e)
    {
        // Unexpected exception: count one failure; fini() is skipped (runtime state untrusted).
        std::cerr << "[FAIL] exception: " << e.what() << std::endl;
        ++g_failed;
        return 2;
    }
    catch (...)
    {
        std::cerr << "[FAIL] unknown exception" << std::endl;
        ++g_failed;
        return 2;
    }

    // Summary line for CI; must agree with the exit code below.
    std::cout << "[SUMMARY] total=" << g_total
              << " failed=" << g_failed
              << " converted=" << g_converted
              << " rejected=" << g_rejected
              << " failedState=" << g_failedState << std::endl;

    return (g_failed == 0) ? 0 : 1;
}
