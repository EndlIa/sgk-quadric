// apps/checks/test_g1_tangent/main.cpp -- G1 evidence for a tangent neighbor face (PDF accuracy target #3).
// Builds a solid whose fillet (cylinder) face is tangent to two planes, then measures the tangent seam G1
// in three states: native analytic -> ToBSpline replacement (baseline) -> pipeline back to analytic.
// Claims g1After <= 0.1 deg (1.7453292519943296e-3 rad) and no degradation. Lengths mm, angles rad.

#include "sgkquadric/pipeline.hpp"
#include "sgkquadric/types.hpp"

#include <Foundation/init.h>
#include <Foundation/RuntimeClass.h>

#include <GeomBase/Dir3D.h>
#include <GeomBase/Number.h>
#include <GeomBase/Point3D.h>
#include <GeomBase/Toler.h>
#include <GeomBase/Ucs3D.h>
#include <GeomBase/UVRange.h>
#include <GeomBase/Vector3D.h>

#include <Geometry/3D/Curve/Line3D.h>
#include <Geometry/3D/Curve/TrimmedCurve3D.h>
#include <Geometry/3D/Geometry3D.h>
#include <Geometry/3D/Surface/BSplineSurface.h>
#include <Geometry/3D/Surface/Plane.h>
#include <Geometry/3D/Surface/Surface.h>
#include <Geometry/SurfaceType.h>

#include <ModelingBase/API.h>
#include <ModelingBase/ModelingRet.h>
#include <ModelingBase/Wire/JoinWiresRet.h>
#include <ModelingPrim/API.h>
#include <ModelingPrim/Extrude/ExtrudeRet.h>
#include <ModelingPrim/Fillet/FilletEdgesOpts.h>
#include <ModelingPrim/Fillet/FilletRet.h>

#include <Topology/Brep/Body.h>
#include <Topology/Brep/Edge.h>
#include <Topology/Brep/Face.h>
#include <Topology/Brep/Wire.h>
#include <Topology/Geom/ModelSurface.h>
#include <Topology/Tools/TopoBuilder.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

using namespace sggk;
using namespace sgkq;

namespace
{

// Non-atomic counters, single-threaded.
int g_total  = 0;
int g_failed = 0;   // >0 -> exit 1

// Three-state G1 readings on the tangent seams [rad]; -1 means the step did not run.
double g_g1Tangent = -1.0;   // state 1: native analytic face
double g_g1Before  = -1.0;   // state 2: after ToBSpline replacement (baseline)
double g_g1After   = -1.0;   // state 3: after pipeline conversion back to analytic

// Thresholds; kG1ThresholdRad = 0.1 deg per spec, matching sgkq::Thresholds defaults.
constexpr double kG1ThresholdRad = 1.7453292519943296e-3;   // rad
constexpr double kAccTol         = 1e-5;                    // mm, maxDist / boundaryCurveHausdorff
constexpr double kTangentEps     = 1e-6;                    // rad, "tangent" cutoff (far below 0.1 deg)

// Geometry [mm]: 100x100 square profile, R20 corner fillet, extruded 100 along +Z.
constexpr double kSide   = 100.0;
constexpr double kRadius = 20.0;
constexpr double kHeight = 100.0;

// Print one assertion and update counters.
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

// Decision enum -> stable text label.
const char* DecisionStr(Decision d)
{
    switch (d)
    {
    case Decision::Converted: return "Converted";
    case Decision::Rejected:  return "Rejected";
    case Decision::Failed:    return "Failed";
    default:                  return "?";
    }
}

// SurfaceType enum -> stable text label.
const char* SurfaceTypeName(SurfaceType t)
{
    switch (t)
    {
    case SurfaceType::Plane:       return "Plane";
    case SurfaceType::Cone:        return "Cone";
    case SurfaceType::Cylinder:    return "Cylinder";
    case SurfaceType::Torus:       return "Torus";
    case SurfaceType::Sphere:      return "Sphere";
    case SurfaceType::RectTrimmed: return "RectTrimmed";
    case SurfaceType::Offset:      return "Offset";
    case SurfaceType::BSpline:     return "BSpline";
    case SurfaceType::Bezier:      return "Bezier";
    case SurfaceType::Extrude:     return "Extrude";
    default:                       return "Unknown";
    }
}

// True for the five analytic surface types (the conversion target set).
bool IsAnalyticType(SurfaceType t)
{
    return t == SurfaceType::Plane || t == SurfaceType::Cylinder || t == SurfaceType::Cone
        || t == SurfaceType::Sphere || t == SurfaceType::Torus;
}

// Format a value in scientific notation with 6 fractional digits (mm or rad).
std::string Fmt(double v)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(6) << v;
    return os.str();
}

// Result of G1 measurement over a set of edges: max normal-angle deviation [rad] plus counts.
struct G1Measurement
{
    double maxVal   = 0.0;
    int    measured = 0;
    int    skipped  = 0;     // not-two-face edges, Default continuity, or exceptions
};

// Call the kernel G1 continuity check per edge and take the max maxVal.
// Same convention as src/pipeline/pipeline.cpp MeasureG1: an edge must have exactly two faces;
// Default continuity means "not measurable".
G1Measurement MeasureEdgesG1(const std::vector<EdgePtr>& edges)
{
    G1Measurement r;

    for (const EdgePtr& e : edges)
    {
        if (!e) { ++r.skipped; continue; }

        FaceList fs;
        try { fs = e->QueryFaces(); }
        catch (...) { ++r.skipped; continue; }

        if (fs.size() != 2) { ++r.skipped; continue; }           // G1 undefined unless exactly 2 faces

        try
        {
            const FacesContinuity fc =                           // kernel G1: 100 samples, DefAngleTol
                e->CalcFacesContinuity(ContinuityType::G1, 100, Toler::DefAngleTol());
            if (fc.continuityType == ContinuityType::Default)    // Default = no face or >2 faces
            {
                ++r.skipped;
                continue;
            }
            if (fc.maxVal > r.maxVal) r.maxVal = fc.maxVal;
            ++r.measured;
        }
        catch (...) { ++r.skipped; continue; }                   // evaluation failed
    }
    return r;
}

// Build a solid with a cylindrical fillet face tangent to two side planes: kernel in-plane edge-edge fillet
// (R=20) at corner (100,0) gives a 90-degree arc tangent to both straight edges at (80,0) and (100,20);
// combine it with 4 straight edges into one closed profile, then extrude along +Z so the arc becomes a
// cylinder tangent to the y=0 and x=100 planes. Fails (why set) if any kernel call returns null/failed.
BodyPtr BuildRoundedSolid(std::string& why)
{
    why.clear();

    // a. Two straight corner edges (standalone, not part of any wire/loop).
    TrimmedCurve3D line1(Point3D(0.0, 0.0, 0.0), Point3D(kSide, 0.0, 0.0));
    TrimmedCurve3D line2(Point3D(kSide, 0.0, 0.0), Point3D(kSide, kSide, 0.0));
    EdgePtr edge1 = TopoBuilder::MakeEdge(line1);
    EdgePtr edge2 = TopoBuilder::MakeEdge(line2);
    if (!edge1 || !edge2) { why = "TopoBuilder::MakeEdge returned null"; return nullptr; }

    // b. In-plane edge-edge fillet; gives an arc tangent to the two straight edges.
    //    fopt = (trimEdge1, trimEdge2, keepEdge1Front, keepEdge2Front), matching the SDK fillet example.
    FilletEdgesOpts fopt(true, true, true, false);
    FilletRetPtr fret = api_fillet_edges(edge1, edge2, kRadius, fopt);
    if (!fret)
    {
        why = "api_fillet_edges returned null (fillet unavailable)";
        return nullptr;
    }
    if (!fret->Succeeded())
    {
        why = "api_fillet_edges failed: " + fret->Status().ErrorMsg();
        return nullptr;
    }
    EdgePtr arc = fret->FilletEdge();
    if (!arc)
    {
        why = "FilletRet::FilletEdge() is null";
        return nullptr;
    }

    // Arc endpoints are the two tangent points: one lies on y=0, the other on x=100.
    Point3D pa, pb;
    try
    {
        pa = arc->StartPoint();
        pb = arc->EndPoint();
    }
    catch (const std::exception& e)
    {
        why = std::string("reading arc endpoints threw: ") + e.what();
        return nullptr;
    }

    Point3D pBottom, pRight;
    if (std::fabs(pa.Y()) < 1e-6)      { pBottom = pa; pRight = pb; }
    else if (std::fabs(pb.Y()) < 1e-6) { pBottom = pb; pRight = pa; }
    else
    {
        why = "arc endpoints not on edges y=0 / x=100 (kernel fillet solved on the other side?)";
        return nullptr;
    }

    // Tangency evidence: each tangent point should be R (=20) from corner (100,0).
    const Point3D corner(kSide, 0.0, 0.0);
    const double dBottom = (pBottom - corner).Length();
    const double dRight  = (pRight - corner).Length();
    std::cout << "[INFO] fillet tangent points bottom=(" << pBottom.X() << "," << pBottom.Y() << "," << pBottom.Z()
              << ") distToCorner=" << dBottom
              << " ; right=(" << pRight.X() << "," << pRight.Y() << "," << pRight.Z()
              << ") distToCorner=" << dRight
              << " ; radius=" << kRadius << std::endl;

    // c. Closed profile: arc + 4 straight edges, endpoints coincident with the tangent points.
    EdgeList loopEdges;                                                             // EdgeList = std::list<EdgePtr>
    loopEdges.push_back(arc);
    loopEdges.push_back(TopoBuilder::MakeEdge(TrimmedCurve3D(Point3D(0.0, 0.0, 0.0), pBottom)));
    loopEdges.push_back(TopoBuilder::MakeEdge(TrimmedCurve3D(pRight, Point3D(kSide, kSide, 0.0))));
    loopEdges.push_back(TopoBuilder::MakeEdge(TrimmedCurve3D(Point3D(kSide, kSide, 0.0), Point3D(0.0, kSide, 0.0))));
    loopEdges.push_back(TopoBuilder::MakeEdge(TrimmedCurve3D(Point3D(0.0, kSide, 0.0), Point3D(0.0, 0.0, 0.0))));
    for (const EdgePtr& e : loopEdges)
    {
        if (!e) { why = "profile edge construction failed (MakeEdge returned null)"; return nullptr; }
    }

    JoinWiresRetPtr jret = api_compose_edges_to_wires(loopEdges);
    if (!jret || !jret->Succeeded())
    {
        why = "api_compose_edges_to_wires failed: " + (jret ? jret->Status().ErrorMsg() : std::string("returned null"));
        return nullptr;
    }
    if (jret->Wires().size() != 1)                                                  // must merge into exactly one wire
    {
        why = "profile did not merge into a single closed Wire (got " + std::to_string(jret->Wires().size()) + " wires)";
        return nullptr;
    }
    WirePtr wire = jret->Wires().front();

    // d. Wire -> planar profile face on z=0 (normal +Z).
    SurfacePtr planeSrf = std::make_shared<Plane>(Ucs3D());
    FacePtr    profile  = api_create_face(wire, planeSrf);
    if (!profile)
    {
        why = "api_create_face returned null (profile not planar / not closed?)";
        return nullptr;
    }

    // e. Extrude along +Z into a solid.
    ExtrudeRetPtr eret = api_extrude_entity(profile, Dir3D(0.0, 0.0, 1.0), 0.0, kHeight);
    if (!eret || !eret->Succeeded() || eret->ResultBodies().empty())
    {
        why = "api_extrude_entity failed: " + (eret ? eret->Status().ErrorMsg() : std::string("returned null"));
        return nullptr;
    }
    return eret->ResultBodies().front();
}

// Return the body's unique Cylinder face (the fillet face); empty if count != 1.
FacePtr FindFilletFace(const BodyPtr& body, int& cylinderCount, std::size_t& faceCount)
{
    cylinderCount = 0;
    faceCount     = 0;
    FacePtr found;

    for (const FacePtr& f : body->QueryFaces())
    {
        ++faceCount;
        const SurfacePtr srf = f ? f->GeomSurface() : SurfacePtr();
        const SurfaceType st = srf ? srf->SurfaceType() : SurfaceType::Unknown;
        if (st == SurfaceType::Cylinder) { ++cylinderCount; found = f; }
    }
    return (cylinderCount == 1) ? found : FacePtr();
}

// Three-state G1 comparison: native tangent -> spline replacement -> pipeline back to analytic.
int Run()
{
    std::cout << std::setprecision(6);

    // 1) Build the tangent-neighbor solid; report face count and per-face types.
    std::string why;
    BodyPtr body = BuildRoundedSolid(why);
    // If the fillet-based construction fails, report why and stop.
    check("step1.buildRoundedSolid", "solid built", body ? "ok" : ("failed: " + why), body != nullptr);
    if (!body) return 0;                                     // setup failure: later steps are meaningless

    int         cylinderCount = 0;
    std::size_t faceCount     = 0;
    FacePtr     filletFace    = FindFilletFace(body, cylinderCount, faceCount);

    std::cout << "[INFO] step1 faceCount=" << faceCount << std::endl;
    {
        std::size_t i = 0;
        for (const FacePtr& f : body->QueryFaces())
        {
            const SurfacePtr srf = f ? f->GeomSurface() : SurfacePtr();
            const SurfaceType st = srf ? srf->SurfaceType() : SurfaceType::Unknown;
            std::cout << "[INFO] step1 face#" << i << " SurfaceType=" << SurfaceTypeName(st) << std::endl;
            ++i;
        }
    }
    // The fillet face is the unique Cylinder face among the body's faces.
    check("step1.filletFaceFound", "exactly 1 Cylinder face",
          std::to_string(cylinderCount) + " Cylinder", filletFace != nullptr);
    if (!filletFace) return 0;

    // 2) Confirm tangency: collect the fillet face's edges; those shared with a plane and lying on a line
    //    are the two tangent seams (cylinder rulings); the rest are 90-degree end arcs.
    std::vector<EdgePtr> seams;      // two tangent seams
    std::vector<EdgePtr> hardEdges;  // two 90-degree end arcs (anchor)

    for (const EdgePtr& e : body->QueryEdges())
    {
        if (!e) continue;
        FaceList fs;
        try { fs = e->QueryFaces(); } catch (...) { continue; }
        if (fs.size() != 2) continue;

        FacePtr other;
        bool    atFillet = false;
        for (const FacePtr& f : fs)
        {
            if (f == filletFace) atFillet = true;
            else                 other    = f;
        }
        if (!atFillet) continue;
        if (!other || !other->GeomSurface()) continue;
        const SurfaceType otherType = other->GeomSurface()->SurfaceType();
        const Curve3DPtr  curve     = e->GeomCurve();
        const bool        isLine    = curve && curve->IsKindOf<Line3D>();

        if (otherType == SurfaceType::Plane && isLine) seams.push_back(e);
        else                                          hardEdges.push_back(e);
    }

    std::cout << "[INFO] step2 seamCount=" << seams.size() << " hardEdgeCount=" << hardEdges.size() << std::endl;
    // Exactly two plane-adjacent line edges must exist; otherwise the construction is wrong.
    check("step2.seamCount", "2", std::to_string(seams.size()), seams.size() == 2);

    // Native-state G1 on both tangent seams.
    const G1Measurement mTangent = MeasureEdgesG1(seams);
    g_g1Tangent = (mTangent.measured > 0) ? mTangent.maxVal : -1.0;
    std::cout << "[REPORT] step2 tangentSeamG1 maxVal=" << Fmt(mTangent.maxVal)
              << " rad measured=" << mTangent.measured << " skipped=" << mTangent.skipped << std::endl;
    // Tangent means both seams are measurable and the angle is < 1e-6 rad (~1/1700 of 0.1 deg).
    check("step2.tangentResidualLt1e-6", "<1e-6 rad",
          Fmt(mTangent.maxVal) + " rad (measured=" + std::to_string(mTangent.measured) + ")",
          mTangent.measured == 2 && mTangent.maxVal < kTangentEps);

    // Anchor: a 90-degree end arc must measure ~pi/2, proving the metric can distinguish hard edges.
    const G1Measurement mHard = MeasureEdgesG1(hardEdges);
    const double         piHalf = PI2 / 4.0;                                         // PI2 = 2*pi
    std::cout << "[REPORT] step2 hardEdgeG1 maxVal=" << Fmt(mHard.maxVal)
              << " rad measured=" << mHard.measured << " skipped=" << mHard.skipped << std::endl;
    check("step2.hardEdgeAnchor", "|maxVal-PI/2|<1e-3 rad (PI/2=" + Fmt(piHalf) + ")",
          Fmt(mHard.maxVal) + " rad", mHard.measured >= 1 && std::fabs(mHard.maxVal - piHalf) < 1e-3);

    // 3) Build the spline input: fillet face -> ToBSpline + FaceAttachSurface.
    const SurfaceType stNative = filletFace->GeomSurface()->SurfaceType();
    std::cout << "[INFO] step3 filletFace SurfaceType(native)=" << SurfaceTypeName(stNative) << std::endl;
    // The native fillet face must be analytic before we convert it.
    check("step3.filletFaceAnalytic", "Cylinder or Torus", SurfaceTypeName(stNative),
          stNative == SurfaceType::Cylinder || stNative == SurfaceType::Torus);

    {
        const SurfacePtr nativeSrf = filletFace->GeomSurface();
        const UVRange    uv        = filletFace->CalcUVBound();
        // Convert the native fillet surface to a spline; must succeed.
        BSplineSurfacePtr bs = nativeSrf->ToBSpline(uv);
        if (!bs)
        {
            check("step3.toBSpline", "non-null", "null", false);
            return 0;
        }
        // Re-attach keeping sense and without copying geometry (false) so the normal convention holds.
        const ModelSurfacePtr ms = TopoBuilder::MakeModelSurface(bs, false);
        TopoBuilder::FaceAttachSurface(filletFace, ms, filletFace->Sense());
    }

    const SurfacePtr srfNow  = filletFace->GeomSurface();
    const bool       isSpline = srfNow && srfNow->IsKindOf<BSplineSurface>();
    // After the swap the face must really hold a spline surface.
    check("step3.isBSplineAfterReplace", "true", isSpline ? "true" : "false", isSpline);
    if (!isSpline) return 0;

    // After the spline swap, the two tangent seams give the pre-conversion baseline g1Before.
    const G1Measurement mBefore = MeasureEdgesG1(seams);
    g_g1Before = mBefore.maxVal;
    std::cout << "[REPORT] step3 g1Before(splineFace) maxVal=" << Fmt(mBefore.maxVal)
              << " rad measured=" << mBefore.measured << " skipped=" << mBefore.skipped << std::endl;
    // Both seams must stay measurable after the swap, else the baseline is meaningless.
    check("step3.seamsMeasurableAfterSpline", "2",
          std::to_string(mBefore.measured), mBefore.measured == 2);

    // 4) Run the pipeline to convert the face back to an analytic surface.
    // Library defaults; the G1 angle threshold is 0.1 deg.
    const Thresholds      th;                                                        // defaults: 1e-5 / 1e-5 / 0.1 deg
    const PipelineOptions opt;
    const FaceReport      rep = ProcessFace(filletFace, th, opt);

    std::cout << "[REPORT] step4 decision=" << DecisionStr(rep.decision)
              << " src=" << ToString(rep.sourceType)
              << " target=" << ToString(rep.params.type)
              << " maxDist=" << Fmt(rep.accuracy.maxDist)
              << " avgDist=" << Fmt(rep.accuracy.avgDist)
              << " boundaryCurveHausdorff=" << Fmt(rep.accuracy.boundaryCurveHausdorff)
              << " faceAllEdgeG1Before=" << Fmt(rep.accuracy.g1Before)
              << " faceAllEdgeG1After=" << Fmt(rep.accuracy.g1After)
              << " err=" << rep.errorCode << std::endl;
    std::cout << "[MSG] step4 " << rep.message << std::endl;
    std::cout << "[NOTE] faceAllEdgeG1* is the max angular deviation over all edges of the face, incl. 2 90-deg end arcs, hence ≈PI/2;"
              << "        this metric only considers the two tangent seams (see step5)" << std::endl;

    // The pipeline must convert this face successfully and store analytic geometry again.
    check("step4.decisionConverted", "Converted", DecisionStr(rep.decision),
          rep.decision == Decision::Converted);

    const SurfacePtr srfAfter  = filletFace->GeomSurface();
    const SurfaceType stAfter  = srfAfter ? srfAfter->SurfaceType() : SurfaceType::Unknown;
    const bool        analytic = IsAnalyticType(stAfter);
    // Verify from the face's own geometry, not merely from the decision.
    check("step4.analyticAfter", "one of {Plane,Cylinder,Cone,Sphere,Torus}",
          SurfaceTypeName(stAfter), analytic);

    // 5) Measure G1 on the two tangent seams after conversion and run the core assertions.
    const G1Measurement mAfter = MeasureEdgesG1(seams);
    g_g1After = mAfter.maxVal;
    std::cout << "[REPORT] step5 g1After(analyticFace) maxVal=" << Fmt(mAfter.maxVal)
              << " rad measured=" << mAfter.measured << " skipped=" << mAfter.skipped << std::endl;

    // Core assertion (1): post-conversion G1 <= 0.1 deg (PDF accuracy target).
    check("step5.g1AfterLE0p1deg", "<=" + Fmt(kG1ThresholdRad) + " rad",
          Fmt(mAfter.maxVal) + " rad (measured=" + std::to_string(mAfter.measured) + ")",
          mAfter.measured == 2 && mAfter.maxVal <= kG1ThresholdRad);

    // Core assertion (2): no degradation -- g1After <= max(g1Before, 0.1 deg).
    const double bound = (std::max)(g_g1Before, kG1ThresholdRad);
    check("step5.g1AfterNotWorse", "<=" + Fmt(bound) + " rad",
          Fmt(mAfter.maxVal) + " rad", mAfter.maxVal <= bound);

    // Core assertions (3)(4): accuracy and boundary Hausdorff.
    check("step5.maxDistLE1e-5", "<=1e-5 mm", Fmt(rep.accuracy.maxDist) + " mm",
          rep.accuracy.maxDist <= kAccTol);
    check("step5.boundaryCurveHausdorffLE1e-5", "<=1e-5 mm", Fmt(rep.accuracy.boundaryCurveHausdorff) + " mm",
          rep.accuracy.boundaryCurveHausdorff <= kAccTol);

    return 0;
}

} // namespace

// Entry point: init runtime -> Run() -> summary/exit code.
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

    // Summary: three-state G1 values [rad] (-1 = step did not run).
    std::cout << "[SUMMARY] total=" << g_total
              << " failed=" << g_failed
              << " g1TangentBefore=" << g_g1Tangent
              << " g1Before=" << g_g1Before
              << " g1After=" << g_g1After
              << " g1ThresholdRad=1.7453292519943296e-3" << std::endl;

    return (g_failed == 0) ? 0 : 1;
}
