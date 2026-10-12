// sgkquadric face pipeline: detect -> fit -> refine -> verify -> write back an analytic quadric.
// Uses the SGK SDK only for type/topology queries (no init()/fini(); lifetime is the caller's).
// Boundary metric follows the PDF target: bidirectional Hausdorff between the input boundary line
// (coedge 3D curve) and the generated boundary line (coedge UV curve mapped to 3D through the face).

#include "sgkquadric/pipeline.hpp"

#include "sgkquadric/detect.hpp"
#include "sgkquadric/fit.hpp"
#include "sgkquadric/math.hpp"

#include <Foundation/RuntimeClass.h>

#include <GeomBase/Dir3D.h>
#include <GeomBase/Interval.h>
#include <GeomBase/Point3D.h>
#include <GeomBase/Toler.h>
#include <GeomBase/Ucs3D.h>
#include <GeomBase/UVRange.h>
#include <GeomBase/Vector3D.h>

#include <Geometry/SurfaceType.h>
#include <Geometry/3D/Geometry3D.h>
#include <Geometry/3D/Surface/BSplineSurface.h>
#include <Geometry/3D/Surface/ConeSurface.h>
#include <Geometry/3D/Surface/CylinderSurface.h>
#include <Geometry/3D/Surface/Plane.h>
#include <Geometry/3D/Surface/SphereSurface.h>
#include <Geometry/3D/Surface/TorusSurface.h>
#include <Geometry/3D/Curve/BoundedCurve3D.h>
#include <Geometry/2D/BoundedCurve2D.h>

#include <Topology/Brep/Coedge.h>
#include <Topology/Brep/Edge.h>
#include <Topology/Brep/Face.h>
#include <Topology/Geom/ModelSurface.h>
#include <Topology/Tools/TopoBuilder.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

using namespace sggk;

namespace sgkq
{
namespace
{

// errorCode values (report traceability only; they do not drive the verdict).
enum : int
{
    kOk                = 0,   // Converted
    kNullFace          = 1,   // null FacePtr
    kAlreadyAnalytic   = 2,   // already analytic, skipped
    kNotBSpline        = 3,   // not a spline surface (Bezier / Extrude / RectTrimmed ...)
    kDetectUnknown     = 4,   // detect returned Unknown
    kFitFailed         = 5,   // FitSurface failed
    kAccuracyReject    = 6,   // accuracy verification failed
    kWriteBackOff      = 7,   // verified but writeBack=false
    kWriteBackThrew    = 8,   // write-back threw
    kInternalThrew     = 9   // pipeline threw
};

// Boundary-metric sampling: 41 iso-parametric points per curve (A and B), endpoints included.
// 41 exceeds the 33-point verification grid, enough to catch mid-band and endpoint deviations;
// the ~2*41^2 pairwise distance work per coedge is negligible.
constexpr int kBoundarySamples = 41;   // iso-parametric sample points per boundary curve

//! @brief SDK SurfaceType -> Canonical (for report.sourceType).
Canonical CanonicalOf(SurfaceType t)
{
    switch (t)
    {
    case SurfaceType::Plane:    return Canonical::Plane;
    case SurfaceType::Cylinder: return Canonical::Cylinder;
    case SurfaceType::Cone:     return Canonical::Cone;
    case SurfaceType::Sphere:   return Canonical::Sphere;
    case SurfaceType::Torus:    return Canonical::Torus;
    default:                    return Canonical::Unknown;
    }
}

//! @brief True if t is already an analytic surface type to keep.
bool IsAnalyticSurfaceType(SurfaceType t)
{
    return t == SurfaceType::Plane || t == SurfaceType::Cylinder || t == SurfaceType::Cone
        || t == SurfaceType::Sphere || t == SurfaceType::Torus;
}

//! @brief True if the Interval is finite and bounded (usable for clipping the spline domain).
bool IsFiniteInterval(const Interval& i)
{
    return i.Type() == IntervalType::FINITE && std::isfinite(i.Min()) && std::isfinite(i.Max());
}

//! @brief Closed-form distance from q to the fitted analytic surface (delegates to math.hpp).
double DistToParams(const SurfaceParams& p, const Point3D& q)
{
    switch (p.type)
    {
    case Canonical::Plane:
        return distPlane(q, p.axisPoint, p.axisDir);
    case Canonical::Cylinder:
        return distCylinder(q, p.axisPoint, p.axisDir, p.radius);
    case Canonical::Cone:
        return distCone(q, p.apex, p.axisDir, p.halfAngle);
    case Canonical::Sphere:
        return distSphere(q, p.axisPoint, p.radius);
    case Canonical::Torus:
        return distTorus(q, p.axisPoint, p.axisDir, p.radius, p.minorRadius);
    default:
        return 0.0;
    }
}

// Boundary metric (PDF wording): bidirectional Hausdorff between the two boundary lines, plus the
// original boundary-to-surface orthogonal distance as an auxiliary quantity.
//   Curve A (input) = the coedge 3D curve (CalcBoundedCurve), 41 iso-parametric points.
//   Curve B (generated) = the coedge UV curve (CalcUVCurve, the pcurve on the face) sampled on its
//     own domain and mapped back to 3D via CalcPoint(u,v) of the surface carrying B: after write-back
//     that is face->GeomSurface() (the new analytic surface), before it the candidate analytic surface.
//   Hausdorff = max( max_{p in A} dist(p,B), max_{q in B} dist(q,A) ).
// The PDF target is <= 1e-5 mm. Discretizing both curves and taking the nearest SAMPLE slightly
// OVERestimates the true Hausdorff (the true nearest point can lie inside a segment), so this is a
// conservative bound. boundaryOnSurface (max distance from the 3D boundary edge to the new surface)
// is kept as an auxiliary check on pcurve correctness.

//! @brief Sampled input boundary line (curve A) of one coedge; empty pts = not measurable.
//! Index-aligned with the coedge list captured in stage 1 (unusable coedges still occupy an index).
struct BoundaryLineA
{
    std::vector<Point3D> pts;
};

//! @brief Bidirectional Hausdorff between two point sets (naive nearest-point search), in mm.
//! @return max( sup_{p in a} d(p,b), sup_{q in b} d(q,a) ); 0 if either side is empty.
//! O(|a|*|b|) (about 2*1681 evaluations for 41 points each) and, being sampled, a harmless
//! overestimate of the true Hausdorff.
double HausdorffSampled(const std::vector<Point3D>& a, const std::vector<Point3D>& b)
{
    double h = 0.0;

    if (a.empty() || b.empty())
    {
        return 0.0;
    }

    for (const Point3D& p : a)
    {
        double nearest = std::numeric_limits<double>::infinity();
        for (const Point3D& q : b)
        {
            const double d = (p - q).Length();
            if (d < nearest) nearest = d;
        }
        if (nearest > h) h = nearest;
    }

    for (const Point3D& q : b)
    {
        double nearest = std::numeric_limits<double>::infinity();
        for (const Point3D& p : a)
        {
            const double d = (q - p).Length();
            if (d < nearest) nearest = d;
        }
        if (nearest > h) h = nearest;
    }

    return h;
}

//! @brief Sample a coedge's generated boundary line (curve B): UV curve -> nSamples points -> CalcPoint.
//! @return 3D points (mm); fewer than 2 means B could not be taken (caller treats it as degenerate).
//! One 3D bounded curve may map to several 2D curves on the surface (Coedge.h:88: a circle through
std::vector<Point3D> SampleBoundaryLineB(const CoedgePtr& c, const SurfacePtr& srf, int nSamples)
{
    std::vector<Point3D> out;

    BoundedCurve2DArray uvCurves;
    try
    {
        uvCurves = c->CalcUVCurve();
    }
    catch (...)
    {
        return out;
    }
    if (uvCurves.empty())
    {
        return out;
    }

    const int k = static_cast<int>(uvCurves.size());
    std::vector<int> perCurve(static_cast<std::size_t>(k), nSamples / k);
    for (int i = 0; i < nSamples % k; ++i) ++perCurve[static_cast<std::size_t>(i)];

    for (int j = 0; j < k; ++j)
    {
        const BoundedCurve2DPtr& cv = uvCurves[static_cast<std::size_t>(j)];
        if (!cv)
        {
            continue;
        }

        Interval dom;
        try
        {
            dom = cv->Domain();
        }
        catch (...)
        {
            continue;
        }

        const double a = dom.Min();
        const double b = dom.Max();
        if (!(b > a) || !std::isfinite(a) || !std::isfinite(b))
        {
            continue;
        }

        const int m = perCurve[static_cast<std::size_t>(j)];
        for (int i = 0; i < m; ++i)
        {
            const double t = (m == 1)
                ? (0.5 * (a + b))
                : (a + (b - a) * (static_cast<double>(i) / static_cast<double>(m - 1)));
            try
            {
                const Point2D uv = cv->CalcPoint(t);
                out.push_back(srf->CalcPoint(uv));
            }
            catch (...)
            {
                continue;
            }
        }
    }

    return out;
}

//! @brief Boundary-metric stage 1: curve A (input boundary line) and boundaryOnSurface (max
//!        orthogonal distance from the 3D boundary to the candidate analytic surface), in mm.
//! Enumerates the face coedges and returns their 3D-curve samples; coedgesOut lets stage 2 reuse the
//! same pointers (write-back does not rebuild coedges, so order and pointers stay valid).
double MeasureBoundaryOnSurface(const FacePtr& face, const SurfaceParams& params,
                                int nSamples, std::vector<BoundaryLineA>& linesA,
                                CoedgeList& coedgesOut, bool& anyMeasured, int& skippedCoedges)
{
    double maxD = 0.0;
    anyMeasured = false;

    CoedgeList coedges;
    try
    {
        coedges = face->QueryCoedges();
    }
    catch (...)
    {
        return 0.0;
    }
    coedgesOut = coedges;

    linesA.reserve(linesA.size() + coedges.size());

    for (const CoedgePtr& c : coedges)
    {
        if (!c)
        {
            ++skippedCoedges;
            linesA.push_back(BoundaryLineA());
            continue;
        }

        BoundedCurve3DPtr crv;
        try
        {
            crv = c->CalcBoundedCurve();
        }
        catch (...)
        {
            ++skippedCoedges;
            linesA.push_back(BoundaryLineA());
            continue;
        }
        if (!crv)
        {
            ++skippedCoedges;
            linesA.push_back(BoundaryLineA());
            continue;
        }

        Interval dom;
        try
        {
            dom = crv->Domain();
        }
        catch (...)
        {
            ++skippedCoedges;
            linesA.push_back(BoundaryLineA());
            continue;
        }

        const double a = dom.Min();
        const double b = dom.Max();
        if (!(b > a) || !std::isfinite(a) || !std::isfinite(b))
        {
            ++skippedCoedges;
            linesA.push_back(BoundaryLineA());
            continue;
        }

        BoundaryLineA line;
        line.pts.reserve(static_cast<std::size_t>(nSamples));

        for (int i = 0; i < nSamples; ++i)
        {
            const double t = a + (b - a) * (static_cast<double>(i) / static_cast<double>(nSamples - 1));
            Point3D q;
            try
            {
                q = crv->CalcPoint(t);
            }
            catch (...)
            {
                continue;
            }
            line.pts.push_back(q);
            const double d = DistToParams(params, q);
            if (d > maxD) maxD = d;
        }

        if (line.pts.size() >= 2) anyMeasured = true;
        linesA.push_back(line);
    }

    return maxD;
}

//! @brief Boundary-metric stage 2: curve B (generated boundary line) and the bidirectional Hausdorff
//!        between the two boundary lines (PDF wording), in mm.
//! Returns +inf if any coedge's curve B cannot be taken (degenerate, never silently ignored).
double MeasureBoundaryCurveHausdorff(const CoedgeList& coedges, const SurfacePtr& uvSurface,
                                     const std::vector<BoundaryLineA>& linesA, int nSamples,
                                     bool& anyMeasured, bool& degenerate, int& skippedCoedges)
{
    double maxH = 0.0;
    anyMeasured = false;
    degenerate  = false;

    if (!uvSurface)
    {
        return std::numeric_limits<double>::infinity();
    }
    if (linesA.size() != coedges.size())
    {
        return std::numeric_limits<double>::infinity();
    }

    std::size_t idx = 0;
    for (const CoedgePtr& c : coedges)
    {
        const BoundaryLineA& line = linesA[idx];
        ++idx;

        if (!c || line.pts.size() < 2)
        {
            ++skippedCoedges;
            continue;
        }

        const std::vector<Point3D> lineB = SampleBoundaryLineB(c, uvSurface, nSamples);
        if (lineB.size() < 2)
        {
            degenerate = true;
            continue;
        }

        const double h = HausdorffSampled(line.pts, lineB);
        if (h > maxH) maxH = h;
        anyMeasured = true;
    }

    if (degenerate)
    {
        return std::numeric_limits<double>::infinity();
    }
    return maxH;
}

// G1: max continuity value over the face's edges; edges with != 2 faces, or returning Default, are skipped.
//! @brief G1 result: max angle deviation (rad) plus measured/skipped counts (explains 'no data').
struct G1Measurement
{
    double maxVal   = 0.0;
    int    measured = 0;
    int    skipped  = 0;
};

G1Measurement MeasureG1(const FacePtr& face)
{
    G1Measurement r;

    std::unordered_set<EdgePtr> edges;
    try
    {
        edges = face->QueryEdges();
    }
    catch (...)
    {
        return r;
    }

    for (const EdgePtr& e : edges)
    {
        if (!e)
        {
            ++r.skipped;
            continue;
        }

        FaceList fs;
        try
        {
            fs = e->QueryFaces();
        }
        catch (...)
        {
            ++r.skipped;
            continue;
        }
        if (fs.size() != 2)
        {
            ++r.skipped;
            continue;
        }

        try
        {
            const FacesContinuity fc =
                e->CalcFacesContinuity(ContinuityType::G1, 100, Toler::DefAngleTol());
            if (fc.continuityType == ContinuityType::Default)
            {
                ++r.skipped;
                continue;
            }
            if (fc.maxVal > r.maxVal) r.maxVal = fc.maxVal;
            ++r.measured;
        }
        catch (...)
        {
            ++r.skipped;
            continue;
        }
    }
    return r;
}

//! @brief Build the SDK analytic surface from fitted params; unknown type returns nullptr.
SurfacePtr MakeAnalytic(const SurfaceParams& p)
{
    switch (p.type)
    {
    case Canonical::Plane:
        return std::make_shared<Plane>(p.axisPoint, Dir3D(p.axisDir));
    case Canonical::Cylinder:
        return std::make_shared<CylinderSurface>(Ucs3D(p.axisPoint, Dir3D(p.axisDir)), p.radius);
    case Canonical::Cone:
        return std::make_shared<ConeSurface>(
            ConeSurface::CreateFromHalfAngle(Ucs3D(p.apex, Dir3D(p.axisDir)), p.halfAngle));
    case Canonical::Sphere:
        return std::make_shared<SphereSurface>(Ucs3D(p.axisPoint, Dir3D(p.axisDir)), p.radius);
    case Canonical::Torus:
        return std::make_shared<TorusSurface>(Ucs3D(p.axisPoint, Dir3D(p.axisDir)),
                                              p.radius, p.minorRadius);
    default:
        return nullptr;
    }
}

//! @brief Formatted verification numbers for the report message.
std::string AccuracyText(const AccuracyReport& a, const Thresholds& th)
{
    std::ostringstream os;
    os << "maxDist=" << std::scientific << std::setprecision(6) << a.maxDist
       << "(th=" << th.maxDist << ")"
       << " avgDist=" << a.avgDist
       << " boundaryCurveHausdorff=" << a.boundaryCurveHausdorff << "(th=" << th.boundaryHausdorff << ")"
       << " boundaryOnSurface=" << a.boundaryOnSurface
       << " g1Before=" << a.g1Before << " g1After=" << a.g1After;
    return os.str();
}

}

//! @brief Face pipeline: detect -> fit -> refine -> verify -> write back.
//! @param face input/write-back target (may be null -> kNullFace); th verification thresholds
//!        (maxDist/boundaryHausdorff, mm; that threshold gates both boundary measurements); opt sampling/refine/writeBack options.
//! @return FaceReport{decision, errorCode, message, params, accuracy}. writeBack=true replaces the
//! face's underlying surface in place (Sense preserved). faceId uses a non-atomic static sequence,
//! so concurrent calls must be serialized by the caller.
FaceReport ProcessFace(const FacePtr& face, const Thresholds& th, const PipelineOptions& opt)
{
    FaceReport report;
    report.decision = Decision::Failed;
    report.errorCode = kNullFace;
    report.message = "null face";

    if (!face)
    {
        return report;
    }

    // faceId: ClassName + per-process sequence number (ProcessFace has no index argument).
    {
        static unsigned long s_seq = 0;
        ++s_seq;
        std::ostringstream os;
        os << (face->ClassName() ? face->ClassName() : "Face") << "#" << s_seq;
        report.faceId = os.str();
    }

    try
    {
        SurfacePtr srf;
        try
        {
            srf = face->GeomSurface();
        }
        catch (const std::exception& e)
        {
            report.decision = Decision::Failed;
            report.errorCode = kInternalThrew;
            report.message = std::string("GeomSurface() threw: ") + e.what();
            return report;
        }
        if (!srf)
        {
            report.decision = Decision::Failed;
            report.errorCode = kInternalThrew;
            report.message = "GeomSurface() returned null";
            return report;
        }

        const SurfaceType st = srf->SurfaceType();
        report.sourceType = CanonicalOf(st);

        // 1) already an analytic surface -> skip.
        if (IsAnalyticSurfaceType(st))
        {
            report.decision = Decision::Rejected;
            report.errorCode = kAlreadyAnalytic;
            report.message = std::string("already analytic (") + ToString(report.sourceType)
                           + "), skipped (no conversion needed)";
            return report;
        }

        // 2) must be a spline surface.
        if (!srf->IsKindOf<BSplineSurface>())
        {
            report.decision = Decision::Rejected;
            report.errorCode = kNotBSpline;
            report.message = std::string("not a spline surface (SurfaceType=")
                           + std::to_string(static_cast<int>(st)) + "), unsupported";
            return report;
        }

        BSplineSurfacePtr bs = sg_cast<BSplineSurface>(srf);
        if (!bs)
        {
            report.decision = Decision::Failed;
            report.errorCode = kInternalThrew;
            report.message = "sg_cast<BSplineSurface> failed";
            return report;
        }

        // 2a) sample domain = CalcUVBound() intersect the spline's own Domain (fall back to the spline domain).
        const Interval du = bs->DomainU();
        const Interval dv = bs->DomainV();

        bool   fallback = false;
        double u0 = du.Min(), u1 = du.Max(), v0 = dv.Min(), v1 = dv.Max();
        try
        {
            const UVRange bound = face->CalcUVBound();
            const Interval bu = bound.IntervalU();
            const Interval bv = bound.IntervalV();
            if (IsFiniteInterval(bu) && IsFiniteInterval(bv))
            {
                u0 = (std::max)(bu.Min(), du.Min());
                u1 = (std::min)(bu.Max(), du.Max());
                v0 = (std::max)(bv.Min(), dv.Min());
                v1 = (std::min)(bv.Max(), dv.Max());
                if (!(u1 > u0) || !(v1 > v0)) fallback = true;
            }
            else
            {
                fallback = true;
            }
        }
        catch (...)
        {
            fallback = true;
        }

        std::string domainNote;
        if (fallback)
        {
            u0 = du.Min(); u1 = du.Max();
            v0 = dv.Min(); v1 = dv.Max();
            domainNote = " [CalcUVBound non-finite or disjoint from spline domain; using spline's own domain]";
        }
        if (!(u1 > u0) || !(v1 > v0))
        {
            report.decision = Decision::Failed;
            report.errorCode = kInternalThrew;
            report.message = "spline parameter domain degenerate; cannot sample" + domainNote;
            return report;
        }

        const UVRange sampleUv(Interval(u0, u1), Interval(v0, v1));

        // 3) detection.
        const std::vector<Sample> samples = SampleSurface(*bs, sampleUv, opt.nu, opt.nv);
        const DetectResult dr = Detect(*bs, sampleUv, th);

        if (dr.candidate == Canonical::Unknown)
        {
            report.decision = Decision::Rejected;
            report.errorCode = kDetectUnknown;
            std::ostringstream os;
            os << "classified Unknown, rejected: residual=" << std::scientific << std::setprecision(6) << dr.residual;
            if (!dr.notes.empty()) os << " note0=" << dr.notes.front();
            os << domainNote;
            report.message = os.str();
            return report;
        }

        std::vector<FitSample> fitSamples;
        fitSamples.reserve(samples.size());
        for (const Sample& s : samples)
        {
            FitSample f;
            f.p = s.p;
            f.n = s.n;
            fitSamples.push_back(f);
        }

        // 4) fit + refine.
        SurfaceParams params;
        std::string   reason;
        if (!FitSurface(dr.candidate, fitSamples, params, reason))
        {
            report.decision = Decision::Failed;
            report.errorCode = kFitFailed;
            report.params = params;
            report.message = std::string("FitSurface(") + ToString(dr.candidate) + ") failed: " + reason;
            return report;
        }

        std::string refineMsg;
        {
            std::string rr;
            SurfaceParams refined = params;
            if (RefineSurface(refined, fitSamples, opt.refineIter, opt.refineTol, rr))
            {
                params = refined;
            }
            else
            {
                refineMsg = " [refine did not converge (closed-form init kept):" + rr + "]";
            }
        }
        report.params = params;

        // 5) accuracy verification on a denser grid.
        const std::vector<Sample> accSamples =
            SampleSurface(*bs, sampleUv, opt.nuAcc, opt.nvAcc);
        std::vector<FitSample> accFit;
        accFit.reserve(accSamples.size());
        for (const Sample& s : accSamples)
        {
            FitSample f;
            f.p = s.p;
            f.n = s.n;
            accFit.push_back(f);
        }

        AccuracyReport acc = MeasureAccuracy(params, accFit, th);

        // 5a) boundary-metric stage 1: curve A + boundaryOnSurface.
        bool  bhMeasured = false;
        int   bhSkipped  = 0;
        std::vector<BoundaryLineA> bhLinesA;
        CoedgeList bhCoedges;
        acc.boundaryOnSurface = MeasureBoundaryOnSurface(face, params, kBoundarySamples, bhLinesA, bhCoedges, bhMeasured, bhSkipped);

        // 5b) candidate analytic surface: write-back target, and the surface used for curve B before write-back.
        SurfacePtr analytic;
        try
        {
            analytic = MakeAnalytic(params);
        }
        catch (const std::exception& e)
        {
            report.decision = Decision::Failed;
            report.errorCode = kWriteBackThrew;
            report.message = std::string("cannot build analytic surface from params (type=") + ToString(params.type) + "): " + e.what();
            return report;
        }
        catch (...)
        {
            report.decision = Decision::Failed;
            report.errorCode = kWriteBackThrew;
            report.message = std::string("cannot build analytic surface from params (type=") + ToString(params.type) + "): unknown exception";
            return report;
        }

        // 6) G1 before write-back (baseline).
        const G1Measurement g1Before = MeasureG1(face);
        acc.g1Before = g1Before.maxVal;

        std::string g1Note = " [G1 before:measured=" + std::to_string(g1Before.measured)
                           + ",skipped=" + std::to_string(g1Before.skipped) + "]";

        // 7) write back (build + attach the analytic surface). This MUST precede measuring curve B:
        // only after the swap are the coedge UV curves (pcurves) the generated boundary lines under
        // the new parameterization. On rejection the original surface is restored.
        SurfacePtr origSurface;
        bool       origSense  = true;
        bool       swapped    = false;

        if (opt.writeBack)
        {
            if (!analytic)
            {
                report.decision = Decision::Failed;
                report.errorCode = kWriteBackThrew;
                report.message = std::string("cannot build analytic surface from params (type=") + ToString(params.type) + ")";
                return report;
            }

            try
            {
                origSurface = face->GeomSurface();
                origSense   = face->Sense();
                const ModelSurfacePtr ms = TopoBuilder::MakeModelSurface(analytic, false);
                TopoBuilder::FaceAttachSurface(face, ms, face->Sense());
                swapped = true;
            }
            catch (const std::exception& e)
            {
                report.decision = Decision::Failed;
                report.errorCode = kWriteBackThrew;
                report.message = std::string("write-back failed (not replaced): ") + e.what()
                               + " " + AccuracyText(acc, th) + g1Note;
                return report;
            }
            catch (...)
            {
                report.decision = Decision::Failed;
                report.errorCode = kWriteBackThrew;
                report.message = "write-back failed (not replaced): unknown exception " + AccuracyText(acc, th) + g1Note;
                return report;
            }
        }

        // 8) boundary-metric stage 2: curve B (generated boundary line) + bidirectional Hausdorff (PDF wording).
        // Written back -> B uses the new analytic surface; not written back -> B uses the candidate
        // surface (its parameterization may differ from the old face, so this is a rough reference).
        bool        bhDegenerate = false;
        int         bhSkippedB   = 0;
        {
            const SurfacePtr uvSurface = swapped ? face->GeomSurface() : analytic;
            bool bhMeasuredB = false;
            const double h = MeasureBoundaryCurveHausdorff(bhCoedges, uvSurface, bhLinesA, kBoundarySamples,
                                                           bhMeasuredB, bhDegenerate, bhSkippedB);
            acc.boundaryCurveHausdorff = bhDegenerate ? std::numeric_limits<double>::infinity() : h;
            if (bhMeasuredB) bhMeasured = true;
        }

        std::string bhNote;
        if (!bhMeasured) bhNote += " [boundary: no coedge measured]";
        if (bhDegenerate) bhNote += " [generated boundary line: no uvCurve; boundaryCurveHausdorff=+inf]";
        if (bhSkipped > 0 || bhSkippedB > 0)
        {
            bhNote += " [skipped shared edges: A=" + std::to_string(bhSkipped) + ", B=" + std::to_string(bhSkippedB) + "]";
        }

        // 9) verdict: in-face max deviation, the two-boundary Hausdorff and boundaryOnSurface must all pass.
        const bool accepted = (acc.maxDist <= th.maxDist)
                           && (acc.boundaryCurveHausdorff <= th.boundaryHausdorff)
                           && (acc.boundaryOnSurface <= th.boundaryHausdorff);
        acc.accepted = accepted;
        report.accuracy = acc;

        if (!accepted)
        {
            if (swapped && origSurface)
            {
                try
                {
                    // cloneGeom=true: original surface still owned by the old ModelSurface, so clone identical geometry.
                    const ModelSurfacePtr back = TopoBuilder::MakeModelSurface(origSurface, true);
                    TopoBuilder::FaceAttachSurface(face, back, origSense);
                }
                catch (...)
                {
                    bhNote += " [restore original surface failed]";
                }
            }
            report.decision = Decision::Rejected;
            report.errorCode = kAccuracyReject;
            report.message = "verification failed; original spline surface kept: " + AccuracyText(acc, th) + g1Note + domainNote + bhNote;
            return report;
        }

        if (!opt.writeBack)
        {
            report.decision = Decision::Rejected;
            report.errorCode = kWriteBackOff;
            report.message = "verification passed but writeBack=false; not written back: " + AccuracyText(acc, th) + g1Note + bhNote;
            return report;
        }

        // 10) G1 after write-back (surface/parameterization changed): re-measure for an A/B comparison.
        const G1Measurement g1After = MeasureG1(face);
        report.accuracy.g1After = g1After.maxVal;
        report.decision = Decision::Converted;
        report.errorCode = kOk;
        {
            std::ostringstream os;
            os << "converted to " << ToString(params.type) << "("
               << AccuracyText(report.accuracy, th) << ")"
               << g1Note
               << " [G1 after:measured=" << g1After.measured
               << ",skipped=" << g1After.skipped << "]"
               << refineMsg << domainNote;
            report.message = os.str();
        }
        return report;
    }
    catch (const std::exception& e)
    {
        report.decision = Decision::Failed;
        report.errorCode = kInternalThrew;
        report.message = std::string("internal exception: ") + e.what();
        return report;
    }
    catch (...)
    {
        report.decision = Decision::Failed;
        report.errorCode = kInternalThrew;
        report.message = "internal unknown exception";
        return report;
    }
}

}
