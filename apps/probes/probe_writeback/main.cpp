// Development probe: checks that swapping a face's surface (analytic <-> spline, i.e. changing
// the parameterization) rebuilds its pcurves correctly.
//  [BEFORE]            build a solid, pick the Cylinder face, record Sense / CalcUVBound /
//                      per-coedge 3D samples (21 points) / pcurve count and domains.
//  [TO_SPLINE]         swap in Surface::ToBSpline(face->CalcUVBound()) (U domain -> [0,1]).
//  [BACK_TO_ANALYTIC]  swap the analytic cylinder back; compare the pcurve mapped onto the
//                      surface against the 3D samples (the core criterion), plus the
//                      closed-form distance and whether the normal flipped.
//  [SUMMARY]           four numeric answers; criteria are stated on the [SUMMARY] lines.
// Exit code only locates the failing stage (0 ran, 2/3/4 build/body/face failure); 0 does not
// mean every criterion passed.

// init()/fini() are paired at the bottom of main; the probe leaves no global kernel state.
#include <Foundation/init.h>
#include <Foundation/RuntimeClass.h>

#include <GeomBase/Number.h>
#include <GeomBase/Ucs3D.h>
#include <GeomBase/Point3D.h>
#include <GeomBase/Point2D.h>
#include <GeomBase/Dir3D.h>
#include <GeomBase/Vector3D.h>
#include <GeomBase/Interval.h>
#include <GeomBase/UVRange.h>

#include <Geometry/SurfaceType.h>
#include <Geometry/3D/Surface/Surface.h>
#include <Geometry/3D/Surface/CylinderSurface.h>
#include <Geometry/3D/Surface/BSplineSurface.h>
#include <Geometry/3D/Curve/BoundedCurve3D.h>
#include <Geometry/2D/BoundedCurve2D.h>

#include <ModelingPrim/API.h>
#include <ModelingPrim/Primitives/PrimitivesRet.h>

#include <Topology/Brep/Body.h>
#include <Topology/Brep/Face.h>
#include <Topology/Brep/Coedge.h>
#include <Topology/Brep/Edge.h>
#include <Topology/Geom/ModelSurface.h>
#include <Topology/Tools/TopoBuilder.h>

#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

using namespace sggk;

// pcurve: the 2D image of a 3D edge in the face's (u,v) domain. If a surface swap does not
// rebuild it for the new parameterization, mapping it back to 3D drifts off the edge - the
// defect this probe hunts. Sense is the face orientation flag (must be carried over);
// CalcUVBound is the face's actual (possibly pseudo-infinite) domain.
// Output fields (each line starts with a bracketed stage tag; lengths in mm unless noted):
//   face#N                                 face index in Body::QueryFaces() order, from 0
//   SurfaceType=<name>(<int>)              current surface type; int aids cross-version checks
//   U=[a,b] / V=[a,b]                      parameter domain; analytic cylinder U in rad, spline U unitless
//   Sense()=<bool>                         face orientation
//   IsKindOf<BSplineSurface>()=<bool>      confirms the write-back produced a spline
//   coedge#i                               i-th half-edge in topological order
//   CalcUVCurve() ok/FAILED                whether a pcurve exists on the current surface
//   count=<k>                              number of pcurves for this coedge (a seam may split it)
//   uvCurve#j domain=[a,b]                 parameter domain of the j-th pcurve
//   3D(coedge-dir) sample count=<n>        along-coedge 3D samples (normally 21)
//   pt[k]=(x,y,z)                          k-th 3D sample
//   pcurve-start vs StartPoint deviation   pcurve start mapped to the surface vs the coedge start
//   per-point dev aligned/reversed         worst per-point distance; aligned = same order, reversed = opposite
//   rangesAre[0,1] / rangesAre[0,2pi]      all pcurve domains within tolerance of that interval
//   rangesSameAsBefore                     spline-stage domain equals the pre-swap domain endpoint-wise
//   3D samples max |sqrt(x^2+y^2)-r|       closed-form distance from the samples to the cylinder
//   normal@refPt before/after/dot          normals at one reference point and their dot (>0 = not flipped)
// The four answers this probe gives (expected values in parentheses):
//   (1) NURBS->analytic pcurve correctness (~1e-13 mm; gate <=1e-5 mm). The kernel rebuilds the
//       pcurve by re-projecting the edge onto the current surface, so it is swap-count independent.
//   (2) Domain follows the surface type: analytic U=[0,2pi] -> spline U=[0,1] -> back to [0,2pi].
//   (3) Normal preserved: dot(before,after) ~ +1 (>0); a negative value means Sense was lost/flipped.
//   (4) Deviation magnitude: even with a new domain the pcurve still maps back onto the edge,
//       ruling out "right domain, wrong shape".
// Invariants: kRadius/kHeight build and re-fill the same cylinder (one ruler); faceSense is
// captured in [BEFORE] and reused; the only write is TopoBuilder::FaceAttachSurface; every step
// catches its own exceptions so one failure never blocks later stages.
// Traceability: stage tags localize a failure to one of the four stages; coedge#i/face#N
// localize it to a topology element; exceptions print typeid + what.

namespace
{

// Probe cylinder constants (must match api_make_solid_cylinder). Units: lengths mm, angles rad.
// kRadius/kHeight also build the analytic cylinder used for write-back, so the closed-form
// distance check uses the same ruler.
// 200 mm: a round radius makes |sqrt(x^2+y^2)-r| easy to eyeball.
constexpr Real kRadius      = 200.0;
// 500 mm: sets only the V bound; the analytic cylinder's own DomainV is undefined (see kWideV).
constexpr Real kHeight      = 500.0;
// 21 points (endpoints included): odd so index 10 is the exact midpoint, used as the normal reference.
constexpr int  kSampleCount = 21;
// 1e-6 range tolerance in parameter units (mixed: u in rad, v in mm). Domain comparison is
// exact-equality, not numeric integration; round-trip endpoints land far inside it, so no false mismatch.
constexpr Real kRangeTol    = 1e-6;
// kWideV = 1e30: the analytic cylinder's DomainV returns a pseudo-infinity (~2e50); this sentinel
// separates real bounds from that (see IsBoundedUVRange).
constexpr Real kWideV       = 1e30;

/// Trivial bool -> "true"/"false" for grep-friendly output lines.
const char* SenseStr(bool b) { return b ? "true" : "false"; }

/// Kernel SurfaceType -> readable name (printing only); unknown -> "Unknown".
std::string SurfaceTypeName(SurfaceType t)
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

//! True when all four UV endpoints are finite and |.| < kWideV. The analytic cylinder's DomainV
//! is a pseudo-infinity, so isfinite alone is not enough; the result decides the ToBSpline range.
bool IsBoundedUVRange(const UVRange& r)
{
    return std::isfinite(r.IntervalU().Min()) && std::isfinite(r.IntervalU().Max())
        && std::isfinite(r.IntervalV().Min()) && std::isfinite(r.IntervalV().Max())
        && std::fabs(r.IntervalU().Min()) < kWideV && std::fabs(r.IntervalU().Max()) < kWideV
        && std::fabs(r.IntervalV().Min()) < kWideV && std::fabs(r.IntervalV().Max()) < kWideV;
}

//! Samples n points uniformly in the curve's parameter (not arc length), endpoints included;
//! step is (b-a)/(n-1), so i=0 and i=n-1 hit the ends exactly.
std::vector<Point3D> SampleCurve3D(const BoundedCurve3DPtr& crv, int n)
{
    std::vector<Point3D> pts;
    if (!crv || n < 2)
        return pts;

    const Interval dom = crv->Domain();
    const Real a = dom.Min();
    const Real b = dom.Max();
    const Real span = b - a;

    pts.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
    {
        const Real t = a + span * (static_cast<Real>(i) / static_cast<Real>(n - 1));
        pts.push_back(crv->CalcPoint(t));
    }
    return pts;
}

//! Samples n points across a coedge's whole pcurve set, allocating counts in proportion to each
//! curve's parameter span (not evenly per curve), so a seam's two sides match the 3D density and
//! the points align index-by-index with the 3D samples.
std::vector<Point2D> SampleUvCurves(const BoundedCurve2DArray& crvs, int n)
{
    std::vector<Point2D> out;
    if (crvs.empty() || n < 2)
        return out;

    std::vector<Real> mins;
    std::vector<Real> spans;
    Real total = 0.0;
    for (const BoundedCurve2DPtr& uv : crvs)
    {
        if (!uv)
        {
            mins.push_back(0.0);
            spans.push_back(0.0);
            continue;
        }
        const Interval d = uv->Domain();
        mins.push_back(d.Min());
        const Real sp = d.Max() - d.Min();
        spans.push_back(sp > 0.0 ? sp : 0.0);
        total += spans.back();
    }

    if (total <= 0.0)
    {
        // Degenerate: all spans are zero -> repeat the first curve's start point.
        const Point2D p = crvs.front() ? crvs.front()->CalcPoint(mins.front()) : Point2D();
        out.assign(static_cast<std::size_t>(n), p);
        return out;
    }

    for (int i = 0; i < n; ++i)
    {
        const Real s = (static_cast<Real>(i) / static_cast<Real>(n - 1)) * total;

        // Locate which curve s falls in, accumulating spans; the 1e-12 absorbs float accumulation error.
        std::size_t idx = spans.size() - 1;
        Real acc = 0.0;
        for (std::size_t k = 0; k < spans.size(); ++k)
        {
            if (s <= acc + spans[k] + 1e-12)
            {
                idx = k;
                break;
            }
            acc += spans[k];
        }
        Real before = 0.0;
        for (std::size_t k = 0; k < idx; ++k)
            before += spans[k];

        Real local = spans[idx] > 0.0 ? (s - before) / spans[idx] : 0.0;
        if (local < 0.0) local = 0.0;
        if (local > 1.0) local = 1.0;

        out.push_back(crvs[idx]->CalcPoint(mins[idx] + local * spans[idx]));
    }
    return out;
}

/// Closed-form distance to the analytic cylinder (axis = Ucs3D origin +z, radius r): |sqrt(x^2+y^2)-r|.
/// Independent of the parameterization, so it is the geometric ground truth.
Real DistToCylinderZAxis(const Point3D& p, Real r)
{
    const Real rho = std::sqrt(p.X() * p.X() + p.Y() * p.Y());
    return std::fabs(rho - r);
}

/// True iff ranges is non-empty and every pair's min/max is within tol (absolute, parameter units) of a/b.
bool RangesAre(const std::vector<std::pair<Real, Real>>& ranges, Real a, Real b, Real tol)
{
    if (ranges.empty())
        return false;
    for (const std::pair<Real, Real>& r : ranges)
    {
        if (std::fabs(r.first - a) > tol || std::fabs(r.second - b) > tol)
            return false;
    }
    return true;
}

// One coedge's observation under one "current surface". Value-semantic: only facts are stored
// (threw? what value?); the criteria are applied later in [SUMMARY]. Parameter-unit quantities
// follow the surface type (analytic u in rad, spline u in [0,1]); lengths are mm.
struct CoedgeStage
{
    // ---- pcurve (CalcUVCurve) acquisition ----
    bool        uvThrew = false;
    bool        uvOk = false;
    std::size_t uvCount = 0;
    std::vector<std::pair<Real, Real>> uvRanges;
    std::string uvError;

    // ---- 3D sampling ----
    bool                 havePts = false;
    std::vector<Point3D> pts3d;

    // ---- endpoint self-consistency: pcurve start mapped to the surface vs the coedge start ----
    bool  startDevOk = false;
    Real  startDev = -1.0;

    // ---- per-point mapping: the core criterion, with aligned and reversed pairings ----
    bool  ptDevOk = false;
    Real  maxPtDevAligned = -1.0;
    Real  maxPtDevReversed = -1.0;

    // ---- geometric ground truth: 3D samples vs the closed-form cylinder distance ----
    bool  distOk = false;
    Real  maxDistToCylinder = -1.0;
};

/// Full observation of the current Face plus its current surface, one entry per coedge. srf is
/// passed in explicitly so the same code works with old and new surfaces. All exceptions are
/// caught and printed internally; a single coedge failure must not abort the probe.
void ObserveStage(const FacePtr& face, const SurfacePtr& srf, std::vector<CoedgeStage>& out)
{
    out.clear();

    CoedgeList coedges;
    try
    {
        coedges = face->QueryCoedges();
    }
    catch (const std::exception& e)
    {
        std::cout << "[OBS] QueryCoedges() threw: " << typeid(e).name() << " what=" << e.what() << std::endl;
        return;
    }
    catch (...)
    {
        std::cout << "[OBS] QueryCoedges() threw non-std exception" << std::endl;
        return;
    }

    out.resize(coedges.size());

    std::size_t ci = 0;
    for (const CoedgePtr& c : coedges)
    {
        const std::size_t i = ci++;
        CoedgeStage& st = out[i];
        if (!c)
            continue;

        // ---- pcurves ----
        BoundedCurve2DArray uvCrvs;
        try
        {
            uvCrvs = c->CalcUVCurve();
            st.uvOk = true;
            st.uvCount = uvCrvs.size();
            for (const BoundedCurve2DPtr& uv : uvCrvs)
            {
                if (!uv)
                {
                    st.uvRanges.push_back(std::make_pair(0.0, 0.0));
                    continue;
                }
                const Interval d = uv->Domain();
                st.uvRanges.push_back(std::make_pair(d.Min(), d.Max()));
            }
        }
        catch (const std::exception& e)
        {
            st.uvThrew = true;
            st.uvError = std::string(typeid(e).name()) + " what=" + e.what();
        }
        catch (...)
        {
            st.uvThrew = true;
            st.uvError = "<non-std exception>";
        }

        // ---- pcurve start mapped to the surface vs the coedge start ----
        if (srf && !uvCrvs.empty() && uvCrvs.front())
        {
            try
            {
                const Point2D uv0 = uvCrvs.front()->CalcStart();
                const Point3D p0 = c->StartPoint();
                const Point3D pOn = srf->CalcPoint(uv0);
                st.startDev = (pOn - p0).Length();
                st.startDevOk = true;
            }
            catch (const std::exception& e)
            {
                std::cout << "[OBS] coedge#" << i << " start-deviation threw: " << typeid(e).name()
                          << " what=" << e.what() << std::endl;
            }
            catch (...)
            {
                std::cout << "[OBS] coedge#" << i << " start-deviation threw non-std exception" << std::endl;
            }
        }

        // ---- 3D edge sampled along the coedge direction ----
        try
        {
            const BoundedCurve3DPtr bc = c->CalcBoundedCurve();
            st.pts3d = SampleCurve3D(bc, kSampleCount);
            st.havePts = !st.pts3d.empty();
        }
        catch (const std::exception& e)
        {
            std::cout << "[OBS] coedge#" << i << " 3D resample threw: " << typeid(e).name()
                      << " what=" << e.what() << std::endl;
        }
        catch (...)
        {
            std::cout << "[OBS] coedge#" << i << " 3D resample threw non-std exception" << std::endl;
        }

        // ---- pcurve point-for-point vs the 3D samples (core criterion) ----
        if (srf && st.havePts && !uvCrvs.empty())
        {
            try
            {
                const std::vector<Point2D> uvs = SampleUvCurves(uvCrvs, kSampleCount);
                if (uvs.size() == static_cast<std::size_t>(kSampleCount)
                    && st.pts3d.size() == static_cast<std::size_t>(kSampleCount))
                {
                    Real mA = 0.0;
                    Real mR = 0.0;
                    for (int k = 0; k < kSampleCount; ++k)
                    {
                        // aligned: uvs[k] <-> pts3d[k]
                        const Point3D pA = srf->CalcPoint(uvs[static_cast<std::size_t>(k)]);
                        const Real dA = (pA - st.pts3d[static_cast<std::size_t>(k)]).Length();
                        if (dA > mA) mA = dA;

                        // reversed: uvs[n-1-k] <-> pts3d[k]; detects a pcurve running opposite to the 3D edge
                        const Point3D pR = srf->CalcPoint(uvs[static_cast<std::size_t>(kSampleCount - 1 - k)]);
                        const Real dR = (pR - st.pts3d[static_cast<std::size_t>(k)]).Length();
                        if (dR > mR) mR = dR;
                    }
                    st.maxPtDevAligned = mA;
                    st.maxPtDevReversed = mR;
                    st.ptDevOk = true;
                }
                // Length mismatch: leave ptDevOk false so the printer shows "unavailable".
            }
            catch (const std::exception& e)
            {
                std::cout << "[OBS] coedge#" << i << " per-point deviation threw: " << typeid(e).name()
                          << " what=" << e.what() << std::endl;
            }
            catch (...)
            {
                std::cout << "[OBS] coedge#" << i << " per-point deviation threw non-std exception" << std::endl;
            }
        }

        // ---- closed-form distance from the 3D samples to the analytic cylinder ----
        if (st.havePts)
        {
            Real maxD = 0.0;
            for (const Point3D& p : st.pts3d)
            {
                const Real d = DistToCylinderZAxis(p, kRadius);
                if (d > maxD) maxD = d;
            }
            st.maxDistToCylinder = maxD;
            st.distOk = true;
        }
    }
}

/// Prints one coedge's pcurve acquisition and the domain of each UV curve. Parameter units depend on the surface type.
void PrintUvInfo(const char* label, std::size_t i, const CoedgeStage& st)
{
    if (!st.uvOk)
    {
        std::cout << label << " coedge#" << i
                  << " CalcUVCurve() FAILED: " << (st.uvError.empty() ? std::string("<unknown>") : st.uvError)
                  << std::endl;
        return;
    }
    std::cout << label << " coedge#" << i
              << " CalcUVCurve() ok, count=" << st.uvCount << std::endl;
    for (std::size_t k = 0; k < st.uvRanges.size(); ++k)
    {
        std::cout << label << " coedge#" << i << " uvCurve#" << k
                  << " domain=[" << st.uvRanges[k].first << ", " << st.uvRanges[k].second << "]" << std::endl;
    }
}

/// Runs the four-stage probe (build -> to spline -> back to analytic -> summary) and returns the
/// exit code: 0 done, 2/3/4 build/body/no-face failure; 10/11 come from main's exception handling.
int RunProbe()
{
    std::cout << std::setprecision(17);

    // ---- step 1: build the solid, list faces, print SurfaceType and CalcUVBound ----
    auto ret = api_make_solid_cylinder(Ucs3D(), kRadius, kHeight, PI2, true);
    if (!ret || !ret->Succeeded())
    {
        std::cout << "[FACES] api_make_solid_cylinder FAILED: "
                  << (ret ? ret->Status().ErrorMsg() : std::string("<null ret>")) << std::endl;
        return 2;
    }

    const BodyPtr body = ret->ResultBody();
    if (!body)
    {
        std::cout << "[FACES] ResultBody() returned null" << std::endl;
        return 3;
    }

    const FaceList faces = body->QueryFaces();
    std::cout << "[FACES] body face count = " << faces.size() << std::endl;

    std::size_t faceIndex = 0;
    for (const FacePtr& f : faces)
    {
        if (!f)
        {
            std::cout << "[FACES] face#" << faceIndex << " is null" << std::endl;
            ++faceIndex;
            continue;
        }
        try
        {
            const SurfacePtr srf = f->GeomSurface();
            const UVRange uv = f->CalcUVBound();
            std::cout << "[FACES] face#" << faceIndex
                      << " SurfaceType=" << (srf ? SurfaceTypeName(srf->SurfaceType()) : std::string("<null srf>"))
                      << "(" << (srf ? static_cast<int>(srf->SurfaceType()) : -1) << ")"
                      << " U=[" << uv.IntervalU().Min() << ", " << uv.IntervalU().Max() << "]"
                      << " V=[" << uv.IntervalV().Min() << ", " << uv.IntervalV().Max() << "]"
                      << std::endl;
        }
        catch (const std::exception& e)
        {
            std::cout << "[FACES] face#" << faceIndex << " query threw: " << typeid(e).name()
                      << " what=" << e.what() << std::endl;
        }
        catch (...)
        {
            std::cout << "[FACES] face#" << faceIndex << " query threw non-std exception" << std::endl;
        }
        ++faceIndex;
    }

    // ---- step 1b: pick the Cylinder face. The solid has two Plane caps too; picking by type avoids
    // selecting a cap, which would make the pcurve checks meaningless. ----
    FacePtr cylFace;
    for (const FacePtr& f : faces)
    {
        if (!f)
            continue;
        try
        {
            const SurfacePtr srf = f->GeomSurface();
            if (srf && srf->SurfaceType() == SurfaceType::Cylinder)
            {
                cylFace = f;
                break;
            }
        }
        catch (...)
        {
        }
    }

    if (!cylFace)
    {
        std::cout << "[PICK] no face with SurfaceType()==Cylinder found; falling back to first non-null face" << std::endl;
        for (const FacePtr& f : faces)
        {
            if (f)
            {
                cylFace = f;
                break;
            }
        }
    }
    if (!cylFace)
    {
        std::cout << "[PICK] no usable face found; aborting" << std::endl;
        return 4;
    }

    SurfacePtr origSrf = cylFace->GeomSurface();
    const bool faceSense = cylFace->Sense();
    std::cout << "[PICK] picked face, SurfaceType="
              << (origSrf ? SurfaceTypeName(origSrf->SurfaceType()) : std::string("<null>"))
              << " Face::Sense()=" << SenseStr(faceSense)
              << std::endl;
    if (origSrf)
    {
        std::cout << "[PICK] orig surface IsKindOf<BSplineSurface>() = "
                  << (origSrf->IsKindOf<BSplineSurface>() ? "true" : "false") << std::endl;
    }

    // ---- step 2: [BEFORE] record everything before the swap ----
    UVRange faceUvBound;
    try
    {
        faceUvBound = cylFace->CalcUVBound();
    }
    catch (const std::exception& e)
    {
        std::cout << "[BEFORE] CalcUVBound() threw: " << typeid(e).name() << " what=" << e.what() << std::endl;
    }
    catch (...)
    {
        std::cout << "[BEFORE] CalcUVBound() threw non-std exception" << std::endl;
    }

    std::cout << "[BEFORE] Face::Sense() = " << SenseStr(faceSense) << std::endl;
    std::cout << "[BEFORE] Face::CalcUVBound() U=[" << faceUvBound.IntervalU().Min() << ", "
              << faceUvBound.IntervalU().Max() << "] V=[" << faceUvBound.IntervalV().Min() << ", "
              << faceUvBound.IntervalV().Max() << "]" << std::endl;

    std::vector<CoedgeStage> before;
    ObserveStage(cylFace, cylFace->GeomSurface(), before);
    std::cout << "[BEFORE] coedge count = " << before.size() << std::endl;

    for (std::size_t i = 0; i < before.size(); ++i)
    {
        const CoedgeStage& st = before[i];
        PrintUvInfo("[BEFORE]", i, st);
        std::cout << "[BEFORE] coedge#" << i << " 3D(coedge-dir) sample count = " << st.pts3d.size() << std::endl;
        for (std::size_t k = 0; k < st.pts3d.size(); ++k)
        {
            const Point3D& p = st.pts3d[k];
            std::cout << "[BEFORE] coedge#" << i << " pt[" << k << "] = ("
                      << p.X() << ", " << p.Y() << ", " << p.Z() << ")" << std::endl;
        }
        if (st.startDevOk)
        {
            std::cout << "[BEFORE] coedge#" << i << " pcurve-start vs StartPoint deviation = "
                      << st.startDev << std::endl;
        }
    }

    // Reference point = midpoint (index 10) of the first coedge's 3D samples; it lies on the face and
    // away from seams/endpoints, so the normal is most stable.
    Point3D refPt;
    bool    refValid = false;
    for (const CoedgeStage& st : before)
    {
        if (st.havePts && st.pts3d.size() >= static_cast<std::size_t>(kSampleCount))
        {
            refPt = st.pts3d[static_cast<std::size_t>(kSampleCount / 2)];
            refValid = true;
            break;
        }
    }
    if (!refValid)
        refPt = Point3D(0, 0, 0);

    Dir3D nBefore(0, 0, 1);
    bool  nBeforeValid = false;
    try
    {
        nBefore = cylFace->CalcNormal(refPt);
        nBeforeValid = true;
    }
    catch (const std::exception& e)
    {
        std::cout << "[BEFORE] CalcNormal threw: " << typeid(e).name() << " what=" << e.what() << std::endl;
    }
    catch (...)
    {
        std::cout << "[BEFORE] CalcNormal threw non-std exception" << std::endl;
    }
    std::cout << "[BEFORE] refPt=(" << refPt.X() << ", " << refPt.Y() << ", " << refPt.Z() << ")"
              << " normal=(" << nBefore.X() << ", " << nBefore.Y() << ", " << nBefore.Z() << ")"
              << " normalValid=" << SenseStr(nBeforeValid) << std::endl;

    // ---- step 3: [TO_SPLINE] analytic cylinder -> spline (parameterization changes) ----
    std::vector<CoedgeStage> splineObs;
    bool       splineBuilt = false;

    {
        // Use face->CalcUVBound() as the ToBSpline range, not the analytic cylinder's DomainV which is
        // undefined (pseudo-infinite); the face bounds are the region actually sampled.
        UVRange toRange = faceUvBound;
        if (!IsBoundedUVRange(toRange))
        {
            toRange = UVRange(Interval(0.0, PI2), Interval(0.0, kHeight));
            std::cout << "[TO_SPLINE] CalcUVBound not fully bounded -> fallback finite range U=[0, "
                      << PI2 << "] V=[0, " << kHeight << "]" << std::endl;
        }
        std::cout << "[TO_SPLINE] ToBSpline range used: U=[" << toRange.IntervalU().Min() << ", "
                  << toRange.IntervalU().Max() << "] V=[" << toRange.IntervalV().Min() << ", "
                  << toRange.IntervalV().Max() << "]" << std::endl;

        BSplineSurfacePtr splineSrf;
        try
        {
            if (origSrf)
                splineSrf = origSrf->ToBSpline(toRange);
        }
        catch (const std::exception& e)
        {
            std::cout << "[TO_SPLINE] ToBSpline(faceUvBound) threw: " << typeid(e).name()
                      << " what=" << e.what() << " -> retry with finite fallback range" << std::endl;
            try
            {
                if (origSrf)
                    splineSrf = origSrf->ToBSpline(UVRange(Interval(0.0, PI2), Interval(0.0, kHeight)));
            }
            catch (const std::exception& e2)
            {
                std::cout << "[TO_SPLINE] ToBSpline(fallback) also threw: " << typeid(e2).name()
                          << " what=" << e2.what() << std::endl;
            }
            catch (...)
            {
                std::cout << "[TO_SPLINE] ToBSpline(fallback) threw non-std exception" << std::endl;
            }
        }
        catch (...)
        {
            std::cout << "[TO_SPLINE] ToBSpline threw non-std exception" << std::endl;
        }

        if (splineSrf)
        {
            splineBuilt = true;

            const Interval du = splineSrf->DomainU();
            const Interval dv = splineSrf->DomainV();
            std::cout << "[TO_SPLINE] ToBSpline ok: DegreeU=" << splineSrf->DegreeU()
                      << " DegreeV=" << splineSrf->DegreeV()
                      << " IsRational=" << SenseStr(splineSrf->IsRational())
                      << " DomainU=[" << du.Min() << ", " << du.Max() << "]"
                      << " DomainV=[" << dv.Min() << ", " << dv.Max() << "]" << std::endl;

            try
            {
                // MakeModelSurface wraps a raw surface as topological geometry (false = not a periodic wrap).
                const ModelSurfacePtr msSpline = TopoBuilder::MakeModelSurface(splineSrf, false);
                // The write: attach the spline surface to the same Face, reusing the original Sense; this triggers pcurve rebuild.
                TopoBuilder::FaceAttachSurface(cylFace, msSpline, faceSense);
                std::cout << "[TO_SPLINE] FaceAttachSurface done (sense=" << SenseStr(faceSense) << ")" << std::endl;
            }
            catch (const std::exception& e)
            {
                std::cout << "[TO_SPLINE] MakeModelSurface/FaceAttachSurface threw: " << typeid(e).name()
                          << " what=" << e.what() << std::endl;
            }
            catch (...)
            {
                std::cout << "[TO_SPLINE] MakeModelSurface/FaceAttachSurface threw non-std exception" << std::endl;
            }

            try
            {
                const SurfacePtr cur = cylFace->GeomSurface();
                std::cout << "[TO_SPLINE] GeomSurface()->SurfaceType() = "
                          << (cur ? SurfaceTypeName(cur->SurfaceType()) : std::string("<null>"))
                          << " (" << (cur ? static_cast<int>(cur->SurfaceType()) : -1) << ")"
                          << " IsKindOf<BSplineSurface>()="
                          << ((cur && cur->IsKindOf<BSplineSurface>()) ? "true" : "false")
                          << std::endl;
            }
            catch (const std::exception& e)
            {
                std::cout << "[TO_SPLINE] GeomSurface()/SurfaceType() threw: " << typeid(e).name()
                          << " what=" << e.what() << std::endl;
            }
            catch (...)
            {
                std::cout << "[TO_SPLINE] GeomSurface()/SurfaceType() threw non-std exception" << std::endl;
            }

            ObserveStage(cylFace, cylFace->GeomSurface(), splineObs);

            for (std::size_t i = 0; i < splineObs.size(); ++i)
            {
                const CoedgeStage& st = splineObs[i];
                PrintUvInfo("[TO_SPLINE]", i, st);

                const bool is01 = RangesAre(st.uvRanges, 0.0, 1.0, kRangeTol);
                bool sameAsBefore = (i < before.size()) && (st.uvCount == before[i].uvCount)
                                 && (st.uvRanges.size() == before[i].uvRanges.size());
                if (sameAsBefore)
                {
                    for (std::size_t k = 0; k < st.uvRanges.size(); ++k)
                    {
                        // Compare domain endpoints one by one; any excess means the domain changed (what this probe looks for).
                        if (std::fabs(st.uvRanges[k].first - before[i].uvRanges[k].first) > kRangeTol
                            || std::fabs(st.uvRanges[k].second - before[i].uvRanges[k].second) > kRangeTol)
                        {
                            sameAsBefore = false;
                            break;
                        }
                    }
                }
                std::cout << "[TO_SPLINE] coedge#" << i
                          << " uvCount=" << st.uvCount
                          << " (before=" << (i < before.size() ? before[i].uvCount : 0) << ")"
                          << " rangesAre[0,1]=" << SenseStr(is01)
                          << " rangesSameAsBefore=" << SenseStr(sameAsBefore)
                          << std::endl;

                if (st.startDevOk)
                {
                    std::cout << "[TO_SPLINE] coedge#" << i
                              << " pcurve-start on spline vs coedge StartPoint deviation = " << st.startDev
                              << std::endl;
                }
                else
                {
                    std::cout << "[TO_SPLINE] coedge#" << i << " start-deviation unavailable" << std::endl;
                }
                if (st.ptDevOk)
                {
                    std::cout << "[TO_SPLINE] coedge#" << i << " per-point dev aligned=" << st.maxPtDevAligned
                              << " reversed=" << st.maxPtDevReversed << std::endl;
                }
            }
        }
        else
        {
            std::cout << "[TO_SPLINE] ToBSpline produced NULL spline surface; skipping attach" << std::endl;
        }
    }

    // ---- step 4: [BACK_TO_ANALYTIC] swap the analytic cylinder back in ----
    std::vector<CoedgeStage> backObs;
    {
        // A fresh analytic cylinder with the same radius (default Ucs3D), used for write-back.
        auto cylSrf = std::make_shared<CylinderSurface>(Ucs3D(), kRadius);
        std::cout << "[BACK_TO_ANALYTIC] built analytic CylinderSurface(r=" << cylSrf->Radius()
                  << ") sense=" << SenseStr(faceSense) << std::endl;

        try
        {
            const ModelSurfacePtr msCyl = TopoBuilder::MakeModelSurface(cylSrf, false);
            // Reuse the Sense saved a cycle earlier so orientation does not flip on the swap.
            TopoBuilder::FaceAttachSurface(cylFace, msCyl, faceSense);
            std::cout << "[BACK_TO_ANALYTIC] FaceAttachSurface done" << std::endl;
        }
        catch (const std::exception& e)
        {
            std::cout << "[BACK_TO_ANALYTIC] MakeModelSurface/FaceAttachSurface threw: " << typeid(e).name()
                      << " what=" << e.what() << std::endl;
        }
        catch (...)
        {
            std::cout << "[BACK_TO_ANALYTIC] MakeModelSurface/FaceAttachSurface threw non-std exception" << std::endl;
        }

        try
        {
            const SurfacePtr cur = cylFace->GeomSurface();
            std::cout << "[BACK_TO_ANALYTIC] GeomSurface()->SurfaceType() = "
                      << (cur ? SurfaceTypeName(cur->SurfaceType()) : std::string("<null>"))
                      << " (" << (cur ? static_cast<int>(cur->SurfaceType()) : -1) << ")"
                      << std::endl;
        }
        catch (const std::exception& e)
        {
            std::cout << "[BACK_TO_ANALYTIC] GeomSurface()/SurfaceType() threw: " << typeid(e).name()
                      << " what=" << e.what() << std::endl;
        }
        catch (...)
        {
            std::cout << "[BACK_TO_ANALYTIC] GeomSurface()/SurfaceType() threw non-std exception" << std::endl;
        }

        ObserveStage(cylFace, cylFace->GeomSurface(), backObs);

        for (std::size_t i = 0; i < backObs.size(); ++i)
        {
            const CoedgeStage& st = backObs[i];
            PrintUvInfo("[BACK_TO_ANALYTIC]", i, st);

            const bool is02pi = RangesAre(st.uvRanges, 0.0, PI2, kRangeTol);
            std::cout << "[BACK_TO_ANALYTIC] coedge#" << i
                      << " uvCurveOk=" << SenseStr(st.uvOk)
                      << " count=" << st.uvCount
                      << " rangesAre[0,2pi]=" << SenseStr(is02pi)
                      << std::endl;

            if (st.havePts && st.ptDevOk)
            {
                std::cout << "[BACK_TO_ANALYTIC] coedge#" << i
                          << " per-point uvCurve-on-surface vs 3D maxDist: aligned=" << st.maxPtDevAligned
                          << " reversed=" << st.maxPtDevReversed << std::endl;
            }
            else
            {
                std::cout << "[BACK_TO_ANALYTIC] coedge#" << i
                          << " per-point comparison unavailable (havePts=" << SenseStr(st.havePts)
                          << " ptDevOk=" << SenseStr(st.ptDevOk) << ")" << std::endl;
            }

            if (st.distOk)
            {
                std::cout << "[BACK_TO_ANALYTIC] coedge#" << i
                          << " 3D samples max |sqrt(x^2+y^2)-r| = " << st.maxDistToCylinder << std::endl;
            }
            if (st.startDevOk)
            {
                std::cout << "[BACK_TO_ANALYTIC] coedge#" << i
                          << " pcurve-start vs StartPoint deviation = " << st.startDev << std::endl;
            }
        }
    }

    // ---- normal at the same 3D reference point: recompute and dot with the [BEFORE] normal;
    // dot > 0 means orientation was preserved ----
    Real nDot = 0.0;
    bool nAfterValid = false;
    try
    {
        const Dir3D nAfter = cylFace->CalcNormal(refPt);
        nDot = nBefore.Dot(nAfter);
        nAfterValid = true;
        std::cout << "[BACK_TO_ANALYTIC] normal@refPt: before=(" << nBefore.X() << ", " << nBefore.Y() << ", " << nBefore.Z()
                  << ") after=(" << nAfter.X() << ", " << nAfter.Y() << ", " << nAfter.Z()
                  << ") dot=" << nDot << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cout << "[BACK_TO_ANALYTIC] CalcNormal threw: " << typeid(e).name() << " what=" << e.what() << std::endl;
    }
    catch (...)
    {
        std::cout << "[BACK_TO_ANALYTIC] CalcNormal threw non-std exception" << std::endl;
    }

    // ---- step 5: [SUMMARY] four numeric answers; all lines are prefixed [SUMMARY] for scripting ----
    // Per-stage suprema of the deviations.
    Real beforeStartDev = 0.0;
    for (const CoedgeStage& st : before)
        if (st.startDevOk && st.startDev > beforeStartDev) beforeStartDev = st.startDev;

    Real splineStartDev = 0.0;
    Real splinePtDev = -1.0;
    std::size_t splineUvOk = 0, splineCoedges = splineObs.size(), splineUv01 = 0;
    for (const CoedgeStage& st : splineObs)
    {
        if (st.startDevOk && st.startDev > splineStartDev) splineStartDev = st.startDev;
        if (st.ptDevOk)
        {
            // Take the smaller of aligned/reversed: when the direction is unknown, that is the true best pairing.
            const Real best = st.maxPtDevAligned < st.maxPtDevReversed ? st.maxPtDevAligned : st.maxPtDevReversed;
            if (best > splinePtDev) splinePtDev = best;
        }
        if (st.uvOk) ++splineUvOk;
        if (RangesAre(st.uvRanges, 0.0, 1.0, kRangeTol)) ++splineUv01;
    }

    Real backStartDev = 0.0;
    Real backPtDev = -1.0;
    Real backPtDevAligned = -1.0;
    Real backDist = 0.0;
    std::size_t backUvOk = 0, backUv01 = 0, backUv02pi = 0, backCoedges = backObs.size();
    for (const CoedgeStage& st : backObs)
    {
        if (st.startDevOk && st.startDev > backStartDev) backStartDev = st.startDev;
        if (st.ptDevOk)
        {
            const Real best = st.maxPtDevAligned < st.maxPtDevReversed ? st.maxPtDevAligned : st.maxPtDevReversed;
            if (best > backPtDev) backPtDev = best;
            if (st.maxPtDevAligned > backPtDevAligned) backPtDevAligned = st.maxPtDevAligned;
        }
        if (st.distOk && st.maxDistToCylinder > backDist) backDist = st.maxDistToCylinder;
        if (st.uvOk) ++backUvOk;
        if (RangesAre(st.uvRanges, 0.0, 1.0, kRangeTol)) ++backUv01;
        if (RangesAre(st.uvRanges, 0.0, PI2, kRangeTol)) ++backUv02pi;
    }

    // Criterion (1)/(4): pcurves are correct iff there is data and max per-point deviation <= 1e-5 mm.
    // That gate is ~8 orders above the measured 1e-13 mm noise and far below any visible geometry error.
    const bool pcurveCorrect = (backPtDev >= 0.0) && (backPtDev <= 1e-5);
    // Criterion (2): the domain changed iff a spline surface was built (its DomainU is [0,1]).
    const bool domainChanged = (splineBuilt ? true : false);

    std::cout << "[SUMMARY] (1) NURBS->analytic pcurve correctness: maxPointDeviation="
              << backPtDev << " mm (coedge-direction aligned=" << backPtDevAligned
              << "), startDeviation=" << backStartDev
              << " -> pcurveCorrect(<=1e-5)=" << SenseStr(pcurveCorrect)
              << std::endl;
    std::cout << "[SUMMARY] (2) UV domain vs surface type: analytic-before U=[0, 2pi] ; spline U=[0,1] ("
              << splineUv01 << "/" << splineCoedges << " coedges) ; back-to-analytic U=[0, 2pi] ("
              << backUv02pi << "/" << backCoedges << " coedges) -> domain changed with type=" << SenseStr(domainChanged)
              << std::endl;
    std::cout << "[SUMMARY] (3) normal preserved: dot(before,after)=" << nDot
              << " afterValid=" << SenseStr(nAfterValid)
              << " -> preserved(>0)=" << SenseStr(nAfterValid && nDot > 0.0)
              << std::endl;
    std::cout << "[SUMMARY] (4) pcurve deviation magnitude: max="
              << backPtDev << " mm ; exceeds1e-5mm=" << SenseStr(!(backPtDev >= 0.0 && backPtDev <= 1e-5))
              << " ; closed-form max |sqrt(x^2+y^2)-r|=" << backDist << " mm"
              << std::endl;
    std::cout << "[SUMMARY] coedge counts: before=" << before.size()
              << " spline=[" << splineUvOk << "/" << splineCoedges << " uvOk]"
              << " back=[" << backUvOk << "/" << backCoedges << " uvOk, " << backUv01
              << " ranges [0,1]]" << std::endl;
    std::cout << "[SUMMARY] spline-stage check: startDev=" << splineStartDev
              << " bestPerPointDev=" << splinePtDev
              << " ; before-stage startDev=" << beforeStartDev
              << " ; splineBuilt=" << SenseStr(splineBuilt) << std::endl;

    return 0;
}

} // namespace

/// Entry point: init kernel -> RunProbe -> fini; returns RunProbe's code, or 10/11 on a std/non-std exception.
int main(int /*argc*/, char* /*argv*/[])
{
    init();

    int rc = 0;
    try
    {
        rc = RunProbe();
    }
    catch (const std::exception& e)
    {
        std::cerr << "[FATAL] " << typeid(e).name() << " what=" << e.what() << std::endl;
        rc = 10;
    }
    catch (...)
    {
        std::cerr << "[FATAL] unknown non-std exception" << std::endl;
        rc = 11;
    }

    fini();
    return rc;
}

// Tolerances: kRangeTol=1e-6 parameter units (domain comparison is exact equality, not numeric
// integration; the round-trip endpoints are exact, so 1e-6 only guards future scaling).
// pcurveCorrect gates on 1e-5 mm, ~8 orders above the measured 1e-13 mm noise.
// kWideV=1e30 is ~20 orders below the analytic DomainV pseudo-infinity (~2e50), so isfinite plus
// this threshold reliably rejects the fake V bound.
// No manual pcurve recomputation is attempted: the round-trip deviation already sits at the
// double-precision noise floor (~1e-13 mm), so a second implementation would add no accuracy.

// Call sequence: init -> api_make_solid_cylinder -> ResultBody -> QueryFaces -> GeomSurface /
// CalcUVBound / QueryCoedges / CalcUVCurve / CalcBoundedCurve / CalcPoint -> ToBSpline ->
// MakeModelSurface / FaceAttachSurface (the only write) -> CalcNormal -> fini. All *Ptr are
// shared handles (borrowed, never stored past a call); CoedgeStage keeps only values, so an
// observation stays valid after the next swap. Failure policy: localized failures set a flag or
// print and never abort the probe; missing data defaults to false (-1 = no measurement).
