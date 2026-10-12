// Development probe: measures whether analytic quadrics (B) are faster than the same
// geometry expressed as BSpline faces (A) for surface/surface intersection, point
// projection and boolean, to test the claim that B is >=50% faster. Each op: warmup 5,
// repeat 50, median microseconds, speedup = A/B; timing wraps only the kernel call.
// B keeps the built analytic faces; A replaces each face with Surface::ToBSpline(uvRange).
// Projection note: GeomProject::PntSrfProject is an unimplemented throwing stub for normal
// projection in SGK 1.7.11, so the projection op uses Surface::CalcNearestPnt instead (see
// the projection-kernel note below).

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <functional>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#include <Foundation/init.h>
#include <Foundation/Common.h>
#include <Foundation/BaseException.h>
#include <Foundation/RuntimeClass.h>

#include <GeomBase/Number.h>
#include <GeomBase/Interval.h>
#include <GeomBase/UVRange.h>
#include <GeomBase/Point2D.h>
#include <GeomBase/Point3D.h>
#include <GeomBase/Dir3D.h>
#include <GeomBase/Ucs3D.h>
#include <GeomBase/Toler.h>

#include <Geometry/3D/Surface/Surface.h>
#include <Geometry/3D/Surface/CylinderSurface.h>
#include <Geometry/3D/Surface/Plane.h>
#include <Geometry/3D/Surface/BSplineSurface.h>
#include <Geometry/Util/SurfaceNearestPntUtil.h>

#include <ModelingPrim/API.h>
#include <Boolean/API.h>

#include <Topology/Entity.h>
#include <Topology/Brep/Topology.h>
#include <Topology/Brep/Body.h>
#include <Topology/Brep/Face.h>
#include <Topology/Geom/ModelSurface.h>
#include <Topology/Tools/TopoBuilder.h>

#include <GeomInt/GeomInt.h>
#include <GeomProject/GeomProject.h>

using namespace sggk;

namespace
{
    constexpr int kWarmup = 5; // kWarmup: warm-up runs
    constexpr int kRepeat = 50; // kRepeat: timed runs

    // timing

    // Warm up kWarmup times, then time kRepeat runs and return the median in microseconds.
    // The timed span covers fn() only; I/O, setup and printing stay outside it.
    double MedianUs(const std::function<void()>& fn)
    {
        for (int i = 0; i < kWarmup; ++i)
            fn();

        std::vector<double> samples;
        samples.reserve(kRepeat);
        for (int i = 0; i < kRepeat; ++i)
        {
            const auto t0 = std::chrono::steady_clock::now();
            fn();
            const auto t1 = std::chrono::steady_clock::now();
            samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
        }
        std::sort(samples.begin(), samples.end());
        return samples[samples.size() / 2];
    }

    // Prints an exception's typeid and what(); for a kernel BaseException it also prints
    // name()/code and the raise site. The [DIAG] prefix separates it from the kernel's own stdout.
    void DiagException(const char* where, const std::exception& e)
    {
        std::printf("[DIAG] %s: type=%s what=%s", where, typeid(e).name(), e.what());
        if (const BaseException* be = dynamic_cast<const BaseException*>(&e))
            std::printf(" name=%s code=0x%08X raised_at=%s(%d) in=%s",
                        be->name(), static_cast<unsigned int>(be->getCode()), be->getFileName(),
                        be->getLineNumber(), be->getFuncName());
        std::printf("\n");
    }

    void DiagUnknownException(const char* where)
    {
        std::printf("[DIAG] %s: unknown non-std::exception\n", where);
    }

    // Counts a body's faces and how many are analytic vs spline.
    struct GeomInfo
    {
        size_t faces = 0;
        size_t analytic = 0;
        size_t spline = 0;
    };

    GeomInfo Inspect(const BodyPtr& body)
    {
        GeomInfo info;
        if (!body)
            return info;
        for (const auto& face : body->QueryFaces())
        {
            ++info.faces;
            SurfacePtr srf = face->GeomSurface();
            if (!srf)
                continue;
            if (sg_cast<BSplineSurface>(srf))
                ++info.spline;
            else
                ++info.analytic;
        }
        return info;
    }

    // Clones src and replaces every face's surface with its BSpline form; converted/total report
    // the replaced and total face counts. face->Sense() is carried over unchanged.
    BodyPtr SplineizeBody(const BodyPtr& src, int& converted, int& total)
    {
        converted = 0;
        total = 0;
        if (!src)
            return nullptr;

        BodyPtr clone = sg_cast<Body>(src->Clone());
        if (!clone)
            return nullptr;

        for (const auto& face : clone->QueryFaces())
        {
            ++total;
            SurfacePtr geom = face->GeomSurface();
            if (!geom)
                continue;

            const UVRange uv = face->CalcUVBound();
            BSplineSurfacePtr bs;
            try
            {
                bs = geom->ToBSpline(uv);
            }
            catch (const std::exception& e)
            {
                DiagException("SplineizeBody ToBSpline", e);
                continue;
            }
            catch (...)
            {
                DiagUnknownException("SplineizeBody ToBSpline");
                continue;
            }
            if (!bs)
                continue;

            ModelSurfacePtr modelSrf = TopoBuilder::MakeModelSurface(bs, true);
            if (!modelSrf)
                continue;

            // Face::Sense(): whether the face and its underlying surface share orientation; preserved here.
            TopoBuilder::FaceAttachSurface(face, modelSrf, face->Sense());
            ++converted;
        }
        return clone;
    }

    struct Case
    {
        std::string op;
        double aUs = 0.0; // A (spline)
        double bUs = 0.0; // B (analytic)
        double speedup = 0.0;
    };

    std::vector<Case> g_cases;

    // Emits one op's timings (A=spline, B=analytic, speedup=A/B); on failure prints N/A, never fabricates.
    void Emit(const std::string& op, bool ok, double aUs, double bUs, const std::string& reason)
    {
        if (!ok)
        {
            std::printf("[TIMING] op=%s N/A reason=%s\n", op.c_str(), reason.c_str());
            std::fprintf(stderr, "%s,NA,NA,NA\n", op.c_str());
            return;
        }
        const double speedup = (bUs > 0.0) ? (aUs / bUs) : 0.0;
        std::printf("[TIMING] op=%s A=%.3f us B=%.3f us speedup(A/B)=%.4f\n",
                    op.c_str(), aUs, bUs, speedup);
        std::fprintf(stderr, "%s,%.3f,%.3f,%.4f\n", op.c_str(), aUs, bUs, speedup);
        g_cases.push_back(Case{op, aUs, bUs, speedup});
    }

    // Surface/surface intersection with a 1e-5 absolute tolerance.
    void RunSrfSrfInt(const Surface& s0, const Surface& s1, const UVRange& u0, const UVRange& u1)
    {
        SrfSrfIntOpts opts;
        opts.tol = Toler(1e-5);
        auto ret = GeomInt::SrfSrfInt(s0, s1, u0, u1, opts);
        (void)ret;
    }

    // Projection kernel: GeomProject::PntSrfProject (the documented API, GeomProject/GeomProject.h:41)
    // is an unimplemented throwing stub for normal projection in SGK 1.7.11 (it raises
    // ProjectionNotImplemented), and its parallel mode likewise rejects BSplineSurface. Use
    // Surface::CalcNearestPnt, which is the same geometric quantity (nearest point along the normal)
    // and is implemented for both analytic and spline faces. SurfaceNearestPntUtil::CalcNearestPoint
    // is the alternative but needs an explicit UV seed.

    enum class ProjKernel
    {
        LegacyProject, // documented API; normal projection is a throwing stub in 1.7.11
        NearestPnt, // the kernel actually used here
        NearestUtil, // alternative, needs an explicit UV seed
    };

    const char* ProjKernelName(ProjKernel k)
    {
        switch (k)
        {
        case ProjKernel::LegacyProject: return "GeomProject::PntSrfProject";
        case ProjKernel::NearestPnt:    return "Surface::CalcNearestPnt";
        case ProjKernel::NearestUtil:   return "SurfaceNearestPntUtil::CalcNearestPoint";
        }
        return "none";
    }

    // Runs one projection kernel over all points. Empty result -> false (no throw); a throwing
    // kernel propagates and is diagnosed by the caller.
    bool ApplyProjKernel(ProjKernel k, const Point3D& pnt, const Surface& srf, const UVRange& uvHint, double& dist)
    {
        switch (k)
        {
        case ProjKernel::LegacyProject:
        {
            // Default options = normal projection.
            const std::vector<PntSrfProjectRet> ret = GeomProject::PntSrfProject(pnt, srf, PntSrfProjectOpts());
            if (ret.empty())
                return false;
            dist = pnt.DistanceTo(ret.front().projPoint);
            return true;
        }
        case ProjKernel::NearestPnt:
        {
            Point2D uv;
            const Point3D nearest = srf.CalcNearestPnt(pnt, uv);
            dist = pnt.DistanceTo(nearest);
            return true;
        }
        case ProjKernel::NearestUtil:
        {
            Point2D uv;
            const Point2D initUV(uvHint.IntervalU().Middle(), uvHint.IntervalV().Middle());
            const SrfNearestPntRet ret =
                SurfaceNearestPntUtil::CalcNearestPoint(pnt, srf, initUV, uv, Toler(1e-6));
            if (ret.type == SrfNearestPntType::IterFailure)
                return false;
            dist = pnt.DistanceTo(ret.pnt);
            return true;
        }
        }
        return false;
    }

    // Dry-runs a kernel over the whole point set; prints [DIAG] on empty result or exception.
    bool ProbeProjKernel(ProjKernel k, const char* label, const std::vector<Point3D>& pnts,
                         const Surface& srf, const UVRange& uvHint, std::string& why)
    {
        char ctx[192];
        std::snprintf(ctx, sizeof(ctx), "%s kernel=%s", label, ProjKernelName(k));
        try
        {
            double dist = 0.0;
            for (const Point3D& p : pnts)
            {
                if (!ApplyProjKernel(k, p, srf, uvHint, dist))
                {
                    why = "empty projection result";
                    std::printf("[DIAG] %s: %s\n", ctx, why.c_str());
                    return false;
                }
            }
            std::printf("[DIAG] %s: ok, points=%zu\n", ctx, pnts.size());
            return true;
        }
        catch (const std::exception& e)
        {
            why = "threw (see [DIAG])";
            DiagException(ctx, e);
            return false;
        }
        catch (...)
        {
            why = "unknown exception (see [DIAG])";
            DiagUnknownException(ctx);
            return false;
        }
    }

    // Times a kernel over the whole point set; returns the median us, or [DIAG] + false on error.
    bool TimeProjKernel(ProjKernel k, const char* label, const std::vector<Point3D>& pnts,
                        const Surface& srf, const UVRange& uvHint, double& outUs, std::string& why)
    {
        char ctx[192];
        std::snprintf(ctx, sizeof(ctx), "%s kernel=%s", label, ProjKernelName(k));
        try
        {
            outUs = MedianUs([&] {
                double dist = 0.0;
                for (const Point3D& p : pnts)
                    ApplyProjKernel(k, p, srf, uvHint, dist);
            });
            return true;
        }
        catch (const std::exception& e)
        {
            why = "threw (see [DIAG])";
            DiagException(ctx, e);
            return false;
        }
        catch (...)
        {
            why = "unknown exception (see [DIAG])";
            DiagUnknownException(ctx);
            return false;
        }
    }

    // Boolean timing: each sample uses its own clone of the inputs (cloning happens outside the
    // timed span), so only api_boolean itself is measured. Returns the number of timed samples.
    bool TimeBoolean(const BodyPtr& targetBase, const BodyPtr& toolBase,
                     double& outMedianUs, int& okTimed, std::string& why)
    {
        if (!targetBase || !toolBase)
        {
            why = "null boolean input body";
            return false;
        }

        const int n = kWarmup + kRepeat;
        std::vector<std::pair<BodyPtr, BodyPtr>> cases;
        cases.reserve(static_cast<size_t>(n));
        try
        {
            for (int i = 0; i < n; ++i)
            {
                BodyPtr tgt = sg_cast<Body>(targetBase->Clone());
                BodyPtr tool = sg_cast<Body>(toolBase->Clone());
                if (!tgt || !tool)
                {
                    why = "clone failed";
                    return false;
                }
                cases.emplace_back(tgt, tool);
            }
        }
        catch (const std::exception& e)
        {
            why = "clone threw";
            DiagException("TimeBoolean clone inputs", e);
            return false;
        }
        catch (...)
        {
            why = "clone threw";
            DiagUnknownException("TimeBoolean clone inputs");
            return false;
        }

        std::vector<double> samples;
        samples.reserve(kRepeat);
        bool diagOnce = false; // print the api_boolean [DIAG] only once to avoid flooding
        for (int i = 0; i < n; ++i)
        {
            BodyPtr tgt = cases[i].first;
            BodyList tools;
            tools.push_back(cases[i].second);
            BooleanOpts opts(BooleanType::SUBTRACTION); // BooleanType::SUBTRACTION

            const auto t0 = std::chrono::steady_clock::now();
            bool ok = false;
            try
            {
                BooleanRetPtr ret = api_boolean(tgt, tools, opts);
                ok = (ret && ret->Succeeded()); // ModelingRet::Succeeded
            }
            catch (const std::exception& e)
            {
                ok = false;
                if (!diagOnce)
                {
                    diagOnce = true;
                    DiagException("TimeBoolean api_boolean", e);
                }
            }
            catch (...)
            {
                ok = false;
                if (!diagOnce)
                {
                    diagOnce = true;
                    DiagUnknownException("TimeBoolean api_boolean");
                }
            }
            const auto t1 = std::chrono::steady_clock::now();

            if (i >= kWarmup && ok)
                samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
        }

        if (samples.empty())
        {
            why = "all boolean runs failed";
            return false;
        }
        std::sort(samples.begin(), samples.end());
        outMedianUs = samples[samples.size() / 2];
        okTimed = static_cast<int>(samples.size());
        return true;
    }
}

int main(int /*argc*/, char* /*argv*/[])
{
    try
    {
        init();

        // CSV header goes to stderr because the kernel pollutes stdout; data rows follow there too.
        std::fprintf(stderr, "op,name_a_us,name_b_us,speedup\n");
        std::printf("[CSV] columns=op,name_a_us,name_b_us,speedup data=stderr\n");

        auto cylRet = api_make_solid_cylinder(Ucs3D(), 200.0, 500.0, PI2, true);
        auto wdgRet = api_make_solid_wedge(
            Ucs3D(Point3D(100.0, 0.0, -250.0), Dir3D(1.0, 0.0, 0.0), Dir3D(0.0, 0.0, 1.0)),
            400.0, 400.0, 500.0);

        if (!cylRet || !cylRet->Succeeded() || !cylRet->ResultBody())
        {
            std::printf("[GEO] fatal: cylinder construction failed\n");
            fini();
            return 2;
        }
        if (!wdgRet || !wdgRet->Succeeded() || !wdgRet->ResultBody())
        {
            std::printf("[GEO] fatal: wedge construction failed\n");
            fini();
            return 2;
        }

        BodyPtr cyl = cylRet->ResultBody();
        BodyPtr wdg = wdgRet->ResultBody();

        // Working geometry: cyl - wdg (the analytic representation B).
        BodyPtr workB;
        {
            BodyPtr target = sg_cast<Body>(cyl->Clone());
            BodyList tools;
            tools.push_back(sg_cast<Body>(wdg->Clone()));
            BooleanOpts opts(BooleanType::SUBTRACTION);
            BooleanRetPtr bret = api_boolean(target, tools, opts);
            if (bret && bret->Succeeded() && !bret->ResultBodies().empty())
                workB = bret->ResultBodies().front();
        }

        const GeomInfo cylInfo = Inspect(cyl);
        const GeomInfo wdgInfo = Inspect(wdg);
        const GeomInfo workInfo = Inspect(workB);
        std::printf("[GEO] cylinder: faces=%zu analytic=%zu spline=%zu\n",
                    cylInfo.faces, cylInfo.analytic, cylInfo.spline);
        std::printf("[GEO] wedge:    faces=%zu analytic=%zu spline=%zu\n",
                    wdgInfo.faces, wdgInfo.analytic, wdgInfo.spline);
        std::printf("[GEO] work(cyl-wedge): faces=%zu analytic=%zu spline=%zu\n",
                    workInfo.faces, workInfo.analytic, workInfo.spline);

        // Two representations: B = analytic (as built), A = per-face ToBSpline replacement.
        int cylConv = 0, cylTotal = 0, wdgConv = 0, wdgTotal = 0, workConv = 0, workTotal = 0;
        BodyPtr cylA = SplineizeBody(cyl, cylConv, cylTotal);
        BodyPtr wdgA = SplineizeBody(wdg, wdgConv, wdgTotal);
        BodyPtr workA = workB ? SplineizeBody(workB, workConv, workTotal) : nullptr;

        std::printf("[GEO] splineize: cylinder %d/%d, wedge %d/%d, work %d/%d faces converted\n",
                    cylConv, cylTotal, wdgConv, wdgTotal, workConv, workTotal);

        // Replacement counts as success only if the work body exists and every face became spline.
        const bool replacementOk =
            (workA != nullptr && workTotal > 0 && workConv == workTotal);

        // Surfaces (and their domains) used for the intersection and projection comparisons.
        SurfacePtr srfB0, srfB1, srfA0, srfA1;
        UVRange uvB0, uvB1, uvA0, uvA1;

        if (replacementOk)
        {
            auto facesB = workB->QueryFaces();
            auto facesA = workA->QueryFaces();
            if (facesB.size() >= 2 && facesA.size() == facesB.size())
            {
                auto itB = facesB.begin();
                auto itA = facesA.begin();
                srfB0 = (*itB)->GeomSurface(); uvB0 = (*itB)->CalcUVBound();
                srfA0 = (*itA)->GeomSurface(); uvA0 = (*itA)->CalcUVBound();
                ++itB; ++itA;
                srfB1 = (*itB)->GeomSurface(); uvB1 = (*itB)->CalcUVBound();
                srfA1 = (*itA)->GeomSurface(); uvA1 = (*itA)->CalcUVBound();
            }
        }

        const bool useReplacement = replacementOk && srfB0 && srfB1 && srfA0 && srfA1;

        // Fallback pair: compare a CylinderSurface (and a Plane) with their own ToBSpline directly.
        // Surface subclasses have no usable copy constructor, so construct them in place.
        // Used when per-face replacement of the work body is unavailable.
        const UVRange fallbackCylUV(Interval(0.0, PI2), Interval(-250.0, 250.0));
        const UVRange fallbackPlaneUV(Interval(-300.0, 300.0), Interval(-300.0, 300.0));

        if (!useReplacement)
        {
            srfB0 = std::make_shared<CylinderSurface>(Ucs3D(), 200.0);
            srfB1 = std::make_shared<Plane>(
                Ucs3D(Point3D(0.0, 0.0, 0.0), Dir3D(1.0, 0.0, 0.0), Dir3D(0.0, 0.0, 1.0)));
            uvB0 = fallbackCylUV;
            uvB1 = fallbackPlaneUV;
            srfA0 = srfB0->ToBSpline(uvB0);
            srfA1 = srfB1->ToBSpline(uvB1);
            uvA0 = uvB0;
            uvA1 = uvB1;

            std::printf("[MODE] pair-of-surfaces (fallback): "
                        "per-face replacement unusable at probe stage "
                        "(work faces converted %d/%d; replacement=%s). "
                        "Compare CylinderSurface vs its ToBSpline directly.\n",
                        workConv, workTotal, replacementOk ? "true" : "false");
        }
        else
        {
            std::printf("[MODE] face-replacement: work geometry faces replaced by "
                        "BSpline (ToBSpline + MakeModelSurface + FaceAttachSurface), %d/%d faces.\n",
                        workConv, workTotal);
        }

        // 3) Operation timing
        // Projection input points: identical set and parameters for A and B. They sit near the
        // cylinder (r=200, z axis) so nearest-point projection is well defined and well conditioned;
        // the first point matches the one used by the previous probe version.
        // Points are the same objects passed to both kernels.
        const std::vector<Point3D> projPnts = {
            Point3D(250.0, 30.0, 120.0),
            Point3D(210.0, 0.0, 0.0),
            Point3D(0.0, 215.0, 80.0),
            Point3D(-205.0, 10.0, -60.0),
            Point3D(12.0, -208.0, 150.0),
            Point3D(180.0, 120.0, -120.0),
            Point3D(-150.0, -150.0, 40.0),
            Point3D(60.0, -160.0, -200.0),
        };

        // 3.1 Surface/surface intersection
        if (srfB0 && srfB1 && srfA0 && srfA1)
        {
            double aUs = 0.0, bUs = 0.0;
            try
            {
                aUs = MedianUs([&] { RunSrfSrfInt(*srfA0, *srfA1, uvA0, uvA1); });
                bUs = MedianUs([&] { RunSrfSrfInt(*srfB0, *srfB1, uvB0, uvB1); });
                Emit("srf_srf_int", true, aUs, bUs, "");
            }
            catch (const std::exception& e)
            {
                DiagException("srf_srf_int", e);
                Emit("srf_srf_int", false, 0.0, 0.0, "intersection threw (see [DIAG])");
            }
            catch (...)
            {
                DiagUnknownException("srf_srf_int");
                Emit("srf_srf_int", false, 0.0, 0.0, "intersection threw (see [DIAG])");
            }
        }
        else
        {
            Emit("srf_srf_int", false, 0.0, 0.0, "surface pair unavailable");
        }

        // 3.2 Point-to-surface projection. The kernel must work for both A and B, so candidates
        // are dry-run in order (untimed); the stub above is expected to fail and fall through to
        // Surface::CalcNearestPnt.
        if (srfB0 && srfA0)
        {
            const ProjKernel candidates[] = { ProjKernel::LegacyProject, ProjKernel::NearestPnt, ProjKernel::NearestUtil };
            ProjKernel kernel = ProjKernel::NearestPnt;
            bool kernelOk = false;
            std::string legacyWhy;

            for (ProjKernel k : candidates)
            {
                std::string whyA, whyB;
                const bool okA = ProbeProjKernel(k, "pt_srf_proj A(BSpline)", projPnts, *srfA0, uvA0, whyA);
                const bool okB = okA && ProbeProjKernel(k, "pt_srf_proj B(analytic)", projPnts, *srfB0, uvB0, whyB);
                if (!okA || !okB)
                {
                    if (k == ProjKernel::LegacyProject)
                        legacyWhy = okA ? ("B(analytic): " + whyB) : ("A(BSpline): " + whyA);
                    continue;
                }
                kernel = k;
                kernelOk = true;
                break;
            }

            if (kernelOk)
            {
                std::printf("[MODE] pt_srf_proj kernel=%s points=%zu (same point set and options for A and B)\n",
                            ProjKernelName(kernel), projPnts.size());
                if (kernel != ProjKernel::LegacyProject && !legacyWhy.empty())
                    std::printf("[MODE] pt_srf_proj: GeomProject::PntSrfProject unusable in this kernel "
                                "version (%s) -> nearest-point projection instead; see [DIAG]\n",
                                legacyWhy.c_str());

                // Sanity check (untimed): the projection distance should agree for A and B (same geometry).
                try
                {
                    double dA = 0.0, dB = 0.0, maxA = 0.0, maxB = 0.0, maxDiff = 0.0;
                    for (const Point3D& p : projPnts)
                    {
                        if (!ApplyProjKernel(kernel, p, *srfA0, uvA0, dA) ||
                            !ApplyProjKernel(kernel, p, *srfB0, uvB0, dB))
                        {
                            std::printf("[DIAG] pt_srf_proj check: empty result, skipped\n");
                            break;
                        }
                        maxA = std::max(maxA, dA);
                        maxB = std::max(maxB, dB);
                        maxDiff = std::max(maxDiff, std::fabs(dA - dB));
                    }
                    std::printf("[DIAG] pt_srf_proj check: kernel=%s points=%zu "
                                "maxDist(A)=%.6f maxDist(B)=%.6f max|distA-distB|=%.6f\n",
                                ProjKernelName(kernel), projPnts.size(), maxA, maxB, maxDiff);
                }
                catch (const std::exception& e)
                {
                    DiagException("pt_srf_proj check", e);
                }
                catch (...)
                {
                    DiagUnknownException("pt_srf_proj check");
                }

                double aUs = 0.0, bUs = 0.0;
                std::string whyA, whyB;
                const bool aOk = TimeProjKernel(kernel, "pt_srf_proj A(BSpline)", projPnts, *srfA0, uvA0, aUs, whyA);
                const bool bOk = aOk && TimeProjKernel(kernel, "pt_srf_proj B(analytic)", projPnts, *srfB0, uvB0, bUs, whyB);
                if (aOk && bOk)
                    Emit("pt_srf_proj", true, aUs, bUs, "");
                else
                    Emit("pt_srf_proj", false, 0.0, 0.0, aOk ? ("B(analytic): " + whyB) : ("A(BSpline): " + whyA));
            }
            else
            {
                std::printf("[MODE] pt_srf_proj: no usable projection kernel for both A and B; see [DIAG]\n");
                Emit("pt_srf_proj", false, 0.0, 0.0, "no usable projection kernel (see [DIAG])");
            }
        }
        else
        {
            Emit("pt_srf_proj", false, 0.0, 0.0, "surface unavailable");
        }

        // 3.3 Boolean (api_boolean destroys its inputs, hence a fresh clone per sample).
        {
            const bool splineInputsOk = (cylA && wdgA && cylConv > 0 && wdgConv > 0);

            double aUs = 0.0;
            int aTimed = 0;
            std::string aWhy;
            bool aOk = false;
            if (splineInputsOk)
            {
                try
                {
                    aOk = TimeBoolean(cylA, wdgA, aUs, aTimed, aWhy);
                }
                catch (const std::exception& e)
                {
                    aOk = false;
                    aWhy = "spline boolean threw (see [DIAG])";
                    DiagException("boolean A(BSpline)", e);
                }
                catch (...)
                {
                    aOk = false;
                    aWhy = "spline boolean threw (see [DIAG])";
                    DiagUnknownException("boolean A(BSpline)");
                }
            }
            else
            {
                aWhy = "spline body conversion failed (cylinder/wedge)";
            }

            double bUs = 0.0;
            int bTimed = 0;
            std::string bWhy;
            bool bOk = false;
            try
            {
                bOk = TimeBoolean(cyl, wdg, bUs, bTimed, bWhy);
            }
            catch (const std::exception& e)
            {
                bOk = false;
                bWhy = "analytic boolean threw (see [DIAG])";
                DiagException("boolean B(analytic)", e);
            }
            catch (...)
            {
                bOk = false;
                bWhy = "analytic boolean threw (see [DIAG])";
                DiagUnknownException("boolean B(analytic)");
            }

            if (aOk && bOk)
            {
                Emit("boolean", true, aUs, bUs, "");
            }
            else
            {
                std::string why;
                if (!aOk) why += "A:" + aWhy + " ";
                if (!bOk) why += "B:" + bWhy + " ";
                Emit("boolean", false, 0.0, 0.0, why);
            }
            (void)aTimed;
            (void)bTimed;
        }

        // 4) Human-readable summary (stdout)
        std::printf("[SUMMARY] warmup=%d repeat=%d (median us)\n", kWarmup, kRepeat);
        for (const auto& c : g_cases)
        {
            const double reduction = (c.aUs > 0.0) ? (1.0 - c.bUs / c.aUs) * 100.0 : 0.0;
            std::printf("[SUMMARY] op=%-12s A=%.3f us  B=%.3f us  speedup=A/B=%.4f  "
                        "B-vs-A reduction=%.1f%%\n",
                        c.op.c_str(), c.aUs, c.bUs, c.speedup, reduction);
        }
        std::printf("[SUMMARY] note: A=BSpline faces, B=analytic quadric faces; "
                    "target: B at least 50%% faster than A\n");

        fini();
        return 0;
    }
    catch (const std::exception& e)
    {
        DiagException("main", e);
        std::printf("[SUMMARY] fatal exception: %s\n", e.what());
        std::fprintf(stderr, "fatal,%s,NA,NA\n", e.what());
        return 1;
    }
    catch (...)
    {
        DiagUnknownException("main");
        std::printf("[SUMMARY] fatal unknown exception\n");
        std::fprintf(stderr, "fatal,unknown,NA,NA\n");
        return 1;
    }
}
