// Performance benchmark: (A) per-corpus-case detect/pipeline timing, (B) paired op timing with
// BSpline (A) vs analytic (B) geometry. Both dump CSV; timing excludes input construction and
// cloning, which happen before t0. Exit code: 0 = both [ASSERT] true, 1 = assertion failure or
// runtime exception, 2 = bad args/IO. Kernel boolean prints "128/128/0" noise on stdout.

#include "sgkquadric/detect.hpp"  // sgkq::Detect
#include "sgkquadric/io.hpp"      // LoadSgt / InputFace
#include "sgkquadric/testkit.hpp" // harness: LoadManifest / ManifestEntry / DefaultDataDir
#include "sgkquadric/pipeline.hpp"// sgkq::ProcessFace
#include "sgkquadric/types.hpp"   // Decision / FaceReport / DetectResult

#include <Foundation/BaseException.h>  // BaseException carries name/code/file:line/function
#include <Foundation/Common.h>
#include <Foundation/RuntimeClass.h>   // sg_cast
#include <Foundation/init.h>           // init()/fini()

#include <GeomBase/Dir3D.h>
#include <GeomBase/Interval.h>
#include <GeomBase/Number.h>   // PI2
#include <GeomBase/Point2D.h>
#include <GeomBase/Point3D.h>
#include <GeomBase/Toler.h>
#include <GeomBase/UVRange.h>
#include <GeomBase/Ucs3D.h>

#include <Geometry/3D/Surface/BSplineSurface.h>
#include <Geometry/3D/Surface/CylinderSurface.h>
#include <Geometry/3D/Surface/Plane.h>
// ToBSpline / CalcNearestPnt
#include <Geometry/3D/Surface/Surface.h>
#include <Geometry/SurfaceType.h>

#include <ModelingBase/ModelingRet.h>
#include <ModelingPrim/API.h>
#include <ModelingPrim/Primitives/PrimitivesRet.h>
#include <Boolean/API.h>
#include <Boolean/BooleanOpts.h>
#include <Boolean/BooleanRet.h>

#include <Topology/Brep/Body.h>
#include <Topology/Brep/Face.h>
#include <Topology/Brep/Topology.h>
#include <Topology/Entity.h>
#include <Topology/Geom/ModelSurface.h>
// MakeModelSurface / MakeFace / FaceAttachSurface
#include <Topology/Tools/TopoBuilder.h>

// GeomInt::SrfSrfInt
#include <GeomInt/GeomInt.h>

#include <algorithm>   // std::sort
#include <chrono>      // steady_clock
#include <cmath>       // std::ceil (P95 rank)
#include <cctype>      // std::isalpha / std::toupper
#include <cstdio>      // std::printf / std::fprintf
#include <cstdlib>     // std::atoi / std::exit
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>  // std::function
#include <iomanip>     // std::fixed / std::setprecision
#include <iostream>
#include <sstream>     // std::ostringstream
#include <string>
#include <system_error>
#include <typeinfo>
#include <utility>     // std::pair
#include <vector>

using namespace sggk;
using namespace sgkq;

namespace
{

// First kWarmup calls are discarded so cold-start costs (page faults, lazy init) stay out of samples.
constexpr int kWarmup = 5; //!< shared by both benchmark parts

// path helpers

//! @brief Normalize WSL `/mnt/<drive>/...` to Windows `<DRIVE>:/...`; anything else is returned
//!        unchanged.
//! @param raw Original path as a byte string (no encoding conversion is attempted).
//! @return The normalized path, or a copy of raw when it is not a `/mnt/<letter>` path.
//! @note The kernel DLL resolves only drive-letter paths, so normalize before any IO. Only the
//!       exact `/mnt/X` form (X an ASCII letter) is handled; no realpath, pure string work, O(len).
std::string NormalizePath(const std::string& raw)
{
    // "/mnt/" prefix plus at least one drive byte; size>5 keeps raw[5] and substr(6) in range.
    if (raw.size() > 5 && raw.compare(0, 5, "/mnt/") == 0)
    {
        const char drive = raw[5];
        // Only an ASCII letter is a drive; this rejects e.g. "/mnt/1foo".
        if (std::isalpha(static_cast<unsigned char>(drive)))
        {
            // A bare "/mnt/X" has no trailing content; fall back to "/" so the result is "X:/".
            const std::string rest = (raw.size() > 6) ? raw.substr(6) : std::string("/");
            std::string       out;
            // Uppercase the drive letter (Windows convention; merges "d:" with "D:").
            out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(drive))));
            // rest already carries its leading '/', so this yields "X:/...".
            out += ':';
            out += rest;
            return out;
        }
    }
    return raw;
}

//! @brief True for a POSIX absolute path ('/' or '\' prefix) or a Windows drive path ("X:").
//! @param p Candidate path.
//! @return true when absolute; false for the empty string.
//! @note Purely syntactic (used by ResolveSgt to decide whether the corpus dir must be prefixed);
//!       it never touches the filesystem.
bool IsAbsolute(const std::string& p)
{
    // POSIX / backslash-UNC prefix, or a "X:" drive letter at offset 1.
    return (!p.empty() && (p[0] == '/' || p[0] == '\\')) || (p.size() >= 2 && p[1] == ':');
}

//! @brief Resolve a manifest sgt path against the corpus dir.
//! @param corpusDir Corpus root (trailing separator optional).
//! @param sgtPath   Path recorded in the manifest; may be absolute or relative.
//! @return sgtPath unchanged when absolute, else corpusDir + '/' + sgtPath.
//! @note No existence check -- a missing file is reported later by LoadSgt via
//!       LoadResult.errorCode.
std::string ResolveSgt(const std::string& corpusDir, const std::string& sgtPath)
{
    if (IsAbsolute(sgtPath))
        return sgtPath;
    std::string base = corpusDir;
    // Add a separator only when absent, so we never produce "dir//file".
    if (!base.empty() && base.back() != '/' && base.back() != '\\')
        base.push_back('/');
    return base + sgtPath;
}

//! @brief Face to benchmark: the first BSpline face, else the first non-null face. This mirrors
//!        run_corpus so both apps see the same input face for a given case.
//! @param faces Faces parsed by LoadSgt.
//! @return Pointer to the SurfacePtr inside faces (borrowed; faces must outlive the pointer);
//!         nullptr when there is no usable face.
//! @note Two passes: pass 1 looks for BSpline (the detect algorithm's intended input), pass 2
//!       degrades to any non-null face so analytic-only corpora still yield timing data.
const SurfacePtr* PickFace(const std::vector<InputFace>& faces)
{
    // Pass 1: prefer BSpline.
    for (const InputFace& f : faces)
    {
        if (f.surface && f.surface->SurfaceType() == SurfaceType::BSpline)
            // Address inside the container; no copy.
            return &f.surface;
    }
    // Pass 2: fall back to the first usable face.
    for (const InputFace& f : faces)
    {
        if (f.surface)
            return &f.surface;
    }
    return nullptr;
}

// stats

//! Benchmark stats, all in microseconds. An empty sample set leaves every field 0 with n==0, so a
//! caller can tell "no samples" from "zero cost" via n.
struct Stats
{
    double median = 0.0; //!< us; preferred over the mean, which long-tail samples bias upward.
    double p95 = 0.0;    //!< us; nearest rank, kept to expose tail risk rather than typical cost.
    double maxUs = 0.0;  //!< us; worst observation, used to spot occasional spikes.
    size_t n = 0;        //!< number of samples (warmup excluded); = repeat unless a run failed.
};

//! @brief median / P95 (nearest rank) / max in us. Empty input -> all zero, n==0.
//! @param v Timing samples in us, taken by value because std::sort reorders them in place.
//! @return Stats; for empty input n==0 and every statistic is 0.
//! @note Median is the upper median v[n/2] (no interpolation). P95 uses nearest rank
//!       idx=ceil(0.95*n) clamped to [1,n], i.e. v[idx-1]; no interpolation, which would invent a
//!       value when n is small.
Stats ComputeStats(std::vector<double> v)
{
    Stats s;
    // The sample count doubles as the "any sample?" test for the caller.
    s.n = v.size();
    if (v.empty())
        return s;
    // Sort the value copy in place so the caller's sequence is untouched.
    std::sort(v.begin(), v.end());
    s.median = v[v.size() / 2];
    size_t idx = static_cast<size_t>(std::ceil(0.95 * static_cast<double>(v.size())));
    // ceil(0.95n) >= 1 for n>=1; guard the float/rounding edge anyway.
    if (idx == 0)
        idx = 1;
    // Clamp against the rounding edge.
    if (idx > v.size())
        idx = v.size();
    // idx is 1-based.
    s.p95 = v[idx - 1];
    s.maxUs = v.back();
    return s;
}

//! @brief Warm up kWarmup times, then time repeat runs of fn(); returns the median in us.
//! @param repeat Number of timed samples (must be > 0; enforced by main's argument check).
//! @param fn     One operation to time; it must have equivalent input size on every call and must
//!               not let its input grow across runs.
//! @return Median of the timed samples in us.
//! @note Warmup runs are neither timed nor collected; they only let the kernel reach steady state.
//!       Only fn() sits between the two steady_clock reads. steady_clock is monotonic, so it is
//!       unaffected by wall-clock jumps.
double MedianUs(int repeat, const std::function<void()>& fn)
{
    // Steady-state warmup: results and time are discarded.
    for (int i = 0; i < kWarmup; ++i)
        fn();
    std::vector<double> samples;
    // Pre-allocate to avoid reallocation jitter inside the timing loop.
    samples.reserve(static_cast<size_t>(repeat));
    for (int i = 0; i < repeat; ++i)
    {
        const auto t0 = std::chrono::steady_clock::now();
        fn(); // the timed region contains this call alone
        const auto t1 = std::chrono::steady_clock::now();
        samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    return ComputeStats(samples).median;
}

// diagnostics

//! @brief Print a std::exception diagnostic, appending BaseException location fields when present.
//! @param where Location tag (e.g. "SplineizeBody ToBSpline") used to find the source in output.
//! @param e     The caught exception (borrowed; this function only reads it).
//! @note Uses printf only, so it is safe to call inside a catch block and cannot throw again.
//!       dynamic_cast detects BaseException, whose name/code/file:line/function fields pinpoint
//!       kernel errors better than what(); non-kernel exceptions print type/what only.
void DiagException(const char* where, const std::exception& e)
{
    std::printf("[DIAG] %s: type=%s what=%s", where, typeid(e).name(), e.what());
    if (const BaseException* be = dynamic_cast<const BaseException*>(&e))
        std::printf(" name=%s code=0x%08X raised_at=%s(%d) in=%s", be->name(),
                    static_cast<unsigned int>(be->getCode()), be->getFileName(), be->getLineNumber(),
                    be->getFuncName());
    std::printf("\n"); // keep later output on its own line
}

//! @brief Diagnostic for a non-std::exception caught by catch(...); only the location is known.
//! @param where Location tag; kept distinct from DiagException output for easier grepping.
void DiagUnknownException(const char* where)
{
    std::printf("[DIAG] %s: unknown non-std::exception\n", where);
}

//! @brief Map Decision to the stable text written to the CSV.
//! @param d Decision value.
//! @return One of "Converted" / "Rejected" / "Failed"; unknown values fall back to "Failed".
//! @note An explicit switch (not a table) so the compiler warns when a new enumerator is added.
const char* DecisionName(Decision d)
{
    switch (d)
    {
    case Decision::Converted: return "Converted";
    case Decision::Rejected:  return "Rejected";
    case Decision::Failed:    return "Failed";
    }
    return "Failed"; // defensive fallback for an out-of-range value
}

//! @brief Control-point count nu*nv of a BSplineSurface (matrix laid out [U][V], matching
//!        gen_corpus::CountCtrl so the CSV column cross-checks against the manifest).
//! @param srf Surface pointer; may be null or non-BSpline.
//! @return nu*nv, or 0 for a non-BSpline/null/empty matrix.
//! @note Read-only; rows are assumed rectangular (guaranteed by the kernel).
int CountCtrl(const SurfacePtr& srf)
{
    BSplineSurfacePtr bs = sg_cast<BSplineSurface>(srf);
    if (!bs)
        // Analytic geometry has no control points.
        return 0;
    // BSplineSurface.h:337
    const Point3DMatrix& cps = bs->ControlPoints();
    const int nu = static_cast<int>(cps.size());
    // Guard cps[0] against an empty matrix.
    const int nv = (nu > 0) ? static_cast<int>(cps[0].size()) : 0;
    return nu * nv;
}

// geometry inspection

//! Snapshot of a Body's face geometry, used to confirm splineization really replaced the faces.
struct GeomInfo
{
    size_t faces = 0;    //!< total faces seen.
    size_t analytic = 0; //!< non-BSpline geometric faces.
    size_t spline = 0;   //!< BSpline faces.
};

//! @brief Walk the Body's faces and classify each geometry as analytic or BSpline.
//! @param body Body pointer; may be null.
//! @return GeomInfo; an all-zero snapshot for a null body.
//! @note A face whose GeomSurface() is null still counts in faces but in neither analytic nor
//!       spline (there is no geometry to classify). Used for [GEO] printing only, not assertions.
GeomInfo Inspect(const BodyPtr& body)
{
    GeomInfo info;
    if (!body)
        // Callers pass workB even when it failed to build.
        return info;
    for (const auto& face : body->QueryFaces()) // Body::QueryFaces [Topology/Brep/Body.h:54]
    {
        ++info.faces;
        SurfacePtr srf = face->GeomSurface(); // Face::GeomSurface [Topology/Brep/Face.h:82]
        if (!srf)
            // Geometryless face: counted in faces only.
            continue;
        if (sg_cast<BSplineSurface>(srf))
            ++info.spline;
        else
            ++info.analytic;
    }
    return info;
}

//! @brief Deep-clone src and replace each face geometry with its BSpline form.
//! @param src       Source body; a null src returns nullptr immediately.
//! @param converted Out: number of faces successfully replaced with BSpline.
//! @param total     Out: number of faces visited (including failed/geometryless ones).
//! @return A deep clone with as many faces splineized as possible; nullptr if cloning failed.
//! @note The clone is owned by the caller and src is never modified.
//!       Each face is converted in its own try/catch: one failing ToBSpline is logged and skipped
//!       while the remaining faces continue. converted==total means the body is 100% splineized,
//!       the precondition for the face-replacement mode in part B.
BodyPtr SplineizeBody(const BodyPtr& src, int& converted, int& total)
{
    converted = 0;
    total = 0;
    if (!src)
        return nullptr;

    BodyPtr clone = sg_cast<Body>(src->Clone()); // Topology::Clone [Topology/Brep/Topology.h:70]
    if (!clone)
        return nullptr;

    for (const auto& face : clone->QueryFaces())
    {
        // Counted whether or not this face converts.
        ++total;
        SurfacePtr geom = face->GeomSurface();
        if (!geom)
            // No geometry: leaves converted < total.
            continue;

        const UVRange uv = face->CalcUVBound(); // Face::CalcUVBound [Face.h:94]
        BSplineSurfacePtr bs;
        try
        {
            bs = geom->ToBSpline(uv); // Surface::ToBSpline [Surface.h:340]
        }
        catch (const std::exception& e)
        {
            // Keep the original geometry for this face.
            DiagException("SplineizeBody ToBSpline", e);
            continue;
        }
        catch (...)
        {
            DiagUnknownException("SplineizeBody ToBSpline");
            continue;
        }
        if (!bs)
            // Implicit failure without an exception.
            continue;

        ModelSurfacePtr modelSrf = TopoBuilder::MakeModelSurface(bs, true); // TopoBuilder.h:88
        if (!modelSrf)
            // Wrap failed (resource/type).
            continue;

        // Keep the original Sense so the face's inside/outside orientation is not flipped.
        TopoBuilder::FaceAttachSurface(face, modelSrf, face->Sense()); // TopoBuilder.h:255 / Face.h:63
        ++converted;
    }
    return clone;
}

// paired operation kernels

//! @brief One surface-surface intersection; the result is discarded, only completion time matters.
//! @param s0,s1 The two surfaces to intersect (read-only).
//! @param u0,u1 Their UV domains -- the search stays inside these domains (no extrapolation).
//! @return Nothing; a call that throws is treated as "operation unavailable".
//! @note opts.tol is a fixed 1e-5 absolute linear tolerance (length, mm), matching the kernel's
//!       default scale. Keeping it constant makes the A and B calls use identical parameters.
void RunSrfSrfInt(const Surface& s0, const Surface& s1, const UVRange& u0, const UVRange& u1)
{
    SrfSrfIntOpts opts;
    // absolute linear tolerance, 1e-5 mm
    opts.tol = Toler(1e-5);
    auto ret = GeomInt::SrfSrfInt(s0, s1, u0, u1, opts);
    // suppress the unused-variable warning
    (void)ret;
}

//! @brief One point projection implemented with CalcNearestPnt instead of the kernel projector.
//! @param pnt Query point in mm.
//! @param srf Target surface.
//! @note GeomProject::PntSrfProject is an unimplemented throwing stub in SGK 1.7.11. For a smooth
//!       surface the nearest point equals the normal projection; probe_speed verified this
//!       equivalence for the corpus, so CalcNearestPnt is used here. Out-of-domain points land on
//!       the domain boundary (no extrapolation). The results are not consumed.
void ApplyNearestPnt(const Point3D& pnt, const Surface& srf)
{
    // uv is required by the interface but not used.
    Point2D uv;
    const Point3D nearest = srf.CalcNearestPnt(pnt, uv);
    // timing only
    (void)nearest;
}

//! @brief Time api_boolean alone; each sample runs on its own pre-made clone, so cloning happens
//!        outside every timed region.
//! @param targetBase Base solid; null makes the call fail immediately.
//! @param toolBase   Tool solid; null makes the call fail immediately.
//! @param repeat     Number of timed samples.
//! @param outMedianUs Out: median over the successful samples, us; untouched on failure.
//! @param okTimed     Out: number of successfully timed samples (failures/throws excluded).
//! @param why         Out: failure reason for [WARN]; not overwritten on success.
//! @return true when at least one sample was timed and outMedianUs was written; false otherwise.
//! @note api_boolean performs the boolean in place and consumes its input bodies, so reusing one
//!       input across samples would leave later samples degenerate/invalid. All kWarmup+repeat
//!       clones are therefore made up front, keeping clone cost out of the timed region. A sample
//!       that throws still closes its t0/t1 interval but is not collected, so a failed run (which
//!       usually returns early) cannot drag the median down. [DIAG] is emitted at most once.
//!       The clones in cases are owned by this function and released on return.
bool TimeBoolean(const BodyPtr& targetBase, const BodyPtr& toolBase, int repeat, double& outMedianUs,
                 int& okTimed, std::string& why)
{
    if (!targetBase || !toolBase)
    {
        why = "null boolean input body";
        return false;
    }

    // One clone pair per sample, warmup included.
    const int n = kWarmup + repeat;
    std::vector<std::pair<BodyPtr, BodyPtr>> cases;
    // Avoid reallocation while cloning.
    cases.reserve(static_cast<size_t>(n));
    try
    {
        for (int i = 0; i < n; ++i)
        {
            BodyPtr tgt = sg_cast<Body>(targetBase->Clone());
            BodyPtr tool = sg_cast<Body>(toolBase->Clone());
            // Any failed clone makes all samples incomparable, so fail the whole call.
            if (!tgt || !tool)
            {
                why = "clone failed";
                return false;
            }
            // Held until the function returns.
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
    samples.reserve(static_cast<size_t>(repeat));
    // Print [DIAG] once, not per failed sample.
    bool diagOnce = false;
    for (int i = 0; i < n; ++i)
    {
        // This sample's own clone.
        BodyPtr tgt = cases[i].first;
        BodyList tools;
        tools.push_back(cases[i].second);
        BooleanOpts opts(BooleanType::SUBTRACTION); // Boolean/BooleanOpts.h:39

        const auto t0 = std::chrono::steady_clock::now();
        bool ok = false;
        try
        {
            BooleanRetPtr ret = api_boolean(tgt, tools, opts);
            // ModelingRet::Succeeded [ModelingRet.h:143]
            ok = (ret && ret->Succeeded());
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

        // Collect only timed (i>=kWarmup) and successful samples; a failed boolean usually returns
        // early and would bias the median downward.
        if (i >= kWarmup && ok)
            samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }

    // No usable sample -> the operation is unavailable.
    if (samples.empty())
    {
        why = "all boolean runs failed";
        return false;
    }
    outMedianUs = ComputeStats(samples).median;
    okTimed = static_cast<int>(samples.size());
    return true;
}

// pair rows

//! One bench_pairs.csv row (A=BSpline representation vs B=analytic representation). When ok is
//! false the whole row degrades to NA.
struct PairRow
{
    std::string op;        //!< stable operation name: srf_srf_int / pt_srf_proj / boolean.
    bool ok = false;       //!< false -> the CSV writes NA.
    double aUs = 0.0;      //!< A (BSpline) median, us.
    double bUs = 0.0;      //!< B (analytic) median, us.
    double speedup = 0.0;  //!< A/B; >1 means analytic is faster (splineization costs more).
    double reduction = 0.0;//!< (1 - B/A)*100; must be >=50 for the assertion.
    std::string reason;    //!< reason when ok=false; printed to stdout [WARN], absent from the CSV.
};

//! @brief Build a pair row from the A/B timings; a failed row is left as NA (never faked).
//! @param op     Operation name (written verbatim).
//! @param ok     Whether a valid measurement exists.
//! @param aUs    A (BSpline) time in us; ignored when ok=false.
//! @param bUs    B (analytic) time in us; ignored when ok=false.
//! @param reason Failure reason; pass an empty string when ok=true.
//! @return The assembled PairRow.
//! @note reduction = (1 - B/A)*100, i.e. how much time B saves relative to A. Division is guarded:
//!       B==0 yields speedup 0 and A==0 yields reduction 0, keeping inf/NaN out of the CSV.
PairRow MakePair(const std::string& op, bool ok, double aUs, double bUs, const std::string& reason)
{
    PairRow r;
    r.op = op;
    r.ok = ok;
    r.reason = reason;
    // Only a valid measurement carries numbers.
    if (ok)
    {
        r.aUs = aUs;
        r.bUs = bUs;
        // Guard against a zero denominator.
        r.speedup = (bUs > 0.0) ? (aUs / bUs) : 0.0;
        r.reduction = (aUs > 0.0) ? (1.0 - bUs / aUs) * 100.0 : 0.0;
    }
    return r;
}

//! @brief Print the command line usage to stderr.
//! @note Also lists the defaults so it is clear what the program reads/writes when no flag is
//!       given; goes to stderr so it never pollutes the stdout data lines.
void PrintUsage()
{
    std::fprintf(stderr,
                 "usage: bench_perf [--corpus <dir>] [--out-dir <dir>] [--repeat N]\n"
                 "       defaults: --corpus <exe_dir>/../data/corpus "
                 "--out-dir <exe_dir>/../data --repeat 50\n");
}

} // namespace

// main
//! @brief Entry point: parse arguments, run benchmark parts A and B, dump the CSVs, then choose
//!        the exit code from the assertions.
//! @param argc,argv Command line; accepts --corpus / --out-dir / --repeat.
//! @return 0 = both [ASSERT] lines true; 1 = assertion failure or runtime exception;
//!         2 = bad arguments / fatal IO.
//! @note init()/fini() are paired on every early-exit path (a missing fini would leak kernel
//!       resources). The outermost try/catch keeps any exception from escaping.
int main(int argc, char* argv[])
{
    // Defaults, derived from the executable location so a zero-argument run is independent of the
    // current working directory: <project root>/data/corpus and <project root>/data.
    std::string corpusDir = sgkq::DefaultDataDir() + "/corpus";
    std::string outDir = sgkq::DefaultDataDir();
    // Default number of timed samples per operation.
    int         repeat = 50;

    // argument parsing
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        // Consume the next argv as an option value; a missing value is a fatal argument error.
        auto next = [&](const char* name) -> const char* {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "[ERR] missing value for %s\n", name);
                PrintUsage();
                // Not yet init()-ed, so no fini() is needed on this path.
                std::exit(2);
            }
            return argv[++i];
        };

        if (a == "--corpus")
            corpusDir = NormalizePath(next("--corpus"));
        else if (a == "--out-dir")
            outDir = NormalizePath(next("--out-dir"));
        else if (a == "--repeat")
            // A non-numeric string yields 0 from atoi, which the guard below rejects.
            repeat = std::atoi(next("--repeat"));
        else
        {
            std::fprintf(stderr, "[ERR] unknown argument: %s\n", a.c_str());
            PrintUsage();
            // Fatal argument error.
            return 2;
        }
    }

    // Normalize the defaults too, so "no flag" and "flag with an equivalent WSL path" agree.
    corpusDir = NormalizePath(corpusDir);
    outDir = NormalizePath(outDir);
    // No meaningful statistics without at least one sample.
    if (repeat <= 0)
    {
        std::fprintf(stderr, "[ERR] --repeat must be > 0 (got %d)\n", repeat);
        // Fatal argument error.
        return 2;
    }

    try
    {
        // Foundation/init.h:39: the kernel runtime must be up before any SDK object is created.
        init();

        {
            // Non-throwing overload so failures can be diagnosed precisely.
            std::error_code ec;
            std::filesystem::create_directories(std::filesystem::path(outDir), ec);
            if (ec)
            {
                std::fprintf(stderr, "[ERR] cannot create out-dir %s: %s\n", outDir.c_str(),
                             ec.message().c_str());
                // Early exit must stay paired with init().
                fini();
                // Fatal argument or IO error.
                return 2;
            }
        }

        // Output/input paths, resolved against the normalized directories.
        const std::string manifestPath = ResolveSgt(corpusDir, "manifest.jsonl");
        const std::string detectCsv = ResolveSgt(outDir, "bench_detect.csv");
        const std::string pairsCsv = ResolveSgt(outDir, "bench_pairs.csv");

        // Startup banner: confirms the parameters actually in effect for this run.
        std::printf("[INFO] corpus=%s out-dir=%s repeat=%d warmup=%d\n", corpusDir.c_str(),
                    outDir.c_str(), repeat, kWarmup);

        // A. detect / pipeline timing over the corpus
        std::vector<ManifestEntry> entries;
        if (!LoadManifest(manifestPath, entries) || entries.empty())
        {
            std::fprintf(stderr, "[ERR] manifest empty or unreadable: %s\n", manifestPath.c_str());
            // Early exit must stay paired with init().
            fini();
            // Fatal argument or IO error.
            return 2;
        }

        // Pooled samples across all cases, used for the global summary lines.
        std::vector<double> pooledDetect;
        std::vector<double> pooledProcess;
        // Number of cases with >=1 valid detect sample; the detect assertion needs this > 0.
        int                 detectCases = 0;

        // Accumulated CSV rows, written in one go at the end.
        std::ostringstream detectCsvBody;

        for (const ManifestEntry& entry : entries)
        {
            const std::string sgtPath = ResolveSgt(corpusDir, entry.sgtPath);

            // Per-case timed samples: detect and pipeline kept separately.
            std::vector<double> dSamples;
            std::vector<double> pSamples;
            // Per-case tracked state, written to this case's CSV row.
            Decision            lastDecision = Decision::Failed; // last pipeline decision
            int                 lastErrorCode = 0;               // last error code (0 = no error)
            int                 ctrlCount = entry.ctrlCount;     // manifest value; measured may override
            std::string         caseNote;                        // non-empty -> [WARN] for this case

            // Warmup iterations + timed iterations.
            const int total = kWarmup + repeat;
            for (int k = 0; k < total; ++k)
            {
                // Reload every iteration: ProcessFace writes back to the .sgt, so reusing a cached
                // LoadResult would feed later samples already-modified input and break the "same
                // distribution across samples" assumption. Loading sits outside the timed region.
                LoadResult lr = LoadSgt(sgtPath);
                // A broken/missing .sgt: record the reason and error code, then skip this case.
                if (!lr.ok)
                {
                    caseNote = "LoadSgt failed code=" + std::to_string(lr.errorCode);
                    lastDecision = Decision::Failed;
                    lastErrorCode = lr.errorCode;
                    // Stop this case; the row becomes NA.
                    break;
                }
                const SurfacePtr* picked = PickFace(lr.faces);
                if (picked == nullptr || !*picked)
                {
                    caseNote = "no usable face";
                    lastDecision = Decision::Failed;
                    // Reuse io's empty-face error code.
                    lastErrorCode = static_cast<int>(LoadErrEmptyFaces);
                    // Stop this case; the row becomes NA.
                    break;
                }
                // Shared ownership, lent by lr.faces (faces outlives srf).
                const SurfacePtr srf = *picked;
                // Geometry is fixed across iterations, so measure the control-point count once.
                if (k == 0)
                {
                    const int cc = CountCtrl(srf);
                    // Only override the manifest value when the measurement is usable.
                    if (cc > 0)
                        ctrlCount = cc;
                }

                // Detect gets the spline face's own domain (U in [0,1] after ToBSpline reparam).
                const UVRange uv(srf->DomainU(), srf->DomainV());

                // pure detect: the timed region wraps Detect alone.
                const auto td0 = std::chrono::steady_clock::now();
                // detect.hpp
                const DetectResult dr = Detect(*srf, uv);
                const auto td1 = std::chrono::steady_clock::now();
                // Warmup samples are dropped.
                if (k >= kWarmup)
                    dSamples.push_back(std::chrono::duration<double, std::micro>(td1 - td0).count());
                // The result is unused here; correctness is checked by other apps.
                (void)dr;

                // full pipeline: ProcessFace takes a topological Face, so wrap the bare surface
                // first -- this construction happens outside the timed region.
                FacePtr face;
                try
                {
                    // TopoBuilder.h:88
                    ModelSurfacePtr ms = TopoBuilder::MakeModelSurface(srf, true);
                    // TopoBuilder.h:248
                    face = TopoBuilder::MakeFace(ms, true);
                }
                catch (const std::exception& e)
                {
                    DiagException("MakeFace", e);
                    caseNote = "MakeFace threw";
                    // Non-fatal: end this case, keeping the detect samples already collected.
                    break;
                }
                catch (...)
                {
                    DiagUnknownException("MakeFace");
                    caseNote = "MakeFace threw";
                    break;
                }
                if (!face)
                {
                    caseNote = "MakeFace returned null";
                    break;
                }

                // timed region wraps ProcessFace alone (decide + fit + verify + write-back).
                const auto tp0 = std::chrono::steady_clock::now();
                // pipeline.hpp
                const FaceReport rep = ProcessFace(face);
                const auto tp1 = std::chrono::steady_clock::now();
                if (k >= kWarmup)
                    pSamples.push_back(std::chrono::duration<double, std::micro>(tp1 - tp0).count());

                // Most recent decision and error code, written to this case's CSV row.
                lastDecision = rep.decision;
                lastErrorCode = rep.errorCode;
            }

            // Per-case stats (median/p95) and the [WARN] for anything that went wrong.
            const Stats ds = ComputeStats(dSamples);
            const Stats ps = ComputeStats(pSamples);
            // Only cases with real detect samples count as "effective".
            if (!dSamples.empty())
                ++detectCases;
            if (!caseNote.empty())
                std::printf("[WARN] case=%s %s\n", entry.id.c_str(), caseNote.c_str());

            // Fold this case's samples into the global pool used for the summary lines.
            for (double v : dSamples)
                pooledDetect.push_back(v);
            for (double v : pSamples)
                pooledProcess.push_back(v);

            // Assemble one CSV row; the column order must match the header exactly, and any empty
            // sample set is written as NA (never a fake 0).
            std::ostringstream row;
            row << entry.id << ',' << (entry.type.empty() ? "Unknown" : entry.type) << ','
                << ctrlCount << ',';
            if (dSamples.empty())
                row << "NA,NA,";
            else
                row << std::fixed << std::setprecision(3) << ds.median << ',' << ds.p95 << ',';
            if (pSamples.empty())
                row << "NA,NA,";
            else
                row << std::fixed << std::setprecision(3) << ps.median << ',' << ps.p95 << ',';
            row << DecisionName(lastDecision) << ',' << lastErrorCode;
            detectCsvBody << row.str() << '\n';
        }

        {
            // Binary mode keeps the bytes reproducible across platforms (no newline translation).
            std::ofstream ofs(detectCsv.c_str(), std::ios::binary);
            if (!ofs.is_open())
            {
                std::fprintf(stderr, "[ERR] cannot open %s for write\n", detectCsv.c_str());
                // Early exit must stay paired with init().
                fini();
                // Fatal argument or IO error.
                return 2;
            }
            // Header field order must match the row assembly above.
            ofs << "id,type,ctrlCount,detect_us_median,detect_us_p95,process_us_median,process_us_p95,"
                   "decision,errorCode\n";
            // The whole body in one write.
            ofs << detectCsvBody.str();
            ofs.flush();
            // A bad stream after flush means the write failed: treat as fatal.
            if (!ofs)
            {
                std::fprintf(stderr, "[ERR] write failed: %s\n", detectCsv.c_str());
                // Early exit must stay paired with init().
                fini();
                // Fatal argument or IO error.
                return 2;
            }
        }
        std::printf("[CSV] wrote %zu data row(s) -> %s\n", entries.size(), detectCsv.c_str());

        // Global summaries over the pooled samples, all in us.
        const Stats dAll = ComputeStats(pooledDetect);
        const Stats pAll = ComputeStats(pooledProcess);
        std::printf("[BENCH] detect median=%.3fus p95=%.3fus max=%.3fus over N=%d,cases=%d\n",
                    dAll.median, dAll.p95, dAll.maxUs, repeat, detectCases);
        std::printf("[BENCH] process median=%.3fus p95=%.3fus max=%.3fus over N=%d,cases=%d\n",
                    pAll.median, pAll.p95, pAll.maxUs, repeat, detectCases);

        // Assertion: at least one valid case and detect median <= 10ms (10000 us).
        const bool detectAssert = (detectCases > 0) && (dAll.median <= 10000.0);
        std::printf("[ASSERT] detect_median_le_10ms=%s\n", detectAssert ? "true" : "false");
        std::printf("[BENCH] detect assert detail: median=%.3fus th=10000us cases=%d\n", dAll.median,
                    detectCases);

        // B. paired ops (projection / intersection / boolean) before vs after splineization
        // One row per operation, in a fixed order.
        std::vector<PairRow> pairRows;

        // Base solids in analytic form. Cylinder: r=200mm, h=500mm, full circumference (PI2), with
        // caps. Wedge: 400x400 base, 500mm tall, placed at (+100,0,-250) (+X 100mm, -Z 250mm) so it
        // actually intersects the cylinder. // ModelingPrim/API.h:49,167
        auto cylRet = api_make_solid_cylinder(Ucs3D(), 200.0, 500.0, PI2, true); // ModelingPrim/API.h:49
        auto wdgRet = api_make_solid_wedge(                                     // ModelingPrim/API.h:167
            Ucs3D(Point3D(100.0, 0.0, -250.0), Dir3D(1.0, 0.0, 0.0), Dir3D(0.0, 0.0, 1.0)), 400.0,
            400.0, 500.0);

        // Analytic cylinder and wedge; each stays null if construction failed.
        BodyPtr cyl;
        BodyPtr wdg;
        if (cylRet && cylRet->Succeeded() && cylRet->ResultBody())
            // PrimitivesRet::ResultBody [PrimitivesRet.h:40]
            cyl = cylRet->ResultBody();
        if (wdgRet && wdgRet->Succeeded() && wdgRet->ResultBody())
            wdg = wdgRet->ResultBody();

        // No base solids -> all three ops are NA, never faked.
        if (!cyl || !wdg)
        {
            std::printf("[WARN] pair bench: base solid construction failed; ops=NA\n");
            pairRows.push_back(MakePair("srf_srf_int", false, 0, 0, "base solid construction failed"));
            pairRows.push_back(MakePair("pt_srf_proj", false, 0, 0, "base solid construction failed"));
            pairRows.push_back(MakePair("boolean", false, 0, 0, "base solid construction failed"));
        }
        else
        {
            // Both base solids should be entirely analytic: the baseline for "what changes".
            const GeomInfo cylInfo = Inspect(cyl);
            const GeomInfo wdgInfo = Inspect(wdg);
            std::printf("[GEO] cylinder: faces=%zu analytic=%zu spline=%zu\n", cylInfo.faces,
                        cylInfo.analytic, cylInfo.spline);
            std::printf("[GEO] wedge:    faces=%zu analytic=%zu spline=%zu\n", wdgInfo.faces,
                        wdgInfo.analytic, wdgInfo.spline);

            // Work geometry B (analytic): cyl - wdg. The one boolean gives faces with trimmed edges,
            // i.e. a "complex" carrier that enlarges the analytic-vs-spline difference.
            BodyPtr workB;
            {
                // boolean consumes its inputs, so clone
                BodyPtr target = sg_cast<Body>(cyl->Clone());
                BodyList tools;
                tools.push_back(sg_cast<Body>(wdg->Clone()));
                BooleanOpts opts(BooleanType::SUBTRACTION);
                try
                {
                    BooleanRetPtr bret = api_boolean(target, tools, opts);
                    if (bret && bret->Succeeded() && !bret->ResultBodies().empty())
                        // Subtraction yields a single body.
                        workB = bret->ResultBodies().front();
                }
                catch (const std::exception& e)
                {
                    // workB stays null -> the fallback mode is used.
                    DiagException("build work geometry (boolean)", e);
                }
                catch (...)
                {
                    DiagUnknownException("build work geometry (boolean)");
                }
            }

            // A = every face BSpline, B = analytic kept as-is.
            int cylConv = 0, cylTotal = 0, wdgConv = 0, wdgTotal = 0, workConv = 0, workTotal = 0;
            BodyPtr cylA = SplineizeBody(cyl, cylConv, cylTotal);
            BodyPtr wdgA = SplineizeBody(wdg, wdgConv, wdgTotal);
            BodyPtr workA = workB ? SplineizeBody(workB, workConv, workTotal) : nullptr;
            std::printf("[GEO] splineize: cylinder %d/%d, wedge %d/%d, work %d/%d faces\n", cylConv,
                        cylTotal, wdgConv, wdgTotal, workConv, workTotal);
            std::printf("[GEO] work(cyl-wedge): faces=%zu analytic=%zu spline=%zu\n",
                        Inspect(workB).faces, Inspect(workB).analytic, Inspect(workB).spline);

            // First two faces of the work geometry: intersection needs a pair, projection one face.
            // The A pair comes from the spline copy and the B pair from the analytic work body.
            SurfacePtr srfB0, srfB1, srfA0, srfA1;
            // Each face's own domain, taken from its own CalcUVBound.
            UVRange uvB0, uvB1, uvA0, uvA1;
            // Face-replacement is only comparable if the work body converted fully (same face
            // count and topology, so faces can be paired by index).
            const bool replacementOk = (workA != nullptr && workTotal > 0 && workConv == workTotal);

            if (replacementOk)
            {
                auto facesB = workB->QueryFaces();
                auto facesA = workA->QueryFaces();
                // Faces can be paired by index.
                if (facesB.size() >= 2 && facesA.size() == facesB.size())
                {
                    auto itB = facesB.begin();
                    auto itA = facesA.begin();
                    srfB0 = (*itB)->GeomSurface();
                    uvB0 = (*itB)->CalcUVBound();
                    srfA0 = (*itA)->GeomSurface();
                    uvA0 = (*itA)->CalcUVBound();
                    ++itB;
                    ++itA;
                    srfB1 = (*itB)->GeomSurface();
                    uvB1 = (*itB)->CalcUVBound();
                    srfA1 = (*itA)->GeomSurface();
                    uvA1 = (*itA)->CalcUVBound();
                }
            }

            // Every one of the four surfaces must be present to use the replacement mode.
            const bool useReplacement = replacementOk && srfB0 && srfB1 && srfA0 && srfA1;

            // Domains for the hand-built fallback surfaces: cylinder u in [0,2pi] rad, v in
            // [-250,250] mm; plane +/-300 mm square.
            const UVRange fallbackCylUV(Interval(0.0, PI2), Interval(-250.0, 250.0));
            const UVRange fallbackPlaneUV(Interval(-300.0, 300.0), Interval(-300.0, 300.0));

            if (!useReplacement)
            {
                // Fallback: compare a bare CylinderSurface and a Plane against their ToBSpline forms,
                // so the "analytic vs spline" contrast still exists when replacement is unusable.
                srfB0 = std::make_shared<CylinderSurface>(Ucs3D(), 200.0);
                // Plane through the origin, normal +X, ref +Z.
                srfB1 = std::make_shared<Plane>(
                    Ucs3D(Point3D(0.0, 0.0, 0.0), Dir3D(1.0, 0.0, 0.0), Dir3D(0.0, 0.0, 1.0)));
                uvB0 = fallbackCylUV;
                uvB1 = fallbackPlaneUV;
                // A side: the analytic surfaces converted to BSpline. Surface::ToBSpline [Surface.h:340]
                srfA0 = srfB0->ToBSpline(uvB0);
                srfA1 = srfB1->ToBSpline(uvB1);
                // A and B share the same domains so the comparison is fair.
                uvA0 = uvB0;
                uvA1 = uvB1;
                std::printf("[MODE] pair bench: fallback pair-of-surfaces (per-face replacement "
                            "unusable; work faces %d/%d).\n",
                            workConv, workTotal);
            }
            else
            {
                std::printf("[MODE] pair bench: face-replacement (work faces -> BSpline), %d/%d faces.\n",
                            workConv, workTotal);
            }

            // Projection points (mm), shared by A and B so input differences cannot skew the timing.
            // The 8 points sit on different sides and heights around the cylinder/plane so the
            // nearest-point search does not collapse to a single trivial optimum.
            const std::vector<Point3D> projPnts = {
                Point3D(250.0, 30.0, 120.0),    Point3D(210.0, 0.0, 0.0),
                Point3D(0.0, 215.0, 80.0),      Point3D(-205.0, 10.0, -60.0),
                Point3D(12.0, -208.0, 150.0),   Point3D(180.0, 120.0, -120.0),
                Point3D(-150.0, -150.0, 40.0),  Point3D(60.0, -160.0, -200.0),
            };

            // 3.1 surface-surface intersection. A uses the spline pair, B the analytic pair; a throw
            // on either side turns the whole row into NA (with a [DIAG] above it).
            if (srfB0 && srfB1 && srfA0 && srfA1)
            {
                try
                {
                    // Time each side separately, same operations and parameters.
                    const double aUs = MedianUs(repeat, [&] { RunSrfSrfInt(*srfA0, *srfA1, uvA0, uvA1); });
                    const double bUs = MedianUs(repeat, [&] { RunSrfSrfInt(*srfB0, *srfB1, uvB0, uvB1); });
                    pairRows.push_back(MakePair("srf_srf_int", true, aUs, bUs, ""));
                }
                catch (const std::exception& e)
                {
                    DiagException("srf_srf_int", e);
                    pairRows.push_back(MakePair("srf_srf_int", false, 0, 0, "intersection threw (see [DIAG])"));
                }
                catch (...)
                {
                    DiagUnknownException("srf_srf_int");
                    pairRows.push_back(MakePair("srf_srf_int", false, 0, 0, "intersection threw (see [DIAG])"));
                }
            }
            else
            {
                pairRows.push_back(MakePair("srf_srf_int", false, 0, 0, "surface pair unavailable"));
            }

            // 3.2 point projection (nearest point). One sample projects all 8 points, so the per
            // sample work is 8 projections; A and B use the same point set.
            if (srfB0 && srfA0)
            {
                try
                {
                    const double aUs = MedianUs(repeat, [&] {
                        for (const Point3D& p : projPnts)
                            ApplyNearestPnt(p, *srfA0);
                    });
                    const double bUs = MedianUs(repeat, [&] {
                        for (const Point3D& p : projPnts)
                            ApplyNearestPnt(p, *srfB0);
                    });
                    pairRows.push_back(MakePair("pt_srf_proj", true, aUs, bUs, ""));
                }
                catch (const std::exception& e)
                {
                    DiagException("pt_srf_proj", e);
                    pairRows.push_back(MakePair("pt_srf_proj", false, 0, 0, "projection threw (see [DIAG])"));
                }
                catch (...)
                {
                    DiagUnknownException("pt_srf_proj");
                    pairRows.push_back(MakePair("pt_srf_proj", false, 0, 0, "projection threw (see [DIAG])"));
                }
            }
            else
            {
                pairRows.push_back(MakePair("pt_srf_proj", false, 0, 0, "surface unavailable"));
            }

            // 3.3 boolean (inputs are consumed, so TimeBoolean clones per sample).
            {
                // A side needs splineized cyl/wdg with at least one converted face each.
                const bool splineInputsOk = (cylA && wdgA && cylConv > 0 && wdgConv > 0);

                double aUs = 0.0;
                int aTimed = 0;
                std::string aWhy;
                bool aOk = false;
                if (splineInputsOk)
                {
                    try
                    {
                        aOk = TimeBoolean(cylA, wdgA, repeat, aUs, aTimed, aWhy);
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
                // B side uses the analytic bodies; TimeBoolean clones internally so cyl/wdg are safe.
                try
                {
                    bOk = TimeBoolean(cyl, wdg, repeat, bUs, bTimed, bWhy);
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
                    pairRows.push_back(MakePair("boolean", true, aUs, bUs, ""));
                }
                else
                {
                    // Join both sides' reasons so it is clear whether A or B failed.
                    std::string why;
                    if (!aOk) why += "A:" + aWhy + " ";
                    if (!bOk) why += "B:" + bWhy + " ";
                    pairRows.push_back(MakePair("boolean", false, 0, 0, why));
                }
                // Successful-sample counts are diagnostic only, not used in CSV or assertions.
                (void)aTimed;
                (void)bTimed;
            }
        }

        {
            std::ofstream ofs(pairsCsv.c_str(), std::ios::binary);
            if (!ofs.is_open())
            {
                std::fprintf(stderr, "[ERR] cannot open %s for write\n", pairsCsv.c_str());
                // Early exit must stay paired with init().
                fini();
                // Fatal argument or IO error.
                return 2;
            }
            ofs << "op,name_a_us,name_b_us,speedup,reduction_pct\n";
            for (const PairRow& r : pairRows)
            {
                if (!r.ok)
                {
                    // Keep the column count consistent (name plus four NA fields).
                    ofs << r.op << ",NA,NA,NA,NA\n";
                }
                else
                {
                    // A/B 3 decimals (us), speedup 4, reduction 1 (%).
                    ofs << r.op << ',' << std::fixed << std::setprecision(3) << r.aUs << ',' << r.bUs
                        << ',' << std::setprecision(4) << r.speedup << ',' << std::setprecision(1)
                        << r.reduction << '\n';
                }
            }
            ofs.flush();
            if (!ofs)
            {
                std::fprintf(stderr, "[ERR] write failed: %s\n", pairsCsv.c_str());
                // Early exit must stay paired with init().
                fini();
                // Fatal argument or IO error.
                return 2;
            }
        }

        // Echo the results to stdout, one full metric line per successful op (failed rows -> [WARN]).
        std::printf("[PAIRS] columns=op,name_a_us,name_b_us,speedup,reduction_pct A=BSpline B=analytic\n");
        for (const PairRow& r : pairRows)
        {
            if (!r.ok)
            {
                std::printf("[WARN] pair op=%s NA reason=%s\n", r.op.c_str(), r.reason.c_str());
            }
            else
            {
                std::printf("[PAIRS] op=%-12s A(BSpline)=%.3fus B(analytic)=%.3fus speedup(A/B)=%.4f "
                            "reduction=%.1f%%\n",
                            r.op.c_str(), r.aUs, r.bUs, r.speedup, r.reduction);
            }
        }

        // Assertion: exactly 3 rows, each ok and with reduction >= 50% (analytic at least halves the
        // elapsed time relative to spline).
        bool pairsAssert = (pairRows.size() == 3);
        for (const PairRow& r : pairRows)
        {
            if (!r.ok || r.reduction < 50.0)
                // Any single failure fails the whole assertion.
                pairsAssert = false;
        }
        std::printf("[ASSERT] each_op_reduction_ge_50pct=%s\n", pairsAssert ? "true" : "false");

        // Exit code contract: 0 when both assertions hold, 1 otherwise.
        const bool allAsserts = detectAssert && pairsAssert;

        // Foundation/init.h:42
        fini();
        return allAsserts ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        // An uncaught exception means runtime state is unreliable; do not call fini() here.
        DiagException("main", e);
        std::fprintf(stderr, "[FAIL] fatal exception: %s\n", e.what());
        // Runtime failure / not up to spec.
        return 1;
    }
    catch (...)
    {
        // Last-resort catch for a non-std::exception escaping main.
        DiagUnknownException("main");
        std::fprintf(stderr, "[FAIL] fatal unknown exception\n");
        // Runtime failure / not up to spec.
        return 1;
    }
}


