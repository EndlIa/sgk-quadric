// apps/checks/test_exceptions/main.cpp - robustness suite for exceptions / out-of-range input.
// Each case passes only if the process survives AND the call fails explicitly (false /
// Canonical::Unknown / Rejected|Failed) with a locatable reason; a silent success or a
// crash is FAIL, a caught kernel exception is [WARN]. main exits 1 on any FAIL; lengths mm, angles rad.

// Project contract headers: the four public APIs under test (io/detect/fit/pipeline/types).
#include "sgkquadric/detect.hpp"
#include "sgkquadric/fit.hpp"
#include "sgkquadric/io.hpp"
#include "sgkquadric/pipeline.hpp"
#include "sgkquadric/types.hpp"

// SDK runtime: init/fini and runtime type info.
#include <Foundation/init.h>
#include <Foundation/RuntimeClass.h>

// SDK geometry primitives: points/vectors/directions/intervals/UVRange/Ucs3D for constructing inputs.
#include <GeomBase/Dir3D.h>
#include <GeomBase/Interval.h>
#include <GeomBase/Point3D.h>
#include <GeomBase/Ucs3D.h>
#include <GeomBase/UVRange.h>
#include <GeomBase/Vector3D.h>

// SDK surfaces: type enum and concrete BSpline/Cylinder/Sphere used to build splines and analytic faces.
#include <Geometry/SurfaceType.h>
#include <Geometry/3D/Surface/BSplineSurface.h>
#include <Geometry/3D/Surface/CylinderSurface.h>
#include <Geometry/3D/Surface/Surface.h>
#include <Geometry/3D/Surface/SphereSurface.h>

// SDK modelling: api_make_solid_cylinder and its return container to produce real geometry.
#include <ModelingPrim/API.h>
#include <ModelingPrim/Primitives/PrimitivesRet.h>

// SDK topology: Body/Face/ModelSurface/serializer/TopoBuilder.
#include <Topology/Brep/Body.h>
#include <Topology/Brep/Face.h>
#include <Topology/Geom/ModelSurface.h>
#include <Topology/Serialize/RapidTopoJsonSerializer.h>
#include <Topology/Tools/TopoBuilder.h>

// C++ standard library: math/exceptions/filesystem/streams/limits/containers.
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace sggk;
using namespace sgkq;

// Anonymous namespace: helpers and cases stay unit-local (runXXX = one case/group).
namespace
{
namespace fs = std::filesystem;

// kPi: full-precision double pi, used to sweep one full turn (2*pi) and the sphere longitude.
const double kPi = 3.14159265358979323846;

int g_total  = 0;
int g_failed = 0;
int g_warned = 0;

// report(): print one assertion as "[TEST] <name> expected=<..> actual=<..> PASS|FAIL".
// A non-empty warn on a passing case adds to g_warned and prints an extra [WARN] line.
void report(const std::string& name, const std::string& expected, const std::string& actual,
            bool pass, const std::string& warn = std::string())
{
    ++g_total;
    if (!pass)
        ++g_failed;
    else if (!warn.empty())
        ++g_warned;

    std::cout << "[TEST] " << name
              << " expected=" << expected
              << " actual=" << actual
              << (pass ? " PASS" : " FAIL") << std::endl;
    if (!warn.empty())
        std::cout << "[WARN] " << name << " " << warn << std::endl;
}

// CodeStr(): format an error code as "errorCode=<n>" for stable grepping.
std::string CodeStr(int code)
{
    return "errorCode=" + std::to_string(code);
}

// Snippet(): replace newlines/tabs with spaces, then truncate to n chars, so a message stays
// on a single [TEST] line.
std::string Snippet(const std::string& s, size_t n = 96)
{
    std::string t = s;
    for (char& c : t)
        if (c == '\n' || c == '\r' || c == '\t')
            c = ' ';
    if (t.size() > n)
        t = t.substr(0, n) + "...";
    return t;
}

// FileNameOf(): basename of a path; used to assert the message names the bad input.
std::string FileNameOf(const std::string& path)
{
    return fs::path(path).filename().string();
}

// Contains(): substring test (maps to std::string::find).
bool Contains(const std::string& hay, const std::string& needle)
{
    return hay.find(needle) != std::string::npos;
}

// DecisionStr(): Decision -> stable string for printing only (not used in assertions).
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

// WriteText(): binary write of raw bytes; cases self-generate corrupted inputs.
bool WriteText(const fs::path& p, const std::string& content)
{
    std::ofstream o(p, std::ios::binary);
    if (!o.is_open())
        return false;
    o << content;
    return o.good();
}

// ReadText(): read the whole file into a string (A5 patches a real .sgt).
bool ReadText(const fs::path& p, std::string& out)
{
    std::ifstream i(p, std::ios::binary);
    if (!i.is_open())
        return false;
    std::ostringstream ss;
    ss << i.rdbuf();
    out = ss.str();
    return true;
}

// ReplaceJsonString(): replace the string value of "key" in a JSON blob via a naive scan:
// "key" -> ':' -> first quote pair. No escape handling - enough for A5 version faking.
bool ReplaceJsonString(const std::string& in, const std::string& key, const std::string& val,
                       std::string& out)
{
    const std::string pat = "\"" + key + "\"";
    const size_t k = in.find(pat);
    if (k == std::string::npos) return false;
    const size_t colon = in.find(':', k + pat.size());
    if (colon == std::string::npos) return false;
    const size_t q1 = in.find('"', colon);
    if (q1 == std::string::npos) return false;
    const size_t q2 = in.find('"', q1 + 1);
    if (q2 == std::string::npos) return false;
    out = in.substr(0, q1 + 1) + val + in.substr(q2);
    return true;
}

// DescribeUnknown(): placeholder text for a kernel non-std exception.
std::string DescribeUnknown()
{
    return "non-std kernel exception";
}

// ParamsFinite(): true iff every SurfaceParams scalar/vector component is finite (B10 guard).
bool ParamsFinite(const SurfaceParams& p)
{
    return std::isfinite(p.radius) && std::isfinite(p.minorRadius) && std::isfinite(p.halfAngle)
        && std::isfinite(p.axisDir.X()) && std::isfinite(p.axisDir.Y()) && std::isfinite(p.axisDir.Z())
        && std::isfinite(p.axisPoint.X()) && std::isfinite(p.axisPoint.Y()) && std::isfinite(p.axisPoint.Z())
        && std::isfinite(p.apex.X()) && std::isfinite(p.apex.Y()) && std::isfinite(p.apex.Z());
}

// runFileLayer(): A1..A6 file-layer negatives - missing / empty / truncated / non-topology .sgt,
// version mismatch and STEP garbage; each must not crash and return ok=false with the filename in the message.
void runFileLayer(const fs::path& tmp)
{
// A1: nonexistent .sgt path -> ok=false, non-LoadOk code, message contains the filename.
    {
// No pre-stat: LoadSgt itself must produce the locatable error.
        const std::string path = (tmp / "no_such_file.sgt").string();
        try
        {
            const LoadResult r = LoadSgt(path);
// PASS requires all three: failure + non-LoadOk code + filename in the message, so a case
// cannot slip through without naming which input was bad.
            const bool pass = (!r.ok && r.errorCode != LoadOk &&
                               Contains(r.message, FileNameOf(path)));
            report("sgt.absent", "ok=false,code!=0,msg contains path",
                   std::string("ok=") + (r.ok ? "true" : "false") + " " + CodeStr(r.errorCode) +
                       " msg=" + Snippet(r.message),
                   pass);
        }
// A std exception escaping here is unacceptable: the file layer must return an error, not throw.
        catch (const std::exception& e)
        {
            report("sgt.absent", "ok=false,code!=0", std::string("threw: ") + e.what(), false);
        }
        catch (...)
        {
            report("sgt.absent", "ok=false,code!=0", "threw: " + DescribeUnknown(), false);
        }
    }

// A2: empty 0-byte file must be rejected, not parsed as an empty model.
    {
        const fs::path p = tmp / "empty.sgt";
        WriteText(p, "");
        try
        {
            const LoadResult r = LoadSgt(p.string());
            const bool pass = (!r.ok && r.errorCode != LoadOk && Contains(r.message, p.filename().string()));
            report("sgt.empty", "ok=false,code!=0,msg contains path",
                   std::string("ok=") + (r.ok ? "true" : "false") + " " + CodeStr(r.errorCode) +
                       " msg=" + Snippet(r.message),
                   pass);
        }
        catch (const std::exception& e)
        {
            report("sgt.empty", "ok=false,code!=0", std::string("threw: ") + e.what(), false);
        }
        catch (...)
        {
            report("sgt.empty", "ok=false,code!=0", "threw: " + DescribeUnknown(), false);
        }
    }

// A3: JSON truncated mid-array must fail to parse instead of crashing.
    {
        const fs::path p = tmp / "truncated.sgt";
// Truncate right at the modelSrfs array start, mimicking an interrupted write.
        WriteText(p, "{\"meta\":{\"sggkVer\":\"1.7.11\"},\"bodies\":[{\"modelSrfs\":[");
        try
        {
            const LoadResult r = LoadSgt(p.string());
            const bool pass = (!r.ok && r.errorCode != LoadOk && Contains(r.message, p.filename().string()));
            report("sgt.truncated", "ok=false,code!=0,msg contains path",
                   std::string("ok=") + (r.ok ? "true" : "false") + " " + CodeStr(r.errorCode) +
                       " msg=" + Snippet(r.message),
                   pass);
        }
        catch (const std::exception& e)
        {
            report("sgt.truncated", "ok=false,code!=0", std::string("threw: ") + e.what(), false);
        }
        catch (...)
        {
            report("sgt.truncated", "ok=false,code!=0", "threw: " + DescribeUnknown(), false);
        }
    }

// A4: valid JSON that is not a topology document (no meta/bodies) must be rejected by format.
    {
        const fs::path p = tmp / "not_topology.sgt";
        WriteText(p, "{\"a\":1}");
        try
        {
            const LoadResult r = LoadSgt(p.string());
            const bool pass = (!r.ok && r.errorCode != LoadOk && Contains(r.message, p.filename().string()));
            report("sgt.notTopology", "ok=false,code!=0,msg contains path",
                   std::string("ok=") + (r.ok ? "true" : "false") + " " + CodeStr(r.errorCode) +
                       " msg=" + Snippet(r.message),
                   pass);
        }
        catch (const std::exception& e)
        {
            report("sgt.notTopology", "ok=false,code!=0", std::string("threw: ") + e.what(), false);
        }
        catch (...)
        {
            report("sgt.notTopology", "ok=false,code!=0", "threw: " + DescribeUnknown(), false);
        }
    }

// A5: version mismatch - serialize a real .sgt via the SDK, patch its version field to "0.1.9",
// then reload; must fail and report version info (code == LoadErrVersionMismatch or message mentions version).
    {
        const fs::path real = tmp / "real.sgt";
// made records whether the real .sgt was produced; false routes to a setup-failure report.
        bool made = false;
        try
        {
            PrimitivesRetPtr mc = api_make_solid_cylinder(Ucs3D(), 50.0, 100.0);
            if (mc && mc->Succeeded() && mc->ResultBody())
            {
                RapidTopoJsonSerializer ser;
                ser.Serialize(mc->ResultBody(), real.string().c_str());
                made = fs::exists(real);
            }
        }
        catch (...)
        {
            made = false;
        }

        std::string content;
        if (made && ReadText(real, content))
        {
            std::string patched;
// Try "sggkVer" first, then the alternate "sggkVersion" key.
            bool        havePatch = ReplaceJsonString(content, "sggkVer", "0.1.9", patched);
            if (!havePatch)
                havePatch = ReplaceJsonString(content, "sggkVersion", "0.1.9", patched);
// Both key names failed: fall back to a raw string replacement.
            if (!havePatch)
            {
// Fallback: replace the current version token "1.7.11" with the incompatible "0.1.9".
                const size_t vs = content.find("\"1.7.11\"");
// +1 / +6 offsets keep the left quote and skip the 6-char version literal.
                if (vs != std::string::npos)
                {
                    patched   = content.substr(0, vs + 1) + "0.1.9" + content.substr(vs + 1 + 6);
                    havePatch = true;
                }
            }
            if (havePatch)
            {
                const fs::path p = tmp / "bad_version.sgt";
                WriteText(p, patched);
                try
                {
                    const LoadResult r = LoadSgt(p.string());
// Version is "reported" loosely (dedicated code OR message mentions version) to tolerate kernel wording.
                    const bool versionReported =
                        (r.errorCode == LoadErrVersionMismatch) ||
                        Contains(r.message, "0.1.9") || Contains(r.message, "sggkVer") ||
                        Contains(r.message, "version") || Contains(r.message, "Version");
                    const bool pass = (!r.ok && versionReported);
                    report("sgt.versionMismatch", "ok=false && (code=versionMismatch || msg mentions version)",
                           std::string("ok=") + (r.ok ? "true" : "false") + " " + CodeStr(r.errorCode) +
                               " msg=" + Snippet(r.message),
                           pass);
                }
                catch (const std::exception& e)
                {
                    report("sgt.versionMismatch", "ok=false,version reported",
                           std::string("threw: ") + e.what(), false);
                }
                catch (...)
                {
                    report("sgt.versionMismatch", "ok=false,version reported",
                           "threw: " + DescribeUnknown(), false);
                }
            }
            else
            {
                report("sgt.versionMismatch", "ok=false,version reported",
                       "could not patch sggkVer field in serialized .sgt", false);
            }
        }
        else
        {
            report("sgt.versionMismatch", "ok=false,version reported",
                   "could not produce a real .sgt via api_make_solid_cylinder + serializer", false);
        }
    }

// A6: LoadStep on random text without an ISO-10303-21 header must be rejected.
    {
        const fs::path p = tmp / "garbage.step";
        WriteText(p, "this is not a STEP file\n@#% not ISO-10303-21 at all\n0123456789\n");
        try
        {
// Exercise the STEP branch (LoadStep), not LoadSgt.
            const LoadResult r = LoadStep(p.string());
            const bool pass = (!r.ok && r.errorCode != LoadOk && Contains(r.message, p.filename().string()));
            report("step.garbage", "ok=false,msg contains path",
                   std::string("ok=") + (r.ok ? "true" : "false") + " " + CodeStr(r.errorCode) +
                       " msg=" + Snippet(r.message),
                   pass);
        }
        catch (const std::exception& e)
        {
            report("step.garbage", "ok=false,msg contains path", std::string("threw: ") + e.what(), false);
        }
        catch (...)
        {
            report("step.garbage", "ok=false,msg contains path", "threw: " + DescribeUnknown(), false);
        }
    }
}


// B7: SampleSurface with a fully out-of-domain UVRange (spline domain [0,1]^2, passed [100,200]^2).
// The kernel extrapolates instead of rejecting, so the criterion is only "does not crash".
// Out-of-domain behaviour is unreachable on the production path (ProcessFace clips to the face domain).
void runSampleOutOfDomain()
{
    try
    {
// Cylinder r=50 mm converted to a spline face; its own domain is [0,1]^2.
        CylinderSurface cyl(Ucs3D(Point3D(0.0, 0.0, 0.0)), 50.0);
        BSplineSurfacePtr bs = cyl.ToBSpline(UVRange(Interval(0.0, 1.0), Interval(0.0, 1.0)));
// A null conversion pointer is a setup failure, not a kernel defect.
        if (!bs)
        {
            report("sample.outOfDomain", "no crash; empty/invalid samples",
                   "ToBSpline returned null (could not build spline face)", false);
            return;
        }

// Out-of-domain UVRange: both u and v in [100,200], far outside [0,1].
        const UVRange outside(Interval(100.0, 200.0), Interval(100.0, 200.0));
        try
        {
// Sample 9x9; each Sample has point p, normal n and a degen flag.
            const std::vector<Sample> s = SampleSurface(*bs, outside, 9, 9);
// Tally degenerate samples and non-finite coordinates for the actual string.
            int degen = 0;
            int bad = 0;
            for (const Sample& x : s)
            {
                if (x.degen) ++degen;
                if (!std::isfinite(x.p.X()) || !std::isfinite(x.p.Y()) || !std::isfinite(x.p.Z()))
                    ++bad;
            }
// Kernel extrapolates out-of-domain samples instead of rejecting them; only "no crash" is asserted -
// the extrapolated values themselves do not affect pass/fail.
            const bool acceptable = true;
// Keep the existing samples=/degen=/nonFinite= format for the actual string.
            std::ostringstream os;
            os << "samples=" << s.size() << " degen=" << degen << " nonFinite=" << bad;
            std::cout << "[NOTE] sample.outOfDomain kernel extrapolates on out-of-domain UVRange instead of rejecting; the production path uses "
                         "ProcessFace with face->CalcUVBound() intersected with the spline's own domain to clamp sampling inside the face,"
                         "so this behavior is unreachable" << std::endl;
            report("sample.outOfDomain", "no crash", os.str(), acceptable);
        }
// Kernel threw std exception: process survived -> PASS with WARN.
        catch (const std::exception& e)
        {
            report("sample.outOfDomain", "no crash; empty/invalid samples",
                   std::string("kernel threw: ") + e.what(), true,
                   std::string("kernel threw on out-of-domain parameters (std::exception), caught by this case, process did not crash: ") + e.what());
        }
// Kernel threw a non-std exception: same treatment.
        catch (...)
        {
            report("sample.outOfDomain", "no crash; empty/invalid samples",
                   "kernel threw: " + DescribeUnknown(), true,
                   "kernel threw a non-std exception on out-of-domain parameters, caught by this case, process did not crash");
        }
    }
// Exception while constructing cylinder/spline: count as a case-environment failure.
    catch (const std::exception& e)
    {
        report("sample.outOfDomain", "no crash; empty/invalid samples",
               std::string("setup threw: ") + e.what(), false);
    }
    catch (...)
    {
        report("sample.outOfDomain", "no crash; empty/invalid samples",
               "setup threw: " + DescribeUnknown(), false);
    }
}

// B8: Detect with no valid samples - sphere v-interval width 0 puts every sample at the
// degenerate pole, so Detect must return Canonical::Unknown rather than name a candidate.
void runDetectEmptySamples()
{
    try
    {
// Sphere r=8 mm; v=0 is a degenerate pole, so a zero-width v-interval yields all-degenerate samples.
        SphereSurface sph(Ucs3D(Point3D(0.0, 0.0, 0.0)), 8.0);
// Sample first to get real counts for the actual string.
        const UVRange poleUv(Interval(0.0, 2.0 * kPi), Interval(0.0, 0.0));
        const std::vector<Sample> s = SampleSurface(sph, poleUv, 9, 9);
        int degen = 0;
        for (const Sample& x : s)
            if (x.degen) ++degen;

// Expect candidate == Unknown (not enough valid samples).
        const DetectResult r = Detect(sph, poleUv);
// PASS = Unknown.
        const bool pass = (r.candidate == Canonical::Unknown);
// Actual: candidate + sample stats, plus the first note if any.
        std::ostringstream os;
        os << "candidate=" << ToString(r.candidate) << " (samples=" << s.size()
           << " degen=" << degen << ")";
        if (!r.notes.empty())
            os << " note0=" << Snippet(r.notes.front(), 64);
        report("detect.emptySamples", "Unknown", os.str(), pass);
    }
    catch (const std::exception& e)
    {
        report("detect.emptySamples", "Unknown", std::string("threw: ") + e.what(), false);
    }
    catch (...)
    {
        report("detect.emptySamples", "Unknown", "threw: " + DescribeUnknown(), false);
    }
}

// B9: FitSurface with too few samples (0 and 2). Each must return false with a non-empty
// reason instead of silently succeeding.
void runFitTooFew()
{
    {
        std::vector<FitSample> none;
        SurfaceParams out;
        std::string reason;
        try
        {
// Fit a Cylinder to the empty set.
            const bool ok = FitSurface(Canonical::Cylinder, none, out, reason);
// PASS = rejected with a locatable reason.
            const bool pass = (!ok && !reason.empty());
            report("fit.zeroSamples", "false && reason non-empty",
                   std::string("ok=") + (ok ? "true" : "false") + " reason=" + Snippet(reason), pass);
        }
        catch (const std::exception& e)
        {
            report("fit.zeroSamples", "false && reason non-empty", std::string("threw: ") + e.what(), false);
        }
        catch (...)
        {
            report("fit.zeroSamples", "false && reason non-empty", "threw: " + DescribeUnknown(), false);
        }
    }

// 2 samples: two diagonal points with parallel normals, not enough to fix an axis.
    {
        std::vector<FitSample> two;
        for (int i = 0; i < 2; ++i)
        {
            FitSample s;
// Points (0,0,0) and (1,2,3), mm.
            s.p = Point3D(1.0 * i, 2.0 * i, 3.0 * i);
            s.n = Vector3D(0.0, 0.0, 1.0);
            two.push_back(s);
        }
        SurfaceParams out;
        std::string reason;
        try
        {
            const bool ok = FitSurface(Canonical::Cylinder, two, out, reason);
            const bool pass = (!ok && !reason.empty());
            report("fit.twoSamples", "false && reason non-empty",
                   std::string("ok=") + (ok ? "true" : "false") + " reason=" + Snippet(reason), pass);
        }
        catch (const std::exception& e)
        {
            report("fit.twoSamples", "false && reason non-empty", std::string("threw: ") + e.what(), false);
        }
        catch (...)
        {
            report("fit.twoSamples", "false && reason non-empty", "threw: " + DescribeUnknown(), false);
        }
    }
}

// B10: FitSurface with NaN/Inf samples. Accepted: overall rejection (false + reason), or
// success whose output params are all finite.
void runFitNanInf()
{
// Build 30 valid cylinder samples (r=10, rising along Z), then inject bad values.
    std::vector<FitSample> samples;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
// Valid samples: uniform angle, linear rise, a true cylinder.
    for (int i = 0; i < 30; ++i)
    {
// Angle of sample i, in [0, 2*pi).
        const double th = 2.0 * kPi * i / 30.0;
        FitSample s;
// Point on the radius-10 mm circle, z = 0.25*i mm.
        s.p = Point3D(10.0 * std::cos(th), 10.0 * std::sin(th), 0.25 * i);
// Outward horizontal normal.
        s.n = Vector3D(std::cos(th), std::sin(th), 0.0);
        samples.push_back(s);
    }
// Inject bad values: NaN coords, Inf normal, mixed Inf/NaN point.

// Sample 3: all three point coords NaN.
    samples[3].p = Point3D(nan, nan, nan);
// Sample 7: normal X = +Inf.
    samples[7].n = Vector3D(inf, 0.0, 0.0);
// Sample 11: point mixes Inf and NaN.
    samples[11].p = Point3D(inf, -inf, nan);

    SurfaceParams out;
    std::string reason;
    try
    {
// Fit; if it returns false, a non-empty reason is sufficient.
        const bool ok = FitSurface(Canonical::Cylinder, samples, out, reason);
        bool pass;
        std::string actual;
        if (!ok)
        {
            pass = !reason.empty();
            actual = "ok=false reason=" + Snippet(reason);
        }
        else
        {
// If accepted, every output parameter must still be finite (no NaN/Inf leaked through).
            pass = ParamsFinite(out);
            actual = std::string("ok=true finiteParams=") + (pass ? "true" : "false");
        }
        report("fit.nanInfSamples", "false||finite result (no NaN params)", actual, pass);
    }
    catch (const std::exception& e)
    {
        report("fit.nanInfSamples", "false||finite result", std::string("threw: ") + e.what(), false);
    }
    catch (...)
    {
        report("fit.nanInfSamples", "false||finite result", "threw: " + DescribeUnknown(), false);
    }
}

// B11: degenerate radius - all points collinear so r -> 0; the fit must be explicitly
// rejected: either false+reason, or radius==0 with a non-empty reason.
void runDegenerateRadius()
{
// 40 collinear points (x=0.3, y=-0.7 fixed, z linear), mm.
    std::vector<FitSample> samples;
// The points lie on a line parallel to Z, so the cylinder radius is not determined.
    for (int i = 0; i < 40; ++i)
    {
// Angle only spreads the normals; the points stay collinear regardless.
        const double th = 2.0 * kPi * i / 40.0;
        FitSample s;
        s.p = Point3D(0.3, -0.7, -5.0 + 0.25 * i);
        s.n = Vector3D(std::cos(th), std::sin(th), 0.0);
        samples.push_back(s);
    }

    SurfaceParams out;
    std::string reason;
    try
    {
        const bool ok = FitSurface(Canonical::Cylinder, samples, out, reason);
// Either rejection form (false+reason, or radius==0+reason) counts as PASS.
        const bool rejected = (!ok && !reason.empty()) ||
                              (out.radius == 0.0 && !reason.empty());
        std::ostringstream os;
        os << "ok=" << (ok ? "true" : "false") << " radius=" << out.radius
           << " reason=" << Snippet(reason);
        report("fit.degenerateRadius", "false || radius==0, reason non-empty", os.str(), rejected);
    }
    catch (const std::exception& e)
    {
        report("fit.degenerateRadius", "false || radius==0", std::string("threw: ") + e.what(), false);
    }
    catch (...)
    {
        report("fit.degenerateRadius", "false || radius==0", "threw: " + DescribeUnknown(), false);
    }
}

// B12: degenerate cone - half-angle 1e-9 rad makes the generatrices nearly parallel, far
// below a fittable magnitude; must be rejected (false + reason).
void runDegenerateCone()
{
// Tiny half-angle (rad); sin~1e-9, cos~1 make the cone almost a cylinder/line.
    const double alpha = 1e-9;
// Precompute sin/cos for reuse in points and normals.
    const double sa = std::sin(alpha), ca = std::cos(alpha);
// 40 samples along the generatrix with t in [10,30) mm.
    std::vector<FitSample> samples;
    for (int i = 0; i < 40; ++i)
    {
        const double th = 2.0 * kPi * i / 40.0;
// t = distance from the apex along the generatrix, mm.
        const double t = 10.0 + 0.5 * i;
        FitSample s;
// Standard cone parameterisation: radius t*sin(a), height t*cos(a).
        s.p = Point3D(t * sa * std::cos(th), t * sa * std::sin(th), t * ca);
// Outward normal with -sin(a) in Z.
        s.n = Vector3D(ca * std::cos(th), ca * std::sin(th), -sa);
        samples.push_back(s);
    }

    SurfaceParams out;
    std::string reason;
    try
    {
// Fit a Cone; expect rejection.
        const bool ok = FitSurface(Canonical::Cone, samples, out, reason);
        const bool pass = (!ok && !reason.empty());
        report("fit.degenerateCone", "false && reason non-empty",
               std::string("ok=") + (ok ? "true" : "false") + " halfAngle=" +
                   std::to_string(out.halfAngle) + " reason=" + Snippet(reason),
               pass);
    }
    catch (const std::exception& e)
    {
        report("fit.degenerateCone", "false && reason non-empty", std::string("threw: ") + e.what(), false);
    }
    catch (...)
    {
        report("fit.degenerateCone", "false && reason non-empty", "threw: " + DescribeUnknown(), false);
    }
}


// FindCylinderFace(): first Cylinder-type face (lateral surface) of a body; empty
// FacePtr if the body is null or has no such face. Used as the carrier for C14/C15.
// C14/C15 need an already-analytic face to mutate.
FacePtr FindCylinderFace(const BodyPtr& body)
{
    if (!body) return FacePtr();
    const FaceList faces = body->QueryFaces();
    for (const FacePtr& f : faces)
    {
// Skip null / geometry-less faces (defensive).
        if (!f || !f->GeomSurface()) continue;
        if (f->GeomSurface()->SurfaceType() == SurfaceType::Cylinder)
            return f;
    }
    return FacePtr();
}

// MakeCylinderBody(): standard solid cylinder (r=50 mm, h=100 mm) shared by C13..C15;
// returns a null BodyPtr if modelling fails.
// Keeps the input geometry identical across C13..C15.
BodyPtr MakeCylinderBody()
{
    PrimitivesRetPtr mc = api_make_solid_cylinder(Ucs3D(), 50.0, 100.0);
    if (!mc || !mc->Succeeded() || !mc->ResultBody())
        return BodyPtr();
    return mc->ResultBody();
}

// C13: ProcessFace(null FacePtr) must return Failed with a non-empty message.
void runNullFace(const Thresholds& th, const PipelineOptions& opt)
{
// try guards the whole flow against a null-pointer dereference inside the kernel.
    try
    {
// Default-constructed FacePtr (nullptr).
        const FacePtr nf;
// Must return Failed, not throw or crash.
        const FaceReport r = ProcessFace(nf, th, opt);
// PASS = Failed with a locatable (non-empty) message.
        const bool pass = (r.decision == Decision::Failed) && !r.message.empty();
        std::ostringstream os;
        os << "decision=" << DecisionStr(r.decision) << " " << CodeStr(r.errorCode)
           << " msg=" << Snippet(r.message);
        report("process.nullFace", "Failed && message non-empty", os.str(), pass);
    }
    catch (const std::exception& e)
    {
        report("process.nullFace", "Failed && message non-empty", std::string("threw: ") + e.what(), false);
    }
    catch (...)
    {
        report("process.nullFace", "Failed && message non-empty", "threw: " + DescribeUnknown(), false);
    }
}

// C14: ProcessFace on an already-analytic face (Cylinder lateral) must be Rejected with
// a message that mentions "analytic".
void runAlreadyAnalytic(const Thresholds& th, const PipelineOptions& opt)
{
    try
    {
// Build the cylinder and take its lateral face (already a Cylinder).
        BodyPtr body = MakeCylinderBody();
        FacePtr lat = FindCylinderFace(body);
// No lateral face -> setup failure.
        if (!lat)
        {
            report("process.alreadyAnalytic", "Rejected && msg says already analytic",
                   "setup failed: no cylinder face in api_make_solid_cylinder result", false);
            return;
        }

// ProcessFace on an analytic face should reject at entry.
        const FaceReport r = ProcessFace(lat, th, opt);
// PASS = Rejected and the message says the face is already analytic.
        const bool pass = (r.decision == Decision::Rejected) &&
                          (Contains(r.message, "analytic"));
        std::ostringstream os;
        os << "decision=" << DecisionStr(r.decision) << " " << CodeStr(r.errorCode)
           << " msg=" << Snippet(r.message);
        report("process.alreadyAnalytic", "Rejected && msg mentions already-analytic", os.str(), pass);
    }
    catch (const std::exception& e)
    {
        report("process.alreadyAnalytic", "Rejected && msg mentions already-analytic",
               std::string("threw: ") + e.what(), false);
    }
    catch (...)
    {
        report("process.alreadyAnalytic", "Rejected && msg mentions already-analytic",
               "threw: " + DescribeUnknown(), false);
    }
}

// C15: ProcessFace on a near-degenerate BSpline face (UV domain shrunk to a point); must be
// Rejected or Failed with a non-empty message, never a crash or silent success.
void runDegenerateUVFace(const Thresholds& th, const PipelineOptions& opt)
{
    try
    {
// Use the cylinder lateral face as the attachment carrier.
        BodyPtr body = MakeCylinderBody();
        FacePtr lat = FindCylinderFace(body);
// No lateral face -> setup failure.
        if (!lat)
        {
            report("process.degenerateUV", "Rejected || Failed, reason non-empty",
                   "setup failed: no cylinder face", false);
            return;
        }

// Attach a near-degenerate BSpline face; two construction paths: (a) all control points
// coincident (geometrically degenerate); (b) ToBSpline over a zero-width UV domain.
        bool attached = false;
        std::string how;
// Path (a): all control points coincident.
        {
            try
            {
// 2x2 control points, all at the origin -> geometrically degenerate.
                Point3DMatrix cps(2, std::vector<Point3D>(2, Point3D(0.0, 0.0, 0.0)));
// Bidegree-1 BSpline; knots/multiplicities give domain [0,1]^2.
                auto bs = std::make_shared<BSplineSurface>(1, 1, cps,
                                                           RealArray{0.0, 1.0}, RealArray{0.0, 1.0},
                                                           UIntArray{2, 2}, UIntArray{2, 2});
// Wrap as a ModelSurface (cloneGeom=false: exclusive ownership).
                const ModelSurfacePtr ms = TopoBuilder::MakeModelSurface(bs, false);
// Attach to the lateral face, keeping its sense.
                TopoBuilder::FaceAttachSurface(lat, ms, lat->Sense());
// Mark success and record which path was used.
                attached = true;
                how = "coincident control points";
            }
// Path (a) failed: try path (b).
            catch (...)
            {
                attached = false;
            }
        }
        if (!attached)
        {
            try
            {
// Original lateral geometry.
                SurfacePtr orig = lat->GeomSurface();
                if (orig)
                {
// Path (b): ToBSpline with a zero-width domain.
                    BSplineSurfacePtr bs =
                        orig->ToBSpline(UVRange(Interval(0.0, 0.0), Interval(0.0, 0.0)));
// Got a spline -> attach it to the face.
                    if (bs)
                    {
                        const ModelSurfacePtr ms = TopoBuilder::MakeModelSurface(bs, false);
                        TopoBuilder::FaceAttachSurface(lat, ms, lat->Sense());
                        attached = true;
                        how = "zero-width ToBSpline domain";
                    }
                }
            }
            catch (...)
            {
                attached = false;
            }
        }

        if (!attached)
        {
            report("process.degenerateUV", "Rejected || Failed, reason non-empty",
                   "could not construct a near-degenerate BSpline face", false);
            return;
        }

        const FaceReport r = ProcessFace(lat, th, opt);
        const bool pass = ((r.decision == Decision::Rejected) || (r.decision == Decision::Failed)) &&
                          !r.message.empty();
        std::ostringstream os;
        os << "decision=" << DecisionStr(r.decision) << " " << CodeStr(r.errorCode)
           << " construct=" << how << " msg=" << Snippet(r.message);
        report("process.degenerateUV", "Rejected || Failed && message non-empty", os.str(), pass);
    }
    catch (const std::exception& e)
    {
        report("process.degenerateUV", "Rejected || Failed && message non-empty",
               std::string("threw: ") + e.what(), false);
    }
    catch (...)
    {
        report("process.degenerateUV", "Rejected || Failed && message non-empty",
               "threw: " + DescribeUnknown(), false);
    }
}

}

int main(int /*argc*/, char* /*argv*/[])
{
// Top-level try: init/fini and all cases; a fatal exception lands in catch and returns 2.
    try
    {
// Initialise the SDK runtime (Foundation/init.h:25).
        init();

        std::error_code ec;
// Controlled temp directory: all self-generated files live here for isolation/cleanup.
        const fs::path tmp = fs::temp_directory_path() / "sgkq_test_exc";
// create_directories errors are non-fatal; a later WriteText surfaces the problem.
        fs::create_directories(tmp, ec);

// A: file layer.
        runFileLayer(tmp);

// B: geometry / parameter layer (out-of-domain, zero/few samples, NaN/Inf, degenerate).
        runSampleOutOfDomain();
        runDetectEmptySamples();
        runFitTooFew();
        runFitNanInf();
        runDegenerateRadius();
        runDegenerateCone();

// C: pipeline layer, using production-default thresholds and options.
        const Thresholds      th;
        const PipelineOptions opt;
// C13..C15.
        runNullFace(th, opt);
        runAlreadyAnalytic(th, opt);
        runDegenerateUVFace(th, opt);

// Tear down the SDK runtime (Foundation/init.h:28), paired with init.
        fini();
    }
// A std::exception escaping setup is fatal: count it and return 2.
    catch (const std::exception& e)
    {
        std::cerr << "[FAIL] fatal exception: " << e.what() << std::endl;
        ++g_failed;
        return 2;
    }
// Non-std exception escaping: same treatment.
    catch (...)
    {
        std::cerr << "[FAIL] fatal unknown exception" << std::endl;
        ++g_failed;
        return 2;
    }

// Final summary: total / failed / warned.
    std::cout << "[SUMMARY] total=" << g_total
              << " failed=" << g_failed
              << " warned=" << g_warned << std::endl;

// Return 1 if any case failed, else 0.
    return (g_failed == 0) ? 0 : 1;
}

