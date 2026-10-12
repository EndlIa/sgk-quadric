// Stability driver: runs every corpus case for --passes rounds (default 2), requiring each
// id's result to be identical across passes, and writes one CSV row per execution.
// The corpus must hold >= --min-unique-cases (default 1000) *distinct* cases; repeated rounds
// never count toward that floor. Defaults derive from the executable location (sgkq::DefaultDataDir):
// --corpus <root>/data/corpus and --out-dir <root>/data. Exit 0 consistent, 1 any
// inconsistent/caseException/loadFailed, 2 bad args/manifest or too few unique cases.

// Project contract headers: io (read), pipeline (process), types (decisions/thresholds).
#include "sgkquadric/io.hpp"
#include "sgkquadric/testkit.hpp"   // harness: LoadManifest / ManifestEntry / DefaultDataDir
#include "sgkquadric/pipeline.hpp"
#include "sgkquadric/types.hpp"

// SDK runtime: init/fini and runtime class info (IsKindOf<...>).
#include <Foundation/RuntimeClass.h>
#include <Foundation/init.h>

// SDK geometry: Surface base class and the SurfaceType enum.
#include <Geometry/3D/Surface/Surface.h>
#include <Geometry/SurfaceType.h>

// SDK topology: Face/ModelSurface and TopoBuilder (wraps the spline surface as a Face).
#include <Topology/Brep/Face.h>
#include <Topology/Geom/ModelSurface.h>
#include <Topology/Tools/TopoBuilder.h>

// C++ standard library only; no third-party dependencies.
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace sggk;
using namespace sgkq;

// Helpers are translation-unit local and stateless, so plain free functions suffice.
namespace
{

// Path helpers: make one command line work in both WSL and native Windows shells by
// normalizing /mnt/<drive>/... to <DRIVE>:/... (same convention as run_corpus and probe_*).
// All are pure string operations, so they never throw and touch no filesystem state.
//! Normalizes WSL/MSYS /mnt/<d>/rest to Win32 <D>:/rest. Only a single-letter drive is
//! recognized (raw.size() > 5); /mnt/wsl/... and every other form is returned unchanged.
//! Pure string operation, no filesystem access.
std::string NormalizePath(const std::string& raw)
{
    if (raw.size() > 5 && raw.compare(0, 5, "/mnt/") == 0)
    {
        const char drive = raw[5];
        if (std::isalpha(static_cast<unsigned char>(drive)))
        {
            const std::string rest = (raw.size() > 6) ? raw.substr(6) : std::string("/");
            std::string       out;
            out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(drive))));
            out += ':';
            out += rest;
            return out;
        }
    }
    return raw;
}

//! True for a POSIX absolute path ('/') or a drive form (':' at index 1); empty -> false.
//! Literal check only, so it is valid for paths that do not exist yet.
bool IsAbsolute(const std::string& p)
{
    return (!p.empty() && p[0] == '/') || (p.size() >= 2 && p[1] == ':');
}

//! Joins dir and leaf with exactly one '/'; a dir already ending in '/' or backslash is not
//! given another, so no "//" appears. An empty dir yields leaf unchanged.
std::string JoinPath(const std::string& dir, const std::string& leaf)
{
    std::string base = dir;
    if (!base.empty() && base.back() != '/' && base.back() != '\\')
        base.push_back('/');
    return base + leaf;
}

//! Roots a relative manifest sgtPath at the corpus dir; an absolute sgtPath is returned
//! as-is. Existence is not checked here (LoadSgt reports it).
std::string ResolveSgt(const std::string& corpusDir, const std::string& sgtPath)
{
    if (IsAbsolute(sgtPath))
        return sgtPath;
    return JoinPath(corpusDir, sgtPath);
}
// Type/decision helpers: map SDK enums to project enums/strings for the CSV and for the
// cross-pass comparison. All are pure mappings (no locale, no time).

//! Maps an SDK SurfaceType to this project's Canonical; non-analytic types (BSpline/Nurbs
//! etc.) map to Unknown. Must stay in sync with pipeline.cpp::CanonicalOf so targetType matches.
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

//! True for the five analytic types (Plane/Cylinder/Cone/Sphere/Torus). Same set as
//! CanonicalOf; done locally to short-circuit already-analytic faces before ProcessFace.
bool IsAnalyticSurfaceType(SurfaceType t)
{
    return t == SurfaceType::Plane || t == SurfaceType::Cylinder || t == SurfaceType::Cone
        || t == SurfaceType::Sphere || t == SurfaceType::Torus;
}

//! Decision -> "Converted"/"Rejected"/"Failed"; unknown -> "Unknown". Must be a pure
//! mapping (no locale/time), else cross-pass string comparison would flag false differences.
const char* DecisionStr(Decision d)
{
    switch (d)
    {
    case Decision::Converted: return "Converted";
    case Decision::Rejected:  return "Rejected";
    case Decision::Failed:    return "Failed";
    default:                  return "Unknown";
    }
    // Unknown and future enum values collapse to Unknown, keeping the string space closed.
}

//! Picks the face to process: the first BSpline face (the converter's target), else the
//! first non-null face. Same strategy as run_corpus::PickFace so both tools choose alike.
const SurfacePtr* PickFace(const std::vector<InputFace>& faces)
{
    for (const InputFace& f : faces)
    // Two-pass scan: prefer BSpline (needs conversion), then fall back to any non-null face.
        if (f.surface && f.surface->SurfaceType() == SurfaceType::BSpline)
            return &f.surface;
    for (const InputFace& f : faces)
        if (f.surface)
            return &f.surface;
    return nullptr;
}

//! Two doubles count as consistent if bitwise equal or |a-b| < 1e-12. NaN vs NaN counts as
//! equal so a pass that produced no value is not flagged unstable. 1e-12 is an absolute length
//! tolerance (mm), tiny against model features (>=1e-3 mm), so it only absorbs last-bit noise.
bool SameDouble(double a, double b)
{
    if (std::isnan(a) && std::isnan(b)) return true;
    if (a == b) return true;
    return std::fabs(a - b) < 1e-12;
}

// One execution record, one field per CSV column.
struct Row
{
    int         round = 0;
    std::string id;
    std::string type;
    int         ctrlCount = 0;
    std::string decision;
    std::string targetType;
    double      maxDist = 0.0; // max normal deviation [mm]
    double      boundaryCurveHausdorff = 0.0; // boundary Hausdorff [mm]
    double      g1Before = 0.0;
    double      g1After = 0.0;
    int         errorCode = 0;
    long long   elapsed_us = 0; // per-case time [us]
// Row assembly lives in main (easy to fill defaults on any branch); AppendRow owns the
// column order, so it has a single source of truth.
};

// Per-id baseline from the first round; only the fields that are compared are kept.
struct Ref
{
    int         round = 0;
    std::string decision;
    std::string targetType;
    int         errorCode = 0;
    double      maxDist = 0.0;
    double      boundaryCurveHausdorff = 0.0;
// Ref is deliberately a subset of Row so the compared fields copy directly.
};

// 17 significant digits is a double's round-trip precision; fewer would lose bits through
// the CSV and defeat the bitwise-equality check. Consume the CSV with a tool, not by eye.
// The CSV is meant for scripts/tools, not eyeballing.
std::string Num(double v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

//! Writes one record in fixed column order (must match the header); numbers go through Num().
//! The line ends with LF (never platform CRLF) so artifacts are byte-identical across platforms.
void AppendRow(std::ostream& os, const Row& r)
{
    os << r.round << ','
       << r.id << ','
       << r.type << ','
       << r.ctrlCount << ','
       << r.decision << ','
       << r.targetType << ','
       << Num(r.maxDist) << ','
       << Num(r.boundaryCurveHausdorff) << ','
       << Num(r.g1Before) << ','
       << Num(r.g1After) << ','
       << r.errorCode << ','
       << r.elapsed_us << '\n';
}

// Usage goes to stderr, keeping stdout parseable.
void PrintUsage()
{
    std::cerr << "usage: batch_stability [--corpus <dir>] [--min-unique-cases N] "
                 "[--passes P] [--rounds P] [--out-dir <dir>]\n";
}

} // namespace

// Entry: parse args -> try { init, load manifest, run rounds, summarize } -> exit code.
// Only argument/init failures return directly; per-case failures are aggregated so the batch
// never aborts on one bad case. No mutable global state (besides SDK init/fini).
// Exit: 0 ok; 1 inconsistent/exception/loadFailed; 2 bad args or init failure.
int main(int argc, char* argv[])
{
    // Default corpus dir matches gen_corpus's output (<root>/data/corpus, derived from the
    // executable location via sgkq::DefaultDataDir), so a zero-arg rerun is cwd-independent.
    std::string corpusDir = sgkq::DefaultDataDir() + "/corpus";
    std::string outDir    = sgkq::DefaultDataDir();
    // passesOpt: execution rounds, default 2; --rounds is a compatible alias for --passes.
    long long   passesOpt = 2;
    // minUniqueOpt: minimum distinct corpus cases, default 1000; rounds never substitute for it.
    long long   minUniqueOpt = 1000;

    // ---- argument parsing ----
    // Value options consume the next argv via next(); a missing value exits 2.
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        // next(): reads the option value (--opt val form only; --opt=val is not accepted); a missing
        // value prints usage and exits 2. Captures i by reference so ++i advances the outer loop.
        auto next = [&](const char* name) -> const char*
        {
            if (i + 1 >= argc)
            {
                std::cerr << "[ERR] missing value for " << name << "\n";
                PrintUsage();
                std::exit(2);
            }
            return argv[++i];
        };

        if (a == "--corpus")
            corpusDir = NormalizePath(next("--corpus"));
        else if (a == "--out-dir")
            outDir = NormalizePath(next("--out-dir"));
        else if (a == "--passes")
            passesOpt = std::atoll(next("--passes"));
        else if (a == "--rounds")
            passesOpt = std::atoll(next("--rounds"));
        else if (a == "--min-unique-cases")
            minUniqueOpt = std::atoll(next("--min-unique-cases"));
        else
        {
            std::cerr << "[ERR] unknown argument: " << a << "\n";
            PrintUsage();
            return 2;
        }
    }

    corpusDir = NormalizePath(corpusDir);
    outDir    = NormalizePath(outDir);

    // Manifest lives in the corpus dir and the CSV in the out dir; both names are fixed.
    const std::string manifestPath = JoinPath(corpusDir, "manifest.jsonl");
    const std::string csvPath      = JoinPath(outDir, "batch_stability.csv");

    int exitCode = 0;

    // Top-level try: fatal exceptions not caught per case reach these catches and return 2.
    // std::exception prints what(); anything else prints a fixed string.
    try
    {
    // init() (Foundation/init.h:25) must be paired with fini() below.
        init();

        // ---- 1) read manifest ----
    // A missing/corrupt or empty manifest is an init error (exit 2): an empty corpus would
    // otherwise look like a clean pass. Entries run in manifest order, keeping the run reproducible.
        std::vector<ManifestEntry> entries;
        if (!LoadManifest(manifestPath, entries))
        {
            std::cerr << "[ERR] failed to load manifest: " << manifestPath << "\n";
            fini();
            return 2;
        }
        if (entries.empty())
        {
            std::cerr << "[ERR] manifest is empty: " << manifestPath << "\n";
            fini();
            return 2;
        }

    // uniqueCases = number of manifest entries (one entry = one distinct case).
        const long long cases = static_cast<long long>(entries.size());

    // Enforce the distinct-case floor before running; never pad with repeated rounds (the
    // requirement is >=1000 distinct faces).
        if (cases < minUniqueOpt)
        {
            std::cerr << "[FAIL] corpus has only " << cases << " unique cases, need >= "
                      << minUniqueOpt
                      << "; regenerate with gen_corpus --cases-per-type ...\n";
            fini();
            return 2;
        }

    // passes rounds; a value < 1 falls back to 1.
        long long passes = passesOpt;
        if (passes < 1)
            passes = 1;

    // CSV is truncated per run; failing to open it is an init error.
        std::ofstream csv(csvPath, std::ios::out | std::ios::trunc);
        if (!csv)
        {
            std::cerr << "[ERR] cannot open CSV for writing: " << csvPath << "\n";
            fini();
            return 2;
        }
    // Header order must match AppendRow exactly or consumers mis-map columns.
        csv << "round,id,type,ctrlCount,decision,targetType,maxDist,boundaryCurveHausdorff,"
               "g1Before,g1After,errorCode,elapsed_us\n";

        std::cout << "[INFO] corpus=" << corpusDir << " manifest=" << manifestPath
                  << " out=" << csvPath << " uniqueCases=" << cases << " passes=" << passes
                  << " totalExecutions=" << (passes * cases) << " minUniqueCases=" << minUniqueOpt
                  << "\n";

    // Keep PipelineOptions/Thresholds at their defaults (nu=nv=9, nuAcc=nvAcc=33, writeBack=true):
    // the defaults are the production path, so this stays aligned with the integration tests.
        const PipelineOptions opts;
        const Thresholds      th{};

    // Thresholds default-constructed (the three detection/fit/accept thresholds).
        std::map<std::string, Ref> refs; // id -> first-round baseline

    // Counters use long long: totalExecutions = passes x uniqueCases (>= 2000) can overflow int.
    // converted/rejected/failedState, loadFailed/skippedAnalytic, inconsistent/caseExceptions.
        long long totalExecutions = 0;
        long long converted = 0, rejected = 0, failedState = 0;
        long long loadFailed = 0, skippedAnalytic = 0;
        long long inconsistent = 0, caseExceptions = 0;

    // steady_clock start for total timing; rt0/ct0 time each round and case.
        const auto t0 = std::chrono::steady_clock::now();

        // ---- 2) run every round and case ----
    // Each round reruns every case independently (see the determinism note above).
        for (long long r = 1; r <= passes; ++r)
        {
            long long roundOk = 0, roundFailed = 0;
            const auto rt0 = std::chrono::steady_clock::now();

        // Case loop in manifest order; roundOk/roundFailed are reset each round.
        // range-for by const ref avoids copies and accidental manifest mutation.
            for (const ManifestEntry& entry : entries)
            {
                ++totalExecutions;

        // Record metadata; round is 1-based for readability. ctrlCount is informational only.
                Row row;
                row.round     = static_cast<int>(r);
                row.id        = entry.id;
                row.type      = entry.type;
                row.ctrlCount = entry.ctrlCount;

                const auto ct0 = std::chrono::steady_clock::now();

                std::string sgtPath;
                try
                {
                    sgtPath = ResolveSgt(corpusDir, entry.sgtPath);
                }
                catch (...)
                {
                // Fall back to the raw manifest sgtPath if resolution throws.
                    sgtPath = entry.sgtPath;
                }

        // Per-case catch: any exception is logged as [CASE-EXCEPTION] and the run continues.
                try
                {
            // Reload every round: ProcessFace with writeBack replaces the face surface, so reusing the
            // same object would feed round N's output into round N+1 (measured behavior).
                    const LoadResult lr = LoadSgt(sgtPath);
            // LoadFailed: keep io.hpp's errorCode and message for traceability.
                    if (!lr.ok)
                    {
                        ++loadFailed;
                        ++roundFailed;
                        row.decision   = "LoadFailed";
                        row.targetType = "Unknown";
                        row.errorCode  = lr.errorCode;
                        std::cout << "[LOADFAIL] id=" << entry.id
                                  << " code=" << lr.errorCode
                                  << " msg=" << lr.message << "\n";
                    }
                    else
                    {
                // PickFace; a null pointer or null surface is treated as an empty face.
                        const SurfacePtr* picked = PickFace(lr.faces);
                    // No usable face -> reuse LoadErrEmptyFaces for a consistent error code.
                        if (picked == nullptr || !*picked)
                        {
                            ++loadFailed;
                            ++roundFailed;
                            row.decision   = "LoadFailed";
                            row.targetType = "Unknown";
                            row.errorCode  = static_cast<int>(LoadErrEmptyFaces);
                            std::cout << "[LOADFAIL] id=" << entry.id
                                      << " code=" << static_cast<int>(LoadErrEmptyFaces)
                                      << " msg=" << sgtPath << ": no usable face\n";
                        }
                        else
                        {
                    // Already analytic -> Skipped, no ProcessFace (errorCode 2 = pipeline's kAlreadyAnalytic).
                            const SurfaceType st = (*picked)->SurfaceType();
                            if (IsAnalyticSurfaceType(st))
                            {
                                ++skippedAnalytic;
                                ++roundOk;
                                row.decision   = "Skipped";
                                row.targetType = ToString(CanonicalOf(st));
                                row.errorCode  = 2;
                            }
                            else
        // Each case sits in its own try so one failure cannot abort the batch.
                            {
                                // LoadSgt yields InputFace (SurfacePtr only); wrap it in a topological Face for ProcessFace.
                                const ModelSurfacePtr ms = TopoBuilder::MakeModelSurface(*picked, false);
                                const FacePtr         face = TopoBuilder::MakeFace(ms, true);

                        // ProcessFace with default options; face is temporary, so writeBack affects only this round.
                                const FaceReport rep = ProcessFace(face, th, opts);

                        // Copy report fields into the row; units pass through unchanged (mm / radians).
                                row.decision          = DecisionStr(rep.decision);
                                row.targetType        = ToString(rep.params.type);
                                row.maxDist           = rep.accuracy.maxDist;
                                row.boundaryCurveHausdorff = rep.accuracy.boundaryCurveHausdorff;
                                row.g1Before          = rep.accuracy.g1Before;
                                row.g1After           = rep.accuracy.g1After;
                                row.errorCode         = rep.errorCode;

                        // Converted/Rejected count as round success; Failed/unknown count against it.
                                switch (rep.decision)
                                {
                                case Decision::Converted: ++converted;   ++roundOk;     break;
                                case Decision::Rejected:  ++rejected;    ++roundOk;     break;
                                default:                  ++failedState; ++roundFailed; break;
                                }
                            }
                        }
                    }
                }
        // std::exception from a case: count and continue, printing id + what().
                catch (const std::exception& e)
                {
                    ++caseExceptions;
                    ++roundFailed;
                    row.decision   = "Exception";
                    row.targetType = "Unknown";
                    row.errorCode  = -1;
                    std::cout << "[CASE-EXCEPTION] id=" << entry.id
                              << " what=" << e.what() << "\n";
                }
        // Non-std exception: same handling with a placeholder description.
                catch (...)
                {
                    ++caseExceptions;
                    ++roundFailed;
                    row.decision   = "Exception";
                    row.targetType = "Unknown";
                    row.errorCode  = -1;
                    std::cout << "[CASE-EXCEPTION] id=" << entry.id
                              << " what=unknown exception\n";
                }

        // Per-case wall time in microseconds, covering path resolution plus processing.
                row.elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::steady_clock::now() - ct0).count();

        // Cross-pass consistency: for each id, decision/targetType/errorCode must be exactly equal,
        // and maxDist/boundaryCurveHausdorff bitwise equal or within 1e-12 mm. The first (or first-seen)
        // round is the baseline; a mismatch is counted and printed but does not stop the batch, so
        // every round still runs and the stability picture stays complete.
                // ---- 3) cross-pass consistency, by id ----
        // Per-id baseline: create on first sight, otherwise compare field by field. std::map (not
        // unordered_map) gives ordered, hash-seed-free iteration, which aids reproducible output.
                auto it = refs.find(entry.id);
                if (it == refs.end())
                {
            // New id: record the six compared fields as the baseline.
                    Ref ref;
                    ref.round             = row.round;
                    ref.decision          = row.decision;
                    ref.targetType        = row.targetType;
                    ref.errorCode         = row.errorCode;
                    ref.maxDist           = row.maxDist;
                    ref.boundaryCurveHausdorff = row.boundaryCurveHausdorff;
                    refs.emplace(entry.id, ref);
                }
                else
                {
            // Compare against the baseline: exact for discrete fields, SameDouble for continuous ones.
                    const Ref& ref = it->second;
                    if (row.decision != ref.decision)
                    {
                        ++inconsistent;
                        // Print id, field, both round numbers and both values for one-shot triage.
                        std::cout << "[INCONSISTENT] id=" << entry.id
                                  << " field=decision round=" << ref.round << " vs " << row.round
                                  << ": " << ref.decision << " vs " << row.decision << "\n";
                    }
                    if (row.targetType != ref.targetType)
                    {
                        // targetType is a Canonical name; a change means the detected/fit target changed.
                        ++inconsistent;
                        std::cout << "[INCONSISTENT] id=" << entry.id
                                  << " field=targetType round=" << ref.round << " vs " << row.round
                                  << ": " << ref.targetType << " vs " << row.targetType << "\n";
                    }
                    if (row.errorCode != ref.errorCode)
                    {
                        // errorCode drift usually tracks non-determinism in fitting/acceptance.
                        ++inconsistent;
                        std::cout << "[INCONSISTENT] id=" << entry.id
                                  << " field=errorCode round=" << ref.round << " vs " << row.round
                                  << ": " << ref.errorCode << " vs " << row.errorCode << "\n";
                    }
                    if (!SameDouble(row.maxDist, ref.maxDist))
                        // Num() prints 17 digits so last-bit differences are visible here.
                    {
                        ++inconsistent;
                        std::cout << "[INCONSISTENT] id=" << entry.id
                                  << " field=maxDist round=" << ref.round << " vs " << row.round
                                  << ": " << Num(ref.maxDist) << " vs " << Num(row.maxDist) << "\n";
                    }
                    if (!SameDouble(row.boundaryCurveHausdorff, ref.boundaryCurveHausdorff))
                        // boundaryCurveHausdorff likewise; both values are printed for cross-checking.
                    {
                        ++inconsistent;
                        std::cout << "[INCONSISTENT] id=" << entry.id
                                  << " field=boundaryCurveHausdorff round=" << ref.round << " vs " << row.round
                                  << ": " << Num(ref.boundaryCurveHausdorff) << " vs "
                                  << Num(row.boundaryCurveHausdorff) << "\n";
                    }
                }

        // Every case writes exactly one CSV row even on failure (Exception, errorCode=-1), and
        // regardless of consistency, so the file can be reviewed after the fact.
                AppendRow(csv, row);
            }

        // Round summary: ok/failed and wall time (ms, 3 decimals); observational only.
            const double roundMs = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - rt0).count();
            std::cout << "[ROUND] r=" << r
                      << " cases=" << cases
                      << " ok=" << roundOk
                      << " failed=" << roundFailed
                      << " elapsed_ms=" << std::fixed << std::setprecision(3) << roundMs
                      << std::defaultfloat << "\n";

        // [MEM] is only a time alignment point for an external memory sampler (this program never
            // measures memory itself).
            std::cout << "[MEM] round=" << r << " note=external-sampled\n";
        }

        // Explicit flush+close so the CSV is complete on disk before any early return/exception.
        csv.flush();
        csv.close();

        const double totalS = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
                                // Empty-face failures reuse LoadErrEmptyFaces to match LoadSgt's error codes.

        // Summary then assertions: the [ASSERT] booleans are derived only from the in-memory
        // counters (O(1), no CSV rescan). One key per line keeps the output machine-parseable.
        // ---- 4) summary and assertions ----
        // Totals plus total wall time in seconds.
        std::cout << "[SUMMARY] uniqueCases=" << cases
                  << " passes=" << passes
                  << " totalExecutions=" << totalExecutions
                  << " converted=" << converted
                  << " rejected=" << rejected
                  << " failedState=" << failedState
                  << " loadFailed=" << loadFailed
                  << " skippedAnalytic=" << skippedAnalytic
                  << " inconsistent=" << inconsistent
                  << " caseExceptions=" << caseExceptions
                  << " totalElapsed_s=" << std::fixed << std::setprecision(3) << totalS
                  << std::defaultfloat << "\n";
                                // Skipped uses errorCode 2, matching pipeline.cpp's kAlreadyAnalytic.

        // [ASSERT] flags: noCrash (no case exceptions), hasInconsistency (a positive marker, not
        // "passed"), uniqueCases_ge_1000 and totalExecutions_ge_1000.
        const bool noCrash           = (caseExceptions == 0); // stability gate
        const bool hasInconsistency  = (inconsistent > 0); // consistency gate
        const bool ge1000            = (totalExecutions >= 1000); // scale (total executions)
        const bool uniqueGe1000      = (cases >= 1000); // >=1000 distinct cases
        std::cout << "[ASSERT] noCrash=" << (noCrash ? "true" : "false")
                  << "/noLeaks:external"
                  << "/hasInconsistency=" << (hasInconsistency ? "true" : "false")
                  << "/uniqueCases_ge_1000=" << (uniqueGe1000 ? "true" : "false")
                  << "/totalExecutions_ge_1000=" << (ge1000 ? "true" : "false")
                  << "\n";

        // Exit 1 if any inconsistent/exception/loadFailed; skippedAnalytic, like converted/rejected,
        // is not a failure.
        exitCode = (inconsistent > 0 || caseExceptions > 0 || loadFailed > 0) ? 1 : 0;

        fini(); // pairs with init()
    }
    catch (const std::exception& e)
    {
        std::cerr << "[FAIL] fatal exception: " << e.what() << "\n";
        return 2;
    }
    catch (...)
    {
        std::cerr << "[FAIL] fatal unknown exception\n";
        return 2;
    }

    return exitCode;
}
