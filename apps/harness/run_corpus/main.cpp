// Detects a corpus of NURBS faces and prints a 6x6 confusion matrix ([CONFUSION]/[SUMMARY]).
// Defaults: --corpus <project root>/data/corpus (root derived from the executable location),
// --manifest <corpus>/manifest.jsonl, nu=nv=9 (81 samples, matching test_detect resolution).
// Paths are normalized so /mnt/d/... and D:/... are equivalent.
// Exit codes: 0 all matched; 1 any mismatch or load failure; 2 bad args / manifest error.
#include "sgkquadric/detect.hpp"
#include "sgkquadric/io.hpp"
#include "sgkquadric/testkit.hpp"   // harness: LoadManifest / ManifestEntry / DefaultDataDir

#include <Foundation/init.h>
#include <GeomBase/Interval.h>
#include <GeomBase/UVRange.h>
#include <Geometry/3D/Surface/Surface.h>
#include <Geometry/SurfaceType.h>

#include <cctype>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace sggk;
using namespace sgkq;

namespace
{

//! Normalizes WSL /mnt/<drive>/rest to <DRIVE>:/rest; the SDK I/O layer expects the X:/ prefix.
//! Only a single-letter drive qualifies (raw.size() > 5); /mnt/wsl/..., already-Windows and
//! relative paths are returned unchanged, so the rest of the code only sees one form.
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

//! True for a Unix root ('/') or a drive form ('X:' at index 1). Literal prefix check only,
//! no filesystem access, so it is valid for paths that do not exist; empty is not absolute.
//! Used to decide whether a manifest sgtPath still needs the corpus dir prepended.
bool IsAbsolute(const std::string& p)
{
    return (!p.empty() && p[0] == '/') || (p.size() >= 2 && p[1] == ':');
}

//! Roots a relative manifest sgtPath at the corpus dir (separator '/', accepted by the
//! Windows API too); an absolute sgtPath is returned as-is. Purely textual: file existence
//! is left to the subsequent LoadSgt, which reports its own error code.
std::string ResolveSgt(const std::string& corpusDir, const std::string& sgtPath)
{
    if (IsAbsolute(sgtPath))
        return sgtPath;
    std::string base = corpusDir;
    if (!base.empty() && base.back() != '/' && base.back() != '\\')
        base.push_back('/');
    return base + sgtPath;
}

//! Maps a Canonical to its confusion-matrix index: Plane/Cylinder/Cone/Sphere/Torus (0..4),
//! everything else including Unknown -> 5. This order must stay in sync with names[] below,
//! or the [CONFUSION] row/column labels shift.
//! Unknown class is index 5.
int CanonIndex(Canonical c)
{
    switch (c)
    {
    case Canonical::Plane:    return 0;
    case Canonical::Cylinder: return 1;
    case Canonical::Cone:     return 2;
    case Canonical::Sphere:   return 3;
    case Canonical::Torus:    return 4;
    default:                  return 5;
    }
}

//! Picks the face to detect: the first BSpline face (corpus targets are always spline faces);
//! if there is none, it falls back to the first non-null face only to avoid an empty run
//! (that fallback face is not necessarily the intended target).
//! Returns a pointer into faces, so the caller must keep faces alive while using it.
const SurfacePtr* PickFace(const std::vector<InputFace>& faces)
{
    for (const InputFace& f : faces)
    {
        if (f.surface && f.surface->SurfaceType() == SurfaceType::BSpline)
            return &f.surface;
    }
    for (const InputFace& f : faces)
    {
        if (f.surface)
            return &f.surface;
    }
    return nullptr;
}

//! Usage goes to stderr so stdout stays machine-parseable ([SUMMARY] etc.).
void PrintUsage()
{
    std::cerr << "usage: run_corpus [--corpus <dir>] [--manifest <path>] [--nu N] [--nv N]\n";
}

}

//! Batch corpus detection entry: read manifest -> per-entry LoadSgt+Detect -> compare to
//! expected -> confusion matrix. nu/nv are resample counts used only for mismatch diagnostics.
//! init()/fini() are paired; the early exits for a bad or empty manifest call fini() before
//! returning 2. Single-threaded, no shared concurrent state.
//! Returns 0 all matched, 1 mismatch or load failure, 2 bad args / manifest error.
int main(int argc, char* argv[])
{
    // Defaults are derived from the executable location (sgkq::DefaultDataDir): <root>/data/corpus,
    // so a zero-argument run works from any cwd; --corpus overrides it. --manifest defaults to
    // <corpus>/manifest.jsonl, and nu=nv=9 (81 samples) keeps test_detect's detection resolution.
    std::string corpusDir = sgkq::DefaultDataDir() + "/corpus";
    std::string manifestPath;
    int         nu = 9;
    int         nv = 9;

    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        // next(): reads the option value (--opt val form only) and prints usage then exits 2 if it
        // is missing; captures i by reference so ++i advances the outer loop.
        // --opt=val is not accepted.
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
        {
            corpusDir = NormalizePath(next("--corpus"));
        }
        else if (a == "--manifest")
        {
            manifestPath = NormalizePath(next("--manifest"));
        }
        else if (a == "--nu")
        {
            nu = std::atoi(next("--nu"));
        }
        else if (a == "--nv")
        {
            nv = std::atoi(next("--nv"));
        }
        else
        {
            std::cerr << "[ERR] unknown argument: " << a << "\n";
            PrintUsage();
            return 2;
        }
    }

    // Normalize the defaults too: --corpus may have overwritten them during parsing.
    corpusDir = NormalizePath(corpusDir);
    if (manifestPath.empty())
        manifestPath = ResolveSgt(corpusDir, "manifest.jsonl");
    else
        manifestPath = NormalizePath(manifestPath);

    int resultCode = 0;
    try
    {
        init();

        // A failed or empty manifest is fatal (fini() first, then exit 2): an empty corpus would
        // otherwise report "0/0 matched", a false pass. Entries follow manifest order, which keeps
        // the run deterministic.
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

        std::cout << "[INFO] corpus=" << corpusDir << " manifest=" << manifestPath
                  << " entries=" << entries.size() << "\n";

        // conf[expected][actual], both indexed by CanonIndex (0..5), so the table shows which
        // classes get confused rather than only a pass/fail count.
        // (Uninitialized = zero counts.)
        int conf[6][6] = {{0}};

        int loadFailed = 0;
        int matched    = 0;
        int total      = 0;

        // Per entry: load the .sgt, pick the target spline face, Detect over that face's own
        // parameter domain, then compare the detected type against the manifest's expected type.
        // A single failure never aborts the loop.
        for (const ManifestEntry& entry : entries)
        {
            ++total;

            const std::string sgtPath = ResolveSgt(corpusDir, entry.sgtPath);
            const LoadResult  lr      = LoadSgt(sgtPath);
            if (!lr.ok)
            {
                ++loadFailed;
                std::cout << "[LOADFAIL] id=" << entry.id
                          << " code=" << lr.errorCode
                          << " msg=" << lr.message << "\n";
                continue;
            }

            const SurfacePtr* picked = PickFace(lr.faces);
            if (picked == nullptr || !*picked)
            {
                ++loadFailed;
                // No usable face: reuse LoadErrEmptyFaces so the code matches the ingest layer's
                // error-code semantics.
                std::cout << "[LOADFAIL] id=" << entry.id
                          << " code=" << static_cast<int>(LoadErrEmptyFaces)
                          << " msg=" << sgtPath << ": no usable face\n";
                continue;
            }
            const SurfacePtr& srf = *picked;

            // Sample/detect over the spline face's own parameter domain (U in [0,1] after ToBSpline),
            // not the trimmed domain: a trim boundary is not a geometric degeneracy.
            const UVRange uv(srf->DomainU(), srf->DomainV());

            const DetectResult r = Detect(*srf, uv);

            const int ei = CanonIndex(entry.expected);
            const int ai = CanonIndex(r.candidate);
            ++conf[ei][ai];

            const bool match = (r.candidate == entry.expected);
            if (match)
            {
                ++matched;
            }
            else
            {
                // On a mismatch, resample once for diagnostics (sample/degenerate counts, UV range) to
                // tell an algorithm error apart from degenerate input geometry or a bad domain.
                const std::vector<Sample> smp = SampleSurface(*srf, uv, nu, nv);
                int                       degen = 0;
                for (const Sample& s : smp)
                {
                    if (s.degen) ++degen;
                }
                const std::string firstNote = r.notes.empty() ? std::string("(no notes)") : r.notes.front();

                std::cout << "[MISMATCH] id=" << entry.id
                          << " expected=" << ToString(entry.expected)
                          << " actual=" << ToString(r.candidate)
                          << " residual=" << std::scientific << std::setprecision(6) << r.residual
                          << " notes=" << firstNote
                          << " [samples=" << smp.size() << " degen=" << degen
                          << " uv=(" << uv.IntervalU().Min() << "," << uv.IntervalU().Max()
                          << ")x(" << uv.IntervalV().Min() << "," << uv.IntervalV().Max() << ")]"
                          << std::defaultfloat << "\n";
            }
        }

        // Confusion table: rows = expected, cols = actual, diagonal = correct. names[] must match
        // the CanonIndex order or the labels shift; the matrix is row-major expected x actual.
        // Diagonal entries are correct decisions.
        const char* names[6] = {"Plane", "Cylinder", "Cone", "Sphere", "Torus", "Unknown"};

        std::cout << "[CONFUSION] rows=expected cols=actual\n";
        std::cout << std::setw(10) << std::left << "expected" << std::right;
        for (int c = 0; c < 6; ++c)
            std::cout << std::setw(10) << names[c];
        std::cout << "\n";
        for (int r = 0; r < 6; ++r)
        {
            std::cout << std::setw(10) << std::left << names[r] << std::right;
            for (int c = 0; c < 6; ++c)
                std::cout << std::setw(10) << conf[r][c];
            std::cout << "\n";
        }

        // accuracy is a percentage; total==0 is unreachable here (early exits above).
        const double accuracy = (total > 0) ? (100.0 * matched / total) : 0.0;
        std::cout << "[SUMMARY] total=" << total
                  << " matched=" << matched
                  << " accuracy=" << std::fixed << std::setprecision(2) << accuracy << "%"
                  << " loadFailed=" << loadFailed << "\n";

        // Exit code: 0 only when nothing failed to load and everything matched; any mismatch or
        // load failure yields 1.
        resultCode = (matched == total && loadFailed == 0) ? 0 : 1;

        fini();
    }
    catch (const std::exception& e)
    {
        // Unexpected exception: exit 2 (distinct from 1 for a mismatch); fini() is intentionally
        // not called on this path.
        std::cerr << "[FAIL] exception: " << e.what() << "\n";
        return 2;
    }
    catch (...)
    {
        std::cerr << "[FAIL] unknown exception\n";
        return 2;
    }

    return resultCode;
}
