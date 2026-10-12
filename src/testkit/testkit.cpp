// src/testkit/testkit.cpp - harness-only code moved out of the library's I/O layer (io.cpp):
// the manifest.jsonl parser, the type-string -> Canonical mapping it uses, and the
// executable-relative default data dir. Built as the sgkqtestkit static library; only the batch
// and benchmark tools link it. Nothing here parses CAD data.
#include "sgkquadric/testkit.hpp"

#if defined(_WIN32)
#include <windows.h>   // GetModuleFileNameA: full path of this executable (see DefaultDataDir)
#elif defined(__linux__)
#include <unistd.h>    // readlink("/proc/self/exe"): full path of this executable
#endif

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace sgkq
{
namespace
{
    //! Maps a type string to Canonical; unrecognized (incl. empty) -> Unknown, which is also how
    //! the manifest marks negative examples. Case-sensitive and order-independent: the literals
    //! must match the corpus generator's emitted type/expected strings. This avoids coupling the
    //! string format to the enum ordering (i.e. it is not the inverse of ToString).
    //! @return the matching enum, or Canonical::Unknown for anything unrecognized.
    Canonical CanonicalFromString(const std::string& s)
    {
        if (s == "Plane")    return Canonical::Plane;
        if (s == "Cylinder") return Canonical::Cylinder;
        if (s == "Cone")     return Canonical::Cone;
        if (s == "Sphere")   return Canonical::Sphere;
        if (s == "Torus")    return Canonical::Torus;
        return Canonical::Unknown;
    }

    // Minimal JSON field extraction, for manifest.jsonl only: each line is a single-line,
    // non-nested object. Nesting, multi-line objects, comments and backslash-u escapes are
    // unsupported. On failure these return false rather than an error code; LoadManifest turns
    // that into a whole-manifest rejection (false + cleared out).
    // Value shapes handled: "..." -> string, [...] -> number array, otherwise -> number.

    //! Returns the offset just after the ':' for "key" (leading whitespace may remain; callers
    //! run SkipWs), or npos. The key is quoted so it cannot match a substring (e.g. "id" vs
    //! "grid"), and the ':' is searched only after the key so a ':' inside a value cannot be
    //! misplaced for the separator. Only the first occurrence of the key is found.
    //! @param line a single JSON line; @param key the bare key name (no quotes).
    //! @return the offset after the ':' (leading whitespace not yet skipped), or npos.
    size_t FindKey(const std::string& line, const std::string& key)
    {
        const std::string pat = "\"" + key + "\"";
        size_t            p   = line.find(pat);
        if (p == std::string::npos)
            return std::string::npos;
        p = line.find(':', p + pat.size());
        if (p == std::string::npos)
            return std::string::npos;
        return p + 1;
    }

    //! Advances p past the four JSON whitespace characters (space/tab/CR/LF); other control
    //! characters are not skipped. Kept separate from FindKey so it can also be called between
    //! array elements.
    //! @param[in,out] p cursor, advanced to the first non-whitespace character or s.size().
    void SkipWs(const std::string& s, size_t& p)
    {
        // Advance until end of string or the first non-whitespace character.
        while (p < s.size() && (s[p] == ' ' || s[p] == '\t' || s[p] == '\r' || s[p] == '\n'))
            ++p;
    }

    //! Extracts a string value that starts and closes with '"' on the same line (a value that
    //! spans lines is rejected). Handles the six common backslash escapes; any other keeps the
    //! character without the backslash (the manifest never emits backslash-u escapes). An empty
    //! value is valid and is distinguished from "key missing" by the return value.
    //! @param out the unescaped value, valid only while the function returns true.
    bool JsonExtractString(const std::string& line, const std::string& key, std::string& out)
    {
        size_t p = FindKey(line, key);
        if (p == std::string::npos)
            return false;
        SkipWs(line, p);
        if (p >= line.size() || line[p] != '"')
            return false;
        ++p;
        std::string v;
        while (p < line.size() && line[p] != '"')
        {
            // An unescaped '"' ends the value; a backslash starts an escape sequence.
            if (line[p] == '\\' && p + 1 < line.size())
            {
                ++p;
                const char c = line[p];
                switch (c)
                {
                case 'n':  v += '\n'; break;
                case 't':  v += '\t'; break;
                case 'r':  v += '\r'; break;
                case '"':  v += '"';  break;
                case '\\': v += '\\'; break;
                case '/':  v += '/';  break;
                default:   v += c;    break;
                }
            }
            else
            {
                v += line[p];
            }
            ++p;
        }
        if (p >= line.size())   // unterminated value -> fail
            return false;
        out = v;
        return true;
    }

    //! Extracts a numeric literal: scans the characters [0-9 + - . e E] then converts with
    //! strtod, so integers, decimals and exponents all work without hand-written parsing; strict
    //! JSON grammar (e.g. rejecting "1.2.3") is not enforced because the generator is well
    //! behaved. The value is always returned as double; integer-semantic fields are rounded by
    //! the caller. A value with no numeric character (e.g. starting with '{' or '"') -> false.
    //! @param out the parsed double, valid only while the function returns true.
    bool JsonExtractNumber(const std::string& line, const std::string& key, double& out)
    {
        size_t p = FindKey(line, key);
        if (p == std::string::npos)
            return false;
        SkipWs(line, p);
        // Scan the maximal run of numeric characters, then let strtod parse it.
        const size_t start = p;
        while (p < line.size())
        {
            const char c = line[p];
            if (std::isdigit(static_cast<unsigned char>(c)) || c == '-' || c == '+' || c == '.' || c == 'e' || c == 'E')
                ++p;
            else
                break;
        }
        if (p == start)   // no numeric char -> value is not a number
            return false;
        out = std::strtod(line.substr(start, p - start).c_str(), nullptr);
        return true;
    }

    //! Extracts a flat numeric array like [u0,u1,v0,v1] (one nesting level only). Clears out on
    //! entry, so a failure leaves it empty, and uses the same character set as JsonExtractNumber
    //! so both agree on what a number is. An empty [] is valid. Used for manifest uvRange, whose
    //! caller checks for at least 4 values.
    //! @param out cleared on entry, then filled with the parsed numbers.
    //! Nested arrays are not supported (one level only).
    bool JsonExtractNumberArray(const std::string& line, const std::string& key, std::vector<double>& out)
    {
        out.clear();
        size_t p = FindKey(line, key);
        if (p == std::string::npos)
            return false;
        SkipWs(line, p);
        if (p >= line.size() || line[p] != '[')
            return false;
        ++p;
        // Elements are separated by optional commas and whitespace.
        for (;;)
        {
            SkipWs(line, p);
            if (p >= line.size())
                return false;
            if (line[p] == ']')
                break;
            const size_t start = p;
            while (p < line.size())
            {
                const char c = line[p];
                if (std::isdigit(static_cast<unsigned char>(c)) || c == '-' || c == '+' || c == '.' || c == 'e' || c == 'E')
                    ++p;
                else
                    break;
            }
            if (p == start)
                return false;
            out.push_back(std::strtod(line.substr(start, p - start).c_str(), nullptr));
            SkipWs(line, p);
            if (p < line.size() && line[p] == ',')
                ++p;
        }
        return true;
    }

    //! Strips leading/trailing whitespace only (isspace), leaving interior whitespace intact;
    //! used so a trailing '\r' on CRLF manifests is not read as line content and cannot leak
    //! into a parsed value.
    //! @return a trimmed copy; the input string is not modified.
    std::string Trim(const std::string& s)
    {
        // Two pointers narrow in from both ends until they meet non-whitespace.
        size_t a = 0, b = s.size();
        while (a < b && std::isspace(static_cast<unsigned char>(s[a])))
            ++a;
        while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1])))
            --b;
        return s.substr(a, b - a);
    }

} // namespace

// LoadManifest - parse manifest.jsonl, one single-line JSON object per non-blank line.
// Key -> ManifestEntry mapping (absent optional fields keep their defaults):
//   "id"        required; a missing or empty id rejects the whole file.
//   "kind"/"type" strings; "sgt" a relative path kept verbatim (never joined or normalized here).
//   "expected"  via CanonicalFromString; missing/unrecognized -> Unknown (the negative marker).
//   "ctrlCount" a double rounded to int; "seed" a double truncated to long long.
//   "uvRange"   the first 4 values; fewer than 4 leaves the default {0,0,0,0}.
// The file is committed atomically: rows accumulate in a local vector and out is assigned only
// on full success, so any bad line leaves out empty (the all-or-nothing contract in testkit.hpp).
// The parser reads only the keys it needs; the manifest also carries a "pose" block (ground-truth
// parameters + the rigid transform) and ctrlTarget/ctrlU/ctrlV, which no program parses today.
// uvRange units follow the surface parameterization (usually radians); getline/ifstream do not
// throw, so no try/catch is needed.
// @param path manifest path; @param out cleared at entry and left empty on failure.
// @return true iff the file opened and every non-blank line yielded a non-empty id.
bool LoadManifest(const std::string& path, std::vector<ManifestEntry>& out)
{
    // Clear up front and accumulate locally so a failure at any line leaves out empty.
    out.clear();
    std::ifstream ifs(path.c_str());
    if (!ifs.is_open())
        return false;

    std::vector<ManifestEntry> tmp;
    std::string                line;
    while (std::getline(ifs, line))
    {
        const std::string trimmed = Trim(line);
        if (trimmed.empty())
            continue;   // skip blank / whitespace-only lines

        // Every string field is optional here except id; a missing one stays empty.
        ManifestEntry e;
        JsonExtractString(trimmed, "id", e.id);
        if (e.id.empty())
            return false; // a line without id rejects the whole manifest

        JsonExtractString(trimmed, "kind", e.kind);
        JsonExtractString(trimmed, "type", e.type);
        JsonExtractString(trimmed, "sgt", e.sgtPath);

        // Missing / unrecognized expected -> Unknown, which is how negatives are marked.
        std::string exp;
        if (JsonExtractString(trimmed, "expected", exp))
            e.expected = CanonicalFromString(exp);

        double num = 0.0;
        if (JsonExtractNumber(trimmed, "ctrlCount", num))
            e.ctrlCount = static_cast<int>(num + (num >= 0.0 ? 0.5 : -0.5));   // round to nearest
        if (JsonExtractNumber(trimmed, "seed", num))
            e.seed = static_cast<long long>(num);   // truncate to integer
        // Need at least 4 values; otherwise uvRange keeps its default {0,0,0,0}.
        std::vector<double> uv;
        if (JsonExtractNumberArray(trimmed, "uvRange", uv) && uv.size() >= 4)
        {
            for (int i = 0; i < 4; ++i)
                e.uvRange[i] = uv[static_cast<size_t>(i)];
        }

        tmp.push_back(std::move(e));   // held locally; committed after the whole file parses
    }

    // Only a stream error (badbit) is a failure; reaching EOF normally means success.
    if (ifs.bad())
        return false;

    out = std::move(tmp);   // commit only after the whole file parsed
    return true;
}

// Derives the data dir from the executable's own location so batch tools run zero-arg from any
// cwd: the build output lives in <root>/build/, so the parent is taken and /data appended;
// otherwise it falls back to the relative "data". argv[0] is not used because on Windows it is
// often just the bare program name, which is not enough to locate the project root. Callers may
// override with --corpus / --out-dir. No existence check is performed here.
// @return the data dir with '/' separators and no trailing separator.
std::string DefaultDataDir()
{
    // Either platform call may fail, leaving exePath empty to trigger the fallback below.
    // Ask the OS for this process's own image path; neither call relies on argv[0].
    std::string exePath;   // empty if the exe path cannot be obtained -> fallback below

#if defined(_WIN32)
    // GetModuleFileNameA returns the full path of the running process image.
    char buf[MAX_PATH] = { 0 };
    const DWORD n = ::GetModuleFileNameA(nullptr, buf, MAX_PATH);
    if (n > 0 && n < MAX_PATH)
    {
        exePath.assign(buf, static_cast<std::size_t>(n));
    }
#elif defined(__linux__)
    // /proc/self/exe is a symlink to this process's own image.
    char buf[4096] = { 0 };
    const ssize_t n = ::readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (n > 0)
    {
        exePath.assign(buf, static_cast<std::size_t>(n));
    }
#endif

    // Derive only when the executable path was actually obtained.
    if (!exePath.empty())
    {
        // Normalize backslashes to '/' so only one separator needs splitting (Windows paths).
        for (std::string::size_type i = 0; i < exePath.size(); ++i)
        {
            if (exePath[i] == '\\')
            {
                exePath[i] = '/';
            }
        }

        // Split at the last separator to get the containing directory.
        const std::string::size_type slash = exePath.find_last_of('/');
        if (slash != std::string::npos)
        {
            // Directory part; the root case (e.g. "/prog") keeps the leading '/'.
            const std::string exeDir = (slash == 0) ? std::string("/") : exePath.substr(0, slash);
            // Build output is <root>/build/, so the parent is <root>, holding data/.
            return exeDir + "/../data";
        }
    }

    // Fallback: relative "data" when the executable location is unknown.
    return "data";
}

} // namespace sgkq
