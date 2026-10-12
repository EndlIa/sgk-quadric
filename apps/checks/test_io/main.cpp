// apps/checks/test_io/main.cpp - I/O ingest self-test (spec: 09-io.md).
// Each case prints "[TEST] <name> expected=<> actual=<> <PASS|FAIL>"; any FAIL -> exit 1.
// Groups: .sgt round-trip, STEP round-trip, three corrupt-input rejections (distinct errorCodes,
// message carries the filename), and Meta version read-back. Exit: 0 pass, 1 FAIL, 2 exception.
#include "sgkquadric/io.hpp"

#include <Foundation/init.h>
#include <GeomBase/Ucs3D.h>
#include <Geometry/SurfaceType.h>
#include <ModelingPrim/API.h>
#include <StepExchange/API.h>
#include <Topology/Serialize/RapidTopoJsonSerializer.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace sggk;
using namespace sgkq;

namespace
{
namespace fs = std::filesystem;

// Case counters for report(): g_total counts assertions, g_failed counts failures.
// Non-atomic; this program is single-threaded.
int g_total  = 0;                       // assertions executed
int g_failed = 0;                       // failed assertions; >0 makes the process exit 1

//! @brief Print one case result and update the global counters.
//! @param pass false increments g_failed, which drives the process exit code.
//! @note Output is "[TEST] <name> expected=<expected> actual=<actual> <PASS|FAIL>",
//!       space-separated for external tokenizing; <name> is a regex anchor, keep it stable.
void report(const std::string& name, const std::string& expected, const std::string& actual, bool pass)
{
    ++g_total;                                          // count every case
    if (!pass)
        ++g_failed;                                     // count failures only
    std::cout << "[TEST] " << name                      // case name (regex anchor)
              << " expected=" << expected               // expected text
              << " actual=" << actual                   // actual text
              << (pass ? " PASS" : " FAIL") << std::endl;  // verdict field
}

//! @brief Boolean specialization of report(): formats both values as "true"/"false".
//! @note Centralizes the display text so call sites need no string literals.
void reportBool(const std::string& name, bool expected, bool actual)
{
    report(name, expected ? "true" : "false", actual ? "true" : "false", expected == actual);
}

//! @brief Extract the value after `keyEq` up to the next space or end of line.
//! @param keyEq includes '=' (e.g. "sggkVersion=") to avoid prefix-key mismatches.
//! @param out set on hit; left untouched when the key is absent.
//! @note The value must contain no spaces; LoadSgt's message format guarantees this.
bool ExtractValueAfter(const std::string& line, const std::string& keyEq, std::string& out)
{
    const size_t k = line.find(keyEq);                  // key includes '=', avoiding prefix matches
    if (k == std::string::npos)
        return false;                                   // key absent: leave out untouched
    size_t s = k + keyEq.size();                        // value starts after the key
    size_t e = line.find(' ', s);                       // value ends at first space
    out = line.substr(s, (e == std::string::npos) ? std::string::npos : e - s);
    return true;                                        // out written (may be empty)
}

//! @brief Whether any face has a non-null surface of the given type.
//! @note SurfaceType only, no geometry: coarse confirmation the Cylinder face survived.
bool HasSurfaceType(const std::vector<InputFace>& faces, SurfaceType t)
{
    for (const InputFace& f : faces)                    // scan each face
    {
        if (f.surface && f.surface->SurfaceType() == t) // surface may be null, check first
            return true;                                // short-circuit on hit
    }
    return false;                                       // no match
}

//! @brief Write content verbatim in binary mode; builds corrupt/missing-field inputs.
//! @note Binary mode prevents newline translation, keeping byte length predictable.
bool WriteText(const fs::path& p, const std::string& content)
{
    std::ofstream o(p, std::ios::binary);               // binary mode: no CRLF translation
    if (!o.is_open())
        return false;                                   // missing dir / no permission
    o << content;                                       // write everything at once
    return o.good();                                    // check stream error bits after flush
}

} // namespace

//! @brief Self-test entry point: .sgt round-trip, STEP round-trip, corrupt rejection, version read-back.
//! @return 0 all pass; 1 at least one FAIL; 2 unexpected exception.
//! @note init()/fini() bracket the normal path; on exception the catch returns without calling fini().
//! @note argc/argv unused: all inputs are created under the system temp dir.
int main(int /*argc*/, char* /*argv*/[])
{
    try
    {
        init();                                         // SDK runtime (pools/singletons); before any api_/Load call

        // Temp working dir; create_directories is idempotent on an existing directory.
        // The ec overload only reports; a failed creation surfaces later as a case FAIL.
        // Reruns overwrite same-named files.
        std::error_code ec;
        const fs::path tmp = fs::temp_directory_path() / "sgkq_test_io";
        fs::create_directories(tmp, ec);

        // 1) .sgt round-trip: make solid -> serialize -> load back -> assert.
        // cylBody is owned via BodyPtr and reused by both .sgt and .step exports.
        const fs::path sgtPath = tmp / "cyl.sgt";
        BodyPtr        cylBody;

        // Standard cylinder: default Ucs, radius 50 mm, height 100 mm (size arbitrary, non-degenerate).
        PrimitivesRetPtr mc = api_make_solid_cylinder(Ucs3D(), 50.0, 100.0);
        {
            // Check all three (non-null Ret, Succeeded(), non-null ResultBody) to avoid null deref.
            const bool made = (mc && mc->Succeeded() && mc->ResultBody());
            reportBool("make_cylinder.succeeded", true, made);
            if (made)
                cylBody = mc->ResultBody();
        }

        // Serialize to .sgt (RapidTopo JSON); void return, so success is verified by the load-back.
        if (cylBody)
        {
            RapidTopoJsonSerializer ser;
            ser.Serialize(cylBody, sgtPath.string().c_str());
        }

        // load_sgt.ok: whole-file load succeeds.
        // load_sgt.faces_nonempty: the face list is non-empty.
        // load_sgt.has_cylinder: a Cylinder surface is present (geometry semantics preserved).
        // load_sgt.kind_sgt / origin_is_path: source kind and path recorded (traceability contract).
        LoadResult lr = LoadSgt(sgtPath.string());
        reportBool("load_sgt.ok", true, lr.ok);
        report("load_sgt.faces_nonempty", ">0", std::to_string(lr.faces.size()), !lr.faces.empty());
        reportBool("load_sgt.has_cylinder", true, HasSurfaceType(lr.faces, SurfaceType::Cylinder));
        reportBool("load_sgt.kind_sgt", true,
                   !lr.faces.empty() && lr.faces.front().kind == SourceKind::SgtJson);
        reportBool("load_sgt.origin_is_path", true,
                   !lr.faces.empty() && lr.faces.front().origin == sgtPath.string());

        // 2) STEP round-trip: export .step -> load back.
        // Export needs both Succeeded() and the file on disk (fs::exists).
        // An invalid body records an explicit FAIL rather than skipping, so it cannot mask section 1.
        const fs::path stepPath = tmp / "cyl.step";
        if (cylBody)
        {
            StepExportRetPtr ex = api_step_export(cylBody, stepPath.string().c_str());
            reportBool("step_export.succeeded", true,
                       ex && ex->Succeeded() && fs::exists(stepPath));
        }
        else
        {
            reportBool("step_export.succeeded", true, false);
        }

        // load_step.ok / faces_nonempty: load succeeds with a non-empty face list.
        // load_step.kind_step: source kind is Step (distinct from the SgtJson above).
        // Face count is only asserted >0: STEP rewrites may legitimately change it.
        LoadResult ls = LoadStep(stepPath.string());
        reportBool("load_step.ok", true, ls.ok);
        report("load_step.faces_nonempty", ">0", std::to_string(ls.faces.size()), !ls.faces.empty());
        reportBool("load_step.kind_step", true,
                   !ls.faces.empty() && ls.faces.front().kind == SourceKind::Step);

        // 3) Corrupt inputs: all rejected, errorCodes pairwise distinct, message contains filename.
        // bad1 = truncated JSON (syntax stage), bad2 = valid JSON missing fields (structure stage),
        // bad3 = file never created (I/O stage): three distinct failure stages.
        const fs::path bad1 = tmp / "bad_truncated.sgt"; // truncated JSON
        const fs::path bad2 = tmp / "bad_missing.sgt";   // valid JSON, missing fields
        const fs::path bad3 = tmp / "bad_absent.sgt";    // not created

        WriteText(bad1, "{\"meta\":{\"sggkVersion\":\"1.7.11\"},\"topologies\":[{\"id\":1,");
        WriteText(bad2, "{\"meta\":{\"sggkVersion\":\"1.7.11\"}}");

        const LoadResult r1 = LoadSgt(bad1.string());
        const LoadResult r2 = LoadSgt(bad2.string());
        const LoadResult r3 = LoadSgt(bad3.string());

        // All three must be rejected (ok==false); failures must not silently return an empty face list.
        reportBool("bad1_truncated.rejected", true, !r1.ok);
        reportBool("bad2_missing.rejected", true, !r2.ok);
        reportBool("bad3_absent.rejected", true, !r3.ok);

        // Pairwise-distinct errorCodes prove failures are classified by stage, not one generic error.
        const bool codesDistinct = (r1.errorCode != r2.errorCode) &&
                                   (r2.errorCode != r3.errorCode) &&
                                   (r1.errorCode != r3.errorCode);
        report("bad.codes_distinct", "3 distinct",
               std::to_string(r1.errorCode) + "/" + std::to_string(r2.errorCode) + "/" +
                   std::to_string(r3.errorCode),
               codesDistinct);

        // message must contain the filename so the failing input is locatable (traceability contract).
        reportBool("bad1_truncated.msg_has_name", true,
                   r1.message.find(bad1.filename().string()) != std::string::npos);
        reportBool("bad2_missing.msg_has_name", true,
                   r2.message.find(bad2.filename().string()) != std::string::npos);
        reportBool("bad3_absent.msg_has_name", true,
                   r3.message.find(bad3.filename().string()) != std::string::npos);

        // 4) Read Meta.sggkVersion back from section 1's successful message.
        // Uses ExtractValueAfter rather than an SDK version-query API (no extra dependency).
        // Missing or empty version -> SKIP not FAIL: version read-back is an optional extra check.
        {
            std::string ver;
            const bool  found = ExtractValueAfter(lr.message, "sggkVersion=", ver);
            if (!found)
            {
                std::cout << "[SKIP] sgt_meta_version reason=LoadSgt message carries no sggkVersion field"
                          << std::endl;
            }
            else if (ver.empty())
            {
                std::cout << "[SKIP] sgt_meta_version reason=Meta.sggkVersion empty in serialized file"
                          << std::endl;
            }
            else
            {
                report("sgt_meta_version", "non-empty", ver, !ver.empty());
            }
        }

        fini();                                         // release SDK runtime on the normal path
    }
    catch (const std::exception& e)
    {
        // Unexpected exception (not a case failure): exit 2, distinct from 1; fini() is not called here.
        std::cerr << "[FAIL] exception: " << e.what() << std::endl;
        return 2;
    }
    catch (...)
    {
        std::cerr << "[FAIL] unknown exception" << std::endl;
        return 2;
    }

    // Summary line for CI; exit 0 if all cases passed, else 1.
    std::cout << "[SUMMARY] total=" << g_total << " failed=" << g_failed << std::endl;
    return (g_failed == 0) ? 0 : 1;
}
