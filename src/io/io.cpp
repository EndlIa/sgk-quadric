// src/io/io.cpp - I/O ingestion: read .sgt / .step model files into InputFace lists.
// Two entry points: LoadSgt (.sgt topology JSON, with a Meta version check) and LoadStep (STEP via
// api_step_import, needs a DATAEXCHANGE license). Every failure sets a specific LoadErrorCode and
// a "LoadSgt:/LoadStep: <path>: <reason>" message; nothing throws and nothing silently returns an
// empty result. No unit conversion is done at this layer.
// Input is a file path (an empty path is treated as "file not found"); output is a LoadResult
// whose faces carry kind (source) and origin (path). The corpus manifest parser and the
// benchmark data dir are harness-side code (testkit.hpp / sgkqtestkit), not part of this layer.
// .sgt runs ParseString (syntax) then DeserializeBody (topology); their exceptions map to
// different error codes, so the two SDK calls stay in separate try blocks. The version check
// is skipped when meta.sggkVersion does not parse (backward compat with old archives).
#include "sgkquadric/io.hpp"

#include <Foundation/SGGK_Version.h>
#include <StepExchange/API.h>
#include <Topology/Brep/Body.h>
#include <Topology/Brep/Face.h>
#include <Topology/Serialize/RapidTopoJsonDeserializer.h>
#include <Topology/Meta.h>

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>

using namespace sggk; // SDK types (Meta / BodyPtr / SurfacePtr / FaceList / StepImportRetPtr)

namespace sgkq
{
namespace
{
    // Internal helpers: file and JSON utilities used only inside this translation unit.
    //! Reads the whole file in binary. Returns false on open failure or badbit; LoadSgt
    //! re-probes with FileOpenable to tell "missing" from "read failed".
    //! Binary mode avoids Windows CRLF->LF translation, so the JSON bytes match the file exactly.
    //! @param out holds the contents only while the function returns true.
    //! @return false on open failure or a badbit stream error; eof/failbit alone are not errors.
    bool ReadWholeFile(const std::string& path, std::string& out)
    {
        std::ifstream ifs(path.c_str(), std::ios::binary);
        if (!ifs.is_open())
            return false;
        std::ostringstream ss;
        ss << ifs.rdbuf();   // one bulk streambuf copy instead of a manual read loop
        if (ifs.bad())       // only badbit is failure; eof/failbit mean normal end of file
            return false;
        out = ss.str();
        return true;
    }

    //! Probes whether the file opens binary/read-only; no bytes are read and the file is not
    //! modified. Uses the same open mode as ReadWholeFile so the two probes agree. Pairs with it
    //! so LoadSgt can map "openable but read failed" to ReadFailed and "not openable" to
    //! FileNotFound, which tells a path problem apart from a storage/permission problem.
    //! @return true iff the binary/read-only open succeeded; no bytes are read.
    bool FileOpenable(const std::string& path)
    {
        std::ifstream ifs(path.c_str(), std::ios::binary);
        return ifs.is_open();
    }

    //! Parses the leading major/minor fields of a version string ("1.7.11" -> 1,7). Only the
    //! first two dot-separated numeric fields count; the patch number is ignored.
    //! Scans numeric fields rather than split+stoi: no allocation and no exceptions on an empty
    //! string or a stray dot. A string with no leading digit (e.g. "v1.7") returns false.
    //! @param s the version string; @param major/@param minor the parsed fields (zero on failure).
    //! @return false when nothing could be parsed; callers then skip the version check instead of
    //! rejecting the file, which keeps old archives (written before versioning) loadable.
    bool ParseMajorMinor(const std::string& s, int& major, int& minor)
    {
        major = minor = 0;
        if (s.empty())
            return false;
        size_t i = 0;
        int    got = 0;
        int    vals[2] = { 0, 0 };
        // Read up to two dot-separated numeric fields; stop at the first non-digit.
        while (i < s.size() && got < 2)
        {
            if (!std::isdigit(static_cast<unsigned char>(s[i])))
                break;
            int v = 0;
            while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i])))
            {
                v = v * 10 + (s[i] - '0');
                ++i;
            }
            vals[got++] = v;
            if (i < s.size() && s[i] == '.')
                ++i;
            else
                break;
        }
        if (got < 1)
            return false;
        major = vals[0];
        minor = (got >= 2) ? vals[1] : 0;
        return true;
    }

    //! Builds the InputFace list from body->QueryFaces(); null faces and faces without a
    //! geometry surface are skipped and counted, never aborting the rest (one bad face must not
    //! drop the others). Order matches QueryFaces(): no sorting, no de-duplication. Ground truth
    //! is not carried in InputFace: the evaluation side reads it from the corpus manifest, so it
    //! cannot reach the algorithm. Both geometry entry points share this helper so InputFace is
    //! assembled in exactly one place (no risk of one path forgetting kind/origin).
    //! @param skipped reset to this call's skip count; passed by reference so the caller can
    //!        splice it into the success/error message.
    //! @return one InputFace per usable face; their SurfacePtr shares the Body's geometry and
    //!         stays valid as long as the Body / result is alive.
    std::vector<InputFace> FacesToInput(const BodyPtr& body, SourceKind kind, const std::string& origin, int& skipped)
    {
        std::vector<InputFace> out;
        skipped = 0;   // reset so the caller cannot read a stale count
        const FaceList fl = body->QueryFaces();
        for (const FacePtr& fp : fl)
        {
            if (!fp)
            {
                ++skipped;
                continue;
            }
            SurfacePtr srf = fp->GeomSurface();   // topology-only faces (no geometry) are unusable
            if (!srf)
            {
                ++skipped;
                continue;
            }
            InputFace in;
            in.surface = srf;
            in.kind    = kind;
            in.origin  = origin;
            out.push_back(std::move(in));
        }
        return out;
    }

} // namespace

// LoadSgt - read a .sgt (Topology/Serialize JSON): read bytes -> ParseString (syntax) ->
// DeserializeBody (topology; fills meta) -> version check -> face extraction.
// Stage order is also the error-code priority (each stage has its own code):
//   A. empty path or unreadable -> FileNotFound; openable but read error -> ReadFailed;
//      an empty file -> JsonParse (an empty file is not valid JSON).
//   B. ParseString threw (truncated / illegal characters) -> JsonParse.
//   C. DeserializeBody threw (valid JSON, bad topology) -> Deserialize; a null Body -> EmptyBody.
//   D. meta.sggkVersion major/minor mismatch -> VersionMismatch; an unparsable version skips the
//      check (backward compat) and the patch number is never compared.
//   E. no usable faces -> EmptyFaces; the message carries total faces and skipped count.
// Distinct codes let callers tell "corrupt file" from "incompatible semantics" instead of a
// generic "import failed". Nothing is converted or unit-transformed here.
// @param path .sgt file path.
// @return on success ok=true with the usable faces and an "ok faces=N skipped=M sggkVersion=...
//         dataVersion=..." message; on failure ok=false with the stage's errorCode and a
//         "LoadSgt: <path>: <reason>" message.
LoadResult LoadSgt(const std::string& path)
{
    // Default failure state; every failure path below sets its own specific error code.
    LoadResult r;
    r.ok        = false;
    r.errorCode = LoadOk;

    if (path.empty())
    {
        // An empty path is treated as "missing"; do not touch the filesystem with an empty name.
        r.errorCode = LoadErrFileNotFound;
        r.message   = "LoadSgt: <empty path>: file not found";
        return r;
    }

    // Re-probe with FileOpenable to split "cannot open" from "openable but read failed".
    std::string content;
    if (!ReadWholeFile(path, content))
    {
        if (FileOpenable(path))
        {
            r.errorCode = LoadErrReadFailed;
            r.message   = "LoadSgt: " + path + ": read failed (I/O error while reading bytes)";
        }
        else
        {
            r.errorCode = LoadErrFileNotFound;
            r.message   = "LoadSgt: " + path + ": file not found or not readable";
        }
        return r;
    }

    if (content.empty())
    {
        // An empty file is not valid JSON (syntax error, not EmptyBody).
        r.errorCode = LoadErrJsonParse;
        r.message   = "LoadSgt: " + path + ": empty file is not valid JSON";
        return r;
    }

    // Syntax layer only: ParseString validates JSON, not topology semantics (checked next).
    // ParseString then DeserializeBody is a fixed two-step SDK API and cannot be merged.
    RapidTopoJsonDeserializer des;
    try
    {
        des.ParseString(content);
    }
    catch (const std::exception& e)
    {
        r.errorCode = LoadErrJsonParse;
        r.message   = "LoadSgt: " + path + ": JSON parse error: " + e.what();
        return r;
    }
    catch (...)
    {
        r.errorCode = LoadErrJsonParse;
        r.message   = "LoadSgt: " + path + ": JSON parse error (unknown exception)";
        return r;
    }

    // Topology layer: valid JSON with bad topology semantics throws here; meta is filled by SDK.
    // des / meta / body are function-local; the returned faces hold shared SurfacePtrs instead.
    Meta    meta;
    BodyPtr body;
    try
    {
        body = des.DeserializeBody(content, meta);
    }
    catch (const std::exception& e)
    {
        r.errorCode = LoadErrDeserialize;
        r.message   = "LoadSgt: " + path + ": deserialize failed: " + e.what();
        return r;
    }
    catch (...)
    {
        r.errorCode = LoadErrDeserialize;
        r.message   = "LoadSgt: " + path + ": deserialize failed (unknown exception)";
        return r;
    }

    if (!body)
    {
        // Deserialization returned no topology payload at all (distinct from EmptyFaces).
        r.errorCode = LoadErrEmptyBody;
        r.message   = "LoadSgt: " + path + ": deserialized body is null (no usable topology payload)";
        return r;
    }

    // Compare major/minor only; an unparsable version string skips the check (backward compat).
    int    maj = 0, min = 0;
    if (ParseMajorMinor(meta.sggkVersion, maj, min))
    {
        if (maj != SGGK_VERSION_MAJOR || min != SGGK_VERSION_MINOR)
        {
            // Strict equality; the message names the version the kernel supports.
            r.errorCode = LoadErrVersionMismatch;
            r.message   = "LoadSgt: " + path + ": unsupported data version sggkVersion=" + meta.sggkVersion +
                          " (supported " + SGGK_VERSION_COMPLETE + ")";
            return r;
        }
    }

    // FacesToInput counts null / geometry-less faces in skipped, reported in the message below.
    int              skipped = 0;
    std::vector<InputFace> faces = FacesToInput(body, SourceKind::SgtJson, path, skipped);
    if (faces.empty())
    {
        // Non-null body but no usable geometry (e.g. a pure topological shell).
        r.errorCode = LoadErrEmptyFaces;
        r.message   = "LoadSgt: " + path + ": body has no usable faces (faces=" +
                      std::to_string(body->QueryFaces().size()) + " skipped=" + std::to_string(skipped) + ")";
        return r;
    }

    // Success: hand over the faces and summarize counts + version for one-line log confirmation.
    r.ok      = true;
    r.faces   = std::move(faces);
    r.message = "ok faces=" + std::to_string(r.faces.size()) +
                " skipped=" + std::to_string(skipped) +
                " sggkVersion=" + meta.sggkVersion +
                " dataVersion=" + std::to_string(meta.dataVersion);
    return r;
}

// LoadStep - STEP import via api_step_import (needs a DATAEXCHANGE license).
// Order: probe the path -> import -> check the result -> extract faces; a bad path is reported
// as FileNotFound rather than as an import failure.
// Failure mapping:
//   B. api_step_import threw, returned null, or Succeeded() is false -> StepImport; the SDK text
//      (Status().ErrorMsg()) is appended to the message, which is the main diagnostic here.
//   C. ResultBody() is null (assembly with no solid, or empty import result) -> EmptyBody.
//   D. no usable faces -> EmptyFaces (message carries total faces and skipped count).
// Unlike LoadSgt there is no Meta version check (STEP carries its own schema version, handled by
// the SDK), and units inside the file are decided by the importer, not converted here.
// @param path .step/.stp file path.
// @return success: ok=true with an "ok faces=N skipped=M" message; failure: ok=false with the
//         stage's errorCode and a "LoadStep: <path>: <reason>" message.
LoadResult LoadStep(const std::string& path)
{
    LoadResult r;
    r.ok        = false;
    r.errorCode = LoadOk;

    if (path.empty())
    {
        r.errorCode = LoadErrFileNotFound;
        r.message   = "LoadStep: <empty path>: file not found";
        return r;
    }
    // Probe before importing so a path problem is not reported as an import failure.
    if (!FileOpenable(path))
    {
        r.errorCode = LoadErrFileNotFound;
        r.message   = "LoadStep: " + path + ": file not found or not readable";
        return r;
    }

    // ret is an SDK smart pointer; failure appears as a null pointer or Succeeded() == false.
    // A missing DATAEXCHANGE license surfaces here as a throw or a failed result.
    StepImportRetPtr ret;
    try
    {
        // StepImportOpts{} uses the importer's defaults.
        ret = api_step_import(path.c_str(), StepImportOpts());
    }
    catch (const std::exception& e)
    {
        r.errorCode = LoadErrStepImport;
        r.message   = "LoadStep: " + path + ": step import threw: " + e.what();
        return r;
    }
    catch (...)
    {
        r.errorCode = LoadErrStepImport;
        r.message   = "LoadStep: " + path + ": step import threw (unknown exception)";
        return r;
    }

    if (!ret)
    {
        // No result means no status text is available; give a generic reason.
        r.errorCode = LoadErrStepImport;
        r.message   = "LoadStep: " + path + ": step import returned null result";
        return r;
    }
    if (!ret->Succeeded())
    {
        // Explicit failure: append the SDK error text, the most actionable diagnostic here.
        r.errorCode = LoadErrStepImport;
        r.message   = "LoadStep: " + path + ": step import failed: " + ret->Status().ErrorMsg();
        return r;
    }

    // A null body is common when the file imports to an assembly with no solid.
    const BodyPtr body = ret->ResultBody();
    if (!body)
    {
        // Null body: assembly without a solid, or an empty import result.
        r.errorCode = LoadErrEmptyBody;
        r.message   = "LoadStep: " + path + ": import produced no body (assembly or empty result)";
        return r;
    }

    int                    skipped = 0;
    std::vector<InputFace> faces   = FacesToInput(body, SourceKind::Step, path, skipped);
    if (faces.empty())
    {
        // Body non-null but no usable geometry; distinct from EmptyBody.
        r.errorCode = LoadErrEmptyFaces;
        r.message   = "LoadStep: " + path + ": body has no usable faces (faces=" +
                      std::to_string(body->QueryFaces().size()) + " skipped=" + std::to_string(skipped) + ")";
        return r;
    }

    // STEP has no Meta, so the success message carries only counts.
    r.ok      = true;
    r.faces   = std::move(faces);
    r.message = "ok faces=" + std::to_string(r.faces.size()) + " skipped=" + std::to_string(skipped);
    return r;
}

} // namespace sgkq

// LoadErrorCode triggers (see io.hpp for the enum). Failure messages are shaped
// "<entry>: <path>: <reason>" (entry in {LoadSgt, LoadStep}) so logs locate the exact file;
// success messages carry the face/skip counts, plus sggkVersion/dataVersion for .sgt.
//   1 FileNotFound    : empty path, or FileOpenable() false.
//   2 ReadFailed      : openable but ReadWholeFile() hit a byte I/O error (.sgt path only).
//   3 JsonParse       : empty content, or ParseString() threw (bad/truncated JSON syntax).
//   4 Deserialize     : DeserializeBody() threw (valid JSON, bad topology semantics).
//   5 VersionMismatch : meta.sggkVersion major/minor differs from the kernel.
//   6 EmptyBody       : DeserializeBody()/ResultBody() yielded a null Body.
//   7 EmptyFaces      : FacesToInput() produced no faces.
//   8 StepImport      : api_step_import threw, returned null, or Succeeded() is false.
//
// Format notes: .sgt is Topology/Serialize JSON (Meta carries sggkVersion/dataVersion/
// originDataVersion; Brep/Body faces expose their geometry via Face::GeomSurface()); .step is
// ISO 10303 imported by StepExchange, whose ResultBody() may be null for assembly-only files.


