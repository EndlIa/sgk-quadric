// sgkquadric/io.hpp - I/O ingestion of real CAD data: read .sgt / .step into a sgkq::InputFace
// list. Downstream detect/fit/write-back modules depend only on InputFace (surface + kind +
// origin), not on the file format. Every failure path sets its own LoadErrorCode and a message
// carrying the origin path; nothing throws and nothing silently returns an empty result. .sgt
// versions are checked against the kernel (a mismatch gives LoadErrVersionMismatch). Shared types
// are in types.hpp; the evaluation harness (corpus manifest, default data dir) lives in testkit.hpp.
#pragma once

#include "sgkquadric/types.hpp"

#include <string>
#include <vector>

namespace sgkq
{
    //! I/O layer error codes (0 = success; each failure class is distinct for reporting/stats).
    enum LoadErrorCode : int
    {
        LoadOk = 0,                 //!< success
        LoadErrFileNotFound = 1,    //!< path missing / cannot open
        LoadErrReadFailed = 2,      //!< opened but reading bytes failed
        LoadErrJsonParse = 3,       //!< JSON syntax error (truncated / illegal characters)
        LoadErrDeserialize = 4,     //!< valid JSON but topology deserialization threw
        LoadErrVersionMismatch = 5, //!< Meta version outside the kernel's supported range
        LoadErrEmptyBody = 6,       //!< deserialization/import produced a null Body
        LoadErrEmptyFaces = 7,      //!< Body has no usable face
        LoadErrStepImport = 8       //!< STEP import failed (threw or kernel error)
    };

    //! Result of one ingestion.
    struct LoadResult
    {
        bool                     ok = false;    //!< success flag
        std::vector<InputFace>   faces;         //!< faces on success (kind/origin set)
        int                      errorCode = LoadOk; //!< a LoadErrorCode; LoadOk when ok == true
        std::string              message;       //!< summary on success; origin-tagged reason on failure
    };

    //! Reads a .sgt (Topology/Serialize JSON) and assembles InputFaces. @return LoadResult; on
    //! failure errorCode identifies the cause and message is "LoadSgt: <path>: <reason>".
    LoadResult LoadSgt(const std::string& path);

    //! Reads a .step/.stp file via the STEP importer (needs a DATAEXCHANGE license). @return
    //! LoadResult; on failure message is "LoadStep: <path>: <reason>".
    LoadResult LoadStep(const std::string& path);

} // namespace sgkq
