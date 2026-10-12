// sgkquadric/testkit.hpp - evaluation-harness support; NOT part of the library's ingestion API.
// manifest.jsonl is this project's own test-corpus format, so its parser and the ground-truth type
// mapping live outside the library target (built as sgkqtestkit), together with the
// executable-relative default data dir used by the batch tools. The library proper only ingests
// real CAD data - see io.hpp: LoadSgt (.sgt) / LoadStep (STEP).
#pragma once

#include "sgkquadric/types.hpp"

#include <string>
#include <vector>

namespace sgkq
{
    //! One case from manifest.jsonl.
    struct ManifestEntry
    {
        std::string id;                 //!< case id
        std::string kind;               //!< "positive" / "negative" etc.
        std::string type;               //!< surface type string written by the generator
        Canonical   expected = Canonical::Unknown; //!< ground-truth type (Unknown = negative)
        std::string sgtPath;            //!< relative .sgt path (kept verbatim)
        int         ctrlCount = 0;      //!< total control points
        double      uvRange[4] = { 0.0, 0.0, 0.0, 0.0 }; //!< [u0,u1,v0,v1]
        long long   seed = 0;           //!< generation random seed
    };

    //! Parses manifest.jsonl (one JSON object per line). @param out is cleared on entry.
    //! @return true only if the file opens and every non-blank line yields an id; otherwise false
    //! with out left empty (all-or-nothing).
    bool LoadManifest(const std::string& path, std::vector<ManifestEntry>& out);

    //! Default data dir derived from the executable location (parent of <root>/build/ plus /data),
    //! so batch tools run zero-arg from any cwd; falls back to the relative "data". Callers may
    //! override with --corpus / --out-dir. @return a '/' -separated path with no trailing separator;
    //! no existence check is performed.
    std::string DefaultDataDir();

} // namespace sgkq
