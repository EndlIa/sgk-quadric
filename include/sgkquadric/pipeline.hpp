// sgkquadric/pipeline.hpp - per-face pipeline: detect -> fit -> verify -> (optional) write-back.
// Chains SampleSurface/Detect (detect.hpp), FitSurface/RefineSurface/MeasureAccuracy (fit.hpp),
// the distXxx closed-form distances (math.hpp) and the FaceReport/Decision/Thresholds contract
// (types.hpp). Implementation in src/pipeline/pipeline.cpp.
#ifndef SGKQUADRIC_PIPELINE_HPP
#define SGKQUADRIC_PIPELINE_HPP

#include <Topology/Brep/Face.h>

#include "sgkquadric/types.hpp"

namespace sgkq
{
// Face.h declares FacePtr (via DECL_SHARED(Face)) in the sggk namespace.
using sggk::FacePtr;

//! Tunable pipeline options (sample density / write-back / refinement).
struct PipelineOptions
{
    int    nu = 9;                 //!< detection / fit sample density (U)
    int    nv = 9;                 //!< detection / fit sample density (V)
    int    nuAcc = 33;             //!< accuracy-check sample density (U)
    int    nvAcc = 33;             //!< accuracy-check sample density (V)
    bool   writeBack = true;       //!< replace the face once verification passes
    int    refineIter = 4;         //!< Gauss-Newton iteration cap
    double refineTol = 1e-12;      //!< refinement convergence threshold
};

//! Processes one face: detect -> fit -> verify -> (optional) write-back; returns a report.
//! An already-analytic surface (Plane/Cylinder/Cone/Sphere/Torus) -> Decision::Rejected with
//! sourceType set; otherwise a BSplineSurface is sampled on face->CalcUVBound() intersect its own
//! domain (Unknown -> Rejected, else fit + refine + verify). Write-back (Decision::Converted)
//! requires maxDist <= th.maxDist and boundaryCurveHausdorff <= th.boundaryHausdorff; a null face or a
//! throwing write-back gives Decision::Failed. @return a fully populated FaceReport.
FaceReport ProcessFace(const FacePtr& face,
                       const Thresholds& th = Thresholds{},
                       const PipelineOptions& opt = PipelineOptions{});

} // namespace sgkq

#endif // SGKQUADRIC_PIPELINE_HPP
