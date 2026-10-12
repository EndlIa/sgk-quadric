// Shared sgkquadric data types (read-only cross-module contract; report changes to the integrator).
#pragma once

#include <Foundation/SGGK_Types.h>
#include <GeomBase/Point2D.h>
#include <GeomBase/Point3D.h>
#include <GeomBase/Vector3D.h>
#include <GeomBase/UVRange.h>
#include <Geometry/3D/Surface/Surface.h>

#include <string>
#include <vector>

namespace sgkq
{
    using sggk::Point2D;
    using sggk::Point3D;
    using sggk::Surface;
    using sggk::SurfacePtr;
    using sggk::UVRange;
    using sggk::Vector3D;

    //! Target canonical surface types; Unknown means unrecognized.
    enum class Canonical
    {
        Unknown = 0,
        Plane,
        Cylinder,
        Cone,
        Sphere,
        Torus
    };

    //! Readable type name for reports/CSV (inline: one definition across translation units).
    inline const char* ToString(Canonical c)
    {
        switch (c)
        {
        case Canonical::Plane:    return "Plane";
        case Canonical::Cylinder: return "Cylinder";
        case Canonical::Cone:     return "Cone";
        case Canonical::Sphere:   return "Sphere";
        case Canonical::Torus:    return "Torus";
        default:                  return "Unknown";
        }
    }

    //! Fitted analytic surface parameters; axisPoint/axisDir: plane = point + normal, cylinder/torus = axis point + dir, cone = apex + dir, sphere = center.
    struct SurfaceParams
    {
        Canonical type = Canonical::Unknown;

        Point3D  axisPoint;
        Vector3D axisDir;

        double radius = 0.0;       // cylinder / sphere radius, or torus major radius R
        double minorRadius = 0.0;  // torus minor radius r
        double halfAngle = 0.0;    // cone half-angle (radians)
        Point3D apex;              // cone vertex
    };

    //! Detection result: candidate + a normalized invariant residual (lower = more credible).
    struct DetectResult
    {
        Canonical candidate = Canonical::Unknown;
        double    residual = 0.0;            // normalized invariant residual
        std::vector<std::string> notes;      // human-readable evidence for the verdict
    };

    //! Accuracy report (mm); boundaryCurveHausdorff = two-way Hausdorff between input and generated boundary curves (sampled, hence a slight overestimate); accepted = all thresholds met.
    struct AccuracyReport
    {
        double maxDist = 0.0;               // max point-to-surface distance
        double avgDist = 0.0;               // mean point-to-surface distance
        double boundaryCurveHausdorff = 0.0; // input vs generated boundary curve, two-way (criterion <= th.boundaryHausdorff)
        double boundaryOnSurface = 0.0;     // auxiliary: max distance from the input boundary edge to the new surface
        double g1Before = 0.0;              // pre-conversion adjacent-face normal angle
        double g1After = 0.0;               // post-conversion
        bool   accepted = false;            // all three thresholds satisfied
    };

    //! Input provenance, for traceable error messages.
    enum class SourceKind
    {
        Synthesized = 0,
        Step,
        SgtJson
    };

    //! One spline face to classify/convert. Ground truth is deliberately not carried here: the
    //! evaluation side reads it from the corpus manifest, so it can never reach the algorithm
    //! (ProcessFace takes a FacePtr, not an InputFace).
    struct InputFace
    {
        SurfacePtr  surface;   // the spline face
        SourceKind  kind = SourceKind::Synthesized;
        std::string origin;    // file path / case id, for error messages
    };

    //! Per-face outcome.
    enum class Decision
    {
        Rejected = 0,   // not a target surface; original face kept
        Converted,      // successfully replaced by an analytic surface
        Failed          // error (bad parameters / exception)
    };

    //! Per-face report row (raw data for the CSV files and the stdout label lines).
    struct FaceReport
    {
        std::string    faceId;
        Canonical      sourceType = Canonical::Unknown;
        Decision       decision = Decision::Rejected;
        SurfaceParams  params;
        AccuracyReport accuracy;
        int            errorCode = 0;
        std::string    message;
    };

    //! Acceptance thresholds.
    struct Thresholds
    {
        double maxDist = 1e-5;              // mm
        double boundaryHausdorff = 1e-5;    // mm
        double g1AngleRad = 0.1 * 3.14159265358979323846 / 180.0;  // 0.1 deg
    };

} // namespace sgkq
