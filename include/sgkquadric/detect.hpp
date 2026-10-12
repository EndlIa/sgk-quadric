// sgkquadric/detect.hpp - surface detection: UV sampler + curvature-signature classification.
// Classification order: plane -> sphere -> torus -> cylinder/cone. Degenerate points (sphere
// poles, cone apex) are skipped. The cone axis must use the smallest covariance eigenvector, not
// the largest second-moment one (the latter is perpendicular to the true axis above 54.7356 deg).
#ifndef SGKQUADRIC_DETECT_HPP
#define SGKQUADRIC_DETECT_HPP

#include <GeomBase/Dir3D.h>
#include <GeomBase/Point2D.h>
#include <GeomBase/Point3D.h>
#include <GeomBase/UVRange.h>
#include <GeomBase/Vector3D.h>
#include <Geometry/3D/Surface/Surface.h>

#include <vector>

// Bring the sggk geometry types into sgkq here (types.hpp declares no using-declarations).
namespace sgkq
{
using sggk::Dir3D;
using sggk::Point2D;
using sggk::Point3D;
using sggk::Surface;
using sggk::SurfacePtr;
using sggk::UVRange;
using sggk::Vector3D;
} // namespace sgkq

#include "sgkquadric/types.hpp"

namespace sgkq
{

//! Curvature signature of one sample point (the minimal unit of classification).
struct Sample
{
    Point3D  p;              //!< 3D point (Surface::CalcPoint / CalcNormal output)
    Vector3D n;              //!< unit normal (Surface::CalcNormal convention, not necessarily outward)
    double   k1 = 0.0;       //!< principal curvature 1 (Surface::CalcCurvature)
    double   k2 = 0.0;       //!< principal curvature 2
    Point2D  uv;             //!< sample parameters
    bool     degen = false;  //!< degenerate point (pole/apex): normal or curvature unusable, not classified
};

//! Samples a nu x nv grid (including the domain boundary) over uv: point / unit normal /
//! principal curvatures per point; nu/nv below 2 are clamped to 2. Pass the spline surface's OWN
//! domain (DomainU()/DomainV(); ToBSpline reparameterizes U to [0,1]), not the source UV.
//! @return row-major samples (U outer, V inner).
std::vector<Sample> SampleSurface(const Surface& srf, const UVRange& uv, int nu, int nv);

//! Classifies by curvature signature (plane -> sphere -> torus -> cylinder/cone); uv must be the
//! surface's own domain, th is unused (accuracy thresholds, kept for the common interface).
//! @return candidate / residual (normalized invariant residual) / notes.
DetectResult Detect(const Surface& srf, const UVRange& uv, const Thresholds& th = Thresholds{});

} // namespace sgkq

#endif // SGKQUADRIC_DETECT_HPP
