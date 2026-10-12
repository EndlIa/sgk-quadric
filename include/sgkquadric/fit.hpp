// sgkquadric/fit.hpp - parameter fitting: per-type closed-form initial guess, optional
// Gauss-Newton refinement, and closed-form distance verification. Uses math.hpp only; no SGGK
// evaluation functions are called; sample points and normals come from the caller (FitSample).
#ifndef SGKQUADRIC_FIT_HPP
#define SGKQUADRIC_FIT_HPP

#include "sgkquadric/math.hpp"
#include "sgkquadric/types.hpp"

#include <string>
#include <vector>

namespace sgkq
{

//! A sample point and its normal (caller-supplied; normally sampled on the spline's own domain).
//! The normal sign (+/- n) does not affect the result (smallest-eigenvector axis; sign-invariant
//! concurrency/collinearity). A zero-length normal still contributes its point to plane/cylinder
//! fits but is dropped entirely for cone/sphere/torus, which need paired point-normal data.
struct FitSample
{
    Point3D  p;   //!< point on the surface
    Vector3D n;   //!< normal at that point (normalized internally)
};

//! Closed-form parameter estimate (initial guess) for `type`; Unknown is rejected.
//! @param out result (out.type == type; only relevant fields set); @param reason failure text,
//!        non-empty on failure (on success it may carry a non-fatal torus-shape note) - confirm
//!        success by the return value, not by reason.
//! @return false on insufficient samples or degeneracy (half-angle->0, radius->0, low normal rank).
//! Axis estimates: cylinder = minVarDirection + fitCircle2D; cone = minVarDirection (smallest
//! eigenvalue; the largest fails above 54.7356 deg); torus = solvePluckerAxis; sphere =
//! solveLineConcurrency; plane = mean normal.
bool FitSurface(Canonical type,
                const std::vector<FitSample>& samples,
                SurfaceParams& out,
                std::string& reason);

//! Gauss-Newton refinement minimizing the sum of squared closed-form orthogonal distances
//! (numerical-difference Jacobian; step halved when the residual rises).
//! @param io in: initial guess; on success: refined params. @param maxIter cap (<= 0 invalid).
//! @param tol stop when max|dparam| <= tol or |d(sum d^2)| <= tol*max(1, sum d^2).
//! @return false if not converged / singular / invalid (io keeps its input); reason = last residual.
bool RefineSurface(SurfaceParams& io,
                   const std::vector<FitSample>& samples,
                   int maxIter,
                   double tol,
                   std::string& reason);

//! Convenience check: max/mean closed-form orthogonal distance from the sample points to the
//! fitted surface (only samples[i].p is used, not the normals); accepted = (maxDist <= th.maxDist).
//! boundaryCurveHausdorff / g1 fields stay 0 here (filled by a later stage). Empty samples or
//! params.type == Unknown give maxDist = avgDist = 0 and accepted = false.
AccuracyReport MeasureAccuracy(const SurfaceParams& params,
                               const std::vector<FitSample>& samples,
                               const Thresholds& th);

} // namespace sgkq

#endif // SGKQUADRIC_FIT_HPP
