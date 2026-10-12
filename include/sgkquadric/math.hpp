// sgkquadric/math.hpp - closed-form numerical core for surface detection/fitting.
// Uses only the SDK geometry types (Point3D/Vector3D), none of the SDK's evaluation functions
// (CalcNormal / CalcCurvature / ToBSpline). All results are closed-form: Jacobi eigendecomposition
// plus pivoted Gaussian elimination, no iterative optimization.
#ifndef SGKQUADRIC_MATH_HPP
#define SGKQUADRIC_MATH_HPP

#include <array>
#include <vector>

#include <GeomBase/Point3D.h>
#include <GeomBase/Vector3D.h>

namespace sgkq
{
using sggk::Point3D;
using sggk::Vector3D;

//! Eigendecomposition of a 3x3 symmetric matrix; eigenvalues ascending (val[k] pairs with vec[k]).
struct SymEig3
{
    double val[3];      //!< ascending eigenvalues
    Vector3D vec[3];    //!< unit eigenvector for val[k]
};

//! Eigendecomposition of a 6x6 symmetric matrix; eigenvalues ascending (val[k] pairs with vec[k]).
struct SymEig6
{
    double val[6];      //!< ascending eigenvalues
    double vec[6][6];   //!< vec[k] is the unit eigenvector for val[k] (6 components)
};

SymEig3 symEig3(const double m[3][3]);

SymEig6 symEig6(const double m[6][6]);

//! Smallest-eigenvector direction of the covariance of v (the constant-projection direction).
//! The mean is subtracted first; the second moment would give a direction perpendicular to the
//! true cone axis once the half-angle exceeds 54.7356 deg. Empty input returns a zero vector.
Vector3D minVarDirection(const std::vector<Vector3D>& v);

//! Solves A x = b by partial-pivot Gaussian elimination (A row-major); false if A is singular/ill-conditioned.
bool solveLinear3(const double A[3][3], const double b[3], double x[3]);

//! Solves the "common point" problem: the point minimizing the sum of squared orthogonal
//! distances to the lines {p_i, n_i}. Normal equations A = sum(I - n n^T), b = sum(I - n n^T) p,
//! solved as A X = b. Used for sphere centers and cone apexes.
//! @param out the common point; @return false if fewer than 3 points, sizes mismatch, or A singular.
bool solveLineConcurrency(const std::vector<Point3D>& p, const std::vector<Vector3D>& n, Point3D& out);

//! Solves the "common axis" (Plucker) problem: the axis crossing every line {p_i, n_i}. Stacks
//! [(p x n)^T | n^T] per line into M, then takes the smallest eigenvector (a, m) of M^T M; the
//! normalized a is the axis direction. Used for cylinder/torus axes (torus normals meet the axis
//! at varying points, so the common-point form does not apply).
//! @param out axisDir the unit axis direction, axisPoint a point on the axis.
//! @return false if fewer than 3 points, sizes mismatch, or the system is degenerate.
bool solvePluckerAxis(const std::vector<Point3D>& p, const std::vector<Vector3D>& n,
                      Vector3D& axisDir, Point3D& axisPoint);

//! 2D circle fit by algebraic least squares (Kasa): fits x^2 + y^2 + Dx + Ey + F = 0, giving
//! center (-D/2, -E/2) and radius sqrt(cx^2 + cy^2 - F), from points (x, y).
//! @return false if fewer than 3 points or degenerate (collinear / singular normal equations).
bool fitCircle2D(const std::vector<std::array<double, 2>>& pts, double& cx, double& cy, double& r);

//! Orthogonal distance from q to the plane through p0 with unit normal n: |(q - p0).n|.
double distPlane(const Point3D& q, const Point3D& p0, const Vector3D& n);

//! Orthogonal distance from q to an infinite cylinder: |rho - r| (rho = distance to the axis).
double distCylinder(const Point3D& q, const Point3D& axisPt, const Vector3D& axisDir, double r);

//! Orthogonal distance from q to a single-napped cone: rho*cos(a) - h*sin(a) on the side; when
//! the foot falls behind the apex, the distance to the apex is used instead.
double distCone(const Point3D& q, const Point3D& apex, const Vector3D& axisDir, double halfAngle);

//! Orthogonal distance from q to a sphere: | |q - C| - r |.
double distSphere(const Point3D& q, const Point3D& center, double r);

//! Orthogonal distance from q to a torus (major R, minor r): | sqrt((rho - R)^2 + h^2) - r |.
double distTorus(const Point3D& q, const Point3D& center, const Vector3D& axisDir, double R, double r);

} // namespace sgkq

#endif // SGKQUADRIC_MATH_HPP
