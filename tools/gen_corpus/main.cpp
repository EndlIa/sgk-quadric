// gen_corpus -- writes a labeled spline-surface corpus (one .sgt + one manifest.jsonl line
// per case). Positives: analytic surface in a local frame -> random pose -> UV trim ->
// Surface::ToBSpline. Negatives: perturbed cylinder / hand-built free-form spline.
// One seeded std::mt19937 drives everything, so a given --seed reproduces the corpus byte for byte.

#include "sgkquadric/io.hpp"
#include "sgkquadric/testkit.hpp"   // harness: DefaultDataDir

#include <Foundation/init.h>
#include <Foundation/SGGK_Types.h>
#include <Foundation/RuntimeClass.h>

#include <GeomBase/Number.h>
#include <GeomBase/Point3D.h>
#include <GeomBase/Vector3D.h>
#include <GeomBase/Dir3D.h>
#include <GeomBase/Ucs3D.h>
#include <GeomBase/Matrix4.h>
#include <GeomBase/Interval.h>
#include <GeomBase/Matrix3.h>
#include <GeomBase/UVRange.h>

#include <Geometry/SurfaceType.h>
#include <Geometry/3D/Surface/Surface.h>
#include <Geometry/3D/Surface/Plane.h>
#include <Geometry/3D/Surface/CylinderSurface.h>
#include <Geometry/3D/Surface/ConeSurface.h>
#include <Geometry/3D/Surface/SphereSurface.h>
#include <Geometry/3D/Surface/TorusSurface.h>
#include <Geometry/3D/Surface/BSplineSurface.h>

#include <Topology/Brep/Body.h>
#include <Topology/Brep/Lump.h>
#include <Topology/Brep/Shell.h>
#include <Topology/Brep/Face.h>
#include <Topology/Geom/ModelSurface.h>
#include <Topology/Tools/TopoBuilder.h>
#include <Topology/Serialize/RapidTopoJsonSerializer.h>
#include <Topology/Serialize/SerializeOpts.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace sggk;

// Notes shared by the code below: analytic surfaces are built in their own local frame
// (axis = Z) and a random pose is applied via Surface::Transform, so manifest truth follows
// from the pose exactly. Control-point targets need exact knot insertion -- a UV trim alone
// is clamped by the kernel's natural segmentation (baselines <= 56). One shared std::mt19937
// means any added/reordered draw shifts every later case.

namespace
{

// Control-point total target range [50,200]: the lower bound keeps enough freedom for a 1/3
// trim, the upper bound caps the downstream least-squares size.
constexpr int kCtrlLo = 50;
constexpr int kCtrlHi = 200;
// Knot-repeat tolerance (absolute half-width in knot parameter) for the multiplicity ctor
// at BSplineSurface.h:75; 1e-9 sits far below the knot spacing yet absorbs round-off.
constexpr double kKnotTol = 1e-9;

//! @brief Rigid + uniform-scale pose: world point = t + s * R * local point, world dir = R * local dir.
//!        Analytic surfaces can only be built in their own frame (axis = Z, centre = origin),
//!        so the pose is applied afterwards; this keeps truth exactly derivable from the pose.
struct Pose
{
    double  s = 1.0;                       // uniform scale, in [0.5,2.0]
    double  R[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };  // rotation matrix (orthonormal, det=+1), row-major
    Point3D t;                             // world translation in mm, in [-40,40]^3
    Vector3D rotAxis;      // rotation axis (unit); recorded for provenance only
    double  rotAngle = 0.0;  // rotation angle in radians, in [0,2*pi]

    //! @brief Local point -> world: t + s*R*p (all values in mm).
    Point3D worldPoint(const Point3D& p) const
    {
        const double x = R[0][0] * p.X() + R[0][1] * p.Y() + R[0][2] * p.Z();
        const double y = R[1][0] * p.X() + R[1][1] * p.Y() + R[1][2] * p.Z();
        const double z = R[2][0] * p.X() + R[2][1] * p.Y() + R[2][2] * p.Z();
        return Point3D(t.X() + s * x, t.Y() + s * y, t.Z() + s * z);
    }

    //! @brief Local Z axis mapped to world (the analytic axis/normal); R is orthonormal, so no renormalisation.
    Vector3D worldZ() const
    {
        return Vector3D(R[0][2], R[1][2], R[2][2]);
    }

    //! @brief Affine matrix whose columns are s*R*e_i and translation t (homogeneous form of worldPoint).
    Matrix4 toMatrix() const
    {
        const Vector3D cx(s * R[0][0], s * R[1][0], s * R[2][0]);
        const Vector3D cy(s * R[0][1], s * R[1][1], s * R[2][1]);
        const Vector3D cz(s * R[0][2], s * R[1][2], s * R[2][2]);
        return Matrix4(cx, cy, cz, Vector3D(t.X(), t.Y(), t.Z()));
    }
};

//! @brief Uniform random double in [a,b).
//! @param rng shared Mersenne-Twister engine (stateful; every call advances it)
//! @param a   lower bound (inclusive)
//! @param b   upper bound (exclusive); requires a < b, else the distribution is undefined
//! @return the sample; fixing the seed and the call sequence is what makes the corpus reproducible
double URand(std::mt19937& rng, double a, double b)
{
    std::uniform_real_distribution<double> d(a, b);
    return d(rng);
}

//! @brief Normalise a WSL-style "/mnt/<drive>/..." argument to Windows "<drive>:/...".
//!        The exe treats "/mnt/d/x" as relative to the current drive, landing it somewhere
//!        unexpected; other path forms are returned unchanged.
//! @param raw path in any style (no IO, existence is not checked)
//! @return the normalised path; the letter+boundary check leaves "/mnt/dir" untouched
std::string NormalizeOutArg(const std::string& raw)
{
    static const std::string kMnt = "/mnt/";
    if (raw.compare(0, kMnt.size(), kMnt) != 0 || raw.size() <= kMnt.size())
        return raw;

    const char drive = raw[kMnt.size()];
    const bool isLetter = (drive >= 'a' && drive <= 'z') || (drive >= 'A' && drive <= 'Z');
    const size_t rest = kMnt.size() + 1;
    const bool isBoundary = (rest == raw.size()) || raw[rest] == '/';  // rejects "/mnt/dir"
    if (!isLetter || !isBoundary)
        return raw;

    std::string out;
    out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(drive))));
    out += ':';
    out += (rest == raw.size()) ? "/" : raw.substr(rest);
    return out;
}

//! @brief Random rigid rotation (uniform axis + angle, Rodrigues) plus uniform scale and translation.
//!        Standard-normal axis components give an isotropic direction on the sphere; uniform
//!        [0,1] components would bias toward the corners. Near-zero vectors are rejected and
//!        resampled (guard cap 100, practically never reached).
//! @param rng engine; the result depends on its current state (see URand)
//! @return a Pose with s in [0.5,2.0], t in [-40,40]^3, angle in [0,2*pi)
Pose MakePose(std::mt19937& rng)
{
    Pose pose;

    std::normal_distribution<double> nd(0.0, 1.0);
    Vector3D axis(0.0, 0.0, 1.0);           // fallback axis if all 100 samples are degenerate
    for (int guard = 0; guard < 100; ++guard)
    {
        Vector3D cand(nd(rng), nd(rng), nd(rng));  // three standard normals -> uniform sphere direction
        if (cand.Length() > 1e-3)                  // reject near-zero vectors before normalising
        {
            cand.Normalize();
            axis = cand;
            break;
        }
    }

    const double angle = URand(rng, 0.0, PI2);  // radians
    const double ax = axis.X(), ay = axis.Y(), az = axis.Z();
    const double c = std::cos(angle), sn = std::sin(angle), k = 1.0 - c;  // Rodrigues coefficients

    // standard Rodrigues matrix for rotation by 'angle' about unit axis (ax,ay,az), row-major
    pose.R[0][0] = c + ax * ax * k;      pose.R[0][1] = ax * ay * k - az * sn; pose.R[0][2] = ax * az * k + ay * sn;
    pose.R[1][0] = ay * ax * k + az * sn; pose.R[1][1] = c + ay * ay * k;      pose.R[1][2] = ay * az * k - ax * sn;
    pose.R[2][0] = az * ax * k - ay * sn; pose.R[2][1] = az * ay * k + ax * sn; pose.R[2][2] = c + az * az * k;

    pose.s = URand(rng, 0.5, 2.0);
    pose.t = Point3D(URand(rng, -40.0, 40.0), URand(rng, -40.0, 40.0), URand(rng, -40.0, 40.0));  // mm
    pose.rotAxis = axis;
    pose.rotAngle = angle;
    return pose;
}

//! @brief Control-point count of a BSplineSurface (matrix is [U][V], rows equal-length).
//! @param bs surface to inspect; may be null
//! @param nu [out] U direction size (row count)
//! @param nv [out] V direction size (columns of the first row)
//! @return nu*nv; 0 when bs is null or the U direction is empty
int CountCtrl(const BSplineSurfacePtr& bs, int& nu, int& nv)
{
    nu = 0;
    nv = 0;
    if (!bs)
        return 0;
    const Point3DMatrix& cps = bs->ControlPoints();  // BSplineSurface.h:337
    nu = static_cast<int>(cps.size());
    if (nu > 0)
        nv = static_cast<int>(cps[0].size());
    return nu * nv;
}

//! @brief Desired control-point total for this case, spread over 50-200 by idx%3:
//!        tier 0 -> 50-62, tier 1 -> 105-125, tier 2 -> 170-190 (perType=8 gives 3/3/2 per tier).
//!        A UV trim alone is clamped by the kernel's natural segmentation, so the tiers can only
//!        be hit together with knot insertion (RefineToTarget decides the achievable total).
//! @param idx position within the type; only idx%3 matters
//! @param rng engine (each tier draws uniformly inside its band)
//! @return the wanted total; whether it is reached is up to RefineToTarget
int TargetCtrlCount(int idx, std::mt19937& rng)
{
    const int tier = ((idx % 3) + 3) % 3;  // normalise to [0,3) for negative idx
    if (tier == 0)
        return static_cast<int>(std::lround(URand(rng, 50.0, 62.0)));
    if (tier == 1)
        return static_cast<int>(std::lround(URand(rng, 105.0, 125.0)));
    return static_cast<int>(std::lround(URand(rng, 170.0, 190.0)));
}

//! @brief Score a candidate range; smaller is better.
//! @param count  control-point total of the transcribed candidate
//! @param target this case's tier target
//! @return |count-target| when count <= kCtrlHi; above kCtrlHi a hard 1e9 penalty is added so no
//!         over-limit range can beat a legal one even by a small margin
long long CandidateScore(int count, int target)
{
    if (count > kCtrlHi)
        return 1000000000LL + static_cast<long long>(count - kCtrlHi);
    return count >= target ? static_cast<long long>(count - target)
                           : static_cast<long long>(target - count);
}

//! @brief Raise the control-point total to target by exact knot insertion (geometry unchanged),
//!        never past kCtrlHi.
//!        One U knot adds nv to the total (one V knot adds nu), so the smaller side is grown
//!        first: bigger per-step gain and the two counts stay balanced. Insert parameters come
//!        from a golden-ratio sequence, giving distinct values well spread over the domain.
//!        Insertions that fail (multiplicity cap, etc.) are caught and retried with a new value.
//! @param bs     surface to densify in place (shape is invariant under knot insertion)
//! @param target desired control-point total
//! @return the achieved total, which may fall short on overshoot protection or repeated failures
int RefineToTarget(const BSplineSurfacePtr& bs, int target)
{
    if (!bs)
        return 0;

    int kU = 0;  // golden-ratio sequence counters, one per direction
    int kV = 0;
    for (int guard = 0; guard < 200000; ++guard)  // hard cap against runaway growth
    {
        int nu = 0, nv = 0;
        const int total = CountCtrl(bs, nu, nv);
        if (total >= target || nu <= 0 || nv <= 0)
            break;  // reached target, or the grid is empty

        const int addU = nv;  // InsertUKnots: nu+1 => total + nv
        const int addV = nu;  // InsertVKnots: nv+1 => total + nu

        bool doU = (nu <= nv);  // grow the smaller side
        if (total + (doU ? addU : addV) > kCtrlHi)
        {
            if (total + (doU ? addV : addU) > kCtrlHi)
                break;          // both directions would overshoot kCtrlHi
            doU = !doU;
        }

        bool done = false;
        for (int attempt = 0; attempt < 16 && !done; ++attempt)  // up to 16 reparameterisations
        {
            const Interval dom = doU ? bs->DomainU() : bs->DomainV();
            const double a = dom.Min(), b = dom.Max();
            if (!(b > a))
                break;  // degenerate domain (zero length)
            const int k = doU ? (kU++) : (kV++);
            // golden-ratio sequence: fractional parts of k*0.618... cover (0,1) without repeats;
            // mapped into [0.03,0.97] to stay off the ends (end clamped multiplicity is full).
            double frac = std::fmod(0.5 + k * 0.6180339887498949, 1.0);
            frac = 0.03 + 0.94 * frac;
            const double param = a + (b - a) * frac;  // absolute knot parameter
            try
            {
                if (doU)
                    bs->InsertUKnots(param, 1);  // BSplineSurface.h:687
                else
                    bs->InsertVKnots(param, 1);  // BSplineSurface.h:690
                done = true;
            }
            catch (...)
            {
                // insertion rejected here (e.g. multiplicity cap): try the next parameter
            }
        }
        if (!done)
            break;
    }

    int fu = 0, fv = 0;
    return CountCtrl(bs, fu, fv);
}

//! @brief Build a clamped knot vector in "with multiplicity" form (repeats written out), as the
//!        ctor at Geometry/3D/Surface/BSplineSurface.h:75 expects. Kernel rule: for n control
//!        points and degree p the vector length is n+p+1, with p+1 repeats at each end and
//!        n-p-1 interior knots splitting [knotMin,knotMax] into n-p equal spans.
//!        Example n=8, p=3: length 12 = [0,0,0,0, 1,2,3,4, 5,5,5,5] over the domain [0,5].
//! @param n       control points in this direction (>= 2); sets the span count n-p
//! @param p       degree in this direction (1 <= p < n)
//! @param knotMin domain start (0 here)
//! @param knotMax domain end (5 here)
//! @return a length-(n+p+1) array; interior knots are equally spaced, so the parametrisation is uniform
RealArray MakeClampedKnots(int n, int p, double knotMin, double knotMax)
{
    RealArray knots;
    knots.reserve(static_cast<size_t>(n + p + 1));

    for (int k = 0; k < p + 1; ++k)
        knots.push_back(knotMin);                   // p+1 repeats at the start (clamped)

    const int spans = n - p;                        // number of spans
    for (int s = 1; s < spans; ++s)                 // interior knots = spans - 1
        knots.push_back(knotMin + (knotMax - knotMin) * (static_cast<double>(s) / spans));

    for (int k = 0; k < p + 1; ++k)
        knots.push_back(knotMax);                   // p+1 repeats at the end (clamped)
    return knots;
}

//! @brief Analytic surface UV domain (two 1-D intervals); a finite fallback is used when the
//!        surface reports an unbounded (periodic / infinite-axis) parameterisation.
struct Domain
{
    Interval u;  // U range (radians for angular surfaces, mm for length-like ones); finite
    Interval v;  // V range; finite
};

//! @brief Finite UV domain of an analytic surface; an unbounded report is replaced per type:
//!        cylinder -> a 2R-long axial segment; cone -> v up to radius ~2, i.e. v = 2/tan(halfAngle)
//!        (guarded by max(1e-3, tan)); anything else -> a fixed [-1,1]. U falls back to [0,2*pi].
//! @param surf          target surface (must be non-null)
//! @param type          surface type; decides the V fallback
//! @param baseCylR      cylinder radius (mm), used as the axial trim length (= 2R)
//! @param coneHalfAngle cone half-angle (rad), used to convert an axial distance to a radius
//! @return a Domain whose two intervals are finite
Domain ResolveDomain(const SurfacePtr& surf, SurfaceType type, double baseCylR, double coneHalfAngle)
{
    Domain d;
    d.u = surf->DomainU();  // Surface.h:119
    d.v = surf->DomainV();  // Surface.h:124

    if (d.u.GetType() != IntervalType::FINITE)  // Interval.h:156 / :36
        d.u = Interval(0.0, PI2);               // fallback: full turn (radians)

    if (d.v.GetType() != IntervalType::FINITE)
    {
        switch (type)
        {
        case SurfaceType::Cylinder:
            d.v = Interval(0.0, 2.0 * baseCylR);
            break;
        case SurfaceType::Cone:
            // radius(v) = v * tan(halfAngle); pick v so the radius is ~2
            d.v = Interval(0.0, 2.0 / std::max(1e-3, std::tan(coneHalfAngle)));
            break;
        default:
            d.v = Interval(-1.0, 1.0);  // plane / sphere / torus: fixed unit fallback
            break;
        }
    }
    return d;
}

//! @brief Sub-interval of length f*L with a random start, always fully inside full.
//! @param full full domain; returned unchanged if Max <= Min
//! @param f    length ratio in (0,1]
//! @param rng  engine; picks the offset so trimmed blocks are not all in the same place
//! @return the sub-interval (same parameter units); the start stays in [a, a+(1-f)L]
Interval TrimInterval(const Interval& full, double f, std::mt19937& rng)
{
    const double a = full.Min(), b = full.Max();
    const double L = b - a;
    if (!(L > 0.0))
        return full;  // degenerate/reversed interval: no trim
    const double len = f * L;
    const double start = a + URand(rng, 0.0, 1.0) * (L - len);  // slack = L-len
    return Interval(start, start + len);
}

//! @brief Assemble the minimal topology (Face->Shell->Lump->Body) and serialize it to path.
//!        Body>Lump>Shell>Face is the minimum the kernel accepts; Face sense is true (outward).
//! @param geomSurf geometric surface (must be non-null); ownership passes to MakeModelSurface
//! @param path      destination .sgt path (already /mnt-normalised by the caller)
//! @note serialization failures are swallowed inside RapidTopoJsonSerializer and not checked here
void SerializeSurface(const SurfacePtr& geomSurf, const std::string& path)
{
    ModelSurfacePtr ms = TopoBuilder::MakeModelSurface(geomSurf, false);  // TopoBuilder.h:88 (2nd arg: non-periodic)
    FacePtr   face  = TopoBuilder::MakeFace(ms, true);                     // TopoBuilder.h:248 (outward normal)
    ShellPtr  shell = TopoBuilder::MakeShell();                            // TopoBuilder.h:289
    TopoBuilder::ShellAddFace(shell, face);                                // TopoBuilder.h:319
    LumpPtr   lump  = TopoBuilder::MakeLump();                             // TopoBuilder.h:333
    TopoBuilder::LumpAddShell(lump, shell);                                // TopoBuilder.h:343
    BodyPtr   body  = TopoBuilder::MakeBody();                             // TopoBuilder.h:347
    TopoBuilder::BodyAddLump(body, lump);                                  // TopoBuilder.h:367

    TopologyPtr topo = body;
    RapidTopoJsonSerializer serializer;
    serializer.Serialize(topo, path.c_str(), true);                        // TopoJsonSerializer.h:40 (3rd: pretty)
}

//! @brief double -> JSON number, 17 significant digits (round-trip lossless).
//! @param x value to serialize; non-finite is replaced by 0.0 (JSON has no NaN/Inf literals)
//! @return a bare decimal string, no quotes
std::string JNum(double x)
{
    if (!std::isfinite(x))
        x = 0.0;  // JSON has no NaN/Inf literals
    std::ostringstream os;
    os << std::setprecision(17) << x;
    return os.str();
}

//! @brief Point -> JSON triple.
//! @param p point (mm)
//! @return "[x,y,z]" (elements via JNum)
std::string JVec(const Point3D& p)
{
    return "[" + JNum(p.X()) + "," + JNum(p.Y()) + "," + JNum(p.Z()) + "]";
}

//! @brief Vector -> JSON triple.
//! @param v direction/axis vector (unit or arbitrary length, by caller convention)
//! @return "[x,y,z]" (elements via JNum)
std::string JVec(const Vector3D& v)
{
    return "[" + JNum(v.X()) + "," + JNum(v.Y()) + "," + JNum(v.Z()) + "]";
}

//! @brief Corpus class label; Unknown marks a negative. Values double as stats[] indices, so
//!        Plane must stay 0 and the order must not change (see stats[6] in main).
enum class Kind
{
    Plane = 0,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    Unknown
};

//! @brief Enum -> manifest type/expected name.
//! @param k type
//! @return a static string literal (Unknown -> "Unknown")
const char* KindName(Kind k)
{
    switch (k)
    {
    case Kind::Plane:    return "Plane";
    case Kind::Cylinder: return "Cylinder";
    case Kind::Cone:     return "Cone";
    case Kind::Sphere:   return "Sphere";
    case Kind::Torus:    return "Torus";
    default:             return "Unknown";
    }
}

//! @brief Enum -> file id prefix (e.g. "plane"/"cyl").
//! @param k type
//! @return the prefix; Unknown -> "neg"
std::string PrefixOf(Kind k)
{
    switch (k)
    {
    case Kind::Plane:    return "plane";
    case Kind::Cylinder: return "cyl";
    case Kind::Cone:     return "cone";
    case Kind::Sphere:   return "sph";
    case Kind::Torus:    return "tor";
    default:             return "neg";
    }
}

//! @brief World-space analytic truth for a positive, written into the manifest pose section.
//!        All values are post-transform (the radius is read back from the surface, not scaled here).
//!        plane: axisPoint+axisDir; cylinder: +radius; cone: apex+axisDir+halfAngle;
//!        sphere: axisPoint(=centre)+radius; torus: axisPoint(=centre)+axisDir+radius(R)+minorRadius(r).
struct Truth
{
    Point3D  axisPoint;    // point on the plane / on the axis / sphere centre / torus centre / apex (mm)
    Vector3D axisDir;      // unit normal / axis (from Pose::worldZ)
    double   radius = 0.0;      // cylinder / sphere radius, or torus major radius R (mm)
    double   minorRadius = 0.0; // torus minor radius r (mm); 0 for other types
    double   halfAngle = 0.0;   // cone half-angle (rad); 0 for other types
    Point3D  apex;              // cone apex (mm); equals axisPoint for other types
};

//! @brief Per-type generation counters for the ending [CORPUS] summary. Raw counts only (no
//!        derived "achievement rate") so warnings cannot be papered over.
struct Stat
{
    int cases = 0;         // cases generated
    int full = 0;          // cases using the whole domain
    int trimmed = 0;       // cases using a trimmed UV range
    int outOfRange = 0;    // final count fell outside [kCtrlLo,kCtrlHi] (warned)
    int belowTarget = 0;   // legal but below the tier target (warned)
    int ctrlMin = 0;       // min control-point count (0 = no sample yet)
    int ctrlMax = 0;       // max control-point count
};

//! @brief Generate one positive: local analytic surface -> random pose -> UV range -> ToBSpline
//!        -> knot-insertion correction -> write .sgt and one manifest line.
//! @param kind    positive type (not Unknown; otherwise returns false)
//! @param idx     index within the type; idx%2 chooses full/trimmed, idx%3 the tier
//! @param id      case id (e.g. "cyl_003"), used as the .sgt name and the manifest key
//! @param rng     shared engine; consumed here more than anywhere else, so any added or reordered
//!                draw shifts every later case
//! @param outDir  directory the .sgt is written to (the caller writes the manifest)
//! @param seed    run seed, recorded verbatim in the manifest
//! @param kindOut [out] "full" or "trimmed"
//! @param err     [out] failure reason, printed by the caller in the [SKIP] line
//! @return true once the .sgt is written and manifestLine is valid; false means the case was skipped
bool GeneratePositive(Kind kind, int idx, const std::string& id, std::mt19937& rng,
                      const std::string& outDir, uint32_t seed, std::string& manifestLine,
                      int& ctrlOut, int& targetOut, std::string& kindOut, std::string& err)
{
    try
    {
        Pose pose = MakePose(rng);

        // radii kept in 10-60 mm so downstream can use one absolute tolerance; cone half-angle in
        // [12,55] deg (smaller nears a cylinder, larger nears a plane, weakening the type signal)
        SurfacePtr geom;
        double baseCylR = 0.0;
        double coneHalf = 0.0;
        double baseMajorR = 0.0;
        double baseMinorR = 0.0;

        switch (kind)
        {
        case Kind::Plane:
            geom = std::make_shared<Plane>(Point3D(0.0, 0.0, 0.0), Dir3D(0.0, 0.0, 1.0));  // Plane.h:42
            break;
        case Kind::Cylinder:
            baseCylR = URand(rng, 10.0, 60.0);
            geom = std::make_shared<CylinderSurface>(Ucs3D(), baseCylR);                  // CylinderSurface.h:33
            break;
        case Kind::Cone:
            // ctor arg is the height at radius 1, i.e. 1/tan(halfAngle)
            coneHalf = URand(rng, 12.0 * PI / 180.0, 55.0 * PI / 180.0);
            geom = std::make_shared<ConeSurface>(Ucs3D(), 1.0 / std::tan(coneHalf));      // ConeSurface.h:34
            break;
        case Kind::Sphere:
            geom = std::make_shared<SphereSurface>(Ucs3D(), URand(rng, 10.0, 60.0));      // SphereSurface.h:33
            break;
        case Kind::Torus:
        {
            baseMajorR = URand(rng, 10.0, 50.0);
            // half doughnuts (r<R, regular) and half apples (r>=R, self-intersecting at the
            // equator with a ridge); decoupled from full/trimmed via (idx/2)%2
            const bool apple = ((idx / 2) % 2) == 0;
            baseMinorR = apple ? baseMajorR * URand(rng, 1.0, 1.15)    // r in [1.00,1.15]R: apple
                               : baseMajorR * URand(rng, 0.35, 0.80);  // r in [0.35,0.80]R: doughnut
            geom = std::make_shared<TorusSurface>(Ucs3D(), baseMajorR, baseMinorR);       // TorusSurface.h:45
            break;
        }
        default:
            err = "unsupported positive kind";
            return false;
        }

        const Matrix4 trsf = pose.toMatrix();
        geom->Transform(trsf);  // Surface.h:365 (in place; surface is now the final world geometry)

        // read radii back from the surface: Transform may normalise/snap them, so the surface's
        // own value is what matches the .sgt. Falls back to the analytic value if sg_cast fails.
        Truth truth;
        truth.axisPoint = pose.worldPoint(Point3D(0.0, 0.0, 0.0));  // world image of the local origin
        truth.axisDir = pose.worldZ();
        truth.apex = truth.axisPoint;
        switch (kind)
        {
        case Kind::Plane:
            break;
        case Kind::Cylinder:
        {
            const CylinderSurfacePtr cy = sg_cast<CylinderSurface>(geom);
            truth.radius = cy ? cy->Radius() : pose.s * baseCylR;  // CylinderSurface.h:69
            break;
        }
        case Kind::Cone:
            truth.halfAngle = coneHalf;  // rigid + uniform scale leaves the half-angle unchanged
            break;
        case Kind::Sphere:
        {
            const SphereSurfacePtr sp = sg_cast<SphereSurface>(geom);
            truth.radius = sp ? sp->Radius() : 0.0;  // SphereSurface.h:58
            break;
        }
        case Kind::Torus:
        {
            const TorusSurfacePtr ts = sg_cast<TorusSurface>(geom);
            truth.radius = ts ? ts->MajorRadius() : pose.s * baseMajorR;        // TorusSurface.h:81
            truth.minorRadius = ts ? ts->MinorRadius() : pose.s * baseMinorR;   // TorusSurface.h:97
            break;
        }
        default:
            break;
        }

        const Domain dom = ResolveDomain(geom, geom->SurfaceType(), baseCylR, coneHalf);

        const bool wantFull = (idx % 2) == 0;  // even: full domain, odd: trimmed block
        kindOut = wantFull ? "full" : "trimmed";

        // full domain -> one candidate; trimmed -> four ratios in [1/3,1/2], best one is picked
        std::vector<std::pair<Interval, Interval>> candidates;
        if (wantFull)
        {
            candidates.emplace_back(dom.u, dom.v);
        }
        else
        {
            const double f0 = URand(rng, 1.0 / 3.0, 0.5);  // first ratio random, rest fixed for coverage
            const double fracs[4] = { f0, 0.5, 0.4, 1.0 / 3.0 };
            for (int c = 0; c < 4; ++c)
                candidates.emplace_back(TrimInterval(dom.u, fracs[c], rng), TrimInterval(dom.v, fracs[c], rng));
        }

        const int target = TargetCtrlCount(idx, rng);
        targetOut = target;

        BSplineSurfacePtr best;                 // best-scoring transcription (gets written)
        Interval bestU = dom.u, bestV = dom.v;  // recorded as the manifest uvRange
        int bestCount = -1;                     // selected control-point count (-1 = none yet)
        long long bestScore = 0;

        for (const auto& cand : candidates)
        {
            BSplineSurfacePtr bs;
            try
            {
                UVRange range(Interval(cand.first.Min(), cand.first.Max()),
                              Interval(cand.second.Min(), cand.second.Max()));  // UVRange.h:34
                bs = geom->ToBSpline(range);                                       // Surface.h:340
            }
            catch (...)
            {
                continue;  // this range failed to transcribe -> try the next candidate
            }
            if (!bs)
                continue;  // kernel returned a null surface -> try the next candidate

            int nu = 0, nv = 0;
            const int count = CountCtrl(bs, nu, nv);
            const long long score = CandidateScore(count, target);
            if (bestCount < 0 || score < bestScore)  // strictly better wins; first round always taken
            {
                best = bs;
                bestU = cand.first;
                bestV = cand.second;
                bestCount = count;
                bestScore = score;
            }
        }

        if (!best)
        {
            err = "ToBSpline returned null for all candidate ranges";
            return false;
        }

        // exact knot insertion (geometry unchanged) raises the count toward the tier target
        const int finalCount = RefineToTarget(best, target);
        int nu = 0, nv = 0;
        CountCtrl(best, nu, nv);  // final U/V split, recorded for cross-checking
        ctrlOut = finalCount;

        const std::string sgtRel = id + ".sgt";            // relative name stored in the manifest
        const std::string sgtAbs = outDir + "/" + sgtRel;
        SerializeSurface(best, sgtAbs);

        // one JSONL line: pose is the truth, expected the label the detector should return; seed is
        // recorded verbatim so the .sgt is reproducible from the fixed rng sequence
        std::ostringstream os;
        os << "{"
           << "\"id\":\"" << id << "\","
           << "\"kind\":\"" << kindOut << "\","
           << "\"type\":\"" << KindName(kind) << "\","
           << "\"ctrlCount\":" << finalCount << ","
           << "\"ctrlTarget\":" << target << ","
           << "\"ctrlU\":" << nu << ","
           << "\"ctrlV\":" << nv << ","
           << "\"uvRange\":[" << JNum(bestU.Min()) << "," << JNum(bestU.Max()) << ","
                            << JNum(bestV.Min()) << "," << JNum(bestV.Max()) << "],"
           << "\"seed\":" << seed << ","
           << "\"pose\":{"
           <<   "\"axisPoint\":" << JVec(truth.axisPoint) << ","
           <<   "\"axisDir\":" << JVec(truth.axisDir) << ","
           <<   "\"radius\":" << JNum(truth.radius) << ","
           <<   "\"minorRadius\":" << JNum(truth.minorRadius) << ","
           <<   "\"halfAngle\":" << JNum(truth.halfAngle) << ","
           <<   "\"apex\":" << JVec(truth.apex) << ","
           <<   "\"scale\":" << JNum(pose.s) << ","
           <<   "\"translation\":" << JVec(pose.t) << ","
           <<   "\"rotAxis\":" << JVec(pose.rotAxis) << ","
           <<   "\"rotAngleRad\":" << JNum(pose.rotAngle)
           << "},"
           << "\"expected\":\"" << KindName(kind) << "\","
           << "\"sgt\":\"" << sgtRel << "\""
           << "}";
        manifestLine = os.str();
        return true;
    }
    catch (const std::exception& e)
    {
        err = std::string("exception: ") + e.what();  // caller prints it in the [SKIP] line
        return false;
    }
    catch (...)
    {
        err = "unknown exception";
        return false;
    }
}

//! @brief Generate one negative: even idx -> a perturbed cylinder, odd idx -> a hand-built
//!        free-form spline. Both look regular but are not exact quadrics, so a detector should
//!        reject them (expected "Unknown"). Draws from the same shared RNG as the positives.
//! @param idx     negative index; even -> perturbed cylinder, odd -> free-form spline
//! @param subType [out] "PerturbedCylinder" or "FreeForm", for the [SKIP] line
//! @param err     [out] failure reason
//! @return true once the .sgt is written and manifestLine is valid; false means the case was skipped
bool GenerateNegative(int idx, const std::string& id, std::mt19937& rng,
                      const std::string& outDir, uint32_t seed, std::string& manifestLine,
                      int& ctrlOut, std::string& subType, std::string& err)
{
    try
    {
        Pose pose = MakePose(rng);
        const Matrix4 trsf = pose.toMatrix();

        SurfacePtr geom;
        if ((idx % 2) == 0)
        {
            // transcribe a cylinder exactly, then jitter every control point: the surface still
            // looks like a cylinder but no longer keeps equal distance to the axis, so a quadric
            // fit exceeds tolerance -> Unknown
            subType = "PerturbedCylinder";
            const double r0 = URand(rng, 20.0, 60.0);
            SurfacePtr base = std::make_shared<CylinderSurface>(Ucs3D(), r0);  // CylinderSurface.h:33
            base->Transform(trsf);                                            // Surface.h:365

            UVRange range(Interval(0.0, PI2), Interval(-r0, r0));          // full turn x axial 2*r0
            BSplineSurfacePtr bs = base->ToBSpline(range);                    // Surface.h:340
            if (!bs)
            {
                err = "negative A: ToBSpline null";
                return false;
            }

            const Point3DMatrix cps = bs->ControlPoints();                    // BSplineSurface.h:337 (copy; write-safe)
            const int cu = static_cast<int>(cps.size());
            const int cv = cu > 0 ? static_cast<int>(cps[0].size()) : 0;      // rows are equal-length
            const double amp = 0.12 * pose.s * r0;                            // ~12% of the radius
            for (int i = 0; i < cu; ++i)
            {
                for (int j = 0; j < cv; ++j)
                {
                    const Point3D p = cps[i][j];
                    // independent uniform jitter per component (mm): enough to break a quadric fit,
                    // small enough not to change the overall scale
                    const double dx = URand(rng, -amp, amp);
                    const double dy = URand(rng, -amp, amp);
                    const double dz = URand(rng, -amp, amp);
                    bs->SetControlPoint(Point3D(p.X() + dx, p.Y() + dy, p.Z() + dz),
                                        static_cast<UInt>(i), static_cast<UInt>(j));  // BSplineSurface.h:282
                }
            }
            geom = bs;
        }
        else
        {
            // free-form spline: z is a sum of two sinusoid sets, so the surface has no quadric
            // structure and a standard fit leaves a large residual -> Unknown
            subType = "FreeForm";
            const int deg = 3;
            const int n = 8;    // 8x8 = 64 control points, inside [50,200]
            Point3DMatrix cps(static_cast<size_t>(n), std::vector<Point3D>(static_cast<size_t>(n)));
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    const double x = 10.0 * i;  // 10 mm grid step
                    const double y = 10.0 * j;
                    // z (mm): amplitudes 6 and 3 give a non-monotone relief
                    const double z = 6.0 * std::sin(1.1 * i + 0.7) * std::cos(0.9 * j + 0.3)
                                   + 3.0 * std::cos(0.6 * i - 0.4 * j);
                    cps[static_cast<size_t>(i)][static_cast<size_t>(j)] = Point3D(x, y, z);
                }
            }

            // knot vector length must be n+p+1 (= 12 for n=8, p=3); a 14-knot array made the
            // kernel infer degree 5 from size-n-1 and reject the degree-3 call. See
            // Geometry/3D/Surface/BSplineSurface.h:75 (multiplicity ctor, tolerance-based).
            const RealArray knotsU = MakeClampedKnots(n, deg, 0.0, 5.0);  // domain [0,5]
            const RealArray knotsV = MakeClampedKnots(n, deg, 0.0, 5.0);

            BSplineSurfacePtr free = std::make_shared<BSplineSurface>(
                static_cast<UInt>(deg), static_cast<UInt>(deg), cps, knotsU, knotsV,
                kKnotTol, kKnotTol);                                          // BSplineSurface.h:75
            Matrix3 uvTrsf;                                                   // UV-space transform filled by Transform
            free->Transform(trsf, uvTrsf);                                    // BSplineSurface.h:476
            geom = free;
        }

        const BSplineSurfacePtr bs = sg_cast<BSplineSurface>(geom);  // both negatives are splines
        int nu = 0, nv = 0;
        ctrlOut = CountCtrl(bs, nu, nv);   // recorded only; negatives are not bound to [50,200]

        const std::string sgtRel = id + ".sgt";
        const std::string sgtAbs = outDir + "/" + sgtRel;
        SerializeSurface(geom, sgtAbs);

        // negatives have no meaningful uvRange (written as [0,0,0,0]); pose is filled for provenance
        std::ostringstream os;
        os << "{"
           << "\"id\":\"" << id << "\","
           << "\"kind\":\"negative\","
           << "\"type\":\"" << subType << "\","
           << "\"ctrlCount\":" << ctrlOut << ","
           << "\"ctrlU\":" << nu << ","
           << "\"ctrlV\":" << nv << ","
           << "\"uvRange\":[0,0,0,0],"
           << "\"seed\":" << seed << ","
           << "\"pose\":{"
           <<   "\"axisPoint\":" << JVec(pose.t) << ","
           <<   "\"axisDir\":" << JVec(pose.worldZ()) << ","
           <<   "\"radius\":0,\"minorRadius\":0,\"halfAngle\":0,"
           <<   "\"apex\":" << JVec(pose.t) << ","
           <<   "\"scale\":" << JNum(pose.s) << ","
           <<   "\"translation\":" << JVec(pose.t) << ","
           <<   "\"rotAxis\":" << JVec(pose.rotAxis) << ","
           <<   "\"rotAngleRad\":" << JNum(pose.rotAngle)
           << "},"
           << "\"expected\":\"Unknown\","
           << "\"sgt\":\"" << sgtRel << "\""
           << "}";
        manifestLine = os.str();
        return true;
    }
    catch (const std::exception& e)
    {
        // include the sub-type so the [SKIP] line identifies the failing variant
        err = (subType.empty() ? std::string("negative") : subType) + ": exception: " + e.what();
        return false;
    }
    catch (...)
    {
        err = (subType.empty() ? std::string("negative") : subType) + ": unknown exception";
        return false;
    }
}

//! @brief Print the command-line usage (the text notes that --out accepts either path style).
void PrintUsage()
{
    std::cout << "usage: gen_corpus [--out <dir>] [--seed <int>] [--cases-per-type <int>]\n";
    std::cout << "  --out  output dir, default <project root>/data/corpus; accepts both D:/x/y and /mnt/d/x/y\n"
                 "         (/mnt/<drive>/... is normalized to <drive>:/... before use)\n";
}

} // namespace

//! @brief Entry point: parse args, generate positives and negatives, write manifest.jsonl, print
//!        the summary. Exit codes: 0 ok, 1 manifest cannot be opened, 2 uncaught exception.
//!        Foundation::init/fini are paired; fini also runs on the manifest early-exit path.
int main(int argc, char** argv)
{
    namespace fs = std::filesystem;

    std::string outDir = sgkq::DefaultDataDir() + "/corpus";  // default; overridable by --out
    uint32_t    seed = 20261008u;         // fixed default so a bare run is reproducible
    int         perType = 8;              // positives per type (5*8 = 40)

    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc)
            outDir = NormalizeOutArg(argv[++i]);  // accepts /mnt/<drive>/... or <drive>:/...
        else if (a == "--seed" && i + 1 < argc)
            seed = static_cast<uint32_t>(std::stoul(argv[++i]));  // invalid input throws, handled below
        else if (a == "--cases-per-type" && i + 1 < argc)
            perType = std::stoi(argv[++i]);
        else if (a == "--help" || a == "-h")
        {
            PrintUsage();
            return 0;
        }
        else
        {
            std::cout << "[WARN] unknown arg: " << a << "\n";
        }
    }
    if (perType < 1)
        perType = 1;  // at least one, to avoid an empty run

    // echo the resolved output path so the caller can confirm where files land
    {
        std::error_code ecAbs;
        const fs::path absOut = fs::absolute(fs::path(outDir), ecAbs);
        std::cout << "[OUT] " << (ecAbs ? outDir : absOut.lexically_normal().string()) << "\n";
        std::cout.flush();  // flush so a piping caller sees it immediately
    }

    try
    {
        init();  // Foundation/init.h:25

        std::error_code ec;
        fs::create_directories(outDir, ec);  // idempotent; a failure surfaces at the first write

        std::mt19937 rng(seed);  // one engine for the whole run: any reorder changes the corpus

        const Kind order[5] = { Kind::Plane, Kind::Cylinder, Kind::Cone, Kind::Sphere, Kind::Torus };
        Stat stats[6];  // indices match Kind values (Unknown = 5)

        std::vector<std::string> lines;               // all manifest lines, written at the end
        int typeCounter[6] = { 0, 0, 0, 0, 0, 0 };    // per-type id counters
        int negCounter = 0;                           // separate negative id counter

        // positives: fixed type order (the order fixes RNG consumption = reproducibility)
        for (int t = 0; t < 5; ++t)
        {
            const Kind kind = order[t];
            for (int i = 0; i < perType; ++i)
            {
                char idbuf[64];  // holds "<prefix>_<3 digits>"
                std::snprintf(idbuf, sizeof(idbuf), "%s_%03d", PrefixOf(kind).c_str(), typeCounter[t]++);
                const std::string id = idbuf;

                std::string line, kindOut, err;
                int ctrl = 0;
                int target = 0;
                if (!GeneratePositive(kind, i, id, rng, outDir, seed, line, ctrl, target, kindOut, err))
                {
                    std::cout << "[SKIP] id=" << id << " reason=" << err << "\n";
                    continue;  // one failed case does not abort the batch
                }

                Stat& st = stats[static_cast<int>(kind)];
                st.cases++;
                if (kindOut == "full")
                    st.full++;
                else
                    st.trimmed++;
                if (st.ctrlMin == 0 || ctrl < st.ctrlMin)  // 0 doubles as the "no sample yet" sentinel
                    st.ctrlMin = ctrl;
                if (ctrl > st.ctrlMax)
                    st.ctrlMax = ctrl;
                if (ctrl < kCtrlLo || ctrl > kCtrlHi)
                {
                    st.outOfRange++;
                    std::cout << "[WARN] id=" << id << " ctrlCount=" << ctrl
                              << " target=" << target
                              << " outside [" << kCtrlLo << "," << kCtrlHi << "]\n";
                }
                else if (ctrl < target)
                {
                    // tier target missed (kernel transcription/insertion limit): warn honestly
                    st.belowTarget++;
                    std::cout << "[WARN] id=" << id << " ctrlCount=" << ctrl
                              << " target=" << target << " below-target\n";
                }

                lines.push_back(line);
            }
        }

        // negatives: at least 10, so the default 40 positives + 10 give 50 lines
        int negCount = perType / 2 + 2;
        if (negCount < 10)
            negCount = 10;
        for (int i = 0; i < negCount; ++i)
        {
            char idbuf[64];
            std::snprintf(idbuf, sizeof(idbuf), "%s_%03d", PrefixOf(Kind::Unknown).c_str(), negCounter++);
            const std::string id = idbuf;

            std::string line, sub, err;
            int ctrl = 0;
            if (!GenerateNegative(i, id, rng, outDir, seed, line, ctrl, sub, err))
            {
                std::cout << "[SKIP] id=" << id << " reason=" << err << "\n";
                continue;
            }

            Stat& st = stats[static_cast<int>(Kind::Unknown)];
            st.cases++;
            if (st.ctrlMin == 0 || ctrl < st.ctrlMin)
                st.ctrlMin = ctrl;
            if (ctrl > st.ctrlMax)
                st.ctrlMax = ctrl;

            lines.push_back(line);
        }

        const std::string manifestPath = outDir + "/manifest.jsonl";
        {
            std::ofstream ofs(manifestPath, std::ios::binary | std::ios::trunc);  // truncate; ASCII + newline
            if (!ofs)
            {
                std::cerr << "[FAIL] cannot open manifest for write: " << manifestPath << "\n";
                fini();  // release the kernel even on this early exit
                return 1;
            }
            for (const std::string& l : lines)
                ofs << l << "\n";  // one JSON object per line (JSONL)
        }

        // per type: control-point range, full/trimmed split, warnings only when non-zero
        for (int t = 0; t < 6; ++t)
        {
            const Stat& st = stats[t];
            std::cout << "[CORPUS] type=" << KindName(static_cast<Kind>(t))
                      << " cases=" << st.cases
                      << " ctrlRange=[" << st.ctrlMin << "," << st.ctrlMax << "]"
                      << " fullDomain=" << st.full
                      << " trimmed=" << st.trimmed;
            if (st.outOfRange > 0)
                std::cout << " outOfRange=" << st.outOfRange;
            if (st.belowTarget > 0)
                std::cout << " belowTarget=" << st.belowTarget;
            std::cout << "\n";
        }
        std::cout << "[CORPUS] total lines=" << lines.size()
                  << " manifest=" << manifestPath << "\n";

        fini();  // Foundation/init.h:28
    }
    catch (const std::exception& e)
    {
        std::cerr << "[FAIL] unhandled exception: " << e.what() << "\n";
        return 2;
    }
    catch (...)
    {
        std::cerr << "[FAIL] unhandled non-std exception\n";
        return 2;
    }
    return 0;
}
