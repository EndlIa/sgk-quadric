// Surface detection: curvature-signature classification (plane, sphere, torus, cylinder/cone).
// Samples the face on a 9x9 UV grid, marks and removes degenerate points (sphere poles, cone apex),
// then classifies. Lengths are mm, angles rad; residuals are per-branch invariant scalars.
// Classification order is fixed (plane -> sphere -> torus -> cylinder/cone) and must not be permuted.
// Scale normalization: diag = AABB diagonal of the non-degenerate samples (mm). Curvature limits are
// kRel/diag (units 1/mm); angle and normalized limits are compared directly, never scaled.
#include "sgkquadric/detect.hpp"
#include "sgkquadric/math.hpp"

#include <GeomBase/Dir3D.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace sgkq
{

namespace
{

// Criteria tolerances. Rel* limits are dimensionless and used as kRel/diag, which for the curvature
// limits has units 1/mm; kNormalAngle is an angle in rad compared directly; kAxisDotTol and
// kRelTorusRes are dimensionless. Dividing all curvature limits by diag makes them scale-invariant.
// kNormalAngle is in rad and is never scaled.
// Consumers: kFlat (curvature zero), kConst (curvature constant), the plane normal test, the c split,
// and the torus residual cap.
const double kRelCurvFlat  = 1e-3;   // curvature counts as zero when |kappa|*diag <= 1e-3
const double kRelCurvConst = 1e-2;   // curvature is 'constant' when the sampled deviation*diag <= 1e-2
const double kNormalAngle  = 1e-4;   // normal-angle tolerance (rad) for the plane test
const double kAxisDotTol   = 1e-3;   // |mean(n.a)| < 1e-3 means cylinder, else cone (the c split)
const double kRelTorusRes  = 1e-3;   // torus profile consistency: upper bound on res/diag

// Classification grid: 9x9 = 81 samples. A coarser 7x7 already resolves plane/sphere/cylinder/cone,
// but the torus axis (Plucker null space) and profile circle fit need enough spread along v, and 9
// samples cover a trimmed arc as short as a quarter turn. Finer grids stop buying accuracy for
// analytic surfaces. Degenerate points (poles/apex) are removed first but still leave far more than
// the 4 samples the signatures require, and a shrunk (clipped) domain stays covered. U and V are
const int kGridU = 9;
const int kGridV = 9;

//! @brief Clamp x into [lo, hi] (typically [-1, 1]).
//! Used before acos: floating point can push a normalized dot product to 1+eps, making acos return
//! NaN. NaN input propagates (both comparisons are false), so callers pass nonzero vectors.
//! @return the clamped value; requires lo <= hi.
double clampd(double x, double lo, double hi)
{
    return x < lo ? lo : (x > hi ? hi : x);
}

//! @brief AABB diagonal length of the sample cloud (model-scale proxy, mm; callers guarantee non-empty).
//! Serves as the global scale factor diag: relative tolerances are multiplied/divided by it, so
//! criteria are invariant to model scaling. An all-points-coincident set returns 0 and Detect falls
//! back to diag = 1.0 to avoid division by zero.
//! O(3N), read-only, no allocation.
double cloudDiagonal(const std::vector<Sample>& s)
{
    Point3D lo = s.front().p;
    Point3D hi = s.front().p;
    for (const Sample& x : s)
    {
        for (int k = 0; k < 3; ++k)
        {
            lo[k] = std::min(lo[k], x.p[k]);
            hi[k] = std::max(hi[k], x.p[k]);
        }
    }
    return (hi - lo).Length();
}

//! @brief Angle between two vectors in radians ([0, pi]); returns 0 if either is the zero vector.
//! The normalized dot product is clamped to [-1, 1] to guard acos against floating-point overflow.
//! Output is rad, the same dimension as kNormalAngle, so the plane test compares them directly.
//! Zero vectors are reported as angle 0, avoiding a 0/0 NaN; detect guarantees nonzero normals.
double angleOf(const Vector3D& a, const Vector3D& b)
{
    const double la = a.Length();
    const double lb = b.Length();
    if (!(la > 0.0) || !(lb > 0.0)) return 0.0;
    return std::acos(clampd(a.Dot(b) / (la * lb), -1.0, 1.0));
}

//! @brief Format a double for notes: scientific notation, 6 significant digits.
//! The fixed format keeps diffs and magnitude checks readable, and 6 digits distinguish pass/reject
//! near a threshold. Units are omitted; callers label them (rad, 1/mm) in the surrounding text.
std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(6) << std::scientific << v;
    return os.str();
}

//! True when torus-criteria evidence should be emitted (env SGKQ_DETECT_TRACE set; static-cached once).
//! The separation statistic (max accepted resNorm vs min rejected resNorm) is a cross-call aggregate,
//! but Detect is a pure function with no global mutable state, so each call prints one machine-
//! readable line that an external run aggregates. Off by default, so production/test output is clean.
//! This is harness diagnostics rather than library output: it writes to std::cout, so a run that sets
//! the variable must redirect stdout (data/torus_evidence_trace.txt is exactly such a captured run).
bool traceEnabled()
{
    static const bool on = (std::getenv("SGKQ_DETECT_TRACE") != nullptr);
    return on;
}

//! @brief Emit one line: [TORUS-EVIDENCE] status=... resNorm=... KsignChange=... R=... r=...
//! status tags: accept; reject-axis (the Plucker axis solve failed); reject-circlefit (the profile
//! (rho,h) circle fit degenerated: too few / collinear / singular); reject-radii (that fit returned
//! illegal major/minor radii); reject-residual (fit succeeded but the normalized residual exceeded
//! kRelTorusRes); reject-unknown-tail (the face never became a torus candidate). resNorm=-1 marks the
//! four tags where no residual exists (reject-residual is the only reject tag that carries one).
//! kSignChange (whether K = k1*k2 changes sign over the set) is informational only. This is
//! diagnostic output: with tracing off Detect behaves byte-for-byte as before.
void traceTorus(const char* status, double resNorm, bool kSignChange, double R, double r)
{
    if (!traceEnabled()) return;
    std::cout << "[TORUS-EVIDENCE] status=" << status
              << " resNorm=" << num(resNorm)
              << " KsignChange=" << (kSignChange ? "true" : "false")
              << " R=" << num(R) << " r=" << num(r) << "\n";
}

}

// Sampler: evaluate the surface on an nu x nv uniform grid (row-major: u outer, v inner), returning
// points with unit normals, principal curvatures and a degenerate flag.
// @param uv sample domain; @param nu, nv sample counts, each raised to at least 2.
// @return exactly nu*nv samples, degenerate ones included (only flagged) so the result stays
//         index-aligned with the grid for debugging; removal happens once, in Detect.
// The grid includes both domain ends: degenerate points (sphere poles, cone apex) and periodic seams
// usually sit on the boundary, so they must be sampled to be flagged and removed.
// Units: point coordinates in model units (mm); normals normalized (dimensionless); k1, k2 in 1/length.
// Degeneracy is flagged if CalcNormal marks it, IsDegeneratePoint holds, the normal is zero/non-finite,
// or a curvature is non-finite. The linear grid assumes a roughly uniform parameterization, which holds
// for the exact analytic-surface inputs used here.
// Ownership: returns a freshly built vector owned by the caller; no reference to srf/uv is kept.
// Concurrency: srf is read-only; sharing one srf across threads is governed by the SDK's contract.
std::vector<Sample> SampleSurface(const Surface& srf, const UVRange& uv, int nu, int nv)
{
    if (nu < 2) nu = 2;
    if (nv < 2) nv = 2;

    const double u0 = uv.IntervalU().Min();
    const double u1 = uv.IntervalU().Max();
    const double v0 = uv.IntervalV().Min();
    const double v1 = uv.IntervalV().Max();

    std::vector<Sample> out;
    out.reserve(static_cast<std::size_t>(nu) * static_cast<std::size_t>(nv));

    for (int i = 0; i < nu; ++i)
    {
        const double u = u0 + (u1 - u0) * (static_cast<double>(i) / static_cast<double>(nu - 1));
        for (int j = 0; j < nv; ++j)
        {
            const double v = v0 + (v1 - v0) * (static_cast<double>(j) / static_cast<double>(nv - 1));

            Sample s;
            s.uv = Point2D(u, v);

            // Surface.h: Dir3D CalcNormal(Real u, Real v, Point3D& thePoint, bool& isDegen);
            // one call yields point + unit normal + degeneracy flag.
            bool isDegen = false;
            const Dir3D nrm = srf.CalcNormal(u, v, s.p, isDegen);
            s.n = nrm;
            s.degen = isDegen;

            // Degeneracy guards (any single hit rejects the sample): CalcNormal failure flag, IsDegeneratePoint, zero/non-finite values
            if (!s.degen && srf.IsDegeneratePoint(s.uv)) s.degen = true;

            const double nl = s.n.Length();
            if (!(nl > 0.0) || !std::isfinite(nl))
            {
                s.degen = true;   // zero / invalid normal => degenerate, skip
            }
            else
            {
                s.n = s.n / nl;
            }

            if (!s.degen)
            {
                // Surface.h: void CalcCurvature(Real u, Real v, Real& k1, Real& k2, Dir3D& d1, Dir3D& d2)
                Dir3D d1, d2;
                srf.CalcCurvature(u, v, s.k1, s.k2, d1, d2);
                if (!std::isfinite(s.k1) || !std::isfinite(s.k2)) s.degen = true;
            }

            out.push_back(s);
        }
    }
    return out;
}

// Detection entry point: sample 9x9, then classify by curvature signature plus evidence checks.
// @param uv sample domain. A range fully outside the face is extrapolated by the SDK, not rejected
//        (see test_exceptions); the production path clips uv to the face, so this is unreachable.
// @param th accuracy-verification thresholds; detection does NOT use them (see the (void)th line).
// @return DetectResult{candidate, residual, notes}; residual is the per-branch invariant scalar and
//         notes is a readable evidence chain naming the deciding step (S1..S4).
// Priority (first match wins; the order must not be permuted):
//   0) sample, drop degenerates; fewer than 4 samples => Unknown.
//   1) plane: both curvatures ~ 0 and the normal constant.
//   2) sphere: k1 ~ k2 ~ k0 != 0, constant between samples.
//   3) torus: S1 exclusion screen + S2 axis + S3 profile + S4 residual.
//   4) cylinder/cone: one curvature is 0, the other nonzero; split by c = mean(n.a).
//   5) otherwise Unknown, reporting the max deviation from developable.
// th is discarded because detection uses only curvature signatures and closed-form residuals;
// Thresholds belongs to the accuracy layer and mixing the two would cross-contaminate them.
// Pure function, no static state; the only I/O is optional tracing (see traceTorus).
// Units: lengths mm; residual angles in rad, curvatures in 1/mm or dimensionless (per branch).
// Idempotent: the same (srf, uv) always yields the same result (no randomness, no cache).
DetectResult Detect(const Surface& srf, const UVRange& uv, const Thresholds& th)
{
    (void)th;

    DetectResult res;

    const std::vector<Sample> all = SampleSurface(srf, uv, kGridU, kGridV);
    std::vector<Sample> s;
    s.reserve(all.size());
    for (const Sample& x : all)
    {
        if (!x.degen) s.push_back(x);
    }

    if (s.size() < 4)
    {
        res.candidate = Canonical::Unknown;
        res.notes.push_back("reject: too few non-degenerate samples (<4), cannot form a curvature signature");
        return res;
    }

    double diag = cloudDiagonal(s);
    if (!(diag > 0.0) || !std::isfinite(diag)) diag = 1.0;

    const double kFlat  = kRelCurvFlat  / diag;
    const double kConst = kRelCurvConst / diag;
    const double n      = static_cast<double>(s.size());

    // 1) plane: both principal curvatures ~ 0 and the normal is constant.
    {
        Vector3D nmean;
        for (const Sample& x : s) nmean += x.n;
        nmean = nmean / n;
        const bool nok = nmean.Normalize();

        bool flat = true;
        double maxAng = 0.0;
        for (const Sample& x : s)
        {
            if (std::fabs(x.k1) >= kFlat || std::fabs(x.k2) >= kFlat) flat = false;
            if (nok) maxAng = std::max(maxAng, angleOf(x.n, nmean));
        }

    // Plane holds iff both curvatures are ~ 0, the mean normal is valid, and the max normal angle is
    // below kNormalAngle; otherwise fall through to the next branch.
        if (flat && nok && maxAng < kNormalAngle)
        {
            res.candidate = Canonical::Plane;
            res.residual  = maxAng;
            res.notes.push_back("plane: all non-degenerate samples have |κ1|,|κ2| < kFlat=" + num(kFlat) +
                                " and max normal deviation " + num(maxAng) + " rad < " + num(kNormalAngle) +
                                " (plane signature)");
            return res;
        }
    }

    // 2) sphere: k1 ~ k2 = k0 != 0, equal to each other and constant between samples.
    {
    // Sphere iff k1~k2 (dmax), each constant between samples (sp1/sp2), and mean |k| nonzero (excludes plane).
        double m1 = 0.0, m2 = 0.0;
        for (const Sample& x : s) { m1 += x.k1; m2 += x.k2; }
        m1 /= n;
        m2 /= n;

        double dmax = 0.0, sp1 = 0.0, sp2 = 0.0;
        for (const Sample& x : s)
        {
            dmax = std::max(dmax, std::fabs(x.k1 - x.k2));
            sp1  = std::max(sp1, std::fabs(x.k1 - m1));
            sp2  = std::max(sp2, std::fabs(x.k2 - m2));
        }
        const double mag = 0.5 * (std::fabs(m1) + std::fabs(m2));

    // residual = diag*(dmax+sp1+sp2) turns the 1/mm deviations into a dimensionless relative value.
        if (dmax < kConst && sp1 < kConst && sp2 < kConst && mag > kFlat)
        {
            res.candidate = Canonical::Sphere;
            res.residual  = diag * (dmax + sp1 + sp2);
            res.notes.push_back("sphere: max|κ1−κ2|=" + num(dmax) + "; max inter-sample deviation of κ1/κ2 " + num(sp1) +
                                "/" + num(sp2) + " both < kConst=" + num(kConst) + "; mean |κ|=" + num(mag) +
                                " > kFlat (sphere signature)");
            return res;
        }
    }

    // 3) torus: exclusion screen (S1) plus axis/profile/residual evidence (S2..S4).
    // S1 (screen, cheap, exclusion only): a face becomes a torus candidate unless it was already
    //    claimed by plane/sphere: every sample has at least one curvature significantly nonzero,
    //    BOTH curvatures have nonzero samples, and it is not sphere-like (k1~k2~const). No positive
    //    torus shape criterion is used here; in particular K = k1*k2 changing sign is NOT required.
    // S2 (axis): the normals meet the axis but the intersection point varies, so only the Plucker
    //    6x6 null space applies (not the concurrency form). solvePluckerAxis failure => Unknown.
    // S3 (profile): in the axis frame h_i=(p_i-C).a, rho_i=||(p_i-C)-h_i*a||; for a true torus the
    //    (rho_i,h_i) lie on the tube-section circle, so fitCircle2D gives center (cx,cy) and radius r.
    //    R = cx is the major radius and r the minor radius. Degenerate fit or non-positive R/r => Unknown.
    // S4 (residual): with H_i = h_i - cy (cy aligns the axial origin to the tube-section center plane),
    //    res = mean|sqrt((rho_i-R)^2 + H_i^2) - r| is the mean closed-form torus distance over the whole
    //    sample set; res/diag > kRelTorusRes => Unknown. S2..S4 must all pass to return Torus.
    //
    // Why K sign change is NOT necessary (it caused the old trimmed-arc false rejects):
    //    torus K = cos(v) / (r(R + r cos v)) changes sign only at v = pi/2, 3pi/2, so a trimmed arc
    //    can lie entirely on one side and never change sign while still being a true torus. K sign
    //    change is sufficient but NOT necessary: report it in notes (KsignChange) but never require it.
    //
    // The zero point of h: the torus distance uses h relative to the TUBE-section center plane, not
    //    to an arbitrary axis point. solvePluckerAxis returns C at the mean axial projection of the
    //    samples, so a full ring has cy ~ 0 and the formula reduces to the plain h form; a trimmed arc
    //    has cy != 0 and MUST use H_i = h_i - cy, or a true torus is rejected.
    //
    // Threshold kRelTorusRes = 1e-3 (normalized by diag, hence scale-free) sits between the machine-
    //    precision residual of an exact circle (~1e-15) and the residual of a two-dimensional (rho,h)
    //    scatter (~5e-2): ~1e12 margin for true tori, ~45x for negatives. Trimming is harmless: u-
    //    trimming repeats one section arc, v-trimming takes a sub-arc of the same circle, leaving the
    //    center and radius unchanged. Apple (R=r) and Lemon (R<r) only change the size relation and
    //    self-intersection, not the 'points lie on a circle' fact, so S4 verifies all three uniformly.
    {
        // S1 screen statistics; K sign change is informational only (does not affect the decision).
        bool notFlat = true;
        int  cnt1 = 0, cnt2 = 0;
        bool kPos = false, kNeg = false;
        double meanMinMag = 0.0;

        for (const Sample& x : s)
        {
            const double a = std::fabs(x.k1);
            const double b = std::fabs(x.k2);
            if (std::max(a, b) <= kFlat) notFlat = false;
            if (a > kFlat) ++cnt1;
            if (b > kFlat) ++cnt2;
            const double K = x.k1 * x.k2;
            if (K > 0.0) kPos = true; else if (K < 0.0) kNeg = true;
            meanMinMag += std::min(a, b);
        }
        meanMinMag /= n;

        // Whether K = k1*k2 changes sign: written to notes/trace only; it never decides acceptance.
        const bool kSignChange = kPos && kNeg;

        // Sphere signature (identical to the branch-2 test), used to exclude sphere-like faces from
        // the torus candidate set (they belong to Sphere, not Torus).
        double sm1 = 0.0, sm2 = 0.0;
        for (const Sample& x : s) { sm1 += x.k1; sm2 += x.k2; }
        sm1 /= n;
        sm2 /= n;
        double sd = 0.0, ss1 = 0.0, ss2 = 0.0;
        for (const Sample& x : s)
        {
            sd  = std::max(sd,  std::fabs(x.k1 - x.k2));
            ss1 = std::max(ss1, std::fabs(x.k1 - sm1));
            ss2 = std::max(ss2, std::fabs(x.k2 - sm2));
        }
        const bool sphereLike = (sd < kConst && ss1 < kConst && ss2 < kConst);

        // S1 screen: (a) notFlat, (b) cnt1>0 && cnt2>0 (both curvatures have nonzero samples),
        // (c) not sphere-like. K sign change is deliberately absent: it is sufficient but not
        // necessary, and only informational (a trimmed arc on one side never changes sign).
        if (notFlat && cnt1 > 0 && cnt2 > 0 && !sphereLike)
        {
            // S2 axis (Plucker 6x6 null space): the normals meet the axis with a varying intersection
            // point, so use the colinearity form, not the concurrency form.
            std::vector<Point3D>  plPts;
            std::vector<Vector3D> plDirs;
            plPts.reserve(s.size());
            plDirs.reserve(s.size());
            for (const Sample& x : s)
            {
                plPts.push_back(x.p);
                plDirs.push_back(x.n);
            }

            Vector3D axisDir;
            Point3D  axisPt;
            // Axis solve failed (degenerate null space / too few points) => Unknown; keep the old invariant value for reference.
            if (!solvePluckerAxis(plPts, plDirs, axisDir, axisPt))
            {
                res.candidate = Canonical::Unknown;
                res.residual  = diag * meanMinMag;
                res.notes.push_back("reject: passed the exclusion pre-screen (S1) into the torus candidate, but the axis solve failed"
                                    "(solvePluckerAxis); diagnostic mean min(|κ1|,|κ2|)·diag=" +
                                    num(res.residual) + ";KsignChange=" +
                                    (kSignChange ? "true" : "false"));
                traceTorus("reject-axis", -1.0, kSignChange, 0.0, 0.0);
                return res;
            }

            // S3 profile: fit a 2D circle to (rho_i, h_i); a true torus lands exactly on the tube-
            // section circle, while a free-form surface is a 2D scatter with a large residual.
            std::vector<std::array<double, 2>> rhoH;
            rhoH.reserve(s.size());
            for (const Sample& x : s)
            {
                const Vector3D d = x.p - axisPt;
                const double   h = d.Dot(axisDir);
                const Vector3D perp = d - axisDir * h;
                rhoH.push_back(std::array<double, 2>{ perp.Length(), h });
            }

            double cx = 0.0, cy = 0.0, rCircle = 0.0;
            // Circle fit degenerate (too few points / collinear / singular) => Unknown.
            if (!fitCircle2D(rhoH, cx, cy, rCircle))
            {
                res.candidate = Canonical::Unknown;
                res.residual  = diag * meanMinMag;
                res.notes.push_back(std::string("reject: passed the exclusion pre-screen (S1) into the torus candidate, but the profile (ρ,h) circle fit degenerated"
                                    "(fitCircle2D; too few points/collinear/singular); KsignChange=") +
                                    (kSignChange ? "true" : "false"));
                traceTorus("reject-circlefit", -1.0, kSignChange, 0.0, 0.0);
                return res;
            }

            const double R = cx;
            const double r = rCircle;
            if (!(R > 0.0) || !(r > 0.0) || !std::isfinite(R) || !std::isfinite(r) ||
                !std::isfinite(cy))
            {
                res.candidate = Canonical::Unknown;
                res.residual  = diag * meanMinMag;
                res.notes.push_back("reject: profile circle fit gave illegal major/minor radii R=" + num(R) + " r=" + num(r) +
                                    "(require R,r positive and finite)");
                traceTorus("reject-radii", -1.0, kSignChange, R, r);
                return res;
            }

            // S4 residual = mean over the whole set of the closed-form torus distance
            // |sqrt((rho-R)^2 + H^2) - r|, with H = h - cy aligned to the tube-section center plane.
            double rsum = 0.0;
            for (const Sample& x : s)
            {
                const Vector3D d     = x.p - axisPt;
                const double   axial = d.Dot(axisDir);
                const double   H     = axial - cy;
                const Vector3D perp  = d - axisDir * axial;
                const double   rho   = perp.Length();
                const double   dRho  = rho - R;
                rsum += std::fabs(std::sqrt(dRho * dRho + H * H) - r);
            }
            const double resNorm = (rsum / n) / diag;

            // Normalized residual over the limit => insufficient profile consistency => Unknown.
            if (!(resNorm <= kRelTorusRes))
            {
                res.candidate = Canonical::Unknown;
                res.residual  = resNorm;
                res.notes.push_back("reject: profile consistency check failed (S4): normalized residual res/diag=" + num(resNorm) +
                                    " > kRelTorusRes=" + num(kRelTorusRes) + "; fitted R=" + num(R) +
                                    " r=" + num(r) + ";KsignChange=" + (kSignChange ? "true" : "false") +
                                    " (torus profile consistency)");
                traceTorus("reject-residual", resNorm, kSignChange, R, r);
                return res;
            }

            // S2 + S3 + S4 all pass => Torus. S1 only excludes, and K sign change does not affect this.
            res.candidate = Canonical::Torus;
            res.residual  = resNorm;
            res.notes.push_back(std::string("torus: exclusion pre-screen passed (both principal curvatures nonzero and non-spherical), axis=Plücker "
                                            "6×6 null space gives R=") + num(R) + ", r=" + num(r) +
                                "; profile consistency residual res/diag=" + num(resNorm) + " ≤ kRelTorusRes=" +
                                num(kRelTorusRes) + "(S2→S3→S4 evidence verified); KsignChange=" +
                                (kSignChange ? "true" : "false") +
                                "(informational, not used in the decision; false means the whole trimmed patch lies on one side of the K sign)"
                                " (torus signature)");
            traceTorus("accept", resNorm, kSignChange, R, r);
            return res;
        }
    }

    // 4) cylinder/cone: one principal curvature is 0; split by c = mean(n.a).
    {
        // Developable signature: min(|k1|, |k2|) is below kFlat everywhere.
        // residual = sqrt(mean((n_i.a - c)^2)) measures the normal spread about a single cone/cylinder angle.
        bool oneZero = true;
        for (const Sample& x : s)
        {
            if (std::min(std::fabs(x.k1), std::fabs(x.k2)) >= kFlat) { oneZero = false; break; }
        }

        if (oneZero)
        {
            std::vector<Vector3D> nn;
            nn.reserve(s.size());
            for (const Sample& x : s) nn.push_back(x.n);

            // Axis = smallest eigenvector of Cov(n) (mean-subtracted). NEVER the largest eigenvector
            // of the second moment: beyond 54.7356 deg it is perpendicular to the true axis.
            const Vector3D a = minVarDirection(nn);

            // a.Length() > 0.5: a true unit vector has length 1; minVarDirection returns the zero
            // vector for degenerate input (the loose bound avoids dithering between valid/invalid).
            if (a.Length() > 0.5)
            {
                double c = 0.0;
                for (const Sample& x : s) c += x.n.Dot(a);
                c /= n;

                double var = 0.0;
                for (const Sample& x : s)
                {
                    const double d = x.n.Dot(a) - c;
                    var += d * d;
                }
                var /= n;

                res.residual = std::sqrt(var);
                // c split: cylinder normals are perpendicular to the axis => c ~ 0; cone normals make angle (90 - alpha) => c = sin(alpha) != 0.
                if (std::fabs(c) < kAxisDotTol)
                {
                    res.candidate = Canonical::Cylinder;
                    res.notes.push_back("cylinder: min(|κ1|,|κ2|) ≡ 0 and mean(n·a)=" + num(c) +
                                        " ≈ 0 (c branch; axis = smallest covariance eigenvector)");
                }
                else
                {
                    res.candidate = Canonical::Cone;
                    res.notes.push_back("cone: min(|κ1|,|κ2|) ≡ 0 and mean(n·a)=" + num(c) +
                                        " = sin α ≠ 0 (c branch; axis = smallest eigenvector, never the second moment)");
                }
                return res;
            }
        }
    }

    // No signature matched: reject, reporting the max normalized deviation from developable.
    res.candidate = Canonical::Unknown;
    double maxMinMag = 0.0;
    for (const Sample& x : s)
    {
        maxMinMag = std::max(maxMinMag, std::min(std::fabs(x.k1), std::fabs(x.k2)));
    }
    res.residual = diag * maxMinMag;
    res.notes.push_back("reject: passed no signature pre-screen (neither a plane/sphere/torus candidate, nor satisfying the cylinder/cone condition 'some principal curvature ≡ 0')"
                        " and matched no evidence check; max min(|κ1|,|κ2|)·diag=" + num(res.residual) +
                        " (no signature match)");
    // No torus branch was entered, so resNorm = -1 (no profile evidence) and KsignChange = false;
    // the status tag distinguishes 'never a candidate' from 'candidate but insufficient evidence'.
    traceTorus("reject-unknown-tail", -1.0, false, 0.0, 0.0);
    return res;
}

}

