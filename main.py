"""gta6-twitch-ccu-forecast-01 entry: release-day Twitch viewer forecaster.

Model (authored for this competition; the kit's starter is used only as a
format reference and as a last-resort fallback):

  For each hour h = 2..24 of release day, forecast five quantiles of that
  hour's concurrent viewers in the game's Twitch category, as a ratio to
  hour 1's print.

  The winning shape is in the case file's own comparables: ~38 past big
  launches' hourly curves. Two facts drive the model:

  1. The release-day curve has a stable shape. Hour 1 (00:00-01:00 UTC)
     only has the regions that unlocked before it (NZ/Australia/Asia,
     Europe from 23:00Z, the UK at 00:00Z). Hours 2-4 dip as Asia sleeps,
     hours 5-6 jump when the US East Coast unlocks (05:00Z), hour 9 when
     the US West Coast unlocks (08:00Z), and the evening (hours 17-22)
     peaks when the US and Europe are both awake and playing.

  2. How far the day ramps above hour 1 depends strongly on hour 1's own
     level. Across the comparables, log(peak/hour-1) correlates about
     -0.94 with log(hour-1 viewers): a launch whose hour 1 is already
     enormous (a huge Asian audience, early access, a simultaneous PC
     unlock) barely ramps, while a launch with a small hour 1 can ramp
     5-12x as the US unlocks. The relation is concave in log-log space --
     steep at low hour-1 levels, flattening at the top (the six
     comparables with hour-1 above 300k all peak at 0.94-1.34x hour 1) --
     so each hour's log-ratio is fitted as a quadratic in log(hour-1
     viewers), not a line. A line extrapolated to GTA VI's mega hour-1
     would wrongly predict the evening below hour 1.

  Concretely, at runtime the program reads the case's comparables,
  keeps the launches with a real release-day curve (hour 1 present,
  at least 20 of the 24 hours present, hour-1 viewers >= 10,000 --
  below that, hour 1 was a pre-launch hour for a fixed-time global
  unlock, and its ratios are meaningless), and for each hour fits
  log(v_h / v_1) as a quadratic in log(v_1) by least squares. The
  forecast median is the fitted value at the case's hour-1 level,
  clamped to the comparables' hour-1 range so nothing extrapolates
  past the data; the quantile spread is the fit's residual standard
  deviation (floored at 0.15 log units), applied with the standard
  normal's quantiles. Quantiles are q_tau = p1 * exp(mu + z_tau * sd).

  The median is not the quadratic alone. A quadratic and a line are
  both reasonable models of the same log-log relation, and averaging
  them hedges the risk that either one's functional form is wrong at
  the case's hour-1 level. But the average is tapered: the line gets
  weight only in the interior of the comparables' hour-1 range, where
  the data is dense and both fits are well supported, and the weight
  fades to zero at both edges. At the edges the quadratic's curvature
  is the structural feature that matters -- at the top end it is what
  keeps the evening at or above hour 1 (the six comparables with
  hour-1 above 300k all peak at 0.94-1.34x hour 1, while a line
  extrapolated there wrongly predicts the evening below hour 1) -- so
  at the edges the forecast is the pure quadratic. The taper is
  w_line = 0.5 * 4 * t * (1 - t), where t is the case's hour-1 level
  as a fraction of the comparables' log hour-1 range: 0 at both ends,
  0.5 (an equal average) at the centre.

  Fallbacks, in order, if the comparables cannot support a fit: the
  empirical median-ratio curve with empirical quantile spread; then the
  house persistence baseline (the starter's formula).

Standard library only, no network, no randomness: two runs on the same
case write byte-identical output.
"""

import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

INPUT_DIR = Path(os.environ.get("MAKEFUN_INPUT_DIR", "/input"))
OUTPUT_DIR = Path(os.environ.get("MAKEFUN_OUTPUT_DIR", "/output"))

LEVELS = [0.05, 0.25, 0.5, 0.75, 0.95]
Z = [-1.645, -0.674, 0, 0.674, 1.645]

# A comparable only counts when its hour-1 print is at least this large.
# Below it, hour 1 was a pre-launch hour (fixed global unlock later in
# the day) and the launch's ratios to it are meaningless.
MIN_H1 = 10000
# A comparable needs at least this many of the 24 release-day hours.
MIN_DAY_HOURS = 20
# The fitted log-ratio is clamped to this band (a ratio of 0.30 to 15).
LOG_RATIO_MIN = -1.2
LOG_RATIO_MAX = 2.7
# The quantile spread is never narrower than this, in log units.
SD_FLOOR = 0.15
# The starter's spread, used only by the persistence fallback.
BASELINE_SIGMA = 0.15


def parse_time(text):
    """Parse an ISO 8601 timestamp (with a Z suffix) to a datetime."""
    return datetime.fromisoformat(text.replace("Z", "+00:00"))


def round2(x):
    """Round to 2 decimals, half up (the starter's convention)."""
    scaled = x * 100
    whole = math.floor(scaled)
    if scaled - whole >= 0.5:
        whole += 1
    value = whole / 100
    return int(value) if value == int(value) else value


def solve3(a, b):
    """Solve a 3x3 linear system by Gaussian elimination with pivoting."""
    n = 3
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-12:
            return None
        m[col], m[piv] = m[piv], m[col]
        for r in range(col + 1, n):
            f = m[r][col] / m[col][col]
            for cc in range(col, n + 1):
                m[r][cc] -= f * m[col][cc]
    x = [0.0] * n
    for i in (2, 1, 0):
        x[i] = (m[i][n] - sum(m[i][j] * x[j] for j in range(i + 1, n))) / m[i][i]
    return x


def fit_quadratic(xs, ys):
    """Least-squares fit of y = c0 + c1*(x-mx) + c2*(x-mx)^2.

    Returns (predict, residual_sd) or None when the fit is degenerate.
    """
    n = len(xs)
    if n < 6:
        return None
    mx = sum(xs) / n
    normal = [[0.0] * 3 for _ in range(3)]
    rhs = [0.0] * 3
    for xi, yi in zip(xs, ys):
        d = xi - mx
        row = [1.0, d, d * d]
        for i in range(3):
            rhs[i] += row[i] * yi
            for j in range(3):
                normal[i][j] += row[i] * row[j]
    coef = solve3(normal, rhs)
    if coef is None:
        return None

    def predict(x):
        d = x - mx
        return coef[0] + coef[1] * d + coef[2] * d * d

    ss = sum((yi - predict(xi)) ** 2 for xi, yi in zip(xs, ys))
    sd = math.sqrt(ss / max(n - 3, 1))
    return predict, sd


def fit_line(xs, ys):
    """Least-squares fit of y = a + b*x. Returns (predict, residual_sd)."""
    n = len(xs)
    if n < 4:
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx <= 1e-12:
        return None
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    a = my - b * mx
    ss = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
    sd = math.sqrt(ss / max(n - 2, 1))
    return (lambda x: a + b * x), sd


def comparable_curves(case):
    """Release-day hourly viewers for each comparable with a real curve.

    Returns a list of (hour_1_viewers, {hour: viewers}) with hours 1..24
    relative to the comparable's day_start.
    """
    curves = []
    for comp in case.get("comparables", []):
        try:
            day_start = parse_time(comp["day_start"])
        except (KeyError, TypeError, ValueError):
            continue
        hourly = {}
        for point in comp.get("hourly", []):
            try:
                t = parse_time(point["from"])
                viewers = point["viewers"]
            except (KeyError, TypeError, ValueError):
                continue
            if viewers is None or viewers <= 0:
                continue
            delta = (t - day_start).total_seconds() / 3600.0
            if 0 <= delta < 24:
                hourly[int(delta) + 1] = viewers
        if 1 in hourly and len(hourly) >= MIN_DAY_HOURS and hourly[1] >= MIN_H1:
            curves.append((hourly[1], hourly))
    return curves


def empirical_quantile(sorted_values, p):
    """Linear-interpolation quantile of an already sorted list."""
    n = len(sorted_values)
    if n == 1:
        return sorted_values[0]
    pos = p * (n - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return sorted_values[lo]
    return sorted_values[lo] * (hi - pos) + sorted_values[hi] * (pos - lo)


def forecast_median_curve(p1, curves, hours):
    """Fallback: empirical median-ratio curve with empirical spread."""
    rows = []
    for h in hours:
        log_ratios = sorted(
            math.log(curve[h] / h1) for h1, curve in curves if h in curve
        )
        if not log_ratios:
            log_ratios = [0.0]
        qs = [empirical_quantile(log_ratios, tau) for tau in LEVELS]
        rows.append(
            {"hour": h, "q": [max(1.0, round2(p1 * math.exp(x))) for x in qs]}
        )
    return {"hours": rows}


def forecast_persistence(p1, hours):
    """Last-resort fallback: the house baseline's persistence forecast."""
    rows = []
    for h in hours:
        qs = [
            round2(p1 * math.exp(z * BASELINE_SIGMA * math.sqrt(h - 1)))
            for z in Z
        ]
        rows.append({"hour": h, "q": [max(1.0, q) for q in qs]})
    return {"hours": rows}


def forecast_case(case, name):
    if case["levels"] != LEVELS:
        raise ValueError(f"{name} asks for levels {case['levels']}, not {LEVELS}")
    first = next((k for k in case["known"] if k["hour"] == 1), None)
    if first is None:
        raise ValueError(f"{name} has no hour-1 print")
    p1 = max(1, first["viewers"])
    hours = [row["hour"] for row in case["forecast"]]

    curves = comparable_curves(case)
    if len(curves) >= 8:
        log_h1s = [math.log(h1) for h1, _ in curves]
        lo_h1, hi_h1 = min(log_h1s), max(log_h1s)
        # The case's hour-1 level, clamped to the comparables' range so
        # the fit never extrapolates past the data it was fitted on.
        lh1 = min(max(math.log(p1), lo_h1), hi_h1)
        # Position of the case's hour-1 level inside the comparables'
        # range, as a fraction: 0 at the bottom, 1 at the top. The
        # linear fit's weight in the median average is tapered by
        # 4*t*(1-t): full (equal) weight at the centre, where the data
        # is dense, fading to zero at both edges, where the quadratic's
        # curvature is the structural feature and the line misbehaves.
        if hi_h1 > lo_h1:
            t = (lh1 - lo_h1) / (hi_h1 - lo_h1)
            w_line = 0.5 * 4.0 * t * (1.0 - t)
        else:
            w_line = 0.0
        rows = []
        ok = True
        for h in hours:
            pts = [
                (math.log(h1), math.log(curve[h] / h1))
                for h1, curve in curves
                if h in curve
            ]
            fitted = fit_quadratic([p[0] for p in pts], [p[1] for p in pts])
            linear = fit_line([p[0] for p in pts], [p[1] for p in pts])
            if fitted is None:
                fitted = linear
                linear = None
            if fitted is None:
                ok = False
                break
            predict, sd = fitted
            mu = predict(lh1)
            if linear is not None and w_line > 0.0:
                lmu = linear[0](lh1)
                mu = (1.0 - w_line) * mu + w_line * lmu
            mu = min(max(mu, LOG_RATIO_MIN), LOG_RATIO_MAX)
            sd = max(sd, SD_FLOOR)
            qs = [mu + z * sd for z in Z]
            rows.append(
                {"hour": h, "q": [max(1.0, round2(p1 * math.exp(x))) for x in qs]}
            )
        if ok:
            return {"hours": rows}

    if curves:
        return forecast_median_curve(p1, curves, hours)
    return forecast_persistence(p1, hours)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in sorted(INPUT_DIR.glob("*.in")):
        case = json.loads(path.read_text(encoding="utf-8"))
        forecast = forecast_case(case, path.name)
        out = OUTPUT_DIR / f"{path.stem}.out"
        out.write_text(
            json.dumps(forecast, separators=(",", ":")) + "\n", encoding="utf-8"
        )


if __name__ == "__main__":
    main()
