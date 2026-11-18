# gta6-ccu-forecast

[![make.fun · gta6-twitch-ccu-forecast-01](https://make.fun/badge/gta6-twitch-ccu-forecast-01.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![rank](https://make.fun/badge/gta6-twitch-ccu-forecast-01/jhw1106/gta6-ccu-forecast/rank.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![reward](https://make.fun/badge/gta6-twitch-ccu-forecast-01/jhw1106/gta6-ccu-forecast/reward.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)

A forecaster for Grand Theft Auto VI's launch-day Twitch viewers, hour by
hour. It reads one case file (the first hour's print plus the case's
comparables: ~38 past big launches' hourly curves) and writes five
quantiles of each of hours 2..24.

## How it forecasts

The release-day curve has a stable shape, driven by the rolling
midnight-local unlock: hour 1 only has the regions that unlocked before
it (New Zealand, Australia, Asia, Europe from 23:00Z, the UK at 00:00Z);
hours 2-4 dip as Asia sleeps; hours 5-6 jump when the US East Coast
unlocks (05:00Z), hour 9 when the US West Coast unlocks (08:00Z); and the
evening (hours 17-22) peaks when the US and Europe are both awake.

How far the day ramps above hour 1 depends strongly on hour 1's own
level: across the comparables, log(peak / hour-1) correlates about -0.94
with log(hour-1 viewers). A launch whose hour 1 is already enormous
barely ramps; one with a small hour 1 can ramp 5-12x. The relation is
concave in log-log space -- steep at low hour-1 levels, flattening at
the top (the six comparables with hour-1 above 300k all peak at
0.94-1.34x hour 1) -- so each hour's log-ratio is fitted as a quadratic
in log(hour-1 viewers), not a line. A line extrapolated to GTA VI's mega
hour-1 would wrongly predict the evening below hour 1.

At runtime the program fits that quadratic per hour on the case's own
comparables (launches with a real release-day curve: hour 1 present, at
least 20 of 24 hours, hour-1 viewers >= 10,000 -- below that, hour 1 was
a pre-launch hour for a fixed-time global unlock), evaluates it at the
case's hour-1 level clamped to the comparables' range, and sets the
quantile spread from the fit's residual standard deviation.

## Results on the practice cases

Locally scored against the four practice launches' real SullyGnome
hourly figures (the kit ships no grader and no answers; the answers were
fetched once from SullyGnome's public chart API for validation only):

| case | this entry | house baseline |
| ---- | ---------: | -------------: |
| Marvel's Wolverine | 0.2699 | 0.7437 |
| Red Dead Redemption 2 | 0.2388 | 0.2868 |
| God of War Ragnarok | 0.1613 | 0.2146 |
| Marvel's Spider-Man 2 | 0.1950 | 0.3033 |
| **mean** | **0.2163** | **0.3871** |

Mean log pinball loss 0.2163 vs the house baseline's 0.3871 (typical miss
about 45% vs 83%). Runtime is about 0.08 s per case; two runs on the
same case write byte-identical output.

## Run it

```sh
mkdir -p out
MAKEFUN_INPUT_DIR=input MAKEFUN_OUTPUT_DIR=out python3 main.py
```

Python 3.13, standard library only, no network.
