# Talking points: call with Gareth Peters, 4 November 2026

These go with slides-2026-11-04 (32 slides). Every number comes from the slides, the paper or the saved experiment results.

## If I only get 30 seconds

We test whether two time series are independent. A test like HSIC makes a promise. If the series really are unrelated, it will wrongly report a link only 5% of the time. That 5% is its false-alarm rate.

The promise rests on proofs that assume the past fades quickly, so that what happened long ago says almost nothing about today. Many real series don't behave that way. They have long memory, which means slow swings that last a long time. We measured what happens to the false-alarm rate on data like that.

We found three things. First, the bigger risk is how the test sets its cut-off, not which test you pick. Second, seasonal patterns cause more false alarms than memory does. Third, kernel tests have a weak spot that ordinary correlation doesn't. Two unrelated series that both become more volatile at the same point in the season look linked to HSIC. Pearson correlation correctly sees nothing.

Our advice is to remove each series' seasonal pattern, both its average level and its spread. Then run HSIC with circular shifts. In our simulations at n = 2000, that kept the false-alarm rate between 3.8% and 7.0%.

## Three things to remember

1. **The cut-off rule matters as much as the test.** Every test compares its statistic with a cut-off, called the critical value. The cut-off has to be estimated from the same data, and there are several ways to do it. Subsampling and wild multipliers learn it from short stretches of the series. A circular shift slides one series past the other and uses the whole path. Take one setting, H = 0.8 and n = 2000. The same HSIC statistic on the same data gives 19.7% false alarms with subsampling, 8.9% with wild multipliers and 5.7% with circular shifts. Short stretches miss the slow swings, so their cut-off comes out too low. Choosing HSIC is only half the decision. You also have to choose the rule.

2. **Kernels see more than correlation does, and that can mislead them.** Pearson correlation only asks whether two series rise and fall together. HSIC asks about any kind of link, including whether they get volatile at the same times. Take two independent series with no memory, whose swings both grow at the same point in each seasonal cycle. HSIC with wild multipliers reports a link 73.3% of the time. Pearson reports one 4.7% of the time, close to the promised 5%. The kernel is right that the swings line up. They line up because both series follow the same calendar, and neither affects the other. That extra sensitivity is the reason to use a kernel. On data ordered in time, it is also how kernels get fooled.

3. **Kernels built to remove memory help, but they don't solve the problem.** The idea was simple. If slow swings cause the false alarms, build a kernel that can't see slow swings. EXP7 tested that idea. It fixed one problem and left three others.

**What we tried**

- The increment kernel (inc1) looks at changes, comparing each value with the one before it. A slow swing barely shows up in day-to-day changes, so most of the memory disappears.
- The Hermite kernels first turn each value into a normal score, then strip out the simplest patterns in those scores. Memory sits mostly in those patterns, so removing them removes most of the memory.
- The Hermite split goes one step further. It still tests the patterns it stripped out, but separately and with a test built for memory, so it doesn't throw that evidence away.

**What got better: the test we ran**

- Setup: two independent fractional Gaussian noise series with the same memory H, each 2,000 points long. We ran 1,000 pairs at each H from 0.5 (no memory) to 0.95 (very strong memory). The series are independent by construction, so every rejection is a false alarm.
- Each pair was tested twice, on the same data with the same random draws. One test used the ordinary Gaussian kernel and the other used the increment kernel (inc1). Any difference comes from the kernel alone.
- For each kernel we tried three ways of setting the cut-off: wild multipliers, subsampling and circular shifts. We counted how often each one rejected at the 5% level.

**What got better: the numbers**

- Wild multipliers. The ordinary kernel goes from 6.4% at H = 0.5 to 19.1% at H = 0.95, almost four times the promised rate. inc1 goes from 7.5% to 6.7%, so memory no longer pushes it up.
- Subsampling. The ordinary kernel goes from 5.2% to 62.4%, so at very strong memory it raises a false alarm on most pairs. inc1 goes from 4.7% to 2.7%.
- Circular shifts. The ordinary kernel goes from 4.5% to 8.1%, so even shifts start to slip at H = 0.95. inc1 stays at 5.0% and 5.2%.
- At H = 0.95 the intervals are tight: 16.8% to 21.7% for the ordinary kernel under wild, and 59.4% to 65.3% under subsampling. With 1,000 pairs, these gaps are not noise.

**What got better: what it means**

- Memory caused the climb. With the same data and the same rules, only the kernel changed, and the climb disappeared. Looking at changes instead of levels removes the slow swings that fooled the cut-off.
- The drop from 19.1% to 6.7% fixes the memory problem, not the wild rule. inc1 already sits at 7.5% with no memory at all. The leftover 2 points above 5% come from the wild rule itself, at every H.
- Subsampling overshoots in the other direction. At 2.7%, inc1 with subsampling is now too cautious, so it may miss real links.
- inc1 with circular shifts is the combination that stays at 5% at both ends, 5.0% with no memory and 5.2% at very strong memory.
- These numbers are for Gaussian data and the HSIC test. On heavy-tailed data and for split MMD the picture is worse, as the next lists show.

**Why the ordinary kernel climbs with memory**

- In EXP7, on 24 independent Gaussian settings, raw HSIC with wild multipliers rose from 7.2% false alarms at H of 0.6 or below to 15.6% at H of 0.8 or above. Subsampling rose from 6.4% to 42.8%. Shifts rose only from 4.7% to 7.1%, and the Harris test did not rise at all.
- First reason: two unrelated series with long memory can drift in step by chance for long stretches, and stronger memory makes those stretches longer. The HSIC statistic picks up these chance alignments and keeps growing with n. At strong memory it grew 4.3 times from n = 500 to n = 8000.
- Second reason: wild multipliers and subsampling estimate the cut-off from short stretches, about 44 points at n = 2000. Those stretches are too short to contain the slow swings, so the cut-off stays low while the statistic grows. The gap widens as memory gets stronger. Shifts and Harris use the whole path, which is why they barely move.

**What did not get better**

- False alarms still sit above 5%. Even with no memory at all, the memory kernels reject 6.4% to 8.3% under wild multipliers. That excess comes from the cut-off rule itself, so changing the kernel can't remove it.
- None of the kernels met the target we wrote down before the run, which was to stay close to 5% under wild multipliers or subsampling.
- The split-sample MMD test still climbs with memory. On heavy-tailed data at H = 0.9, four of the memory kernels reject between 8.0% and 17.4%. One version of the Hermite split (hsplit2 with block cut-offs) did hold near 5% in simulation, but nothing proves it.

**What each kernel can no longer see**

Removing memory also removes real signal that looks like memory.
- inc1 misses two series that move together slowly. It detects a slow common factor 8.6% of the time, where the ordinary kernel detects it every time. It also misses a flat moving average, at 4.3% against 71.1%.
- The first Hermite kernel misses plain linear dependence. It detects it 6.0% of the time, against 58.0% for the ordinary kernel.

**A simpler tool did better on false alarms**

- The Harris (2020) shift test was published before we started. It keeps the ordinary kernel and changes only the cut-off rule.
- It held false alarms between 1.6% and 3.8% in all 44 stationary cases we ran. No memory kernel beat it on every type of real link.
- It has a cost. With no memory it detects less than circular shifts do, trailing by more than 5 points on five of six alternatives at H = 0.5.

**What this means in practice (a proposal, not yet our procedure)**

- What we tested: inc1 combined with circular shifts held its size in simulation. It gave 5.0% at H = 0.5, 5.2% at H = 0.95 and 5.8% under slow drift. On the real-data surrogate null in EXP8 it gave 4.0% to 6.2% across the 16 cells. inc1 with the Harris test gave 1.4% to 3.1% there.
- What we have not done: the paper's recommended procedure still uses the ordinary kernel with shifts, and inc1 is not part of it. Whether to add it is Ask 3 for Gareth. No proof covers inc1 with shifts.
- The trade-off: adding inc1 buys protection against memory and slow drift, at the cost of missing slow common factors and flat moving averages. Anyone who uses it should say so.

## Slide 1: Title

- Skip it. One sentence on the aim: where kernel tests break on long memory, why, and what to do instead.

## Slide 2: Research problem and questions

- We test whether two series are independent, using one observed path. Either series can have long memory. Its correlations fade so slowly that they never add up to a finite total.
- H, the Hurst index, measures memory. At 0.5 there's none. Close to 1 it's strong. Our baseline process is fractional Gaussian noise (fGn).
- The rule estimates the critical value from that same path. The critical value is the cut-off the statistic has to beat. If the null is true, a calibrated rule rejects 5% of the time.
- Four questions. Does the theory cover these data? How wrong is each rule? What causes the error? What should someone actually do?

## Slide 3: Why the 5% guarantee is in doubt

The chain on this slide is not ours. It is the proof that the published kernel tests use to justify their 5% false-alarm rate. The two papers our HSIC rules come from are the main ones:

- Chwialkowski and Gretton (ICML 2014), "A Kernel Independence Test for Random Processes". This is the HSIC test with circular shifts.
- Chwialkowski, Sejdinovic and Gretton (NeurIPS 2014), "A Wild Bootstrap for Degenerate Kernel Tests". This is the wild-multiplier rule for HSIC and MMD.

Newer papers follow the same pattern. Betken, Dehling and Kroll (2024) do this for a distance-covariance bootstrap. Ghoshal and Shao (2026) do it for self-normalised RKHS tests. Diz-Castro, Febrero-Bande and González-Manteiga (2026) do it for an HSIC wild bootstrap.

Each of these proofs runs through the same four steps:

1. **Dependence fades fast.** The data must be "mixing" at a fast enough rate. Mixing means that values far apart in time are close to unrelated, and the rate says how quickly that happens. The papers use different versions (alpha-, beta- or tau-mixing), but all need the dependence to die out quickly enough to add up to a finite total.
2. **The statistic settles down.** Given step 1, the HSIC statistic (times n) has a known limiting distribution when the series are independent. This is the "null limit".
3. **The rule finds the cut-off.** Given step 2, the rule (wild multipliers, shifts, subsampling) estimates the 95th percentile of that limit from the data.
4. **The 5% holds.** Given step 3, the test rejects a true null 5% of the time as n grows.

Where our data break it:

- Fractional Gaussian noise with H above 1/2 fails step 1. Its correlations fade like a power of the lag, too slowly to add up. Our own check shows the wild-bootstrap paper's rate conditions fail for every H above 1/2.
- Steps 2 to 4 are therefore not shown to be wrong. Their proofs simply start from step 1, so for these data nothing supports them.
- No paper gives the null limit of Gaussian-kernel HSIC for independent fGn. Without it, nobody can work out the false-alarm rate on paper. That is why we measure it by simulation.
- For a seasonal mean with a random phase, step 1 also fails, and there we could go further. Proposition 1 proves that the false-alarm rate goes to 100% for Pearson, wild multipliers and subsampling.
- None of these papers made an error. Their assumptions simply rule out long-memory data, mostly without saying so. The source audit on slide 11 goes through each condition.

## Slide 4: Three kinds of statistic, each needing a critical value

- Split MMD compares the two halves of one series. HSIC and distance covariance compare two series lined up in time. Pearson is the linear benchmark.
- The dependence problem lands on the critical-value rule, so we treat the rule as a separate factor.
- An oracle that knows the true null shows how much each rule costs.

## Slide 5: Long memory reduces the information in a path

- The gap between the two half-sample means has variance (4 - 2^(2H)) m^(2H-2). That's exact, not an approximation (Proposition 2).
- At H = 0.5 the gap shrinks like 1/m. At H = 0.8 it shrinks only like m^-0.4.
- So 2000 points of a long-memory series tell you much less than 2000 independent points would. The power cap on slide 19 comes straight from this.

## Slide 6: Three critical-value rules

- Subsampling uses windows of about the square root of n points, which is 44 at n = 2000. The wild rule weights lags on the same scale.
- Long memory moves slowly, so a 44-point window misses most of the swing. The critical value comes out too low, and the test rejects too often.
- That gives a prediction we can check. Rules that look at short stretches should fail. A rule that uses the whole path should hold up better.

## Slide 7: The circular shift

- Wrap Y around by some offset and recompute the statistic. Y keeps its own memory. The only thing that moves is how it lines up with X.
- Doing this over many offsets gives a reference distribution built from the full series. We leave out offsets below 0.1n.
- The original result (Chwialkowski and Gretton, 2014) assumes mixing, and fGn isn't mixing. So shifts don't have a proof here either.

## Slide 8: How we test each question

- Theory: read the conditions in the primary sources and check fGn against each one.
- Calibration and attribution: Monte Carlo with a known null. We switch on one process feature at a time and cross it with H.
- Procedure: compare power with every rule tuned to the same true size. A test can hold its size just by never rejecting, and this guards against that.
- The stock durations are only an outside check. We can't verify their null, so their rates are descriptive.

## Slide 9: Monte Carlo design

- We choose the generator, so the null is true by construction. Every rejection is a false alarm.
- With 1,000 replications, a 5% rate has a 95% Wilson interval of roughly ±1.4 points. That covers counting error and nothing else.
- For EXP2 we wrote down the hypotheses and generator checks before the main run.

## Slide 10: The experiments

- EXP1 maps calibration across H and n. EXP2 crosses seven generators with H.
- EXP3 measures power with each rule tuned to a true 5%. EXP4 tests seasonal adjustment.
- EXP5 and EXP6 check larger n and other bandwidths. The stock data are AAPL, MSFT and NVDA trade durations from January to June 2024.

## Slide 11: Q1, the published conditions exclude fGn with H above 1/2

- Every source we audited assumes either i.i.d. data or mixing at a summable rate. fGn with H above 1/2 meets neither.
- Mostly the exclusion is implied rather than stated. Our own derivation puts tau(r) at order r^(2H-2) or larger, which breaks the wild bootstrap's rate condition.
- Betken and Dehling partly cover distance covariance. No source gives a Gaussian-kernel HSIC limit, and none deals with seasonality.
- The theory tells you to worry but not how much. From the papers alone you can't tell whether the error is 6% or 60%.

## Slide 12: Q2, the rule changes the size

- Same generator, same n. Only the rule changes, so the gaps between the lines come from the rule.
- Some of the excess has nothing to do with memory. The wild rule is already at 6.5% at H = 0.5.
- "HSIC rejected" isn't a complete result. You have to say which rule.

## Slide 13: Q2, higher Hermite rank delays the distortion

- Expand a feature in Hermite polynomials. Its Hermite rank is the order of the first term that isn't zero. Rank m gives covariance of order rho(k)^m.
- Rank one is affected at every H above 1/2. Rank two balances at H = 3/4. Products of independent rank-one terms balance at H_X + H_Y = 3/2.
- This explains the shape of the curves. It doesn't give a critical value, and it doesn't pin down a proved boundary.

## Slide 14: Q3, seasonality distorts size more than memory

- We switch on one feature at a time and cross it with memory. Heavy tails barely matter.
- HSIC with shifts stays between 3.8% and 7.0% across these cells.
- The threshold AR also adds short-range persistence, so its excess isn't pure nonlinearity. The seasonal settings are illustrations. We didn't calibrate them to markets.
- If you're worried about long memory, you may be worrying about the wrong thing. In our settings seasonality is the bigger risk.

## Slide 15: Q3, a seasonal scale distorts kernel tests unless they use shifts

- The two series are independent, each has its own random phase, and there's no memory.
- HSIC wild, HSIC subsampling and dCov subsampling reject between 46% and 73%. HSIC with shifts rejects 5.0%, and Pearson 4.7%.
- We measured this. We haven't proved it.
- Only kernel tests fall for this one. Pearson doesn't react, so the usual linear checks won't warn you.

## Slide 16: Q3, a seasonal mean distorts every rule except shifts

- One path shows exactly one phase gap between the two seasonal cycles. The sample correlation settles on a random non-zero value fixed by that gap.
- With P = 50, 46 of the 50 possible gaps give a limit clearly away from zero. That predicts 92% rejection, and Pearson rejects 91.4%.
- Proposition 1 proves the size goes to 1 for Pearson, wild and subsampling. Shifts re-randomise the gap and stay at 5.8%. For that part we only have a sketch.
- This is Yule's nonsense correlation from 1926, and we now have a proof for the kernel rules too. More data doesn't help, because a single path never shows you a second gap.

## Slide 17: Q4, shifts keep their power at a true 5%

- Tune every rule to a true 5%. Shifts then come within 4.9 points of the best other HSIC rule in all 18 cells.
- Raw power is misleading. HSIC subsampling shows 86.4% power at its nominal cut-off, and only 14.6% once it's tuned to a true 5%.
- The tuning uses the known null, so it's a way of checking the rules. You can't do it on real data.
- So shifts don't get their good size by refusing to reject.

## Slide 18: Q4, adjust the seasonal scale as well as the mean

- Kernel tests react to scale, so removing only the mean barely helps. Box 1 has 73.8% after mean-only adjustment and 7.2% after adjusting the scale too.
- A period estimated from the data did as well as the true period.
- Adjusting costs a little size under memory. HSIC with shifts goes from 5.3% to 6.8% on Gaussian pairs at H = 0.8.
- The habit of removing a seasonal mean isn't enough for kernel tests. Under a shared fixed phase, shifts still reject 28% to 34%, so adjust first.

## Slide 19: Q4, memory caps power against a mean shift

- At n = 2000, a shift of 0.1 is 2.24 null standard deviations at H = 0.5. At H = 0.8 it's only 0.40.
- At H = 0.8 no location-invariant 5% test can get above 12% power. MMD reaches 6.4% to 7.5%.
- The loss comes from the data, and a better rule can't win it back. A weak split-MMD result here doesn't mean the method failed.

## Slide 20: Formal results, what we prove and what we do not

- Proposition 1: with a random-phase seasonal mean, the correlation has a random non-zero limit. Any rule whose critical value is of smaller order than the seasonal term has size going to 1.
- Proposition 2: the half-mean variance is exact for every m, and the power cap holds for every location-invariant test.
- The asymptotic size of the shift rule is not yet proved. The argument depends on an unreviewed preprint (Düker and Zoubouloglou), and we don't treat it as established.
- We did prove one thing directly: inc1's statistic stays bounded at every H below 1, at a fixed bandwidth. The seasonal-scale failure is only measured.

## Slide 21: Q4, what is and is not proved

- At H = 0.8 the shift rates are 5.7%, 6.3% and 8.3% at n = 2000, 4000 and 8000. The intervals overlap, and the n = 8000 cell has only 144 replications.
- The shift rate might creep up with n. These runs can't show that, and they can't rule it out.
- Shift-rule size under long memory is still not yet proved.
- A window covering a quarter of a cycle can't be bin-adjusted, so the raw kernel has no fix for slow drift. In EXP7, inc1 with shifts gives 5.8% in that setting. That's simulation only.

## Slide 22: Stock trade durations

- Taking out the time-of-day pattern moves no cell by more than 0.6 points.
- About 44% of windows still have a variance trend, where 5% would be expected with no trend. Pairs with trends reject more often under subsampling. The shift rule doesn't go up with them.
- These are descriptive placebo rates. They are not false-alarm rates, because we can't verify independence and the same days show up in many pairs.
- The real data show the slow-drift pattern that the simulations flagged as the rule's weak spot.

## Slide 23: A kernel that forgets memory

- The increment kernel (inc1) compares each point with the one before it. Working with changes throws away the slow swings that long memory is made of.
- The Hermite split takes out the low polynomial layers of the normal scores, which is where the memory sits. It then tests those layers separately.
- For HSIC the climb with H stops. At H = 0.95, wild falls from 19.1% to 6.7% and subsampling from 62.4% to 2.7%.
- The wild excess doesn't go away. Even with no memory, the memory kernels reject 6.4% to 8.3% under wild, so the rule's own error is still there.

## Slide 24: What it fixes and what it does not

- Split MMD still climbs. On heavy tails at H = 0.9, datadiff rejects 17.4%, sym 11.6%, herm1 9.4% and hsplit1 8.0%.
- hsplit2 with block critical values held its size on both data types in simulation. No proof covers it.
- Every kernel has a blind spot. inc1 misses a flat moving average (4.3% against 71.1%) and a slow common factor (8.6% against 100%). herm1 misses linear dependence (6.0% against 58.0%).
- On the raw kernel, the Harris non-circular shift rejects 1.6% to 3.8% in all 44 stationary nulls. No memory kernel beats it everywhere.
- How I'd pitch inc1: it cures the memory climb in HSIC when you pair it with shifts or Harris. On its own it doesn't fix the size.

## Slide 25: Real stock data with a known null (EXP8)

**The problem it solves**

- On real stock data we never knew whether two series were truly independent. So the paper could only report placebo rates. Those describe how the tests behave, but they aren't false-alarm rates.
- EXP8 builds a null we know is true, out of real market data.

**What we did**

- The data are Nasdaq trade durations for AAPL, MSFT and NVDA, January to June 2024. We used the same windows as the paper's stock section. Each window stays inside one day, and paired days are at least 20 trading sessions apart.
- In each pair, X is a real window from one stock. Y is a shuffled copy of a real window from another stock. The copy keeps the window's exact values and, roughly, its spectrum, so it keeps the spread and the memory. Its timing is randomised, so it can't line up with X. Every rejection is therefore a false alarm.
- That gave 1,000 pairs per setting, over four stock pairings and two window lengths (500 and 2,000). We ran 8,000 pairs with a plain shuffled copy and 8,000 with a copy that also keeps the window's slow drift. We also ran 4,000 real-against-real placebo pairs for comparison.
- Before touching the real data, we checked the method on simulated fGn at H = 0.8. A shuffled copy gave the same rates as a fresh independent series, with shifts at 5.5% against 6.0%.

**What we found (n = 2000, four pairings pooled)**

- Circular shifts held: 4.9% with the ordinary kernel and 4.7% with inc1.
- Wild multipliers gave 6.3% (interval 5.6% to 7.1%) and subsampling 7.6% (6.8% to 8.5%). Both intervals sit above 5%.
- inc1 with subsampling gave 4.2%. The Harris test gave 2.3%, which is safe but cautious.
- The real-against-real placebo rates are close to these: 6.6% against 6.3% for wild, and 8.8% against 7.6% for subsampling.

**What it means**

- These are the project's first measured false-alarm rates on market data. The simulation result carries over: shifts stay at 5%, and the short-window rules run high.
- The placebo looks like the known null, so most of the excess in the paper's stock table is probably the rules misfiring. Real links between the stocks would add to it, and we see little room for them. The Y windows differ between the two arms, so this comparison is descriptive.
- Memory in these windows is mild (H around 0.61), so the excess is 1 to 3 points. That is far below the 19% we saw in strong-memory simulations.

## Slide 26: EXP8, what held, what failed, what it changes

- We wrote five predictions down before the run. Three failed, and we report them as failed.
- What held. The windows match the paper's saved windows exactly, so the placebo arm uses the same data as the paper's table. Shifts gave 4.1% to 5.1% in each pairing. Harris stayed at or below 3.2% in all 16 cells.
- Failed: the memory trend. Wild and subsampling are above 5%, but the rate didn't rise with each window's H. The windows' H values are too close together (spread 0.05) to show a trend.
- Failed: drift. Keeping each copy's slow drift added about one point (interval −0.3 to 1.9 for wild). So drift on its own explains at most a small part of the stock excess.
- Failed: inc1 under wild. The rate fell from 6.6% to 5.7%, but p = 0.11 isn't significant. Two other results did come through for inc1. The trend-linked gap shrank from 3.1 to 0.8 points. Under subsampling, inc1 fell from 8.8% to 4.3%, a result we didn't pre-register.
- Limits. The copy keeps values, spectrum and drift, but not volatility clustering. Days recur across pairs, so the pairs aren't fully independent of each other.
- For the paper: the stock section can now report measured false-alarm rates, not only placebo fractions.

## Slide 27: Answers to the four research questions

- Theory: the assumptions rule these data out by implication, and they say nothing about the size.
- Calibration: the rule matters as much as the statistic. On one cell it's 19.7%, 8.9% or 5.7% depending on the rule.
- Attribution: seasonality distorts size more than memory, and heavy tails barely matter.
- Procedure: adjust, then run HSIC with shifts. Slow drift is still a problem.

## Slide 28: The procedure for testing independence

- Step 1: remove the seasonal mean and the seasonal scale from each series. Step 2: run HSIC with circular shifts. Step 3: if a window covers less than one cycle, treat a rejection as weak evidence.
- The 3.8% to 7.0% range leaves out two cases: slow drift (10.1%) and a shared fixed phase (28% to 34%).
- For split MMD, no raw-kernel rule held size under memory. The hsplit2 block rule is a partial fix. It held in simulation, with no proof.
- If someone remembers one slide, it should be this one. Each step comes from one experiment, so each can be checked on its own.

## Slide 29: What was already known, and what is new

- Already known: i.i.d. assumptions, mixing assumptions, the H_X + H_Y = 3/2 line, oversized dCov subsampling, and Yule's nonsense correlation.
- New here: sizes for each rule side by side, process features crossed with memory, seasonal scale as a failure specific to kernels, Propositions 1 and 2, and a decision rule whose power we checked.
- Say plainly that the first list isn't ours.

## Slide 30: Closest prior work, Betken and Dehling

- They already have the null, the 3/2 line and oversized dCov subsampling.
- We add more rules (wild, shifts), more statistics (HSIC, split MMD), process features crossed with H, and the two propositions.
- Our dCov subsampling size is 10.6% at H = 0.7 and n = 2000. Theirs is 11.2% to 17.5% at H = 0.7.
- A referee will ask how we differ from them. Have this slide ready.

## Slide 31: Questions for the call (1 of 2), the proofs

- Ask 1: is a proof of shift-rule size realistic, and which case should come first? Right now it's not yet proved, and the asymptotic argument depends on an unreviewed preprint.
- The n = 8000 rate has an interval from 4.8% to 14.0%, so the finite-sample size is open as well.
- Ask 2: will Gareth read the proofs appendix? So far only a numerical script has checked the proofs.
- Point him to three steps: the ergodic lemma with periodic weights, the HSIC limit under median bandwidths, and the invariance reduction behind the power bound.

## Slide 32: Questions for the call (2 of 2), method and venue

- Ask 3: do we adopt inc1 with shifts as the slow-drift fix, or state slow drift as a limitation? Keep in mind inc1 misses a slow common factor.
- Ask 4: is a data-driven bandwidth experiment worth the time? Shifts held size at all four fixed multiples, and no multiple wins on power.
- Ask 5: which venue? The draft runs 25 pages, and AISTATS means cutting to about 8. The JFec special issue deadline has passed.

## Questions he might ask

1. "Is circular-shift size proved under long memory?" No, not yet. There's an argument, but it depends on an unreviewed preprint. For now the size is something we measured.
2. "Does the shift rate rise with n under strong memory?" It might. At H = 0.8 it's 5.7%, 6.3% and 8.3% at n = 2000, 4000 and 8000, but the intervals overlap, and the n = 8000 cell has only 144 replications.
3. "How general are the seasonal results?" They use one period, one amplitude and a uniform random phase. Think of them as illustrations. They aren't calibrated to market data.
4. "Do the stock results show false alarms?" No. We can't verify independence, and the same days recur across pairs. The rejection fractions are descriptive placebo rates.
5. "Why not just use the Harris shift test?" It has a finite-sample size bound, and it held 1.6% to 3.8% in every stationary null. But it gives up power when there's no memory. At H = 0.5 it trailed circular shifts by more than 5 points on five of six alternatives.
