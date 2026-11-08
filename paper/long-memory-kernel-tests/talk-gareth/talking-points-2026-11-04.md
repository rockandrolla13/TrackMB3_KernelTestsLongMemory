# Talking points: call with Gareth Peters, 4 November 2026

Companion to slides-2026-11-04 (30 slides). Glance at the bold line under each slide. Numbers come from the slides or the paper only.

## The point of the research in 30 seconds

Kernel tests like HSIC promise a 5% false-alarm rate, but the proofs assume memory dies out fast. Long-memory data break that assumption, and we measure how badly. The answer depends more on the critical-value rule than on the statistic, and seasonality does more damage than memory. A seasonal scale shared by two independent series fools kernel tests even when Pearson correlation is fine. Our advice is to remove the seasonal mean and scale, then run HSIC with circular shifts, which held 3.8% to 7.0% in our simulations at n = 2000.

## Three things to remember

1. **The rule that sets the critical value matters as much as the test.** On one cell (H = 0.8, n = 2000), HSIC rejects 19.7% with subsampling, 8.9% with wild multipliers and 5.7% with circular shifts. A user who picks a test but not a rule has not fixed the false-alarm rate.
2. **Kernel tests have an extra exposure that Pearson does not.** A seasonal scale shared by two independent series drives HSIC wild to 73.3% with no memory at all, while Pearson stays at 4.7%. The kernel sees that both series swing widely at the same times. This is the reason to use a kernel, and on ordered data it is also the risk.
3. **Changing the kernel moves the problem. It does not remove it.** In EXP7, kernels built to remove memory stopped the HSIC climb with H. But none passed the pre-registered size rule under wild or subsampling, and each is blind to some dependence. The Harris (2020) non-circular shift test already held 1.6% to 3.8% in every stationary null.

## Slide 1: Title

- Skip. Say the aim: where kernel tests fail on long memory, why, and what to do.

## Slide 2: Research problem and questions

- We test independence of two series from one observed path. Each series may have long memory, meaning correlations fade so slowly they never add up to a finite total.
- The Hurst index H measures memory. H = 0.5 means none, values near 1 mean strong memory. Our baseline is fractional Gaussian noise (fGn).
- A rule estimates the critical value (the cut-off the statistic must pass to reject) from that same path. The rule is calibrated if it rejects a true null 5% of the time.
- Four questions: does the theory cover these data, how wrong is each rule, what causes the error, and what should a user do.

## Slide 3: Why the 5% guarantee is in doubt

- The 5% guarantee rests on four proof steps, each assuming the one before.
- Step 1 assumes distant parts of a series are nearly unrelated (mixing). fGn with H above 1/2 is not mixing.
- So steps 2 to 4 are not wrong, just unproved. No one knows the null limit of HSIC for independent fGn, so we measure the false-alarm rate.
- **Why it matters:** the failure is not a bug in one paper. The whole chain has no foundation for these data, so simulation is the only evidence available.

## Slide 4: Three kinds of statistic, each needing a critical value

- Split MMD compares the two halves of one series. HSIC and distance covariance compare two aligned series. Pearson is the linear reference.
- The critical-value rule is the part that has to handle the dependence. We treat it as a factor in its own right.
- An oracle that uses the known null shows what each rule costs.

## Slide 5: Long memory reduces the information in a path

- The gap between two half-sample means has variance (4 - 2^(2H)) m^(2H-2). This is exact (Proposition 2).
- At H = 0.5 it shrinks like 1/m. At H = 0.8 it shrinks only like m^-0.4.
- **Why it matters:** a long-memory series of 2000 points carries far less information than 2000 independent points. The power cap on slide 19 follows from this.

## Slide 6: Three critical-value rules

- Subsampling uses windows of about the square root of n points, which is 44 points at n = 2000. The wild rule weights lags on the same scale.
- Long memory moves slowly, so a short window misses most of the movement. The critical value then comes out too low, and the test over-rejects.
- **Why it matters:** this predicts which rules fail. Rules that look at short stretches should fail, and a rule that uses the full path should do better.

## Slide 7: The circular shift

- We wrap Y around by an offset and recompute the statistic. Y keeps its own memory. Only its alignment with X moves.
- Many offsets give the reference distribution from the full series. We exclude offsets below 0.1n.
- The original result (Chwialkowski and Gretton, 2014) assumes mixing, which fGn lacks. So shifts have no proof here either.

## Slide 8: How we test each question

- Theory: read the conditions in the primary sources and check fGn against each one.
- Calibration and attribution: Monte Carlo under a known null, with one process feature switched on at a time and crossed with H.
- Procedure: compare power at matched size, because a test can hold size by never rejecting.
- Stock durations act as an external check only. Their null is unverified, so the rates are descriptive.

## Slide 9: Monte Carlo design

- We set the generator, so the null is true by construction. Any rejection is a false alarm.
- With 1,000 replications, a 5% rate has a 95% Wilson interval of about plus or minus 1.4 points. That is counting error only.
- For EXP2, a design record fixed the hypotheses and generator checks before the main run.

## Slide 10: The experiments

- EXP1 maps calibration across H and n. EXP2 crosses seven generators with H.
- EXP3 measures power with each rule tuned to a true 5%. EXP4 tests seasonal adjustment.
- EXP5 and EXP6 check larger n and bandwidths. Stocks are AAPL, MSFT and NVDA trade durations, January to June 2024.

## Slide 11: Q1, the published conditions exclude fGn with H above 1/2

- Every audited source assumes i.i.d. data or mixing at a summable rate. fGn with H above 1/2 fails both.
- The exclusion is mostly by implication. Our own derivation gives tau(r) at least of order r^(2H-2), which breaks the wild bootstrap rate condition.
- Betken and Dehling cover distance covariance in part. No source gives a Gaussian-kernel HSIC limit. None treats seasonality.
- **Why it matters:** the theory warns, but it gives no size for the failure. A user cannot tell from the papers whether the error is 6% or 60%.

## Slide 12: Q2, the rule changes the size

- Same generator and same n, only the rule differs, so the rule accounts for the gaps between the lines.
- Part of the excess is not memory. The wild rule already sits at 6.5% at H = 0.5.
- **Why it matters:** reporting "HSIC rejected" without the rule is not enough. The rule is a first-class choice.

## Slide 13: Q2, higher Hermite rank delays the distortion

- Hermite rank is the order of the first non-zero term when a feature is expanded in Hermite polynomials. Rank m gives covariance of order rho(k)^m.
- Rank one is affected for every H above 1/2. Rank two balances at H = 3/4. Independent rank-one products balance at H_X + H_Y = 3/2.
- The calculation explains the shape of the curves. It gives no critical value and does not locate a proved boundary.

## Slide 14: Q3, seasonality distorts size more than memory

- We switch on one feature at a time and cross it with memory. Heavy tails matter little.
- HSIC shift stays between 3.8% and 7.0% across these cells.
- The threshold AR adds short-range persistence, so its excess is not purely nonlinearity. The seasonal settings are illustrative, not calibrated to markets.
- **Why it matters:** a user worried about long memory may be looking at the wrong problem. Seasonality is the bigger risk in our settings.

## Slide 15: Q3, a seasonal scale distorts kernel tests unless they use shifts

- The series are independent, each with its own random phase, and no memory.
- HSIC wild, HSIC subsampling and dCov subsampling reject 46% to 73%. HSIC shift rejects 5.0%. Pearson rejects 4.7%.
- This is measured, with no proof yet.
- **Why it matters:** this failure is specific to kernel tests. Pearson ignores it, so standard linear diagnostics will not warn you.

## Slide 16: Q3, a seasonal mean distorts every rule except shifts

- One path contains one phase gap between the two seasonal cycles. The sample correlation converges to a random non-zero limit set by that gap.
- With P = 50, 46 of 50 gaps give a limit clearly away from zero. That predicts 92%, and Pearson rejects 91.4%.
- Proposition 1 proves size tends to 1 for Pearson, wild and subsampling. Shifts re-randomize the gap and stay at 5.8%, but we can only sketch why.
- **Why it matters:** this is Yule's nonsense correlation (1926), now proved for kernel rules too. More data does not help, because one path never shows a second gap.

## Slide 17: Q4, shifts keep their power at a true 5%

- With every rule tuned to a true 5%, shifts stay within 4.9 points of the best other HSIC rule in all 18 cells.
- Raw power misleads. HSIC subsampling shows 86.4% at its nominal cut-off and only 14.6% at a true 5%.
- The tuning uses the known null, so it checks the rules. It is not a method for real data.
- **Why it matters:** the good size of shifts is not bought by never rejecting.

## Slide 18: Q4, adjust the seasonal scale as well as the mean

- Kernel tests respond to the scale, so removing only the mean barely helps. Box 1 gives 73.8% after mean-only adjustment, against 7.2% after scale adjustment.
- A period estimated from the data worked as well as the true period.
- Adjustment costs some size under memory. HSIC shift goes from 5.3% to 6.8% on Gaussian pairs at H = 0.8.
- **Why it matters:** the usual habit of removing a seasonal mean is not enough for kernel tests. Under a shared fixed phase, shifts still reject 28% to 34%, so adjust first.

## Slide 19: Q4, memory caps power against a mean shift

- At n = 2000, a shift of 0.1 is 2.24 null standard deviations at H = 0.5, and only 0.40 at H = 0.8.
- No location-invariant 5% test can pass 12% power at H = 0.8. MMD reaches 6.4% to 7.5%.
- **Why it matters:** this loss is in the data, not in the kernel. No better rule can recover it, so a weak split-MMD result is not a method failure.

## Slide 20: Formal results, what we prove and what we do not

- Proposition 1: under a random-phase seasonal mean, the correlation has a random non-zero limit. Any rule with a critical value of smaller order than the seasonal term has size tending to 1.
- Proposition 2: the half-mean variance is exact for every m, and the power cap holds for every location-invariant test.
- The asymptotic size of the shift rule is **not yet proved**. The argument rests on an unreviewed preprint (Düker and Zoubouloglou). We do not treat it as established.
- Proved directly: inc1's statistic stays bounded at every H below 1, at a fixed bandwidth. The seasonal-scale failure is only measured.

## Slide 21: Q4, what is and is not proved

- Shift rates at H = 0.8 are 5.7%, 6.3% and 8.3% at n = 2000, 4000 and 8000. The intervals overlap, and n = 8000 has only 144 replications.
- So the shift rate may drift up with n. These runs do not establish a drift, and they do not rule one out.
- Shift-rule size under long memory remains **not yet proved**.
- A quarter-cycle window cannot be bin-adjusted, so the raw kernel has no correction for slow drift. In EXP7, inc1 with shifts gives 5.8% there, in simulation only.

## Slide 22: Stock trade durations

- Removing the time-of-day pattern moves no cell by more than 0.6 points.
- About 44% of windows still show a variance trend, against 5% expected with no trend. Pairs with trends reject more under subsampling. The shift rule does not rise in step.
- These are descriptive placebo rates, **not false-alarm rates**. We cannot verify independence, and days recur across pairs.
- **Why it matters:** real data show the slow-drift pattern that simulation flagged as the rule's weak point.

## Slide 23: A kernel that forgets memory

- The increment kernel (inc1) compares each point with the one before it. Changes forget the slow swings that make long memory.
- The Hermite split removes the low polynomial layers of the normal scores, where memory sits, and tests those layers on their own.
- For HSIC the climb with H stops. At H = 0.95, wild falls from 19.1% to 6.7%, and subsampling from 62.4% to 2.7%.
- The wild excess stays. Memory kernels reject 6.4% to 8.3% under wild even with no memory, so the rule's own error remains.

## Slide 24: What it fixes and what it does not

- Split MMD still climbs. On heavy tails at H = 0.9, datadiff rejects 17.4%, sym 11.6%, herm1 9.4% and hsplit1 8.0%.
- hsplit2 with block critical values held size on both data types in simulation. No proof covers it.
- Each kernel has a blind spot. inc1 misses a flat moving average (4.3% against 71.1%) and a slow common factor (8.6% against 100%). herm1 misses linear dependence (6.0% against 58.0%).
- The Harris non-circular shift on the raw kernel rejects 1.6% to 3.8% in all 44 stationary nulls. No memory kernel beats it everywhere.
- **Why it matters:** present inc1 as a cure for the memory climb in HSIC, paired with shifts or Harris. It is not a fix for size.

## Slide 25: Answers to the four research questions

- Theory: the assumptions exclude these data by implication and give no size.
- Calibration: the rule matters as much as the statistic (19.7%, 8.9%, 5.7% on one cell).
- Attribution: seasonality distorts size more than memory. Heavy tails matter little.
- Procedure: adjust, then HSIC with shifts. Slow drift remains.

## Slide 26: The procedure for testing independence

- Step 1: remove the seasonal mean and the seasonal scale from each series. Step 2: HSIC with circular shifts. Step 3: treat a rejection as weak if a window covers under one cycle.
- The 3.8% to 7.0% range excludes slow drift (10.1%) and a shared fixed phase (28% to 34%).
- For split MMD, no raw-kernel rule held size under memory. The hsplit2 block rule is one partial fix, held in simulation without a proof.
- **Why it matters:** this is the one slide a practitioner takes home. Each step comes from one experiment, so each can be checked.

## Slide 27: What was already known, and what is new

- Known: i.i.d. assumptions, mixing assumptions, the H_X + H_Y = 3/2 line, oversized dCov subsampling, and Yule's nonsense correlation.
- New: rule-by-rule sizes side by side, features crossed with memory, seasonal scale as a kernel-specific failure, Propositions 1 and 2, and a decision rule with power checked.
- Be explicit that we do not claim the known items.

## Slide 28: Closest prior work, Betken and Dehling

- They already have the null, the 3/2 line and oversized dCov subsampling.
- We add more rules (wild, shifts), more statistics (HSIC, split MMD), process features crossed with H, and the two propositions.
- Our dCov subsampling size is 10.6% at H = 0.7 and n = 2000, against their 11.2% to 17.5% at H = 0.7.
- **Why it matters:** any referee will ask how we differ from them. This slide is the answer.

## Slide 29: Questions for the call (1 of 2), the proofs

- Ask 1: is a proof of shift-rule size realistic, and for which case first? Today it is **not yet proved**. The asymptotic argument rests on an unreviewed preprint.
- The n = 8000 rate has an interval of 4.8% to 14.0%, so the finite-sample size is open too.
- Ask 2: will Gareth read the proofs appendix? Only a numerical script has checked the proofs so far.
- Point him to three steps: the ergodic lemma with periodic weights, the HSIC limit under median bandwidths, and the invariance reduction behind the power bound.

## Slide 30: Questions for the call (2 of 2), method and venue

- Ask 3: adopt inc1 with shifts as the slow-drift correction, or state slow drift as a limit? Remember inc1 misses a slow common factor.
- Ask 4: is a data-driven bandwidth experiment worth the time? Shifts held size at all four fixed multiples, and no multiple wins on power.
- Ask 5: which venue? The draft has 25 pages, and AISTATS means cutting to about 8. The JFec special issue date has passed.

## Likely questions and short answers

1. **"Is circular-shift size proved under long memory?"** No, it is not yet proved. An argument exists, but it rests on an unreviewed preprint. Today the size is a measured result only.
2. **"Does the shift rate rise with n under strong memory?"** It may. At H = 0.8 the rates are 5.7%, 6.3% and 8.3% at n = 2000, 4000 and 8000, but the intervals overlap. The n = 8000 cell has only 144 replications.
3. **"How general are the seasonal results?"** They use one period, one amplitude and a uniform random phase. They are illustrative and not calibrated to market data.
4. **"Do the stock results show false alarms?"** No. Independence is unverified and days recur across pairs. The rejection fractions are descriptive placebo rates only.
5. **"Why not just use the Harris shift test?"** It has a finite-sample size bound and held 1.6% to 3.8% in every stationary null. But it lost power without memory, trailing circular shifts by more than 5 points on five of six alternatives at H = 0.5.

EXP8 (real-data surrogate null): results to be added.
