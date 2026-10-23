# Paired multiplier control design

The candidate comparison was defined before reading the 1,000-replication control result. The earlier 100-replication run is a pilot only.

- **Null:** two independently generated `H=0.5` Gaussian fGn processes.
- **Sample sizes:** 500, 1,000, 2,000; 1,000 independent pairs at each size.
- **Tests:** Gaussian-kernel HSIC with the repository's median bandwidth and distance covariance with the repository's double-centred distances.
- **Only experimental change:** use the same Markov Rademacher draws either after empirical column centering (current rule) or raw (candidate). Both use 500 draws, length `floor(sqrt(n))`, the same observed statistic, and p-value `(1+#(bootstrap >= observed))/501`.
- **Randomness:** one deterministic seed per `(n, replication)`; each pair and its multiplier draws are shared by the two calibrations. `run_check.py` records the script, core-code and output hashes.
- **Descriptive advancement rule:** the raw variant advances to a matched E1-window comparison only if all six raw control rates lie within 3–7%, no raw cell has a larger absolute distance from 5% than its centred counterpart, and the maximum absolute distance across cells is smaller. This is a practical pilot gate, not a theorem or a guarantee at other memory values.
- **If the gate fails:** keep the existing E1 rates as placebo observations, report the control failure, and remove bootstrap-based interpretation of duration dynamics from the paper headline. Do not subtract a baseline percentage from E1 rates.

The two multiplier implementations were checked to be numerically identical for the centred case at all three sample sizes under a shared random seed. The 100-pair pilot had wide Wilson intervals and was not used to choose a winner.
