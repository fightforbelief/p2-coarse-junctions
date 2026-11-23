# External trajectory reanalysis

Set before computing transfer scores on either measured record. The original 84 grid cases and sampling experiments remain frozen.

## Inputs and geometry

Use the passive calibration window of the LUMICKS C-Trap passive calibration file (Zenodo version 7729824, publisher MD5 cdc5f55039e4ff52d9c78eebf11f33c5), bead 1 high-frequency x and y channels if both are present and synchronous. Use the one-dimensional harmonic trace x in Argun et al., DeepCalib Example 1b at Git revision bea01cc09413d9a19c5e23aad9e594fe879fea61 (SHA-256 992deadeaf290153c45d7ff2bbb974d520a32587b1fa4500dbedf18d16ea2059). The harmonic trap uses active thermophoretic feedback. A reversible stochastic surrogate is a model class for the measured coordinates. Neither experiment provides a known microscopic state count. Exclude driven rotational examples. Report the channel limitation if a source only provides one axis.

Use the first fifth of each acquisition window to set each coordinate's mean and standard deviation. Score only the remaining four fifths, without thinning or detrending. The hub is the interval or square with standardized coordinates of absolute value at most 0.75. A port requires a coordinate beyond either face and, in two dimensions, transverse absolute coordinate at most 0.375. All remaining points form an exterior label. In one dimension the two tails are the ports.

## Endpoints and controls

Evaluate half-lags of 1, 4, 16, 64, 256 and 1024 samples and retain every result. Use the average of forward and backward joint masses as the input to the reversible model comparison. Retain unaveraged mass matrices and the maximum orientation mismatch. Estimate hub mass from the pooled endpoints of the one-step pairs. Record the all-pair two-step score, exact and occupancy-relaxed floors, and the single-state-per-label Chapman–Kolmogorov prediction. Report positive gaps exceeding 1e-12 as plug-in comparisons. They have no independent-binomial coverage interpretation.

Repeat separately in eight contiguous held-out blocks, discarding pairs crossing a block boundary. A time-shuffled held-out control with fixed NumPy seed 20261123 retains all one-time occupancies. Order is randomized once per record. Eight block results describe sensitivity to time within each record. They are not independent trials.

## Single-state trajectory test

Use consecutive excursions after the first visit to the hub on the held-out path at the original frame rate. Discard incomplete excursions and those whose first or last outside label is not a port. Form disjoint successive pairs of retained excursions. Count matches between the returning port of the first and the departing port of the second. Under a homogeneous single-state hub the completed excursions regenerate. Given the collection of completed excursions, window length and boundary information, their order is exchangeable. Randomize the even excursion marks conditional on both margins. Use 9,999 permutations and the plus-one right-tail p-value. Apply Holm's procedure over the two actual recordings. For time-shuffled controls report descriptive p-values. Report the number of pairs, match fraction and null expected fraction; the test is only for one hub state.

## Addendum before outcome access

Use the full path as a separate comparator for multi-state inference. After entry into the hub directly from port i, score the first subsequent outside frame: one if it is in port i, zero otherwise, including the unmonitored exterior. An unfinished terminal visit supplies no outcome. Include every completed eligible visit without pairing. For each of the six ports across the two recordings (or four if the passive record has a single axis), test whether any hub state returns to that port with probability above one half. Use the equally weighted mixture of Bernoulli likelihood ratios with null 0.5 and alternatives 0.55, 0.65, 0.75 and 0.90; use the final e-value and threshold 6/0.05 for the two actual recordings. Report the number of port-specific rejections as a lower bound, with minimum one. Apply the same calculation descriptively to shuffled controls. This test uses entrance and departure events, not just endpoint marginals. It was specified after the HDF5 schema inspection and before running any transfer or excursion calculation.
