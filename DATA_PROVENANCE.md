# Data provenance

The structural test graph is a processed source-target edge list from the public
SFO/P3R model in the UCF-RISES/OutageMap repository. That repository contains
OpenDSS files and an `edgeList.csv` for P3R and documents its SMART-DS-based
power-network workflow.

Source repository:
https://github.com/UCF-RISES/OutageMap

Exact source file used for the processed edge list:
`P3R/edgeList.csv`

Pinned Git commit inspected for this package:
`7f989a1f424e384dcd22d0ab10a76d4b5fce90ee`

Git blob SHA for the source `P3R/edgeList.csv`:
`2946e67ebd871ac80e7c215b94e270ac8db7f32d`

The file `data/smartds_sfo_p3r_edge_pairs.csv` retains only the `source` and
`target` columns used by the deletion-radius experiment. It has 766 directed
edges over 767 nodes.

The underlying SMART-DS dataset is the NREL/OEDI dataset titled
"SMART-DS Synthetic Electrical Network Data OpenDSS Models for SFO, GSO, and AUS"
(Palmintier et al., 2020). SMART-DS is a realistic synthetic power-distribution
benchmark, not real utility operating data.

OEDI landing page/search record:
https://data.openei.org/submissions/2981

Important scope statement: this P3R topology is from the SMART-DS San Francisco
data family. It is not claimed to be the exact 10,090-node subnetwork used in
Xu and Moghaddass (IISE Transactions, 2023).
