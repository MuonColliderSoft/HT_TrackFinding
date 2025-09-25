[HT_toyPrototype](https://github.com/MuonColliderSoft/HT_TrackFinding/tree/main/HT_toyPrototype/README.md)

[HT_FullSimulation](https://github.com/MuonColliderSoft/HT_TrackFinding/tree/main/FullSimulation/README.md)

To check out only the FullSimulation package:
```
git clone --filter=blob:none --no-checkout git@github.com:MuonColliderSoft/HT_TrackFinding.git
cd HT_TrackFinding
git sparse-checkout init --cone
git sparse-checkout set HT_FullSimulation
git checkout
```
