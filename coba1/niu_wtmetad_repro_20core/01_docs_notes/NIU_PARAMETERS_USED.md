# Niu-style parameters used in this workspace

| Item | Niu paper/SI | Workspace setting |
|---|---:|---:|
| Target phase | beta-cristobalite | beta-cristobalite |
| Main CV | XRD `{111}` intensity | Si-only Debye-style XRD `{111}` |
| System size for Fig. 7 | 1536 atoms / 512 Si | 1536 atoms / 512 Si |
| Temperature | 2300 K | 2300 K |
| WTMETAD bias factor | 100 | 100 |
| Hill width | 5 CV units | 5 CV units |
| Hill height | 40 kJ/mol | 9.560 kcal/mol |
| Hill deposition | every 1 ps | every 4000 steps at 0.25 fs |
| Snapshot interval | 40 ps | dump stride 160000 at 0.25 fs |
| Local entropy rm | 0.75 nm | 7.5 Å |
| Local entropy ra | 0.45 nm | 4.5 Å |
| Local entropy sigma | 0.05 nm | 0.5 Å |

Q for `{111}`:

d111 = a / sqrt(3), with a = 7.15 Å.
Q111 = 2*pi/d111 ≈ 1.522 Å^-1.
