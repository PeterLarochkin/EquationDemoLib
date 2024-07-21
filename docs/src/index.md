# EquationDemoLib.jl
## Experiments
This library compares the compilation and evaluation speeds of two identical programs. Results are shown below.
```bash
\# Juila program
n:0, diff:0.0306972956
n:1000, diff:0.0029081786
n:2000, diff:0.0009298924
n:3000, diff:0.0003381923
n:4000, diff:0.0001284034
n:5000, diff:0.0000494494
n:6000, diff:0.0000191328
n:7000, diff:0.0000074143
n:8000, diff:0.0000028746
n:9000, diff:0.0000011147
total count:9117
diff:0.0006959054
Elapsed time: 5.566 seconds

```

```bash
\# C program
n:0, diff:0.0306972956
n:1000, diff:0.0029081786
n:2000, diff:0.0009298924
n:3000, diff:0.0003381923
n:4000, diff:0.0001284034
n:5000, diff:0.0000494494
n:6000, diff:0.0000191328
n:7000, diff:0.0000074143
n:8000, diff:0.0000028746
n:9000, diff:0.0000011147
total count: 9117
diff:0.0006959054
Elapsed time: 5.85 seconds
```



## Documentation

```@docs
q
```
```@docs
u
```
```@docs
F
```
```@docs
psi
```
```@docs
apply_A
```
```@docs
get_B
```
```@docs
scalar_product
```
```@docs
proximity_search
```
