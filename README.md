# Sampling from CMP
Considering

$$
  \lambda \sim Unif(0.5, 5) \quad \text{ and } \quad \nu \sim Unif(0.5, 2)
$$

### Comparison methods
|                   | rejection RoU | rejection_benson | exact  |
|-------------------|-----------|------------------|--------|
| Time(s)           | 0.0278    | 0.3576           | 0.0700 |
| Chi2 p-val        | 0.1864    | 0.7500           | 1.0000 |
| \|Mean - TrueMean\| | 0.2233    | 4.0772           | 0.0309 |


### Comparison with $\nu = 1$

|                   | rejection RoU | rejection_benson | exact  | poisson (rpois)  |
|-------------------|-----------|------------------|--------|--------|
| Time(s)           | 0.0303    | 0.0139           | 0.0724 | 0.0003 |
| Chi2 p-val        | 0.2808    | 0.9500           | 1.0000 | 0.9500 |
| \|Mean - TrueMean\| | 0.1772    | 0.0311           | 0.0324 | 0.0335 |


## PMF comparation
Run `multiplots.R` for

<img width="1291" height="785" alt="image" src="https://github.com/user-attachments/assets/c2a2aede-601b-4dbb-8feb-ee64792fe4a1" />
