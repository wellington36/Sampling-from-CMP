# Sampling from CMP
First we sample

$$
\lambda_i \sim \mathrm{Unif}\left(\tfrac{1}{2}, 20 \right), \quad
\nu_i \sim \mathrm{Unif}\left(\tfrac{1}{2}, 5 \right)
$$

for $i = 1, 2, \cdots, 200$, then we sample two sets of variables. For each method we sample

$$
X_{i,j} \sim \mathrm{CMP}(\lambda_i, \nu_i), \quad
Y_{i,j} \sim \mathrm{CMP}(\lambda_i, 1)
$$

for $i = 1, 2, \cdots, 200$ and $j = 1, 2, \cdots, 10000$.

We show the results of $X_{i,j}$ for the mean time in seconds, the mean chi-square test, and the mean error:

|                   | rejection RoU | rejection_benson | exact  |
|-------------------|-----------|------------------|--------|
| Time(s)           | 0.0593    | 0.0305           | 0.1963 |
| Chi2 p-val        | 0.1884    | 0.9949           | 0.9959 |
| \|Mean - TrueMean\| | 0.0578    | 0.0110           | 0.0135 |

And we show the equivalent result but for $Y_{i,j}$, comparing also with the `rpois` function from the R package (since $\nu = 1$ corresponds to the Poisson distribution):

|                   | rejection RoU | rejection_benson | exact  | poisson (rpois)  |
|-------------------|-----------|------------------|--------|--------|
| Time(s)           | 0.0732    | 0.0182           | 0.2506 | 0.0005 |
| Chi2 p-val        | 0.7182    | 0.9700           | 0.9816 | 0.9678 |
| \|Mean - TrueMean\| | 0.0434    | 0.0216           | 0.0236 | 0.0230 |


## PMF comparation
Run `multiplots.R` for

<img width="1291" height="785" alt="image" src="https://github.com/user-attachments/assets/c2a2aede-601b-4dbb-8feb-ee64792fe4a1" />
