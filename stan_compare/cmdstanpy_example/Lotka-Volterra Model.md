# Lotka-Volterra Model

The Lotka-Volterra equations are a pair of first-order, non-linear, differential equations frequently used to describe the dynamics of biological systems in which two species interact, predator and prey.

## Equations

The equations are:

$$
\begin{align*}
\frac{dx}{dt} &= \alpha x - \beta xy \\
\frac{dy}{dt} &= \delta xy - \gamma y
\end{align*}
$$

where:
- \( x \) is the number of prey,
- \( y \) is the number of predators,
- \( alpha \) : The growth rate of the prey population.
- \( beta \)  : The death rate of prey per predator.
- \( gamma \) : The growth rate of the predator per prey.
- \( delta \) : The rate at which predators die.

## Python Implementation

```python
def ode_system(t, z, alpha, beta, delta, gamma):
    x, y = z
    dxdt = alpha * x - beta * x * y
    dydt = gamma * x * y - delta * y
    return [dxdt, dydt]

```

## Stan Implemenatation

```
array[] real ODE_sys(real time,  array[] real x, array[] real parms, array[] real rdata, array[] int idata) {
     // the array of parameters invloved in ode equations
     real alpha = parms[1];
     real beta = parms[2];
     real delta = parms[3];
     real gamma = parms[4];

     // the system of ODEs
     array[2] real dxdt;
     dxdt[1] = alpha * x[1] - beta * x[1] * x[2];
     dxdt[2] = gamma * x[1] * x[2]  - deltas * x[2];

     return dxdt;
   }
```