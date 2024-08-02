# Two Compartment Model

This two compartment system is an uncoupled pair of first-order, non-linear, differential equations used to describe the dynamics of biological systems in which two species DO NOT interact.

The equations are:

$$
\begin{align*}
\frac{dx}{dt} &= \alpha x - \beta x \\
\frac{dy}{dt} &= \delta y - \gamma y
\end{align*}
$$

where:
- \( x \) is the number of entities in Compartment 1,
- \( y \) is the number of entities in Compartment 2,
- \( alpha \) : The rate of increase in x.
- \( beta \)  : The rate of loss of x.
- \( gamma \) : The rate of increase in y.
- \( delta \) : The rate at which y deacreses.


## Python Implementation

To run the model, execute the following command in your terminal:

```bash
python Two_compartment_model.py