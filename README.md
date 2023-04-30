# Two-body problem
Simulation of the two-body problem with inital values / parameters:
- **Body 1:** mass $m_1$, 2D space vector $\vec{r_1}$, 2D velocity vector $\vec{v_1}$
- **Body 2:** mass $m_2$, 2D space vector $\vec{r_2}$, 2D velocity vector $\vec{v_2}$
- **Algorithm:** stepsize $h$

## 1. Run simulation and plot trajectories
```
make
```
> write coordinates to *build/name.csv*

> plots of trajectories in *build/name.pdf*

for each time-stepsize **h** and for both nummerically solving methods **euler** and **verlet**


## 2. Create an animation of the trajectories
```
make animation
```
> movies in videos/build/name.mp4


Hint: all file names include the method **euler/verlet** and the choosen stepsize **h**
