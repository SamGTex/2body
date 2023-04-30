# Two-body problem
Simulation of the two-body problem, described by Newton's equation of motion and Newton's law of gravitation.

Using two numerical methods to solve the differential equation:
> euler algorithm

> verlet algorithm

With inital values / parameters:
- **Body 1:** mass $m_1$, 2D space vector $\vec{r_1}$, 2D velocity vector $\vec{v_1}$
- **Body 2:** mass $m_2$, 2D space vector $\vec{r_2}$, 2D velocity vector $\vec{v_2}$
- **Algorithm:** stepsize $h$

Attention: Currently you must enter the named parameters in main.cpp!
> (Comming soon: config file with inital values)

## 1. Run simulation and plot trajectories
```
make
```
- write coordinates to *build/name.csv*

- plots of trajectories in *build/name.pdf*

> for each time-stepsize **h** and for both nummerically solving methods **euler** and **verlet**


## 2. Create an animation of the trajectories
```
make animation
```
- movies in videos/build/name.mp4


Hint: all file names contain the used method **euler/verlet** with chosen stepsize **h**

## Demo
Time Interval: [0, 100]

Too small stepsize h=1.0
https://user-images.githubusercontent.com/55881123/235359591-abadd867-2ad6-46d2-b170-4b4369e06764.mp4

Adequate stepsize h=0.01
https://user-images.githubusercontent.com/55881123/235359713-9f502876-1003-444e-9146-5ef892bdd4e0.mp4


## Required Packages/Tools
Install python libraries with Conda:
```
conda install numpy matplotlib pandas
```

Install ffmpeg with package manager (so save matplotlib animations as .mp4)
```
sudo apt-get install ffmpeg
```
