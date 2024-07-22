# Trajectories Near a Black Hole
Computational project to model the trajectory of a negligable mass particle in the vicinity of a black hole. The project is written in Python using the [SciPy](https://scipy.org/), [NumPy](https://numpy.org/), and [Matplotlib](https://matplotlib.org/) libraries. 

Close to a massive, non-rotating body, the curvature of spacetime is given by the Schwarzchild metric, which is the simplest non-flat solution to Einstein's field equation. The metric can be specified through the length element $\text{d}s$, given here in spherical polar coordinates $t$, $r$, $\theta$, $\phi$ and geometrised units $r_\text{S}$, $c$ as:

$$\text{d}s^2 = \left(1-\frac{1}{r}\right)\text{d}t^2-\left(1-\frac{1}{r}\right)^{-1}\text{d}r^2-r^2\text{d}\theta^2 -r^2\sin^2\theta \\,\text{d}\phi^2.$$

The Schwarzchild radius, $r_\text{S}$, and the speed of light in a vacuum, $c$, constitute length and velocity scales. In these units and the Schwarzchild metric considering trajectories in the plane with $\theta=\frac{\pi}{2}$, the geodesic equations may be expressed as a set of six first order, coupled ODEs in terms of the four velocity components $U^t$, $U^r$, $U^\phi$:

$$\begin{split}         &\dfrac{\text{d}t}{\text{d}\tau} = U^t, \\
         &\dfrac{\text{d}r}{\text{d}\tau} =  U^r, \\
         &\dfrac{\text{d}\phi}{\text{d}\tau} = U^\phi,\\
         &\dfrac{\text{d}U^t}{\text{d}\tau} = -\dfrac{1}{r^2\left(1 - \frac{1}{r}\right)}U^r U^t,\\
         &\dfrac{\text{d}U^r}{\text{d}\tau} = -\dfrac{1}{2r^2}\left(1 - \frac{1}{r}\right)(U^t)^2+\dfrac{1}{2r^2\left(1 - \frac{1}{r}\right)}(U^r)^2+(r - 1)(U^\phi)^2,\\
    &\dfrac{\text{d}U^\phi}{\text{d}\tau} = -\dfrac{2}{r}U^rU^\phi.\end{split}$$

The program solves these equations from a given set of initial coordinates $t_0$, $r_0$ and $\phi_0$ and velocities $U^r_0$ and $U^\phi_0$ using the SciPy [``solve_ivp``](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#solve-ivp) function. To do it employs an implicit Runge-Kutta method. It is then able to plot the trajectory of the particle in the plane.



    
![1706624608850](https://github.com/TobyWoodcock/Black_hole_project/assets/160004842/9a9a7207-4b05-4116-a750-547f86026ba6)

![BK_orbit-ezgif com-speed](https://github.com/TobyWoodcock/Black_hole_project/assets/160004842/6820cbee-4e94-47f1-88ac-f8c3400444dc)

