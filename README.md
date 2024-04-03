# ParallelBoid

ParallelBoid is a parallelization effort of the boids algorithm developed by Craig Reynolds in 1987, which can be seen [here](https://www.cs.toronto.edu/~dt/siggraph97-course/cwr87/). The algorithm will be implemented in serial in C++, then adapted to both OpenMP and CUDA.

![ParallelBoid](Pics/ParallelBoid.png "ParallelBoid")

## Build
'make boidgen'

## Set up Parameters
- Create "parameters.ini" file in root
- Populate it with this format:

```
[boids]
nboids=3
time=20
[weights]
w_collision=0.3
w_alignment=0.4
w_centering=0.3
```

## Run
'./build/boids'