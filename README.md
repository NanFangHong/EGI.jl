# Introduction

Suppose you need to recover the shape of a convex polyhedron, but all you know is the **normal** and **area** of each surface. Can you do that?

The answer is **yes**.

This is called **Minkowski problem** raised 100 years ago. Its result deeply influenced the development of modern geometry. It was proved that any convex object can be uniquely determined by its surface normal and surface area. In another words, one may store a convex 3-D object only by its surface normal and surface area. Such representation is called **Extended Gauss Image (EGI)**, which is broadly used for 3-D object alignment and recognition.

![1](https://github.com/NanFangHong/EGI.jl/blob/master/assets/Picture1.png)

# EGI.jl

![5](https://github.com/NanFangHong/EGI.jl/blob/master/assets/Animation.gif)

```
using Pkg
Pkg.add.(["LinearAlgebra", "Polyhedra", "QHull", "Makie"])
Pkg.clone("https://github.com/NanFangHong/EGI.jl")
```


The package **EGI.jl** can let you play with Minkowski problem effortlessly. 

```
using EGI
```

1. It defines a type called **ExtendedGaussImage**. You can either input your own EGI, e.g. 

```
faces = 6
img = ExtendedGaussImage((Float64.([1 -1 0 0 0 0;
    0 0 1 -1 0 0;
    0 0 0 0 1 -1]),Float64.([1,1,1,1,1,1][:,:])))
```

or randomly generated a valid EGI by,  

```
faces = 6 # suppose we want 6-faced polyhedron 
p, img = randPolyhedron(faces)
```
Now you have a polyhedron **p** and its EGI **img**


2. You can visualize **p** by, 

```
RenderPolyhedron(p)
```


3. You can have access to surface area, surface normal, volume by, 

```
Surface(p)
Volume(p)
Normal(p)
```

4. You can reconstruct 3-D shape from Extended Gauss Image 

```
k = 0.5 # Stepsize, if warning "attempt to access *-element Array{Float64,1} at index [*+1]" , try smaller k
ϵ = 1e-11 # Error tolerance 
qArray, KL = Little(k, ϵ, img) # Get a sequence of polyhedra q that approaches polyhedron p. Also get a sequence of error between those candidate polyhedra q and true polyhedron p. 
RenderPolyhedron(qArray[end]) 
```

5. You may wonder how James Little's algorithm works. I will show it later. You can make an animation by saving images of this converging process by, 

```
SavePolyhedronArray(qArray)
```

# Expected Output

```
faces = 8
```
8

```
p, img = randPolyhedron(faces)
```

(HalfSpace([-0.611479, -0.294878, 1.33617], 1.0) ∩ HalfSpace([-0.583018, 0.791774, 0.721503], 1.0) ∩ HalfSpace([-1.73855, -0.781748, -1.16332], 1.0) ∩ HalfSpace([0.748682, 1.29045, 1.34069], 1.0) ∩ HalfSpace([3.33446, -0.955458, 2.56791], 1.0) ∩ HalfSpace([0.472213, 0.531693, -1.67085], 1.0) ∩ HalfSpace([1.04939, 2.39937, -1.47103], 1.0) ∩ HalfSpace([-0.314485, -0.968992, -0.991615], 1.0) : convexhull([-0.222029, 0.134775, -0.618359], [0.619998, -0.141498, -0.468301], [0.270455, 0.41707, 0.193414], [-0.175158, -1.16343, 0.183981], [0.0396359, -0.334823, -0.693841], [0.604504, -0.596353, -0.617423], [-0.302647, -1.1395, 0.358431], [-0.1905, 0.161454, 0.696861], [-0.931669, 0.23603, 0.374133], [-0.432029, 0.368863, 0.632102], [-0.840786, 0.727945, -0.0922543], [-0.395061, 0.743413, 0.250945]), ExtendedGaussImage([-0.407998 -0.478047 … 0.349369 -0.221208; -0.196752 0.649217 … 0.798808 -0.681587; 0.891531 0.591598 … -0.489741 -0.6975], [0.621574; 0.338138; … ; 1.01718; 0.375791]))

```
RenderPolyhedron(p)
```
![2](https://github.com/NanFangHong/EGI.jl/blob/master/assets/GroundTruth.png)

```
k = 0.5 
ϵ = 1e-11
qArray, KL = Little(k, ϵ, img)
```

p2 ./ p1 = [1.15971; 1.48769; 0.733839; 1.12192; 0.874397; 1.60854; 0.824566; 1.89001]
KL = 0.04892604631423933
p2 ./ p1 = [1.04057; 1.31113; 0.866845; 1.09485; 0.94101; 1.37476; 0.874362; 1.44819]
KL = 0.015377621330277508
... 
p2 ./ p1 = [0.999993; 1.00001; 0.999999; 1.00001; 0.999999; 0.999996; 0.999997; 1.00001]
KL = 1.2129134430645393e-11
p2 ./ p1 = [0.999995; 1.00001; 0.999999; 1.00001; 1.0; 0.999997; 0.999998; 1.00001]
KL = 6.859450877437155e-12

(QHull.Polyhedron[HalfSpace([-0.756371, -0.364751, 1.65277], 1.0) ∩ HalfSpace([-0.873056, 1.18566, 1.08043], 1.0) ∩ HalfSpace([-2.00878, -0.90326, -1.34415], 1.0) ∩ ...)

```
RenderPolyhedron(qArray[1])
```
![3](https://github.com/NanFangHong/EGI.jl/blob/master/assets/Initial.png)

```
RenderPolyhedron(qArray[end])
```
![4](https://github.com/NanFangHong/EGI.jl/blob/master/assets/End.png)

# James Little's Algorithm 

![6](https://github.com/NanFangHong/EGI.jl/blob/master/assets/Little1.png)
![7](https://github.com/NanFangHong/EGI.jl/blob/master/assets/Little2.png)
![8](https://github.com/NanFangHong/EGI.jl/blob/master/assets/Little3.png)



# Reference

1. Gu, Xianfeng, et al. "Variational principles for Minkowski type problems, discrete optimal transport, and discrete Monge-Ampere equations." arXiv preprint arXiv:1302.5472 (2013).

2. Little, James J. "An Iterative Method for Reconstructing Convex Polyhedra From External Guassian Images." AAAI. 1983.
