# TightBinding.jl
This can construct the tight-binding model and calculate energies in Julia 1.0. 


This can 
1. construct the Hamiltonian as a functional of a momentum k.
2. plot the band structure.
3. show the crystal structure.
4. plot the band structure of the finite-width system with one surface or boundary.

There is the sample jupyter notebook. 

## Install

```
add https://github.com/cometscome/TightBinding.jl
```

## sample
Here is a Graphene case

```julia
using TightBinding
#Primitive vectors
a1 = [sqrt(3)/2,1/2]
a2= [0,1]
#set lattice
la = set_Lattice(2,[a1,a2])
#add atoms
add_atoms!(la,[1/3,1/3])
add_atoms!(la,[2/3,2/3])
```
Then we added two atoms (atom 1 and atom 2). 
We can see the possible hoppings.

```julia
show_neighbors(la)
```

Output is

```
Possible hoppings
(1,1), x:-1//1, y:-1//1
(1,2), x:-2//3, y:-2//3
(2,2), x:-1//1, y:-1//1
(1,1), x:-1//1, y:0//1
(1,2), x:-2//3, y:1//3
(2,2), x:-1//1, y:0//1
(1,1), x:-1//1, y:1//1
(1,2), x:-2//3, y:4//3
(2,2), x:-1//1, y:1//1
(1,1), x:0//1, y:-1//1
(1,2), x:1//3, y:-2//3
(2,2), x:0//1, y:-1//1
(1,1), x:0//1, y:0//1
(1,2), x:1//3, y:1//3
(2,2), x:0//1, y:0//1
(1,1), x:0//1, y:1//1
(1,2), x:1//3, y:4//3
(2,2), x:0//1, y:1//1
(1,1), x:1//1, y:-1//1
(1,2), x:4//3, y:-2//3
(2,2), x:1//1, y:-1//1
(1,1), x:1//1, y:0//1
(1,2), x:4//3, y:1//3
(2,2), x:1//1, y:0//1
(1,1), x:1//1, y:1//1
(2,2), x:1//1, y:1//1
```

If you want to construct the Graphene, you choose hoppings from atom 1 to atom 2: 

```julia
#construct hoppings
t = 1.0
add_hoppings!(la,-t,1,2,[1/3,1/3])
add_hoppings!(la,-t,1,2,[-2/3,1/3])
add_hoppings!(la,-t,1,2,[1/3,-2/3])
```


```julia
#show the lattice structure
plot_lattice_2d(la)
```

![68747470733a2f2f71696974612d696d6167652d73746f72652e73332e616d617a6f6e6177732e636f6d2f302f3234363131332f37346633306563662d356137632d643337362d393235642d6561663563343634376362632e706e67](https://user-images.githubusercontent.com/21115243/46902071-aa0da680-cef9-11e8-9d4b-3cfa41633dc9.png)


```julia
# Density of states
nk = 100 #numer ob meshes. nk^d meshes are used. d is a dimension.
plot_DOS(la, nk) 
```

![68747470733a2f2f71696974612d696d6167652d73746f72652e73332e616d617a6f6e6177732e636f6d2f302f3234363131332f39343635643263312d643466332d333634372d363036652d3836626263313462313530622e706e67](https://user-images.githubusercontent.com/21115243/46902081-cc072900-cef9-11e8-8e22-908f91b132a8.png)

```julia
#show the band structure
klines = set_Klines()
kmin = [0,0]
kmax = [2π/sqrt(3),0]
あdd_Kpoints!(klines,kmin,kmax,"G","K")

kmin = [2π/sqrt(3),0]
kmax = [2π/sqrt(3),2π/3]
add_Kpoints!(klines,kmin,kmax,"K","M")

kmin = [2π/sqrt(3),2π/3]
kmax = [0,0]
add_Kpoints!(klines,kmin,kmax,"M","G")
calc_band_plot(klines,la)
```

![68747470733a2f2f71696974612d696d6167652d73746f72652e73332e616d617a6f6e6177732e636f6d2f302f3234363131332f32616530653833392d633239642d333166332d336533332d3136343164376431636230382e706e67](https://user-images.githubusercontent.com/21115243/46902092-f22cc900-cef9-11e8-85be-948a0e7d3dae.png)


### Graphene nano-ribbon

```julia
#We have already constructed atoms and hoppings.
#We add the line to plot
klines = set_Klines()
kmin = [-π]
kmax = [π]
add_Kpoints!(klines,kmin,kmax,"-pi","pi")
```

```julia
#We consider the periodic boundary condition along the primitive vector
direction = 1
#Periodic boundary condition
calc_band_plot_finite(klines,la,direction,periodic=true)
```

![68747470733a2f2f71696974612d696d6167652d73746f72652e73332e616d617a6f6e6177732e636f6d2f302f3234363131332f66323033656365632d393835322d303931612d336332382d3662633463386138356666312e706e67](https://user-images.githubusercontent.com/21115243/46902101-315b1a00-cefa-11e8-9e41-01ffe464a3b1.png)



```julia
#We introduce the surface perpendicular to the premitive vector
direction = 1
#Open boundary condition
calc_band_plot_finite(klines,la,direction,periodic=false)
```
![68747470733a2f2f71696974612d696d6167652d73746f72652e73332e616d617a6f6e6177732e636f6d2f302f3234363131332f36313038363162632d316538302d343364632d303064322d3035643237663865383435652e706e67](https://user-images.githubusercontent.com/21115243/46902102-34eea100-cefa-11e8-8abf-9216a3163ac4.png)

This software is released under the MIT License, see LICENSE.
