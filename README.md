# TightBinding.jl
This can construct the tight-binding model and calculate energies. 
This is not completed, yet.

This can 
1. construct the Hamiltonian as a functional of a momentum k.
2. plot the band structure.
3. show the crystal structure.
4. plot the band structure of the finite-width system with one surface or boundary.

There is the sample jupyter notebook. 

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
show_neighbers(la)
```

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

```julia
#show the band structure
        klines = set_Klines()
        kmin = [0,0]
        kmax = [2π/sqrt(3),0]
        add_Kpoints!(klines,kmin,kmax,"G","K")

        kmin = [2π/sqrt(3),0]
        kmax = [2π/sqrt(3),2π/3]
        add_Kpoints!(klines,kmin,kmax,"K","M")

        kmin = [2π/sqrt(3),2π/3]
        kmax = [0,0]
        add_Kpoints!(klines,kmin,kmax,"M","G")
        calc_band_plot(klines,la)
```

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

```julia
#We introduce the surface perpendicular to the premitive vector
direction = 1
#Open boundary condition
calc_band_plot_finite(klines,la,direction,periodic=false)
```


This software is released under the MIT License, see LICENSE.
