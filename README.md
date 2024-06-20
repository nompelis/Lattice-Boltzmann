# Lattice-Boltzmann learning collab.

Never stop learning...

This is a very simple setup that will constantly be getting updated with
more functionality to bring this to a functional sover that we can benchmark.
Presently, it uses a 3D notional box (Cartesian grid of cells) that has the
so-called "D3Q19" Lattice-Botlzmann velocity discretization. It uses the
following convention:

```
          y |                                        back                     |
            |_______                              |--17---|                   |
           /|      /|               mid           |   |   |                   |
          / |     / |                8---3---7   14---6---13                  |
         /__|____/  |                | \ | / |    |   |   |    |---15---|     |
         |  |___ |_ |______ x        2---0---1    |--18---|    |    |   |     |
         | /     | /                 | / | \ |                12----5---11    |
         |/      |/                 10---4---9                 |    |   |     |
         |_______/                                             |---16---|     |
      z /                                                        front        |
```

______________________________________________________________________________
`` References ``

[1] Alessandro De Rosis & Christophe Coreixas, "Multiphysics flow simulations
    using D3Q19 lattice Boltzmann methods based on central moments," Physics
    of Fluids 32, 117101 (2020)

(To be continued)
