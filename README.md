# Source Code for: It Is Rocket Science Part III: How to get to the Moon

> Date: April 2026  
> Author: Louis Finegan  
> Email: louis.m.finegan@gmail.com  
> Website: [louis-finegan.github.io](https://louis-finegan.github.io)

This repository contains the source code used to generate the animations of the orbital manoeuvre simulations featured in the blog post:

>[`It Is Rocket Science Part III: How to get to the Moon`](https://louis-finegan.github.io/2026/04/13/Orbital-Manoeuvres.html)

The blog post includes animations of:

- A Hohmann transfer orbit 

<video width="600" controls>
  <source src="figures/satellite_orbit.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>

- A Trans-Lunar Injection

<video width="600" controls>
  <source src="figures/lunar_transfer_animation.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>

## Requirments

- Python 3.10+ recommended
- Dependencies listed in `requirements.txt`

Install dependencies with:

```bash
pip install -r requirements.txt
```

## Usage

This project generates orbital simulation data using C code located in `src/`, which is then visualised using Python (`matplotlib`) to produce animated trajectory plots.

The code was developed using gcc 11.4.0 on Ubuntu 22.04.2.

### Hohmann transfer orbit: run

```bash
gcc src/TwoBody.c -o build/sim_TB -lm
./build/sim_TB
python Orbitplot.py
```

### Trans-Lunar Injection: run

```bash
gcc src/TransLunarInjection.c -o build/sim_TLI -lm
./build/sim_TLI
python TLIplot.py
```

## License

This project is released under the MIT License. 