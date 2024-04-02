# mscm

Multi-species count models for Liesel.

For more info on multi-species count models, see Hannes Riebl, Jonas Glatthorn and Thomas Kneib (2023), "A Structured Additive Multi-Species Count Model for Assessing the Relation Between Site Conditions and Species Diversity", in: Hannes Riebl, "Semi-Parametric Distributional Regression in Forestry and Ecology: Software, Models and Applications", chapter C (appendix), pp. 101--124, https://doi.org/10.53846/goediss-10051.

This repository contains the code to replicate the application (`paper/application`) and the simulation study (`paper/simulation-study`) in the manuscript. After following the instructions below, run the Quarto notebooks or the `run.sh` shell scripts for parallelized execution using GNU Parallel.

## Installation

This project uses [conda](https://docs.conda.io/en/latest/miniconda.html) and
[pdm](https://pdm.fming.dev) for dependency management. To install `mscm`, run:

```
conda env create -f environment.yml
conda activate mscm
Rscript -e "remotes::install_github('liesel-devs/rliesel')"

# optional: install the tex gyre heros font
wget "https://www.gust.org.pl/projects/e-foundry/tex-gyre/heros/qhv2.004otf.zip"
mkdir -p ~/.local/share/fonts
unzip qhv2.004otf.zip -d ~/.local/share/fonts
rm qhv2.004otf.zip

pdm install --dev
```

To update the dependencies in the lock file, run:

```
pdm update --dev
```

To install exactly the same versions as in the lock file, run:

```
pdm sync --dev
```

## Development

- `pdm run docs`: Serves the docs with [pdoc](https://pdoc.dev).
- `pdm run lint`: Runs the [pre-commit](https://pre-commit.com) hooks.
- `pdm run test`: Runs [pytest](https://pytest.org).
