# Jets + Pythia

[![Github Codespace](https://img.shields.io/badge/open-GH_Codespaces-blue?logo=github)](https://codespaces.new/aprozo/pythia-jets-simple-tutorial?quickstart=1)

This is a self-contained tutorial for simple generating Pythia event, running [FastJet](https://fastjet.fr/) over it and analyzing jets

---

## How to start:

For running on Github Codespaces submit an application for Free [Github Education](https://github.com/education) benefits.
Then click on [Github Codespace button](https://codespaces.new/aprozo/pythia-jets-simple-tutorial?quickstart=1)
to start a container (predefined software environment).

## Running Steps

Format for executable arguments : `pTHatMin pTHatMax|inf [nEvents=50000]`
For example, generating events with 10 < ptHat< 20 Gev/c for 1000 events would be:
```bash
make
./makeTree 10 20 1000
```

Parameters can be tuned in `makeTree.cc`

```cpp
  const double jetRadius = 0.4;
  const double jetPtMin = 3.0;
  // particle parameters
  const double particlePtMin = 0.15;
  const double particleEtaMax = 1.5;
```

After that run for analysis of jet tree which will fill some basic histograms based on that file:

```bash
root -l -b -q anaTree.cpp+
```

## Remark: Running on your own laptop

In case you want to enter and run on your own laptop:

- You need to install either [Docker engine](https://docs.docker.com/get-started/get-docker/) on MacOS or [Apptainer (singularity)](https://apptainer.org/docs/admin/main/installation.html)

  For simplier Apptainer (singularity) installation on Linux:

```bash
sudo apt update
sudo apt install -y software-properties-common
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
```

- And then run commands:

```bash
git clone https://github.com/aprozo/pythia-jets-simple-tutorial.git
cd pythia-jets-simple-tutorial
apptainer pull rivet-pythia.sif docker://hepstore/rivet-pythia:main
apptainer exec rivet-pythia.sif bash
```
