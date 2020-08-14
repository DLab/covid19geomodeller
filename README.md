# COVID SEIR Model - Dlab

## Install Dependencies

Clone this repository in your computer

Open a terminal and go to repository's directory.

Execute:

`chmod +x install.sh`

`./install.sh`



### For a more robust simulation install scikit-odes (Advanced)

#### Scikits-odes:

https://scikits-odes.readthedocs.io/en/latest/installation.html





## For parameter optimization (beta):

### Pygmo:

​	 https://esa.github.io/pygmo2/install.html


## Running the app from a Docker container

### Full version (with scikits.odes and pygmo)

1. Build docker image: `docker build -t cv19gm:0.3 .`
2. Run docker container (command for Linux based distributions): `docker run -it --rm --name cv19gm -e DISPLAY -v "$HOME/.Xauthority:/root/.Xauthority" --net=host cv19gm:0.3 /bin/bash`

### Lite version (without scikits.odes and pygmo)

1. Build docker image: `docker build -t cv19gm-lite:0.1 -f Dockerfile.lite .`
2. Run docker container (command for Linux based distributions): `docker run -it --rm --name cv19gm -e DISPLAY -v "$HOME/.Xauthority:/root/.Xauthority" --net=host cv19gm-lite:0.1 /bin/bash`


# Repo Structure
Home/
├── Examples/
├── Tests/
├── Documents/
├── Data/ (En duda)
├── SRC/
   ├── SEIR/
   ├── SEIRHVD/
   ├── SEIRStar/
   ├── SEIR/
   ├── utils/
      ├── plots.py
      ├── utils.py       
Readme.md
Install.sh
requirements.txt
Docker
Licence 
etc

 
