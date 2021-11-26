# CV19GM Python Library - Dlab

## Install Dependencies

Clone this repository in your computer

Open a terminal and go to repository's directory.

Execute:

`python setup.py install`



### Pygmo:

​	 https://esa.github.io/pygmo2/install.html


## Running the app from a Docker container

1. Build docker image: `docker build -t cv19gm:0.3 .`
2. Run docker container (command for Linux based distributions): `docker run -it --rm --name cv19gm -e DISPLAY -v "$HOME/.Xauthority:/root/.Xauthority" --net=host cv19gm:0.3 /bin/bash`


# Repo Structure
```
Home/  
├── Examples/  
├── Tests/  
├── Documents/  
├── Data/ (En duda)  
├── SRC/  
   ├── SEIR/  
   ├── SEIRHVD/  
   ├── SEIRStar/  
   ├── SEIRHDVStar/  
   ├── utils/  
      ├── plots.py  
      ├── utils.py         
Readme.md  
Install.sh  
requirements.txt  
Docker  
Licence   
etc.  
```

