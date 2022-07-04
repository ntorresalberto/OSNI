# OSNI
This repository for the Optimized Spectral Numerical Integration (OSNI) Library.


# compiling inside docker

install [docker](https://docs.docker.com/engine/install/ubuntu/) (and follow the [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/)), then install [dogi](https://github.com/ntorresalberto/dogi).


```
git clone https://github.com/aGotelli/OSNI.git
cd OSNI
dogi run ubuntu       # this spawns a docker container and mounts the repository
dogi                  # this command will confirm if you're inside a container
./installdeps.sh      # containers are barebones by default, you need to install dependencies
```

You can now do whatever you want inside the container without the risk of polluting your host system (like `sudo make install`).

You can also open new terminals inside the same container with `dogi exec`.

# testing the installed library

```
cd OSNI/test
mkdir build
cmake ..
make -j
```
