## Create and run the docker image for modrunner

#### Build image
```
docker build -t modrunnerLite .  
```

### Note you probably need a fair few resources to run this
I get errors for my data if I don't increase the memory allocated to docker.

https://docs.docker.com/config/containers/resource_constraints/#:~:text=By%20default%2C%20if%20an%20out,the%20%2Dm%2F%2D%2Dmemory%20option.

```
docker run -v /Users/X/test_data/mouse:/data/mouse modrunner minimap2 --c map --f data/mouse/mouse_docker.csv --genome
```
# Convert geneome to transcriptome
```
 docker run -v /Users/ariane/Documents/code/TomboRunner/test_data/mouse:/data/mouse modrunner minimap2 --c gents  --f data/mouse/mouse_docker.csv 
 ```
# Map reads to reference
docker run -v /Users/ariane/Documents/code/TomboRunner/test_data/mouse:/data/mouse modrunner minimap2 --c map --f data/mouse/mouse_docker.csv

### Eligos2
docker run -v /Users/ariane/Documents/code/TomboRunner/test_data/mouse:/data/mouse modrunner minimap2 --c tobed --f data/mouse/mouse_docker.csv 
docker run -v /Users/ariane/Documents/code/TomboRunner/test_data/mouse:/data/mouse modrunner eligos2 --c model --f data/mouse/mouse_docker.csv --sample1 wt

### Run via CLI 
This is the normal and recommended use, for running things on a HPC or server.

#### Step 1:
Basecall: 

Quality control step
#### Step 2:


#### Run jupyter to debug
```
docker run -it -p -m "8g" 8083:8083 modrunner 
```

Then once on the container, go:
```
cd examples
jupyter notebook --ip 0.0.0.0 --port 8083 --no-browser --allow-root
```
Now you can access it through the local browser on http://localhost:8888 :)

#### Run with data and notebooks you have hosted locally

```
docker run -it -p 8083:8083 --memory="8g" -v {full path to where your data is located}:/examples modrunner 
```

For example, for me the command looks like:
```
docker run -it -p 8083:8083 --memory="8g"  -v /Users/ariane/Documents/code/TomboRunner/examples:/examples modrunner 
```

#### Save built image
```
docker save modrunner | gzip > modrunner.tar.gz
```

#### Load saved image
```
docker load < modrunner.tar.gz
```

#### Remove image
```
docker rmi modrunner
```

#### Remove any running containers
```
docker ps -a
```
Remove the desired containers via their ID with `docker rm ID`


#### Convert to singularity 

Given we run our code on Tinnaroo, and docker isn't allowed, we have to wrap the docker image in a singularity container...

See below.
```
docker build . -t "modrunner:latest"
```

```
mkdir -p /tmp/test

```

```
docker run -v /var/run/docker.sock:/var/run/docker.sock \
  -v /Users/ariane/Documents/code/TomboRunner/singularity:/output \
  --privileged -t --rm \
  quay.io/singularity/docker2singularity \
  modrunner:latest
  ```

https://github.com/singularityhub/docker2singularity
```
Tips for making Docker images compatible with Singularity
Define all environmental variables using the ENV instruction set. Do not rely on .bashrc, .profile, etc.
Define an ENTRYPOINT instruction set pointing to the command line interface to your pipeline
Do not define CMD - rely only on ENTRYPOINT
You can interactively test the software inside the container by overriding the ENTRYPOINT docker run -i -t --entrypoint /bin/bash bids/example
Do not rely on being able to write anywhere other than the home folder and /scratch. Make sure your container runs with the --read-only --tmpfs /run --tmpfs /tmp parameters (this emulates the read-only behavior of Singularity)
Don’t rely on having elevated user permissions
Don’t use the USER instruction set
```