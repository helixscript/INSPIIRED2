
#### Download the INSPIIRED Docker image.
```wget https://bushmanlab.org/export/inspiired2_latest.tar.gz```

#### Load the image.
```docker load -i inspiired2_latest.tar.gz```

#### Prepare your run.
Place your sequencing data and processing script (run.sh) in a directory, eg. ~/workspace, and include the directory in the Docker call.  
  
```docker run --rm --shm-size=30g -v ~/workspace:/workspace -w /workspace inspiired2 bash run.sh```

#### Processing with databasing.  

To include databasing features in your processing, additional parameters need to be passed to Docker.
Databasing is only available for the buildFragments module and is implimented by including the --dbConfigFile --dbConfigID flags.   
  
```docker run --rm --shm-size=30g --user $(id -u):$(id -g) -v /media/md0/data/inspiired:/data -v ~/workspace:/workspace -w /workspace inspiired2 bash run.sh```
