
### Download the INSPIIRED Docker image.
wget https://bushmanlab.org/export/inspiired2_latest.tar.gz

### Load the image.
docker load -i inspiired2_latest.tar.gz

### Place your sequencing data and processing script (run.sh) in a directory, eg. /home/everett/workspace, and include the directory in the Docker call.
docker run --rm --shm-size=30g -v /home/everett/workspace:/workspace -w /workspace inspiired2 bash run.sh
