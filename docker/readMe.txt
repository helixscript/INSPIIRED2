# Run one dir level above INSPIIRED code base

docker build -t inspiired2:latest .
docker save inspiired2:latest | gzip > inspiired2.tar.gz
