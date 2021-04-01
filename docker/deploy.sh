#!/bin/bash

echo "./deploy.sh $*" > redeploy.sh
chmod +x redeploy.sh

existing=$(docker ps -aq -f name=rmsd)
if [ -n "$existing" ]; then
    echo "removing existing container"
    docker rm -f $existing
fi

docker run -d \
--name rmsd \
--restart unless-stopped \
-e ARGS="$*" \
rmsd
