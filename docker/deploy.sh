if [ "$(docker ps -aq -f name=rmsd)" != "" ]; then
    # cleanup
    echo "removing exited container"
    docker rm -f rmsd
fi

docker run -d \
--name rmsd \
--restart unless-stopped \
-e ARGS="$*" \
rmsd
