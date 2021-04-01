# Nanome - RMSD

A Nanome Plugin to run an RMSD calculation on multiple structures, and auto-align them.

## Dependencies

[Docker](https://docs.docker.com/get-docker/)

## Usage

To run RMSD in a Docker container:

```sh
$ cd docker
$ ./build.sh
$ ./deploy.sh -a <plugin_server_address> [optional args]
```

## Development

To run RMSD with autoreload:

```sh
$ python3 -m pip install -r requirements.txt
$ python3 run.py -r -a <plugin_server_address> [optional args]
```

## License

MIT
