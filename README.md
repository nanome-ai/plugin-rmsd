# Nanome - RMSD Plugin

This python plugin allows for users to select multiple molecules, run an RMSD calculation, and auto-align the structures.

Plugin supported on Windows, Linux, and Mac for Nanome v1.14

### Installation

```sh
$ pip install nanome-rmsd
```

### Usage

To start the plugin:

```sh
$ nanome-rmsd -a <plugin_server_address>
```

### Docker Usage

To run in a Docker container:

```sh
$ cd docker
$ ./build.sh
$ ./deploy.sh -a <plugin_server_address>
```

### License

MIT
