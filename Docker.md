# Docker

## Introduction

Unlike a virtual machine, [Docker](https://www.docker.com/) does not include a separate operating system. Rather, it exploits Linux kernel resource isolation (CPU, memory, block I/O, network) and namespaces to allow independent "containers" to run within a single Linux instance.

A Docker **container** is a basic version of Linux.

A Docker **image** is software that is loaded into the container and then run.

When an image is run, Docker checks to see that you have a copy of the image. If not, it is downloaded from [Docker Hub](https://hub.docker.com/). Docker then loads the image into a container and runs it.

A Dockerfile can be used to script the building of images. It can specify:

* Existing Docker images to extend and other environment settings.
* Software that makes up the image.
* Software to copy into the image.
* Commands to run in the container.

Windows users cannot escape the need for virtual machines. [Install Docker for Windows](https://docs.docker.com/windows/step_one/) comments that:

> Because the Docker Engine daemon uses Linux-specific kernel features,
> you can't run Docker Engine natively in Windows. Instead, you must use
> the Docker Machine command, docker-machine, to create and attach to a
> small Linux VM on your machine. This VM hosts Docker Engine for you on
> your Windows system.

Installing Docker for Windows installs Oracle VirtualBox for use as the virtual machine provider.

For more information see:

* [Docker](https://en.wikipedia.org/wiki/Docker_(software)) on Wikipedia.
* [Docker](https://www.docker.com/) web site
* [Docker Docs](https://docs.docker.com/)

---

## Install Docker on Ubuntu

[Get Docker CE for Ubuntu](https://docs.docker.com/install/linux/docker-ce/ubuntu/):

```console
$ sudo apt-get remove docker docker-engine docker.io
$ sudo apt-get update
$ sudo apt-get install -y \
  apt-transport-https \
  ca-certificates \
  curl \
  software-properties-common
$ curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
$ sudo apt-key fingerprint 0EBFCD88
pub   rsa4096 2017-02-22 [SCEA]
      9DC8 5822 9FC7 DD38 854A  E2D8 8D81 803C 0EBF CD88
uid           [ unknown] Docker Release (CE deb) <docker@docker.com>
sub   rsa4096 2017-02-22 [S]
$ sudo add-apt-repository \
  "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) \
  stable"
$ sudo apt-get update
$ sudo apt-get install -y docker-ce
```

Check:

```console
$ sudo docker version
Client: Docker Engine - Community
 Version:           19.03.7
 API version:       1.40
 Go version:        go1.12.17
 Git commit:        7141c199a2
 Built:             Wed Mar  4 01:22:36 2020
 OS/Arch:           linux/amd64
 Experimental:      false

Server: Docker Engine - Community
 Engine:
  Version:          19.03.7
  API version:      1.40 (minimum version 1.12)
  Go version:       go1.12.17
  Git commit:       7141c199a2
  Built:            Wed Mar  4 01:21:08 2020
  OS/Arch:          linux/amd64
  Experimental:     false
 containerd:
  Version:          1.2.13
  GitCommit:        7ad184331fa3e55e52b890ea95e65ba581ae3429
 runc:
  Version:          1.0.0-rc10
  GitCommit:        dc9208a3303feef5b3839f4323d9beb36df0a9dd
 docker-init:
  Version:          0.18.0
  GitCommit:        fec3683
$ sudo docker run hello-world
Unable to find image 'hello-world:latest' locally
latest: Pulling from library/hello-world
1b930d010525: Pull complete 
Digest: sha256:fc6a51919cfeb2e6763f62b6d9e8815acbf7cd2e476ea353743570610737b752
Status: Downloaded newer image for hello-world:latest
...
Hello from Docker!
...
```

Configure Docker to start on boot:

```console
$ sudo systemctl enable docker
$ sudo systemctl status docker
? docker.service - Docker Application Container Engine
   Loaded: loaded (/lib/systemd/system/docker.service; enabled; vendor preset: e
   Active: active (running) since Thu 2020-03-05 05:53:56 PST; 4min 11s ago
     Docs: https://docs.docker.com
 Main PID: 64998 (dockerd)
    Tasks: 13
```

[Post-installation steps for Linux](https://docs.docker.com/install/linux/linux-postinstall/):

```console
$ sudo groupadd docker
$ sudo usermod -aG docker ubuntu # Assuming user is called ubuntu
```

Log out and in again, or restart machine/VM:

```console
$ groups
... docker ...
$ docker run hello-world
...
Hello from Docker!
...
```
