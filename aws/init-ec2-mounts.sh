#!/bin/bash
## Mounts ec2 instances using sshfs
## Requires:
##	sshfs

# Set the mount directory
EC2_MOUNTDIR="${HOME}/EC2"

# set the instance name
EC2_INSTANCE_NAME=cellranger

mkdir -p "${EC2_MOUNTDIR}/${EC2_INSTANCE_NAME}"

# undo any previous mounts
# umount "${EC2_MOUNTDIR}/${EC2_INSTANCE_NAME}"
sshfs cellranger:/home/ubuntu "${EC2_MOUNTDIR}/${EC2_INSTANCE_NAME}"

