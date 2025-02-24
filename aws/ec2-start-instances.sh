#!/bin/sh
## Given an ec2 instance id, start the instance and print the public IP address

if [ $# -ne 1 ]; then
	echo "Usage: $0 <instance-id>"
	exit 1
fi

aws ec2 start-instances --instance-ids $1
