#!/bin/sh
## Lists amazon EC2 instance and displays 3 columns; name, id; and ip (if available)
aws ec2 describe-instances --query 'Reservations[*].Instances[*].[Tags[?Key==`Name`].Value | [0], InstanceId, PublicIpAddress]' --output text
