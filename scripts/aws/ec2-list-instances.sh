#!/bin/bash
## Lists amazon EC2 instance and displays 3 columns; name, id; and ip (if available)

aws ec2 describe-instances \
  --query 'Reservations[*].Instances[*].[Tags[?Key==`Name`].Value | [0], InstanceId, PublicDnsName]' \
  --output json \
  | jq -r '.[][] | "\(.[0]//"")\t\(.[1]//"")\t\(.[2]//"")"' \
  | column -t -s $'\t'
