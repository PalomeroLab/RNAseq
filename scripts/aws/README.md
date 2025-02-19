# AWS scripts

## EC2

Query EC2 instances with the script

```sh
./ec2-list-instances.sh
```

 Profiles

You can work with two accounts by creating two profiles on the aws command line.
It will prompt you for your AWS Access Key ID, AWS Secret Access Key and desired
region, so have them ready.

> [Source](https://stackoverflow.com/a/34246053)

```sh
aws configure --profile palomerolab
aws configure --profile ziweicollaborator
```

Pass the `--profile` flag to the `aws` command to use the desired profile.

```sh
aws s3 ls --profile ziweicollaborator ziwei-collaborator-upload
touch grass && aws s3 cp grass s3://ziwei-collaborator-upload --profile ziweicollaborator
```

