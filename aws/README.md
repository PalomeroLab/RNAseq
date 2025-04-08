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

## Copy data from an S3 bucket to another account

[Source.](https://docs.aws.amazon.com/prescriptive-guidance/latest/patterns/copy-data-from-an-s3-bucket-to-another-account-and-region-by-using-the-aws-cli.html)

First check permissions on the source bucket.

```sh
aws sts assume-role \
 --role-arn "arn:aws:iam::925204023671:user/ziwei-collaborator-ryan:role/S3MigrationRole" \
 --role-session-name AWSCLI-Session

aws iam list-attached-user-policies --user-name ziwei-collaborator-ryan --profile aml

aws s3 cp s3://source-bucket-name/path/to/file s3://destination-bucket-name/path/to/file --profile source_profile --source-region source_region --profile destination_profile --destination-region destination_region

aws s3 sync s3://labaf-missionbio-tapestri/AML_singlecell_atlas_JB/

aws s3 cp s3://lab-aaf-ngs-data-archive/RNAseq/10X_Chromium/AML_singlecell_atlas_JB/README.txt \
s3://ziwei-collaborator-upload --profile aml
```
