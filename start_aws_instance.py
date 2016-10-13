import boto3

ec2 = boto3.resource('ec2')

rc = ec2.create_instances(
    ImageId = 'ami-17178000',
    MinCount = 1,
    MaxCount = 1,
    KeyName = 'nb_mit',
    InstanceType = 'c4.4xlarge',
    NetworkInterfaces=[{ "DeviceIndex": 0, 'SubnetId': 'subnet-61c9c716', 'Groups' : ['sg-9d7fe4e5'], "AssociatePublicIpAddress": True }] 
)

for instance in rc:
   instance.wait_until_running()
   instance.load()
   print(instance.public_ip_address)
  
