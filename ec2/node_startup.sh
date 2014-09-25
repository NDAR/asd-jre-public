#!/bin/bash

# NODES=`qconf -sel | grep -v "master"`
# for NODE in $NODES; do
# qsub -l h=$NODE /data/asd-jre/node_startup.sh
# done


# install RAID
export DEBIAN_FRONTEND=noninteractive

apt-get -y --force-yes install mdadm

# unmount
# On EC2, /mnt is usually already on
mounted=$(mount | grep /mnt | wc -l)
if [ $mounted -eq 1 ]; then
    umount -l /mnt
    # remove /mnt from fstab
    sed -i '/\/mnt/d' /etc/fstab
fi

# Create a RAID0 volume, don't ask for permission (--run)
mdadm --create /dev/md0 --level 0 --raid-devices=2 --run /dev/xvdb1 /dev/xvdc1 
sleep 5
 
# Wait for volume to be finished (RAID0 does not take long...)
while [[ `sudo mdadm --detail /dev/md0 | grep 'Rebuild Status'` != '' ]]; do
  sleep 10
done

mkfs.ext4 /dev/md0
# remount!
mount /dev/md0 /mnt


# SET UP S3 ACCESS
cd /usr/bin/s3cmd/
tar -xzf v1.5.0-alpha3.tar.gz
cd s3cmd-1.5.0-alpha3
python setup.py install

# copy over all the access files
cp /data/config/.s3cfg* /root/

# install the qplot software
cp /data/bin/qplot /usr/local/bin

# install shutdown script
cp /data/asd-jre/ec2-shutdown.sh /etc/init.d/ec2-shutdown
sudo update-rc.d ec2-shutdown start 10 0 6 .


cp -R /data/asd-jre /var/tmp
cd /var/tmp/asd-jre
python setup.py install

# Install ggplot2 into R
Rscript --vanilla -e 'install.packages("ggplot2", repos="http://cran.r-project.org")'

# Install perl scripts needed for fingerprint
apt-get install -q -y --force-yes cpanminus
cpanm Parallel::ForkManager
export SAMTOOLS=/data/bin/samtools-0.1.19
cpanm Bio::DB::Sam --force


# install gnu parallel
cd /data/bin/parallel-20131222/
make install
