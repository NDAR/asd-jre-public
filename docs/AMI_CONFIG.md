


## AMI Image host setup
Starting with this image:
`ami-765b3e1f us-east-1 starcluster-base-ubuntu-12.04-x86_64`

```shell
BASE_AMI="ami-765b3e1f"
NODE_TYPE="m1.xlarge"
starcluster start -o -s 1 -b 0.50 -i $NODE_TYPE -n $BASE_AMI imagehost

starcluster listclusters --show-ssh-status imagehost

starcluster sshmaster imagehost
```


## AMI Configuration


### Python/iPython versions

```shell

$ python --version
Python 2.7.3

$ ipython --version
0.13.1
```

### Python packages


```shell
easy_install -U distribute

pip install ruffus
pip install flask
pip install flask-restful
pip install redis
pip install requests
pip install pymongo
pip install mongokit
pip install rq
pip install fabric
pip install pysam
pip install -U pandas
# see below for pytables installation
```
Testing python package versions

```python
packages = ["numpy", "scipy", "pandas", "ruffus", "flask", 
            "socket", "redis", "pymongo", "mongokit", "rq",
            "pysam", "requests"]
imported = {}
versions = {}

for p in packages:
    try:
        imported[p] = __import__(p)
    except ImportError:
        print "%s\tERROR: could not be imported" % p
        continue
    try:
        versions[p] = imported[p].__version__
    except AttributeError:
        versions[p] = "No version found"
    print "%s\t%s" % (p, versions[p])
```
Output:

```
OpenBLAS : Your OS does not support AVX instructions. OpenBLAS is using Nehalem kernels as a fallback, which may give poorer performance.
numpy   1.7.1
scipy   0.9.0
pandas  0.12.0
ruffus  2.2
flask   0.10.1
socket  No version found
redis   2.8.0
pymongo No version found
mongokit    0.9.0
rq  0.3.11
pysam   0.7.6
requests    2.0.0
```

### Upgrading to java 7
```
apt-get update
apt-get install -q -y openjdk-7-jre
mv /usr/bin/java /usr/bin/java6
ln -s /usr/lib/jvm/java-7-openjdk-amd64/jre/bin/java /usr/bin/java7
ln -s /usr/lib/jvm/java-7-openjdk-amd64/jre/bin/java /usr/bin/java
```

### Install s3cmd
```
mkdir /usr/bin/s3cmd
cd /usr/bin/s3cmd
wget --no-check-certificate https://github.com/s3tools/s3cmd/archive/v1.5.0-alpha3.tar.gz
tar -xzf v1.5.0-alpha3.tar.gz
cd /usr/bin/s3cmd/s3cmd-1.5.0-alpha3/
python setup.py install
cd ~
s3cmd --configure
```
**Note** that this configuration file is saved to '/root/.s3cfg' (contains ACCESS KEYS!)

### Databases
#### Redis
```
mkdir /usr/bin/redis
cd /usr/bin/redis
wget http://redis.googlecode.com/files/redis-2.6.14.tar.gz
tar xzf redis-2.6.14.tar.gz
cd redis-2.6.14
make
```
To start Redis run:

```shell
/usr/bin/redis/redis-2.6.14/src/redis-server
```

And to use the CLI:

```shell
/usr/bin/redis/redis-2.6.14/src/redis-cli
```

#### Mongodb
```shell
sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 7F0CEB10
echo 'deb http://downloads-distro.mongodb.org/repo/ubuntu-upstart dist 10gen' | sudo tee /etc/apt/sources.list.d/10gen.list
sudo apt-get update
sudo apt-get -q -y install mongodb-10gen
sudo service mongodb stop
```

Must configure file at `/etc/mongodb.conf` to point to db data directory etc. Default configuration set to /data mount point. **Must** set bind_ip **and** port as well to keep only localhost connections!
```ini
dbpath = /data/db/mongodb/
logpath = /data/log/mongodb/mongodb.log

port = 27017
bind_ip = 127.0.0.1
```

And then, the folders must be set to the mongodb/mongodb user/group 
```shell
chown -R mongodb:mongodb /data/db/mongodb/
chown -R mongodb:mongodb /data/log/mongodb/
```

To start and stop the service:
```shell
sudo service mongodb start
sudo service mongodb stop
sudo service mongodb restart
```

### Genomics specific software

#### BWA
```shell
mkdir /usr/bin/bwa
curl -L http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.5a.tar.bz2/download -o /usr/bin/bwa/bwa-0.7.5a.tar.bz2
cd /usr/bin/bwa
tar -jxf bwa-0.7.5a.tar.bz2
cd bwa-0.7.5a
make
ln -s /usr/bin/bwa/bwa-0.7.5a/bwa /usr/local/bin/bwa
```

#### GATK (saved nightly on S3)
```shell
mkdir /usr/bin/gatk
cd /usr/bin/gatk
wget https://s3.amazonaws.com/asdjre/GenomeAnalysisTK-2.7-2.tar.bz2
tar -jxf GenomeAnalysisTK-2.7-2.tar.bz2
ln -s /usr/bin/gatk/GenomeAnalysisTK-2.7-2-g6bda569 /usr/local/bin/gatk
# check version
java -d64 -Xmx10g -jar /usr/local/bin/gatk/GenomeAnalysisTK.jar -version
```

#### Picard
```shell
mkdir /usr/bin/picard
curl -L http://sourceforge.net/projects/picard/files/latest/download -o /usr/bin/picard/picard-tools-latest.zip
cd /usr/bin/picard/
unzip picard-tools-latest.zip
ln -s /usr/bin/picard/picard-tools-1.100 /usr/local/bin/picard
```

#### Samtools
```shell
apt-get install -q -y samtools
```

#### Bedtools
```shell
mkdir /usr/bin/bedtools
curl http://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz -o /usr/bin/bedtools/BEDTools.v2.17.0.tar.gz
cd /usr/bin/bedtools
tar -zxvf BEDTools.*.tar.gz
cd bedtools-2.17.0/
make
cd /usr/local/bin
ln -s /usr/bin/bedtools/bedtools-2.17.0/bin/* .
```

#### vcftools
```shell
mkdir /usr/bin/vcftools
curl -L http://sourceforge.net/projects/vcftools/files/vcftools_0.1.11.tar.gz/download -o /usr/bin/vcftools/vcftools_0.1.11.tar.gz
cd /usr/bin/vcftools
tar xzf vcftools_0.1.11.tar.gz
cd vcftools_0.1.11/
make
ln -s /usr/bin/vcftools/vcftools_0.1.11/bin /usr/local/bin/vcftools
```

#### qplot (statistics on bam files)
Note that this downloads a ton of files. Much better to just use executable file and reference files in S3.
```shell
mkdir /usr/local/bin/qplot
wget http://www.sph.umich.edu/csg/zhanxw/software/qplot/qplot.20130627.tar.gz
```

#### Freebayes
```
git clone --recursive git://github.com/ekg/freebayes.git
cd freebayes/
make
```

### Other software

#### R
```shell
apt-get install -q -y r-base-core
```
#### R packages
```shell
Rscript --vanilla -e 'install.packages("ggplot2", repos="http://cran.r-project.org")'
```

#### PyTables
```shell
sudo apt-get install -q -y liblzo2-dev
sudo apt-get install -q -y libhdf5-serial-dev
sudo pip install numexpr
sudo pip install -e git+https://github.com/PyTables/PyTables.git@v.2.4.0#egg=tables
```


#### Perl stuff for Sanders' QC fingerprinting script
```shell
# perl package manager
apt-get install cpanminus
cpanm Parallel::ForkManager
# install and build samtools with special flags
cd ~
wget http://downloads.sourceforge.net/project/samtools/samtools/0.1.19/samtools-0.1.19.tar.bz2
tar -xjf samtools-0.1.19.tar.bz2
cd samtools-0.1.19
make clean
make CXXFLAGS=-fPIC CFLAGS=-fPIC CPPFLAGS=-fPIC

#  Run on all nodes
apt-get install -q -y cpanminus
cpanm Parallel::ForkManager
export SAMTOOLS=/data/bin/samtools-0.1.19
cpanm Bio::DB::Sam
```

### Monitoring and server software

#### collectd
```shell
mkdir /usr/bin/collectd
cd /usr/bin/collectd
wget http://collectd.org/files/collectd-5.4.0.tar.bz2
tar jxf collectd-*.tar.bz2
cd collectd-*
./configure
make all install
```

#### iotop
```shell
mkdir /usr/bin/iotop
cd /usr/bin/iotop
wget http://guichaz.free.fr/iotop/files/iotop-0.6.tar.gz .
tar -xzf iotop-0.6.tar.gz
cd iotop-0.6
python setup.py install
```

## Save AMI to S3

### Remove ~/.s3cfg file with AWS keys!
```shell
rm ~/.s3cfg
```

### Use StarCluster s3image to save
*Note*: currently requires use of development branch-- fix to be included in 0.94.1
```shell
INSTANCE="i-6a51df17"
AMI_NAME="asdjre_base"
starcluster s3image $INSTANCE $AMI_NAME asdjre_ami
``` 


## RAID

RAID setup

```shell
apt-get install -q -y mdadm

# remove /mnt from fstab
sed -i '/\/mnt/d' /etc/fstab
# unmount
umount /mnt
# create raid
mdadm --create /dev/md0 --level 0 --raid-devices 2 /dev/xvdb1 /dev/xvdc1 
#wait briefly
sleep 5
mkfs.ext4 /dev/md0 
# use these to immediately init the FS: -E lazy_itable_init=0,lazy_journal_init=0
# remount!
mount /dev/md0 /mnt
```
