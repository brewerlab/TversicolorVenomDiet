#Tetragnatha versicolor population-level diet and venom gland analyses

Population-level analysis of both sexes trophic level (stable isotopes), differential expression venom glands, and molecular evolution of venom transcripts.

##Transcriptome Assemblies

###Assemble individual venom gland transcriptomes
~/Oyster_River_Protocol/oyster.mk main MEM=150 CPU=24 READ1=reads/CA421ANXX_s3_1_GSLv3-7_73_SL214002.fastq.gz READ2=reads/CA421ANXX_s3_2_GSLv3-7_73_SL214002.fastq.gz RUNOUT=SL214002
~/Oyster_River_Protocol/oyster.mk main MEM=150 CPU=24 READ1=reads/CA421ANXX_s3_1_GSLv3-7_74_SL214003.fastq.gz READ2=reads/CA421ANXX_s3_2_GSLv3-7_74_SL214003.fastq.gz RUNOUT=SL214003
~/Oyster_River_Protocol/oyster.mk main MEM=150 CPU=24 READ1=reads/CA421ANXX_s3_1_GSLv3-7_75_SL214004.fastq.gz READ2=reads/CA421ANXX_s3_2_GSLv3-7_75_SL214004.fastq.gz RUNOUT=SL214004
~/Oyster_River_Protocol/oyster.mk main MEM=150 CPU=24 READ1=reads/CA421ANXX_s3_1_GSLv3-7_76_SL214005.fastq.gz READ2=reads/CA421ANXX_s3_2_GSLv3-7_76_SL214005.fastq.gz RUNOUT=SL214005
~/Oyster_River_Protocol/oyster.mk main MEM=150 CPU=24 READ1=reads/CA421ANXX_s3_1_GSLv3-7_77_SL214006.fastq.gz READ2=reads/CA421ANXX_s3_2_GSLv3-7_77_SL214006.fastq.gz RUNOUT=SL214006
~/Oyster_River_Protocol/oyster.mk main MEM=150 CPU=24 READ1=reads/CA421ANXX_s3_1_GSLv3-7_78_SL214007.fastq.gz READ2=reads/CA421ANXX_s3_2_GSLv3-7_78_SL214007.fastq.gz RUNOUT=SL214007
~/Oyster_River_Protocol/oyster.mk main MEM=150 CPU=24 READ1=reads/CA421ANXX_s3_1_GSLv3-7_79_SL214008.fastq.gz READ2=reads/CA421ANXX_s3_2_GSLv3-7_79_SL214008.fastq.gz RUNOUT=SL214008
~/Oyster_River_Protocol/oyster.mk main MEM=150 CPU=24 READ1=reads/CA421ANXX_s3_1_GSLv3-7_80_SL214009.fastq.gz READ2=reads/CA421ANXX_s3_2_GSLv3-7_80_SL214009.fastq.gz RUNOUT=SL214009
~/Oyster_River_Protocol/oyster.mk main MEM=150 CPU=24 READ1=reads/CA421ANXX_s3_1_GSLv3-7_81_SL214010.fastq.gz READ2=reads/CA421ANXX_s3_2_GSLv3-7_81_SL214010.fastq.gz RUNOUT=SL214010
~/Oyster_River_Protocol/oyster.mk main MEM=150 CPU=24 READ1=reads/CA421ANXX_s3_1_GSLv3-7_82_SL214011.fastq.gz READ2=reads/CA421ANXX_s3_2_GSLv3-7_82_SL214011.fastq.gz RUNOUT=SL214011
~/Oyster_River_Protocol/oyster.mk main MEM=150 CPU=24 READ1=reads/CA421ANXX_s3_1_GSLv3-7_83_SL214012.fastq.gz READ2=reads/CA421ANXX_s3_2_GSLv3-7_83_SL214012.fastq.gz RUNOUT=SL214012
~/Oyster_River_Protocol/oyster.mk main MEM=150 CPU=24 READ1=reads/CA421ANXX_s3_1_GSLv3-7_84_SL214013.fastq.gz READ2=reads/CA421ANXX_s3_2_GSLv3-7_84_SL214013.fastq.gz RUNOUT=SL214013

###Assemble combined representatives each combination of sex and population

~/Oyster_River_Protocol/oyster.mk main MEM=150 CPU=24 READ1=reads/TversicolorGlands_MF_004-007-008-012_1.fastq.gz READ2=reads/TversicolorGlands_MF_004-007-008-012_2.fastq.gz RUNOUT=TversicolorGlandsMF

##Read pseudomapping and abundance estimates

####Map reads from each individual to the combined assembly

```
salmon index -t ../ORP_assemblies/assemblies/TversicolorGlandsMF.orthomerged.fasta -i TversicolorGlandsMF.orthomerged_salmon_index
for fn in ../ORP_assemblies/rcorr/SL2140{02..13}.TRIM; do samp=`basename ${fn}`; echo "Processing sample ${samp}"; salmon quant --seqBias --gcBias -i TversicolorGlandsMF.orthomerged_salmon_index -l A -1 ${fn}_1P.cor.fq -2 ${fn}_2P.cor.fq -p 8 -o quants/${samp}_quant; done
```

##Annotate combined assembly

###Install DAMMIT and dependencies

```
sudo apt-get update

sudo apt-get -y install python3-dev hmmer unzip \
    infernal ncbi-blast+ liburi-escape-xs-perl emboss liburi-perl \
    build-essential libsm6 libxrender1 libfontconfig1 \
    parallel libx11-dev python3-venv last-align transdecoder
    
python3.5 -m venv ~/py3
. ~/py3/bin/activate
pip install -U pip

pip install -r <(curl https://raw.githubusercontent.com/camillescott/shmlast/master/requirements.txt)
pip install --upgrade pip
pip install shmlast

cd
(wget -O - pi.dk/3 || curl pi.dk/3/ || fetch -o - http://pi.dk/3) | bash
sudo cp $HOME/bin/parallel /usr/bin/parallel

cd
git clone https://gitlab.com/ezlab/busco.git
pushd busco && python setup.py install && popd

export PATH=$HOME/busco/scripts:$PATH
echo 'export PATH=$HOME/busco/scripts:$PATH' >> $HOME/.bashrc

pip install https://github.com/camillescott/dammit/archive/refactor/1.0.zip
```

###Install DAMMIT databases and annotate combined assembly

```
mkdir /media/brewerlab/RAID/dammit_dbs
#cd /media/brewerlab/RAID/dammit_dbs
dammit databases --install --database-dir /media/brewerlab/RAID/dammit_dbs --full --busco-group arthropoda
dammit annotate TversicolorGlandsMF.orthomerged.fasta --busco-group arthropoda --n_threads 12 --database-dir /media/brewerlab/RAID/dammit_dbs --full
```
