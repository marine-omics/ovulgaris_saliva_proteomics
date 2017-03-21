# Proteogenomics Analysis of O. vulgaris salivary glands and saliva

To run the analysis you must first ensure that the required hpc data exists as follows

`../msdata`
`../trinity`

Once these are in place you can run `01_link_files.sh` to create symbolic links to required input files in this directory

```bash
	./01_link_files.sh
```

With input files in place it is now possible to run the `rake` script which automates all of the proteogenomics steps including

	- File format conversion (with msconvert)
	- 6-frame and decoy database construction
	- Tandem MS database search with `xtandem` and `msgfplus`
	- Peptide and Protein inference with `PeptideProphet`, `iProphet` and `ProteinProphet`
	- Conversion of `ProteinProphet` results to gff

In order for the rake script to run successfully you will need to install some dependencies

The TPP

```bash
	brew install trans_proteomic_pipeline
```

A working ruby and rubygems installation.  I recommend using rvm to manage this

```bash
	gpg --keyserver hkp://keys.gnupg.net --recv-keys 409B6B1796C275462A1703113804BB82D39DC0E3
	\curl -sSL https://get.rvm.io | bash -s stable --ruby
```

The `protk` ruby gem

```bash
	brew install libxml2
	gem install protk -- --with-xml2-config=/usr/local/Cellar/libxml2/2.9.4_2/bin/xml2-config
```

The MSGFPlus search engine

```bash
	mkdir ~/bin
	cd ~/bin
	wget http://proteomics.ucsd.edu/Software/MSGFPlus/MSGFPlus.20140630.zip
	unzip MSGFPlus.20140630.zip
	chmod u+x MSGFPlus.jar
	echo "PATH=${PATH}:${HOME}/bin" >> ~/.bash_profile
```

To run the script simply type `rake`

```bash
	rake
```

This script will take a long time to run (days) and will result in two useful `gff3` files.  Both files contain coordinates of peptides in reference to the Trinity.fasta transcriptome.  They can be used to visualise the locations of peptides against transcripts and also to find genuine novel coding sequences not otherwise predicted by trinotate.

	- `novel_pg.gff3` contains peptides identified from 6-frame translations
	- `known_pg.gff3` contains peptides that match to coding sequences predicted by transdecoder as part of trinotate
