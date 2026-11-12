#rule minimap2_index:

#rule ngmlr_index:


rule get_vep_cache:
	output:
		directory("resources/vep/cache"),
	params:
		species="homo_sapiens_merged",
		build="GRCh38",
		release="113",
	log:
		"logs/vep_cache.log",
	cache: "omit-software"  # save space and time with between workflow caching (see docs)
	wrapper:
		"v3.3.5/bio/vep/cache"

rule download_vep_plugins:
	output:
		directory("resources/vep/plugins")
	params:
		release=113
	wrapper:
		"v3.3.5/bio/vep/plugins"

#rule install_jaffal:
#    output: 
#        ""
#        """
#        mkdir -p repo/ && \
#        cd resources/ && \
#        wget -c https://github.com/Oshlack/JAFFA/releases/download/version-2.3/JAFFA-version-2.3.tar.gz -O resources/JAFFA-version-2.3.tar.gz && \
#        cd resources/ && \
#        tar -zxvf JAFFA-version-2.3.tar.gz && \
#        cd JAFFA-version-2.3 && \
##        wget -v https://figshare.com/ndownloader/files/25410494 -O JAFFA_REFERENCE_FILES_HG38_GENCODE22.V2.tar.gz && \
#        tar -zxvf JAFFA_REFERENCE_FILES_HG38_GENCODE22.V2.tar.gz && \
#        mv JAFFA_REFERENCE_FILES_HG38_GENCODE22.V2.tar.gz 
#        """

#rule install_MetaRNN:
#    output:
#        """
#
#        """

rule install_hmmcopy_utils:
	output:
		"resources/hmmcopy_utils/bin/readCounter"
	threads:
		1
	params:
		"https://github.com/shahcompbio/hmmcopy_utils.git"
	shell:
		"""
		cd {input} && \
		git clone {params} && \
		cd hmmcopy_utils/ && \
		cmake . && \
		make && \
		cd -
		"""

rule install_ichorCNA:
	output:
		"resources/ichorCNA-0.2.0/scripts/runIchorCNA.R"
	threads:
		1
	shell:
		"""
		wget https://github.com/broadinstitute/ichorCNA/archive/refs/tags/v0.2.0.tar.gz -O resources/ichorCNA.tar.gz && \
		tar -xvzf resources/ichorCNA.tar.gz
		"""