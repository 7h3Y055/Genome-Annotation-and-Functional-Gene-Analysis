all: setup_conda setup_env setup_db wheat_model

setup_conda:
	@if ! command -v conda > /dev/null 2>&1; then \
		echo "[+] Installing Miniconda..."; \
		wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh; \
		bash miniconda.sh -b -p $$HOME/miniconda; \
		rm miniconda.sh; \
		export PATH="$$HOME/miniconda/bin:$$PATH"; \
		$$HOME/miniconda/bin/conda init; \
		echo "[+] Miniconda installed."; \
	else \
		echo "[*] Conda already installed."; \
	fi

setup_env:
	@echo "[+] Creating and setting up bioenv Conda environment..."
	@export PATH="$$HOME/miniconda/bin:$$PATH"; \
	conda create -n bioenv -y python=3.10; \
	conda install -n bioenv -y -c bioconda augustus gffread blast biopython matplotlib numpy wget

setup_db:
	@echo "[+] Downloading UniProt DB and creating BLAST DB..."
	@conda run -n bioenv wget -N https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
	@conda run -n bioenv gunzip -f uniprot_sprot.fasta.gz
	@conda run -n bioenv makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot_sprot_db

wheat_model:
	@echo "[+] Setting up wheat model for Augustus..."
	@AUG_CFG=$$(conda run -n bioenv augustus --print-cfg=CFG); \
	SP_DIR=$$AUG_CFG/species; \
	if [ ! -d "$$SP_DIR/wheat" ]; then \
		git clone --depth 1 https://github.com/Gaius-Augustus/Augustus.git tmp_aug_repo; \
		cp -r tmp_aug_repo/config/species/wheat $$SP_DIR/wheat; \
		rm -rf tmp_aug_repo; \
		echo "[+] Wheat model installed."; \
	else \
		echo "[*] Wheat model already present."; \
	fi

clean:
	rm -f uniprot_sprot.fasta*
	rm -rf uniprot_sprot_db tmp_aug_repo


# source /opt/anaconda/bin/activate