# Phicore


### Installation 
**Pre-requisities**
- anaconda or miniconda already installed 
- python3

**Setting up the environment**

    ```
    conda create -n phicore 
    conda activate phicore
    #install trnascan-se
    conda install -c bioconda trnascan-se
    #install phanotate dev version
    git clone -b develop https://github.com/deprekate/PHANOTATE.git
    cd PHANOTATE/
    python3 -m venv env
    source env/bin/activate
    python setup.py install
    pip3 install numpy
    ```

**Installing phicore**

    ```
    git clone https://github.com/npbhavya/Phicore.git
    ```

**Deactivating environment**

    ```
    #deactivate the virtual environment 
    deactivate 
    #deactivate the conda environment 
    conda deactivate
    ```

**Activating environment if already installed**

    ```
    conda activate phicore
    cd PHANOTATE
    source env/bin/activate
    ```

