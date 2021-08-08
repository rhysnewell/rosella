#!/bin/bash -e

PREFIX=$CONDA_PREFIX
echo $PREFIX

# Build statically linked binary with Rust
C_INCLUDE_PATH=$CONDA_PREFIX/include \
LIBRARY_PATH=$CONDA_PREFIX/lib \
LD_LIBRARY_PATH=$CONDA_PREFIX/lib \
LIBCLANG_PATH=$CONDA_PREFIX/lib/libclang.so \
RUST_BACKTRACE=1 cargo build --release

# Install the binaries
C_INCLUDE_PATH=$PREFIX/include \
LIBRARY_PATH=$PREFIX/lib \
LD_LIBRARY_PATH=$CONDA_PREFIX/lib \
LIBCLANG_PATH=$PREFIX/lib/libclang.so \
RUST_BACKTRACE=1 cargo install --force --root $PREFIX

# Install flight
cd flight/ && pip install . && cd ../

mkdir -p $CONDA_PREFIX/etc/conda/activate.d
touch $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
touch $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh


echo "#!/bin/bash" >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
echo "export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
echo "export OLD_LIBRARY_PATH=$LIBRARY_PATH" >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
echo "export LIBRARY_PATH=$CONDA_PREFIX/lib:${LIBRARY_PATH}" >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
echo "export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:${LD_LIBRARY_PATH}" >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh

echo "#!/bin/bash" >> $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
echo "export LIBRARY_PATH=${OLD_LIBRARY_PATH}" >> $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
echo "export LD_LIBRARY_PATH=${OLD_LD_LIBRARY_PATH}" >> $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
echo "unset OLD_LIBRARY_PATH" >> $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
echo "unset OLD_LD_LIBRARY_PATH" >> $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
# Install spacegraphcats
#cd spacegraphcats/ && pip install . && cd ../

# move Rscript and python
# cp src/bin/snp_density_plots.R $CONDA_PREFIX/bin/
#cp src/bin/cluster.py $CONDA_PREFIX/bin/
