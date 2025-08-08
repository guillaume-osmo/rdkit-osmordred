set -e

echo "Building RDKit from source"

# Clean previous build
rm -rf build
mkdir build
cd build

# Run CMake configuration
cmake .. \
    -DCMAKE_BUILD_TYPE=Debug \
    -DRDK_INSTALL_INTREE=ON \
    -DRDK_INSTALL_STATIC_LIBS=OFF \
    -DRDK_BUILD_CPP_TESTS=ON \
    -DRDK_BUILD_WRAPPER_Chem=ON \
    -DRDK_BUILD_PYTHON_WRAPPERS=ON \
    -DRDK_BUILD_COORDGEN_SUPPORT=ON \
    -DRDK_BUILD_MAEPARSER_SUPPORT=ON \
    -DRDK_OPTIMIZE_POPCNT=ON \
    -DRDK_BUILD_TEST_GZIP=ON \
    -DRDK_BUILD_FREESASA_SUPPORT=ON \
    -DRDK_BUILD_AVALON_SUPPORT=ON \
    -DRDK_BUILD_INCHI_SUPPORT=ON \
    -DRDK_BUILD_YAEHMOP_SUPPORT=ON \
    -DRDK_BUILD_XYZ2MOL_SUPPORT=ON \
    -DRDK_BUILD_CAIRO_SUPPORT=ON \
    -DRDK_BUILD_QT_SUPPORT=OFF \
    -DRDK_BUILD_SWIG_WRAPPERS=OFF \
    -DRDK_SWIG_STATIC=OFF \
    -DRDK_BUILD_THREADSAFE_SSS=ON \
    -DRDK_TEST_MULTITHREADED=ON \
    -DRDK_BUILD_CFFI_LIB=ON \
    -DRDK_BUILD_OSMORDRED_SUPPORT=ON \
    -DPYTHON_EXECUTABLE=/opt/conda/envs/rdkit_build/bin/python \
    -DCMAKE_MAKE_PROGRAM=/usr/bin/make \
    -DCMAKE_C_COMPILER=/usr/bin/gcc \
    -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
    -DPYTHON3_NUMPY_INCLUDE_PATH=/opt/conda/envs/rdkit_build/lib/python3.11/site-packages/numpy/core/include \
    -DCMAKE_CXX_FLAGS="-D_GLIBCXX_USE_CXX11_ABI=1"

# Build
make -j$(nproc)

# Install
make install

# os.environ['RDBASE'] = '/work'
# os.environ['RDK_BUILD'] = '/work/build'
# os.environ['PYTHONPATH'] = f"{os.environ['RDBASE']}:{os.environ['RDK_BUILD']}{os.environ.get('PYTHONPATH', '')}"
# os.environ['LD_LIBRARY_PATH'] = f"{os.environ['RDK_BUILD']}/lib:{os.environ.get('CONDA_PREFIX', '')}/lib:{os.environ.get('LD_LIBRARY_PATH', '')}"  # Linux
